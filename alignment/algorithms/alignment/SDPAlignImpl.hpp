#ifndef _BLASR_SDP_ALIGN_IMPL_HPP_
#define _BLASR_SDP_ALIGN_IMPL_HPP_

#include <cstddef>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <ostream>
// pbdata
#include "../../../pbdata/Types.h"
#include "../../../pbdata/defs.h"
#include "../../../pbdata/utils.hpp"
#include "../../../pbdata/Enumerations.h"
#include "../../../pbdata/DNASequence.hpp"
#include "../../../pbdata/matrix/FlatMatrix.hpp"

#include "../../datastructures/alignment/Alignment.hpp"
#include "GraphPaper.hpp"
#include "AlignmentUtils.hpp"
#include "SWAlign.hpp"
#include "sdp/SDPFragment.hpp"
#include "sdp/SparseDynamicProgramming.hpp"
#include "SDPAlign.hpp"

template<typename T_QuerySequence, typename T_TargetSequence, typename T_ScoreFn>
int SDPAlign(T_QuerySequence &query, T_TargetSequence &target,
        T_ScoreFn &scoreFn, int wordSize, 
        int sdpIns, int sdpDel, float indelRate,
        blasr::Alignment &alignment, 
        AlignmentType alignType,
        bool detailedAlignment,
        bool extendFrontByLocalAlignment,
        DNALength noRecurseUnder,
        bool fastSDP,
        unsigned int minFragmentsToUseGraphPaper) {
    /*
       Since SDP Align uses a large list of buffers, but none are
       provided with this mechanism of calling SDPAlign, allocate the
       buffers on the stack.
       */
    std::vector<Fragment> fragmentSet, prefixFragmentSet, suffixFragmentSet;
    TupleList<PositionDNATuple> targetTupleList;
    TupleList<PositionDNATuple> targetPrefixTupleList;
    TupleList<PositionDNATuple> targetSuffixTupleList;
    std::vector<int> maxFragmentChain;

    return SDPAlign(query, target,
            scoreFn, wordSize, 
            sdpIns, sdpDel, indelRate,
            alignment, 
            fragmentSet, prefixFragmentSet, suffixFragmentSet, 
            targetTupleList, targetPrefixTupleList, targetSuffixTupleList,
            maxFragmentChain,
            alignType,
            detailedAlignment,
            extendFrontByLocalAlignment,
            noRecurseUnder,
            fastSDP,
            minFragmentsToUseGraphPaper);
}

template<typename T_QuerySequence, typename T_TargetSequence, typename T_ScoreFn, typename T_BufferCache>
int SDPAlign(T_QuerySequence &query, T_TargetSequence &target,
        T_ScoreFn &scoreFn, int wordSize, 
        int sdpIns, int sdpDel, float indelRate,
        blasr::Alignment &alignment, 
        T_BufferCache &buffers,
        AlignmentType alignType,
        bool detailedAlignment,
        bool extendFrontByLocalAlignment, 
        DNALength noRecurseUnder,
        bool fastSDP,
        unsigned int minFragmentsToUseGraphPaper) {

    return SDPAlign(query, target, scoreFn, wordSize, 
            sdpIns, sdpDel, indelRate,
            alignment,  
            buffers.sdpFragmentSet,
            buffers.sdpPrefixFragmentSet,
            buffers.sdpSuffixFragmentSet,
            buffers.sdpCachedTargetTupleList,
            buffers.sdpCachedTargetPrefixTupleList,
            buffers.sdpCachedTargetSuffixTupleList,
            buffers.sdpCachedMaxFragmentChain,
            alignType, detailedAlignment, 
            extendFrontByLocalAlignment, noRecurseUnder,
            fastSDP, minFragmentsToUseGraphPaper);
}

template<typename T_QuerySequence, typename T_TargetSequence, typename T_ScoreFn, typename T_TupleList>
int SDPAlign(T_QuerySequence &query, T_TargetSequence &target,
        T_ScoreFn &scoreFn,
        int wordSize, 
        int sdpIns, int sdpDel, float indelRate,
        blasr::Alignment &alignment, 
        std::vector<Fragment> &fragmentSet,
        std::vector<Fragment> &prefixFragmentSet,
        std::vector<Fragment> &suffixFragmentSet,
        T_TupleList &targetTupleList,
        T_TupleList &targetPrefixTupleList,
        T_TupleList &targetSuffixTupleList,
        std::vector<int> &maxFragmentChain,
        // A few optinal parameters, should delete that last one.
        AlignmentType alignType,
        bool detailedAlignment,
        bool extendFrontByLocalAlignment, 
        DNALength noRecurseUnder,
        bool fastSDP,
        unsigned int minFragmentsToUseGraphPaper) {
    // minFragmentsToUseGraphPaper: minimum number of fragments to 
    // use Graph Paper for speed up.

    fragmentSet.clear();
    prefixFragmentSet.clear();
    suffixFragmentSet.clear();
    targetTupleList.clear();
    targetPrefixTupleList.clear();
    targetSuffixTupleList.clear();
    maxFragmentChain.clear();

    //
    //	Collect a set of matching fragments between query and target.
    //  Since this function is an inner-loop for alignment, anything to
    //	speed it up will help.  One way to speed it up is to re-use the
    //	vectors that contain the sdp matches. 
    //

    TupleMetrics tm, tmSmall;
    tm.Initialize(wordSize);

    int smallWordSize = (wordSize < SDP_DETAILED_WORD_SIZE ? wordSize : SDP_DETAILED_WORD_SIZE);
    tmSmall.Initialize(smallWordSize);

    //
    // Partition the read into a prefix, middle, and suffix.  The prefix
    // and suffix are matched using a smaller word size allowing for
    // higher sensitivity at the ends of reads, which are more likely to
    // be misaligned.
    //
    int prefixLength, middleLength, suffixLength, suffixPos; // prefix pos is 0
    prefixLength = min(target.length, (DNALength) SDP_PREFIX_LENGTH);
    suffixLength = min(target.length - prefixLength, (DNALength) SDP_SUFFIX_LENGTH);
    middleLength = target.length - prefixLength - suffixLength;

    DNASequence prefix, middle, suffix;
    DNASequence qPrefix, qMiddle, qSuffix;
    DNALength pos = 0;
    prefix.seq = &target.seq[pos];
    prefix.length = prefixLength;
    pos += prefixLength;

    // Align the entire query against the entire target to get alignments 
    // in the middle.
    middle.seq = &target.seq[0];
    middle.length = target.length;

    pos += middleLength;
    suffixPos = pos; // while = prefixLength + middleLength
    suffix.seq = &target.seq[suffixPos];
    suffix.length = suffixLength;

    int qPrefixLength, qMiddleLength, qSuffixLength, qMiddlePos, qSuffixPos; // prefix pos is 0
    qPrefixLength = min(query.length, (DNALength) SDP_PREFIX_LENGTH);
    qSuffixLength = min(query.length - qPrefixLength, (DNALength) SDP_SUFFIX_LENGTH);
    qMiddleLength = query.length - qPrefixLength - qSuffixLength;

    pos = 0;
    qPrefix.seq = &query.seq[pos];
    qPrefix.length = qPrefixLength;

    // Align the entire query against the entire target to get alignments 
    // in the middle.
    qMiddle.seq = &query.seq[0];
    qMiddle.length = query.length;
    qMiddlePos = pos += qPrefixLength;
    (void)(qMiddlePos);

    qSuffixPos = pos += qMiddleLength; // = qPrefixLength + qMiddleLength
    qSuffix.seq = &query.seq[qSuffixPos];
    qSuffix.length = qSuffixLength;

	fragmentSet.clear();
    SequenceToTupleList(prefix, tmSmall, targetPrefixTupleList);
    SequenceToTupleList(suffix, tmSmall, targetSuffixTupleList);
    SequenceToTupleList(middle, tm, targetTupleList);

    targetPrefixTupleList.Sort();
    targetSuffixTupleList.Sort();
    targetTupleList.Sort();


    //
    // Store in fragmentSet the tuples that match between the target
    // and query.
    //
    StoreMatchingPositions(qPrefix, tmSmall, targetPrefixTupleList, prefixFragmentSet);
    StoreMatchingPositions(qSuffix, tmSmall, targetSuffixTupleList, suffixFragmentSet);
    StoreMatchingPositions(qMiddle, tm, targetTupleList, fragmentSet); 

    // 
    // The method to store matching positions is not weight aware.
    // Store the weight here.
    //
    VectorIndex f;

    for (f = 0; f < suffixFragmentSet.size(); f++) {
        (suffixFragmentSet)[f].weight = tm.tupleSize;
        (suffixFragmentSet)[f].length = tmSmall.tupleSize;
    }
    for (f = 0; f < prefixFragmentSet.size(); f++) {
        (prefixFragmentSet)[f].weight = tm.tupleSize;
        (prefixFragmentSet)[f].length = tmSmall.tupleSize;
    }
    for (f = 0; f < fragmentSet.size(); f++) {
        (fragmentSet)[f].weight = tm.tupleSize;
        (fragmentSet)[f].length = tm.tupleSize;
    }

    //
    // Since different partitions of the read are matched, the locations
    // of the matches do not have the correct position because of the
    // offsets.  Fix that here.
    for (f = 0; f < suffixFragmentSet.size(); f++) {
        (suffixFragmentSet)[f].x += qSuffixPos;
        (suffixFragmentSet)[f].y += suffixPos;
    }

    //
    // Collect all fragments into one.
    //
    fragmentSet.insert(fragmentSet.begin(), prefixFragmentSet.begin(), prefixFragmentSet.end());
    fragmentSet.insert(fragmentSet.end(), suffixFragmentSet.begin(), suffixFragmentSet.end());

    FlatMatrix2D<int> graphScoreMat;
    FlatMatrix2D<Arrow> graphPathMat;
    FlatMatrix2D<int> graphBins;
    if (fragmentSet.size() > minFragmentsToUseGraphPaper and fastSDP) {
        int nCol = 50;
        vector<bool> onOptPath(fragmentSet.size(), false);
        GraphPaper<Fragment>(fragmentSet, nCol, nCol,
                                      graphBins, graphScoreMat, 
                                      graphPathMat, onOptPath);
        RemoveOffOpt(fragmentSet, onOptPath);
    } 
    graphScoreMat.Clear();
    graphPathMat.Clear();
    graphBins.Clear();

    //
    // Because there are fragments from multiple overlapping regions, remove
    // any fragments that have the same starting coordinate as a
    // previous fragment.
    //
    std::sort(fragmentSet.begin(), fragmentSet.end(), LexicographicFragmentSort<Fragment>());
    f = 0;
    int fCur = 0;
    while (f + 1 <= fragmentSet.size()) {
        fragmentSet[fCur] = fragmentSet[f];
        while (f < fragmentSet.size() and fragmentSet[fCur].x == fragmentSet[f].x and fragmentSet[fCur].y == fragmentSet[f].y) {
            f++;
        }
        fCur++;
    }
    fragmentSet.resize(fCur);
  
	if (fragmentSet.size() == 0) {
		//
		// This requires at least one seeded tuple to begin an alignment.
		//
		return 0;
	}

    //
    // Find the longest chain of anchors.
    //
    SDPLongestCommonSubsequence(query.length, fragmentSet, tm.tupleSize, sdpIns, sdpDel, scoreFn.scoreMatrix[0][0], maxFragmentChain, alignType);

    //
    // Now turn the max fragment chain into a real alignment.
    //
    int startF;
    blasr::Alignment chainAlignment;
    alignment.qPos = 0;
    alignment.tPos = 0;
    Block block;
    std::vector<int> fragScoreMat;
    std::vector<Arrow> fragPathMat;

    //
    // Patch the sdp fragments into an alignment, possibly breaking the
    // alignment if the gap between two fragments is too large.
    //
    for (f = 0; f < maxFragmentChain.size(); f++ ){ 
        startF = f;
        // Condense contiguous stretches.
        while(f < maxFragmentChain.size()  - 1 and 
              fragmentSet[maxFragmentChain[f]].x == fragmentSet[maxFragmentChain[f+1]].x - 1 and
              fragmentSet[maxFragmentChain[f]].y == fragmentSet[maxFragmentChain[f+1]].y - 1) {
            f++;
        }

        block.qPos = fragmentSet[maxFragmentChain[startF]].x;
        block.tPos = fragmentSet[maxFragmentChain[startF]].y;

        // Compute the block length as the difference between the starting
        // point of the current block and ending point of the last
        // overlapping block.  This was previously calculated by adding
        // the number of merged blocks - 1 to the first block length.
        // When the block lengths are heterogenous, this does not work,
        // and it would be possible to have a block with a length that
        // extends past the end of a sequence.  By taking the length as
        // the difference here, it ensures this will not happen.
        //
        block.length = fragmentSet[maxFragmentChain[f]].x + fragmentSet[maxFragmentChain[f]].length - fragmentSet[maxFragmentChain[startF]].x;
		chainAlignment.blocks.push_back(block);
	}

	//
	// It may be possible that in regions of low similarity, spurious matches fit into the LCS.  
	// Assume that indels cause the matches to diverge from the diagonal on a random walk.  If they 
	// walk more than 3 standard deviations away from the diagonal, they are probably spurious. 
	//
	unsigned int b;
	chainAlignment.qPos = 0;
	chainAlignment.tPos = 0;

    for (b = 0; b < chainAlignment.size()-1; b++){ 
        if (chainAlignment.blocks[b].qPos + chainAlignment.blocks[b].length > chainAlignment.blocks[b+1].qPos) {
            chainAlignment.blocks[b].length = (chainAlignment.blocks[b+1].qPos - chainAlignment.blocks[b].qPos);
        }
        if (chainAlignment.blocks[b].tPos + chainAlignment.blocks[b].length > chainAlignment.blocks[b+1].tPos) {
            chainAlignment.blocks[b].length = (chainAlignment.blocks[b+1].tPos - chainAlignment.blocks[b].tPos);
        }
        // the min indel rate between the two chain blocks is the difference in diagonals between the two sequences.
        int curDiag, nextDiag, diffDiag;
        curDiag = chainAlignment.blocks[b].tPos - chainAlignment.blocks[b].qPos;
        nextDiag = chainAlignment.blocks[b+1].tPos - chainAlignment.blocks[b+1].qPos;
        diffDiag = std::abs(curDiag - nextDiag);

        //
        // It is expected that the deviation is at least 1, so discount for this
        //
        diffDiag--;
        // compare the alignment distances.  
    }

    std::vector<bool> blockIsGood;
    blockIsGood.resize(chainAlignment.size());
    fill(blockIsGood.begin(), blockIsGood.end(), true);

    //
    // The hack that allows anchors of different lengths at the front
    // and end of alignments (to increase sensitivity at the ends of
    // sequences) has the side effect that there may be blocks that have
    // zero length.  This shouldn't happen, so to balance this out
    // remove blocks that have zero length.
    //
    for (b = 0; b < chainAlignment.size(); b++){ 
        if (chainAlignment.blocks[b].length == 0) {
            blockIsGood[b] = false;
        }
    }
    for (b = 1; b < chainAlignment.size()-1; b++){ 
        // the min indel rate between the two chain blocks is the difference in diagonals between the two sequences.
        int prevDiag = abs(((int)chainAlignment.blocks[b].tPos -   (int)chainAlignment.blocks[b].qPos)  -
                ((int)chainAlignment.blocks[b-1].tPos - (int)chainAlignment.blocks[b-1].qPos));

        int prevDist = std::min(chainAlignment.blocks[b].tPos - chainAlignment.blocks[b-1].tPos,
                chainAlignment.blocks[b].qPos - chainAlignment.blocks[b-1].qPos);

        int nextDiag = abs(((int)chainAlignment.blocks[b+1].tPos - (int)chainAlignment.blocks[b+1].qPos)  -
                ((int)chainAlignment.blocks[b].tPos -   (int)chainAlignment.blocks[b].qPos));

        int nextDist = std::min(chainAlignment.blocks[b+1].tPos - chainAlignment.blocks[b].tPos,
                chainAlignment.blocks[b+1].qPos - chainAlignment.blocks[b].qPos);

        if (prevDist * indelRate < prevDiag and nextDist * indelRate < nextDiag) {
            blockIsGood[b] = false;
        }
    }

    for (b = chainAlignment.size(); b > 0; b--) {
        if (blockIsGood[b-1] == false) {
            chainAlignment.blocks.erase(chainAlignment.blocks.begin() + b-1);
        }
    }

    if (chainAlignment.blocks.size() > 0) {
        T_QuerySequence  qFragment;
        T_TargetSequence tFragment;
        blasr::Alignment fragAlignment;
        unsigned int fb;
        if (alignType == Global) {
            //
            // For Global alignment, refine the alignment from the beginnings of the
            // sequences to the start of the first block.
            //
            if (chainAlignment.blocks[0].qPos > 0 and
                chainAlignment.blocks[0].tPos > 0) {

                qFragment.seq = &query.seq[0];
                qFragment.length = chainAlignment.blocks[0].qPos;
                tFragment.seq = &target.seq[0];
                tFragment.length = chainAlignment.blocks[0].tPos;
                for (fb = 0; fb < alignment.blocks.size(); fb++) {
                    alignment.blocks.push_back(fragAlignment.blocks[b]);
                }
            }
        }
        else if (alignType == Local) {
            // Perform a front-anchored alignment to extend the alignment to
            // the beginning of the read.
            if (chainAlignment.blocks[0].qPos > 0 and 
                    chainAlignment.blocks[0].tPos > 0) {
                qFragment.seq = (Nucleotide*) &query.seq[0];
                qFragment.length = chainAlignment.blocks[0].qPos;

                tFragment.seq = (Nucleotide*) &target.seq[0];
                tFragment.length = chainAlignment.blocks[0].tPos;
                blasr::Alignment frontAlignment;
                // Currently, there might be some space between the beginning
                // of the alignment and the beginning of the read.  Run an
                // EndAnchored alignment that allows free gaps to the start of
                // where the alignment begins, but normal, ungapped alignment
                // otherwise. 
                if (extendFrontByLocalAlignment) {
                    if (noRecurseUnder == 0 or qFragment.length * tFragment.length < noRecurseUnder) {
                    SWAlign(qFragment, tFragment, fragScoreMat, fragPathMat, frontAlignment, scoreFn, EndAnchored);
                    } else {
                        // cout << "running recursive sdp alignment. " << endl;
                        vector<int> recurseFragmentChain;
                        SDPAlign(qFragment, tFragment, scoreFn,
                                std::max(wordSize/2, 5),
                                sdpIns, sdpDel,  indelRate,
                                frontAlignment,
                                fragmentSet,
                                prefixFragmentSet,
                                suffixFragmentSet,
                                targetTupleList,
                                targetPrefixTupleList,
                                targetSuffixTupleList,
                                recurseFragmentChain,
                                alignType, detailedAlignment, 
                                extendFrontByLocalAlignment, 0);
                    }

                    unsigned int anchorBlock;
                    for (anchorBlock = 0; anchorBlock < frontAlignment.blocks.size(); anchorBlock++) {
                        //
                        // The front alignment needs to be transformed to the
                        // coordinate offsets that the chain alignment is in.  This
                        // is an alignment starting at position 0 in the target and
                        // query.  Currently, the front alignment is offset into the
                        // sequences by frontAlignment.[q/t]Pos.
                        //
                        frontAlignment.blocks[anchorBlock].tPos += frontAlignment.tPos;
                        frontAlignment.blocks[anchorBlock].qPos += frontAlignment.qPos;
                        alignment.blocks.push_back(frontAlignment.blocks[anchorBlock]);
                    }
                }
            }
        }

		// 
		// The chain alignment blocks are not complete blocks, so they
		// must be appended to the true alignment and then patched up.
		//
        for (b = 0; b < chainAlignment.size() - 1; b++) {
            alignment.blocks.push_back(chainAlignment.blocks[b]);

            //
            // Do a detaied smith-waterman alignment between blocks, if this
            // is specified.  
            fragAlignment.Clear();
            qFragment.Free();
            qFragment.ReferenceSubstring(query, chainAlignment.blocks[b].qPos + chainAlignment.blocks[b].length);
            qFragment.length = chainAlignment.blocks[b+1].qPos - 
                (chainAlignment.blocks[b].qPos + chainAlignment.blocks[b].length);

            tFragment.seq    = &(target.seq[chainAlignment.blocks[b].tPos + chainAlignment.blocks[b].length]);
            tFragment.length = (chainAlignment.blocks[b+1].tPos - 
                    (chainAlignment.blocks[b].tPos + chainAlignment.blocks[b].length));

            if (qFragment.length > 0 and 
                    tFragment.length > 0 and
                    detailedAlignment == true) {

                if (noRecurseUnder == 0 or qFragment.length * tFragment.length < noRecurseUnder) {
                    SWAlign(qFragment, tFragment, fragScoreMat, fragPathMat, fragAlignment, scoreFn, Global);
                }
                else {
                    //          cout << "running recursive sdp alignment on " << qFragment.length * tFragment.length << endl;
                    std::vector<int> recurseFragmentChain;
                    SDPAlign(qFragment, tFragment, scoreFn,
                            std::max(wordSize/2, 5),
                            sdpIns, sdpDel,  indelRate,
                            fragAlignment,
                            fragmentSet,
                            prefixFragmentSet,
                            suffixFragmentSet,
                            targetTupleList,
                            targetPrefixTupleList,
                            targetSuffixTupleList,
                            recurseFragmentChain,
                            alignType, detailedAlignment, 0, 0);
                }
                fragAlignment.qPos = 0;
                fragAlignment.tPos = 0;

                int qOffset = chainAlignment.blocks[b].qPos + chainAlignment.blocks[b].length;
                int tOffset = chainAlignment.blocks[b].tPos + chainAlignment.blocks[b].length;

                for (fb = 0; fb < fragAlignment.blocks.size(); fb++) {
                    fragAlignment.blocks[fb].qPos += qOffset;
                    fragAlignment.blocks[fb].tPos += tOffset;
                    alignment.blocks.push_back(fragAlignment.blocks[fb]);
                }
            }
        }
        int lastBlock = chainAlignment.blocks.size() - 1;
        if (alignType == Global or alignType == Local) {
            if (chainAlignment.size() > 0) {
                // Add the last block.
                alignment.blocks.push_back(chainAlignment.blocks[lastBlock]);
                if (alignType == Global and detailedAlignment == true) {
                    //
                    // When doing a global alignment, the sequence from the end of
                    // the last block of the query should be aligned to the end of 
                    // the text.
                    //
                    qFragment.Free();
                    qFragment.ReferenceSubstring(query, chainAlignment.blocks[lastBlock].qPos + chainAlignment.blocks[lastBlock].length,
                            query.length -  
                            (chainAlignment.blocks[lastBlock].qPos + chainAlignment.blocks[lastBlock].length));

                    tFragment.seq    = &(target.seq[chainAlignment.blocks[lastBlock].tPos + chainAlignment.blocks[lastBlock].length]);
                    tFragment.length = (target.length - 
                            (chainAlignment.blocks[lastBlock].tPos + chainAlignment.blocks[lastBlock].length));
                    if (qFragment.length > 0 and
                            tFragment.length > 0 ) {

                        if (extendFrontByLocalAlignment) {

                            fragAlignment.Clear();
                            if (qFragment.length * tFragment.length > 10000) {
                                //              cout << "Cautin: slow alignment crossing! " << qFragment.length  << " " << tFragment.length << endl;
                            }

                            if (noRecurseUnder == 0 or qFragment.length * tFragment.length < noRecurseUnder) {
                                SWAlign(qFragment, tFragment, fragScoreMat, fragPathMat, fragAlignment, scoreFn, EndAnchored);
                            }
                            else {
                                std::vector<int> recurseFragmentChain;
                                SDPAlign(qFragment, tFragment, scoreFn,
                                        std::max(wordSize/2, 5),
                                        sdpIns, sdpDel,  indelRate,
                                        fragAlignment,
                                        fragmentSet,
                                        prefixFragmentSet,
                                        suffixFragmentSet,
                                        targetTupleList,
                                        targetPrefixTupleList,
                                        targetSuffixTupleList,
                                        recurseFragmentChain,
                                        alignType, detailedAlignment, extendFrontByLocalAlignment, 0);
                            }



                            int qOffset = chainAlignment.blocks[lastBlock].qPos + chainAlignment.blocks[lastBlock].length;
                            int tOffset = chainAlignment.blocks[lastBlock].tPos + chainAlignment.blocks[lastBlock].length;
                            unsigned int fb;
                            for (fb = 0; fb < fragAlignment.size(); fb++) { 
                                fragAlignment.blocks[fb].qPos += qOffset;
                                fragAlignment.blocks[fb].tPos += tOffset;
                                alignment.blocks.push_back(fragAlignment.blocks[fb]);
                            }
                        }
                    }
                }
            }
        }
    }

    if (alignType == Local) {
        alignment.tPos = alignment.blocks[0].tPos;
        alignment.qPos = alignment.blocks[0].qPos;
        VectorIndex b;
        for (b = 0; b < alignment.blocks.size(); b++) { 
            alignment.blocks[b].qPos -= alignment.qPos;
            alignment.blocks[b].tPos -= alignment.tPos;
        }
    }
    int alignmentScore;
    /*	ofstream queryOut("query.fasta");
        FASTASequence tmp;
        ((DNASequence&)tmp).Copy(query);
        tmp.CopyTitle("query");
        tmp.PrintSeq(queryOut);
        queryOut.close();
        ofstream targetOut("target.fasta");
        ((DNASequence&)tmp).Copy(target);
        tmp.CopyTitle("target");
        tmp.PrintSeq(targetOut);
        targetOut.close();
        */
    alignmentScore = ComputeAlignmentScore(alignment, query, target, scoreFn);
    return alignmentScore;
}

#endif //_BLASR_SDP_ALIGN_IMPL_HPP_

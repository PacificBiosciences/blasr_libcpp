#include "SAMPrinter.hpp"
#include <algorithm> //reverse

using namespace SAMOutput; 

void SAMOutput::BuildFlag(T_AlignmentCandidate &alignment, AlignmentContext &context, uint16_t &flag) {

    /*
     *  Much of the flags are commented out for now since they do not
     *  generate PICARD compliant SAM.  This needs to be worked out. 
     */


    //
    // Without supporting strobe, assume 1 segment per template.
    flag = 0;
    /*
       if (context.nSubreads > 1) {
       flag |= MULTI_SEGMENTS;
       }*/

    //    if (context.AllSubreadsAligned() and context.nSubreads > 1) {
    //      flag |= ALL_SEGMENTS_ALIGNED;
    //    }

    if (alignment.tStrand == 1) {
        flag |= SEQ_REVERSED;
    }
    /*
       if (context.hasNextSubreadPos == false and context.nSubreads > 1) {
       flag |= NEXT_SEGMENT_UNMAPPED;
       }
       if (context.nextSubreadDir == 1) {
       flag |= SEQ_NEXT_REVERSED;
       }
       if (context.IsFirst() and context.nSubreads > 1) {
       flag |= FIRST_SEGMENT;
       }
       */
       else if (context.nSubreads > 1) {
           /*
            * Remember, if you're not first, you're last.
            */
           //      flag |= LAST_SEGMENT;
       }
       if (context.isPrimary == false) {
           flag |= SECONDARY_ALIGNMENT;
       }
}

void SAMOutput::CreateDNAString(DNASequence &seq, 
        DNASequence &clippedSeq,
        //
        // Trimming is used for both hard non-clipping
        // so it is called trim instead of clip.
        //
        int trimFront,
        int trimEnd) {
    assert(seq.length - trimEnd >= trimFront);

    clippedSeq.seq    = &seq.seq[trimFront];
    clippedSeq.length = seq.length - trimEnd - trimFront;
}

void SAMOutput::AddGaps(T_AlignmentCandidate &alignment, int gapIndex,
        std::vector<int> &opSize, std::vector<char> &opChar) {
    int g;
    for (g = 0; g < alignment.gaps[gapIndex].size(); g++) {
        if (alignment.gaps[gapIndex][g].seq == blasr::Gap::Query) {
            opSize.push_back(alignment.gaps[gapIndex][g].length);
            opChar.push_back('D');
        }
        else if (alignment.gaps[gapIndex][g].seq == blasr::Gap::Target) {
            opSize.push_back(alignment.gaps[gapIndex][g].length);
            opChar.push_back('I');
        }
    }
}

void SAMOutput::AddMatchBlockCigarOps(DNASequence & qSeq, DNASequence & tSeq, 
        blasr::Block & b, DNALength & qSeqPos, DNALength & tSeqPos,
        std::vector<int> & opSize, std::vector<char> & opChar) {
    DNALength qPos = qSeqPos + b.qPos, tPos = tSeqPos + b.tPos, n = 0;
    bool started = false, prevSeqMatch = false;
    for(DNALength i = 0; i < b.length; i++) {
        bool curSeqMatch = (qSeq[qPos + i] == tSeq[tPos + i]);
        if (started) {
            if (curSeqMatch == prevSeqMatch) opSize[opSize.size()-1]++;
            else {
                opSize.push_back(1);
                opChar.push_back(curSeqMatch?'=':'X');
            }
        } else {
            started = true;
            opSize.push_back(1);
            opChar.push_back(curSeqMatch?'=':'X');
        }
        prevSeqMatch = curSeqMatch;
    }
}

void SAMOutput::MergeAdjacentIndels(std::vector<int> &opSize, 
                                    std::vector<char> &opChar,
                                    const char mismatchChar) {
    assert(opSize.size() == opChar.size());
    int i, j;
    for (i = 0, j = 1; i < opSize.size() and j < opSize.size(); j++) {
        const int  ni = opSize[i], nj = opSize[j];
        const char ci = opChar[i], cj = opChar[j];
        if (ci == cj) {
            opSize[i] = ni + nj; // merge i and j to i
        } else {
            if ((ci == 'I' and cj == 'D') or (ci == 'D' and cj == 'I')) {
                opSize[i] = std::min(ni, nj);
                opChar[i] = mismatchChar; // merge 'I' and 'D' to Mismatch
                if (i != 0 and i != opSize.size() and opChar[i] == opChar[i - 1]) {
                    opSize[i - 1] += opSize[i]; // merge with i - 1?
                    i--;
                }
                if (ni != nj) {
                    i++; // the remaining 'I' or 'D'
                    opSize[i] = std::abs(ni - nj);
                    opChar[i] = (ni > nj) ? ci: cj;
                }
            } else  {
                i++; // move forward.
                opSize[i] = nj;
                opChar[i] = cj;
            }
        }
    }
    opSize.erase(opSize.begin() + i + 1, opSize.end());
    opChar.erase(opChar.begin() + i + 1, opChar.end());
}

void SAMOutput::CreateNoClippingCigarOps(T_AlignmentCandidate &alignment, 
        std::vector<int> &opSize, std::vector<char> &opChar,
        bool cigarUseSeqMatch, const bool allowAdjacentIndels) {
    //
    // Create the cigar string for the aligned region of a read,
    // excluding the clipping.
    //
    int b;
    // Each block creates a match NM (N=block.length)
    int nBlocks = alignment.blocks.size();
    int nGaps   = alignment.gaps.size();
    opSize.clear();
    opChar.clear();
    //
    // Add gaps at the beginning of the alignment.
    //
    if (nGaps > 0) {
        AddGaps(alignment, 0, opSize, opChar);
    }
    for (b = 0; b < nBlocks; b++) {
        //
        // Determine the gap before printing the match, since it is
        // possible that the qurey and target are gapped at the same
        // time, which merges into a mismatch.
        //
        int qGap=0, tGap=0, commonGap=0;
        int matchLength = alignment.blocks[b].length;
        if (nGaps == 0) {
            if (b + 1 < nBlocks) {
                qGap = alignment.blocks[b+1].qPos - alignment.blocks[b].qPos - alignment.blocks[b].length;
                tGap = alignment.blocks[b+1].tPos - alignment.blocks[b].tPos - alignment.blocks[b].length;
                int commonGap;
                commonGap = abs(qGap - tGap);
                qGap -= commonGap;
                tGap -= commonGap;
                matchLength += commonGap;
                if (cigarUseSeqMatch) {
                    AddMatchBlockCigarOps(alignment.qAlignedSeq, 
                                          alignment.tAlignedSeq, 
                                          alignment.blocks[b], 
                                          alignment.qPos, alignment.tPos,
                                          opSize, opChar);
                } else {
                    opSize.push_back(matchLength);
                    opChar.push_back('M');
                }
                assert((qGap > 0 and tGap == 0) or (qGap == 0 and tGap > 0));
                if (qGap > 0) {
                    opSize.push_back(qGap);
                    opChar.push_back('I');
                }
                if (tGap > 0) {
                    opSize.push_back(tGap);
                    opChar.push_back('D');
                }
            }
        }
        else {
            if (cigarUseSeqMatch) {
                AddMatchBlockCigarOps(alignment.qAlignedSeq, 
                                      alignment.tAlignedSeq, 
                                      alignment.blocks[b], 
                                      alignment.qPos, alignment.tPos,
                                      opSize, opChar);
            } else {
                opSize.push_back(matchLength);
                opChar.push_back('M');
            }
            int gapIndex = b+1;
            AddGaps(alignment, gapIndex, opSize, opChar);
        }
    }
    if (alignment.tStrand == 1) {
        std::reverse(opSize.begin(), opSize.end());
        std::reverse(opChar.begin(), opChar.end());
    }

    if (not allowAdjacentIndels) {
        MergeAdjacentIndels(opSize, opChar, (cigarUseSeqMatch?'X':'M'));
    }
}

void SAMOutput::CigarOpsToString(std::vector<int> &opSize, std::vector<char> &opChar, 
        std::string &cigarString) {
    std::stringstream sstrm;
    int i, nElem;
    for (i = 0, nElem = opSize.size(); i < nElem; i++) {
        sstrm << opSize[i] << opChar[i];
    }
    cigarString = sstrm.str();
}


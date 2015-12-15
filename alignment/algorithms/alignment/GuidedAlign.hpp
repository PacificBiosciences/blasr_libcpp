#ifndef _BLASR_GUIDE_ALIGNMENT_HPP_
#define _BLASR_GUIDE_ALIGNMENT_HPP_

#include <vector>
#include <cmath>
#include <limits.h>
#include <ostream>
#include "Types.h"
#include "defs.h"
#include "NucConversion.hpp"
#include "DNASequence.hpp"
#include "sdp/SDPFragment.hpp"
#include "matrix/FlatMatrix.hpp"
#include "datastructures/alignment/Alignment.hpp"
#include "datastructures/anchoring/MatchPos.hpp"
#include "matrix/Matrix.hpp"
#include "tuples/TupleList.hpp"
#include "tuples/TupleMetrics.hpp"
#include "tuples/DNATuple.hpp"
#include "utils/LogUtils.hpp"
#include "utils/PhredUtils.hpp"
#include "qvs/QualityValue.hpp"
#include "qvs/QualityValueVector.hpp"
#include "AlignmentUtils.hpp"
#include "DistanceMatrixScoreFunction.hpp"
#include "algorithms/alignment/SDPAlign.hpp"

#define LOWEST_LOG_VALUE  -700

#define MAX_BAND_SIZE 250

//FIXME: Change data type of target pos to GenomeLength in order to support
//       > 4G genome.

class GuideRow {
public:
    int q, t;
    int tPre, tPost;
    unsigned int matrixOffset; // Where the center (q) is in the score
    // and path matrices.
    int GetRowLength(); 
};

typedef std::vector<GuideRow> Guide;

class GetBufferIndexFunctor {
// row + rowSeqOffset - 1 == seqRow
// so, rowInSeq = -1, rowSeqOffset = 0 -> row = 0
public:
    int seqRowOffset;
    int guideSize;
    int operator()(Guide &guide, int seqRow, int seqCol, int &index); 
};

int ComputeMatrixNElem(Guide &guide);

void StoreMatrixOffsets(Guide &guide); 

float QVToLogPScale(char qv);

void QVToLogPScale(QualityValueVector<QualityValue> &qualVect, size_t phredVectLength, std::vector<float> &lnVect); 

int AlignmentToGuide(blasr::Alignment &alignment, Guide &guide, int bandSize); 

template<typename QSequence, typename TSequence, typename T_ScoreFn>
int GuidedAlign(QSequence &origQSeq, TSequence &origTSeq,  blasr::Alignment &guideAlignment,
        T_ScoreFn &scoreFn,
        int bandSize,
        blasr::Alignment &alignment,
        std::vector<int>    &scoreMat,
        std::vector<Arrow>  &pathMat,
        std::vector<double> &probMat,
        std::vector<double> &optPathProbMat,
        std::vector<float>  &lnSubPValueVect,
        std::vector<float>  &lnInsPValueVect,
        std::vector<float>  &lnDelPValueVect,
        std::vector<float>  &lnMatchPValueVect,
        AlignmentType alignType=Global, 
        bool computeProb=false) {

    Guide guide;
    AlignmentToGuide(guideAlignment, guide, bandSize);
    StoreMatrixOffsets(guide);
    ComputeMatrixNElem(guide);

    //
    // Make a copy of the sequences that is guaranteed to be in 3-bit format for faster alignment.
    // (mabybe eventually reuse the qseq and tseq memory)
    //
    QSequence qSeq;
    TSequence tSeq;
    qSeq.Assign(origQSeq);
    tSeq.Assign(origTSeq);

    unsigned int matrixNElem = ComputeMatrixNElem(guide);
    assert(matrixNElem >= 0);
    StoreMatrixOffsets(guide);

    /*
     * The following code is useful to produce images of the dp-matrix.
     * Make sure the sequenences are less than 5kb each though.

     Matrix<float> probMatrix;
     probMatrix.Resize(qSeq.length, tSeq.length);
     probMatrix.Initialize(0);
     ofstream matrixOut;
     stringstream matrixOutNameStrm;
     matrixOutNameStrm << "probMatrix_"<< runIndex << ".dat";
     matrixOut.open(matrixOutNameStrm.str().c_str());
     */

    if (computeProb) { 
        //
        // Convert phred scale to proper ln for faster manipulation later on.
        //
        QVToLogPScale(qSeq.substitutionQV, qSeq.length, lnSubPValueVect);
        QVToLogPScale(qSeq.insertionQV,    qSeq.length, lnInsPValueVect);
        QVToLogPScale(qSeq.deletionQV,     qSeq.length, lnDelPValueVect);
        if (lnMatchPValueVect.size() < qSeq.length) {
            lnMatchPValueVect.resize(qSeq.length);
        }
        //
        // Normalize probability vectors so that the probability of transition from each cell is 1.
        //
        unsigned int i;
        for (i = 0; i < qSeq.length; i++) {
            float subSum  = LogSumOfTwo(lnSubPValueVect[i], QVToLogPScale(scoreFn.substitutionPrior)); // prior on substitution rate
            float denominator = LogSumOfThree(lnDelPValueVect[i], lnInsPValueVect[i], subSum);
            lnDelPValueVect[i] = lnDelPValueVect[i] - denominator;
            lnSubPValueVect[i]   = lnSubPValueVect[i] - denominator;
            lnInsPValueVect[i]   = lnInsPValueVect[i] - denominator;
            lnMatchPValueVect[i] = subSum - denominator;
        }

    }



    // 
    // Make sure the alignments can fit in the reused buffers.
    //

    if (scoreMat.size() < matrixNElem) {
        scoreMat.resize(matrixNElem);
        pathMat.resize(matrixNElem);
        fill(scoreMat.begin(), scoreMat.end(), 0);
        fill(pathMat.begin(), pathMat.end(), NoArrow);
    }
    if (computeProb) {
        if (probMat.size() < matrixNElem) {
            probMat.resize(matrixNElem);
            optPathProbMat.resize(matrixNElem);
        }
    }

    // 
    // Initialze matrices.  Only initialize up to matrixNElem rather
    // than matrix.size() because the matrix.size() unnecessary space
    // may be allocated.
    //

    std::fill(scoreMat.begin(), scoreMat.begin() + matrixNElem, 0);
    std::fill(pathMat.begin(), pathMat.begin() + matrixNElem, NoArrow);
    if (computeProb) {
        std::fill(probMat.begin(), probMat.begin() + matrixNElem, 0);	
        std::fill(optPathProbMat.begin(), optPathProbMat.begin() + matrixNElem, 0);
    }
    //
    // Initialize boundary conditions.
    //
    int q, t;
    int bufferIndex;
    // start alignemnt at the beginning of the guide, and align to the
    // end of the guide.
    if (guide.size() == 0) {
        qSeq.Free();
        tSeq.Free();
        return 0;
    }
    int qStart = guide[1].q;
    int tStart = guide[1].t;
    int qEnd   = guide[guide.size()-1].q+1;
    int tEnd   = guide[guide.size()-1].t+1;

    GetBufferIndexFunctor GetBufferIndex;
    GetBufferIndex.seqRowOffset = qStart;
    GetBufferIndex.guideSize    = guide.size();
    int indicesAreValid, delIndexIsValid, insIndexIsValid;
    bufferIndex = -1;
    indicesAreValid = GetBufferIndex(guide, qStart-1, tStart-1, bufferIndex);
    assert(indicesAreValid);
    scoreMat[bufferIndex] = 0;
    pathMat[bufferIndex]  = NoArrow;
    int matchIndex, insIndex, delIndex, curIndex;

    //
    // Initialize deletion row.
    //
    if (computeProb) {
        probMat[0] = optPathProbMat[0] = 0;
    }
    for (t = tStart; t < tStart + guide[0].tPost; t++) {
        curIndex=-1;
        indicesAreValid = GetBufferIndex(guide, qStart-1, t, curIndex);
        if (indicesAreValid == 0 ) {
            std::cout << "QSeq" << std::endl;
            (static_cast<DNASequence*>(&origQSeq))->PrintSeq(std::cout);
            std::cout << "TSeq" << std::endl;
            (static_cast<DNASequence*>(&origTSeq))->PrintSeq(std::cout);
            assert(0);
        }
        delIndex = -1;
        delIndexIsValid = GetBufferIndex(guide, qStart-1, t-1, delIndex);

        if (delIndexIsValid) {
            if (alignType == Global) {
                scoreMat[curIndex] = scoreMat[delIndex] + scoreFn.del;
            }
            else if (alignType == Local) {
                scoreMat[curIndex] = 0;
            }
            pathMat[curIndex] = Left;
            if (computeProb) {
                if (qSeq.qual.Empty() == false) {
                    optPathProbMat[curIndex] = probMat[curIndex] = probMat[delIndex] + QVToLogPScale(scoreFn.globalDeletionPrior);
                }
            }
        }
    }

    // 
    // Initialize stripe along the top of the grid.
    //
    for (q = qStart ; q < qStart + bandSize and q < qEnd; q++) {
        insIndex = -1;
        insIndexIsValid = GetBufferIndex(guide, q-1, tStart-1, // diagonal from t-start
                insIndex);
        curIndex = -1;
        indicesAreValid = GetBufferIndex(guide, q, tStart-1, curIndex);

        if (insIndexIsValid and indicesAreValid) {
            assert(insIndex >= 0);
            assert(curIndex >= 0);
            if (alignType == Global) {
                scoreMat[curIndex] = scoreMat[insIndex] + scoreFn.ins;
            }
            else {
                scoreMat[curIndex] = 0;
            }
            pathMat[curIndex] = Up;
            if (computeProb) {
                if (qSeq.qual.Empty() == false) {
                    optPathProbMat[curIndex] = probMat[curIndex] = probMat[insIndex] + QVToLogPScale(scoreFn.Insertion(tSeq,(DNALength) 0, qSeq, (DNALength)q));
                }
            }
        }
    }

    int matchScore, insScore, delScore;

    for (q = qStart; q < qEnd; q++) {
        int qi = q - qStart + 1;
        int tp = guide[qi].t;
        curIndex = matchIndex = insIndex = delIndex = -1;

        //
        // Do some work that will help define when matchIndex and insIndex
        // may be used.  Once delIndex is computed once, it is valid for
        // all t positions.
        //
        int prevRowTEnd = -1;
        if ( qi > 0) { 
            //
            // Define the boundaries of the column which may access previously
            // computed cells with a match.
            //
            prevRowTEnd = guide[qi-1].t + guide[qi-1].tPost;
        }    

        for (t = tp - guide[qi].tPre ; t < guide[qi].t + guide[qi].tPost +1; t++) {


            if (q < qStart + bandSize and t == tp - guide[qi].tPre - 1) {
                // On the boundary condition, don't access the 1st element;
                t++;
                continue;
            }
            // Make sure the index is not past the end of the sequence.
            if (t < -1) continue;
            if (t >= tEnd) continue;

            //
            // No cells are available to use for insertion cost
            // computation. 
            //
            if (t > prevRowTEnd) {
                insIndex = -1;
            }
            if (t > prevRowTEnd + 1) {
                matchIndex = -1;
            }

            //
            // Find the indices in the buffer.  Since the rows are of
            // different sizes, one can't just use offsets from the buffer
            // index. 
            //

            if (GetBufferIndex(guide, q-1,t-1, matchIndex)) {
                assert(matchIndex >= 0);
                matchScore = scoreMat[matchIndex] + scoreFn.Match(tSeq, t, qSeq, q);
            }
            else {
                matchScore = INF_INT;
            }

            if (GetBufferIndex(guide, q-1, t, insIndex)) {
                assert(insIndex >= 0);
                insScore = scoreMat[insIndex] + scoreFn.Insertion(tSeq,(DNALength) t, qSeq, (DNALength)q);
            }
            else {
                insScore = INF_INT;
            }
            if (GetBufferIndex(guide, q, t-1, delIndex)) {
                assert(delIndex >= 0);
                delScore = scoreMat[delIndex] + scoreFn.Deletion(tSeq, (DNALength) t, qSeq, (DNALength)q);
            }
            else {
                delScore = INF_INT;
            }

            int minScore = MIN(matchScore, MIN(insScore, delScore));
            int result   = GetBufferIndex(guide, q, t, curIndex);
            // This should only loop over valid cells.
            assert(result);
            assert(curIndex >= 0);
            scoreMat[curIndex] = minScore;
            if (minScore == INF_INT) {
                pathMat[curIndex] = NoArrow;
                if (computeProb) {
                    probMat[curIndex] = 1;
                    optPathProbMat[curIndex] = 0;
                }
            }
            else {
                assert(result == 1);
                if (minScore == matchScore) {
                    pathMat[curIndex] = Diagonal;
                }
                else if (minScore == delScore) {
                    pathMat[curIndex] = Left;
                }
                else {
                    pathMat[curIndex] = Up;
                }

                float pMisMatch, pIns, pDel;
                // Assign these to anything over 1 to signal they are not assigned.
                pMisMatch = 2;
                pIns      = 2;
                pDel      = 2;
                if (computeProb) {
                    if (matchScore != INF_INT) {
                        pMisMatch = QVToLogPScale(scoreFn.NormalizedMatch(tSeq, t, qSeq, q));
                    }
                    if (insScore != INF_INT) {
                        pIns = QVToLogPScale(scoreFn.NormalizedInsertion(tSeq, t, qSeq, q));
                    }
                    if (delScore != INF_INT) {
                        pDel = QVToLogPScale(scoreFn.NormalizedDeletion(tSeq, t, qSeq, q));
                    }

                    if (qSeq.qual.Empty() == false) {
                        if (matchScore != INF_INT and delScore != INF_INT and insScore != INF_INT) {
                            probMat[curIndex] = LogSumOfThree(probMat[matchIndex] + pMisMatch,
                                    probMat[delIndex] + pDel,
                                    probMat[insIndex] + pIns);
                        }
                        else if (matchScore != INF_INT and delScore != INF_INT) {
                            probMat[curIndex] = LogSumOfTwo(probMat[matchIndex] + pMisMatch,
                                    probMat[delIndex] + pDel);
                        }
                        else if (matchScore != INF_INT and insScore != INF_INT) {
                            probMat[curIndex] = LogSumOfTwo(probMat[matchIndex] + pMisMatch,
                                    probMat[insIndex] + pIns);					
                        }
                        else if (insScore != INF_INT and delScore != INF_INT) {
                            probMat[curIndex] = LogSumOfTwo(probMat[delIndex] + pDel,
                                    probMat[insIndex] + pIns);
                        }
                        else if (matchScore != INF_INT) {
                            probMat[curIndex] = probMat[matchIndex] + pMisMatch;
                        }
                        else if (delScore != INF_INT) {
                            probMat[curIndex] = probMat[delIndex] + pDel;
                        }
                        else if (insScore != INF_INT) {
                            probMat[curIndex] = probMat[insIndex] + pIns;
                        }
                        //
                        // Not normalizing probabilities, but using value as if it
                        // was a probability later on, so cap at 0 (= log 1).
                        //
                        if (probMat[curIndex] > 0) {
                            probMat[curIndex] = 0;
                        }
                        assert(!std::isnan(probMat[curIndex]));
                    }
                }
            }
        }
    }		
    // Ok, for now just trace back from qend/tend
    q = qEnd-1;
    t = tEnd-1;
    std::vector<Arrow>  optAlignment;
    int bufferIndexIsValid;
    while(q >= qStart or t >= tStart) {
        bufferIndex = -1;
        bufferIndexIsValid = GetBufferIndex(guide, q, t, bufferIndex);
        assert(bufferIndexIsValid);
        assert(bufferIndex >= 0);
        Arrow arrow;
        arrow = pathMat[bufferIndex];
        if (arrow == NoArrow) {
            tSeq.ToAscii();
            qSeq.ToAscii();
            unsigned int gi;
            for (gi = 0; gi < guide.size(); gi++) {
                std::cout << guide[gi].q << " " << guide[gi].t << " " << guide[gi].tPre << " " << guide[gi].tPost << std::endl;
            }

            std::cout << "qseq: "<< std::endl;
            (static_cast<DNASequence*>(&qSeq))->PrintSeq(std::cout);
            std::cout << "tseq: "<< std::endl;
            (static_cast<DNASequence*>(&tSeq))->PrintSeq(std::cout);
            std::cout << "ERROR, this path has gone awry at " << q << " " << t << " !" << std::endl;
            exit(1);
        }
        optAlignment.push_back(arrow);
        if (arrow == Diagonal) {
            q--;
            t--;
        }
        else if (arrow == Up) {
            q--;
        }
        else if (arrow == Left) {
            t--;
        }
    }

    alignment.nCells = ComputeMatrixNElem(guide);
    std::reverse(optAlignment.begin(), optAlignment.end());
    alignment.qPos = qStart;
    alignment.tPos = tStart;
    alignment.ArrowPathToAlignment(optAlignment);
    RemoveAlignmentPrefixGaps(alignment);
    int lastIndex = 0;
    tSeq.Free();
    qSeq.Free();
    lastIndex = -1;
    if (GetBufferIndex(guide, qEnd - 1, tEnd - 1, lastIndex)) {
        alignment.score = scoreMat[lastIndex];
        if (computeProb) {
            alignment.probScore = probMat[lastIndex];
        }
        return scoreMat[lastIndex];
    }
    else {
        return 0;
    }
}

template<typename QSequence, typename TSequence, typename T_ScoreFn> //, typename T_BufferCache>
int GuidedAlign(QSequence &origQSeq, TSequence &origTSeq,
        T_ScoreFn &scoreFn,
        int bandSize,
        int sdpIns, int sdpDel, float sdpIndelRate,
        blasr::Alignment &alignment,
        AlignmentType alignType=Global,
        bool computeProb = false,
        int sdpTupleSize= 8) {
    blasr::Alignment sdpAlignment;

    int alignScore = SDPAlign(origQSeq, origTSeq,
            scoreFn, sdpTupleSize, 
            sdpIns, sdpDel, sdpIndelRate,
            sdpAlignment); //, Local, false, false);

    int b;
    for (b = 0; b < sdpAlignment.blocks.size(); b++) {
        sdpAlignment.blocks[b].qPos += sdpAlignment.qPos;
        sdpAlignment.blocks[b].tPos += sdpAlignment.tPos;
    }
    sdpAlignment.tPos = 0;
    sdpAlignment.qPos = 0;

    return GuidedAlign(origQSeq, origTSeq, sdpAlignment, scoreFn, bandSize, alignment,
            // fill in optional parameters
            alignType, computeProb);

}


//
// Use Case: No guide yet exists for this alignment, but using
// buffers.  Run SDP alignment first, then refine on the guide.
//

template<typename QSequence, typename TSequence, typename T_ScoreFn, typename T_BufferCache>
int GuidedAlign(QSequence &origQSeq, TSequence &origTSeq, 
        T_ScoreFn &scoreFn,
        int bandSize,
        int sdpIns, int sdpDel, float sdpIndelRate,
        T_BufferCache &buffers,
        blasr::Alignment &alignment, 
        AlignmentType alignType=Global,
        bool computeProb = false,
        int sdpTupleSize= 8) {

    blasr::Alignment sdpAlignment;

    int alignScore = SDPAlign(origQSeq, origTSeq,
            scoreFn, sdpTupleSize, 
            sdpIns, sdpDel, sdpIndelRate,
            sdpAlignment, buffers, Local, false, false);

    int b;
    for (b = 0; b < sdpAlignment.blocks.size(); b++) {
        sdpAlignment.blocks[b].qPos += sdpAlignment.qPos;
        sdpAlignment.blocks[b].tPos += sdpAlignment.tPos;
    }
    sdpAlignment.tPos = 0;
    sdpAlignment.qPos = 0;

    return GuidedAlign(origQSeq, origTSeq, sdpAlignment, scoreFn, bandSize, buffers, alignment,
            // fill in optional parameters
            alignType, computeProb);
}

//
// Use case, guide exists, using buffers
//
template<typename QSequence, typename TSequence, typename T_ScoreFn, typename T_BufferCache>
int GuidedAlign(QSequence &origQSeq, TSequence &origTSeq,  blasr::Alignment &guideAlignment,
        T_ScoreFn &scoreFn,
        int bandSize,
        T_BufferCache &buffers,
        blasr::Alignment &alignment, 
        AlignmentType alignType=Global, 
        bool computeProb=false) {
    return GuidedAlign(origQSeq, origTSeq, guideAlignment, scoreFn, bandSize, alignment, 
            buffers.scoreMat,
            buffers.pathMat,
            buffers.probMat,
            buffers.optPathProbMat,
            buffers.lnSubPValueMat,
            buffers.lnInsPValueMat,
            buffers.lnDelPValueMat,
            buffers.lnMatchPValueMat, alignType, computeProb);
}

//
// Missing the use case for guide does not exist, and not using
// buffers.  This is just a very long function declaration, so it's
// not worth writing until it is needed.
//


//
// Use case, guide exists, but not using buffers.                   
//
template<typename QSequence, typename TSequence, typename T_ScoreFn>
int GuidedAlign(QSequence &origQSeq, TSequence &origTSeq,  blasr::Alignment &guideAlignment,
        T_ScoreFn &scoreFn,
        int bandSize,
        blasr::Alignment &alignment, 
        AlignmentType alignType=Global, 
        bool computeProb=false) {

    //  Make synonyms for members of the buffers class for easier typing.
    std::vector<int>    scoreMat;
    std::vector<Arrow>  pathMat;
    std::vector<double> probMat;
    std::vector<double> optPathProbMat;
    std::vector<float>  lnSubPValueVect;
    std::vector<float>  lnInsPValueVect;
    std::vector<float>  lnDelPValueVect;
    std::vector<float>  lnMatchPValueVect;

    return GuidedAlign(origQSeq, origTSeq, guideAlignment, scoreFn, bandSize, alignment,
            scoreMat,
            pathMat,
            probMat,
            optPathProbMat,
            lnSubPValueVect,
            lnInsPValueVect,
            lnDelPValueVect,
            lnMatchPValueVect, alignType, computeProb);
}

#endif // _BLASR_GUIDE_ALIGNMENT_HPP_

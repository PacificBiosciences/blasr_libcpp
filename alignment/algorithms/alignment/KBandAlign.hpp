#ifndef _BLASR_K_BAND_ALIGN_HPP_
#define _BLASR_K_BAND_ALIGN_HPP_

#include <algorithm>
#include <vector>
#include <limits.h>
// pbdata
#include "../../../pbdata/defs.h"
#include "../../../pbdata/NucConversion.hpp"
#include "../../../pbdata/matrix/FlatMatrix.hpp"

#include "AlignmentUtils.hpp"
#include "../../datastructures/alignment/Alignment.hpp"
#include "../../statistics/StatUtils.hpp"

class DefaultGuide {
 public:
	DefaultGuide();
	int operator()(int i);
};

template<typename T_QuerySequence, typename T_TargetSequence, typename T_Alignment, typename T_ScoreFn>
int KBandAlign(T_QuerySequence &pqSeq, T_TargetSequence &ptSeq,
							 int matchMat[5][5], int ins, int del, DNALength k,
							 std::vector<int> &scoreMat,
							 std::vector<Arrow> & pathMat,
							 T_Alignment &alignment,
							 T_ScoreFn &scoreFn,
							 AlignmentType alignType=Global,
               bool samplePaths=false) {
	return KBandAlign(pqSeq, ptSeq, matchMat, ins, del, k, scoreMat, pathMat, alignment, alignType, scoreFn, samplePaths);
}

void SetKBoundedLengths(DNALength tLength, DNALength qLength, DNALength k, DNALength &tLen, DNALength &qLen);

template<typename T_Sequence>
void CreateThreeBitSequence(T_Sequence &origSeq, T_Sequence &threeBitSeq) {
	//
	// Make a copy of the sequences that is guaranteed to be in 3-bit format.
	// This is 2 bits for A,C,T,and G, and an extra bit to signal masked sequence.
	//

	ResizeSequence(threeBitSeq, origSeq.length);
	
	VectorIndex i;
	for (i = 0; i < origSeq.length; i++ ) {
		threeBitSeq.seq[i] = ThreeBit[origSeq.seq[i]];
	}
}


template<typename T_QuerySequence, typename T_TargetSequence, typename T_Alignment, typename T_ScoreFn>
int KBandAlign(T_QuerySequence &qSeq, T_TargetSequence &tSeq,
							 int matchMat[5][5], int ins, int del, int k,
							 std::vector<int>   &scoreMat,
							 std::vector<Arrow> &pathMat,
							 T_Alignment   &alignment, 
							 AlignmentType alignType, 
							 T_ScoreFn &scoreFn, bool samplePaths=false) {
    (void)(matchMat);

	DNALength qLen, tLen;	
	SetKBoundedLengths(tSeq.length, qSeq.length, k, tLen, qLen);

	//
	//
	// Allow for length up to diaonal + k + 1 for boundary.
	// 
	// Allow for width:
	//   diagonal (1)
	//   up to k insertions (k)
	//   up to k deletions  (k)
	//   boundary on left side of matrix (1)
	// 
  DNALength nCols = 2*k + 1;
	DNALength totalMatSize = (qLen + 1) * nCols;
	alignment.nCells = totalMatSize;
	if (scoreMat.size() < totalMatSize) {
		scoreMat.resize(totalMatSize);
		pathMat.resize(totalMatSize);
	}

	// 
	// Initialze matrices
	//
	std::fill(scoreMat.begin(), scoreMat.begin() + totalMatSize, 0);
	std::fill(pathMat.begin(),  pathMat.begin()  + totalMatSize, NoArrow);
	
	//
	// Initialize the boundaries of the score and path matrices.
	//
	int q, t;

	for (q = 1; q <= k && q < qLen + 1; q++) {
		scoreMat[rc2index(q, k - q, nCols)] = q * ins;
		pathMat[rc2index(q, k - q , nCols)] = Up;
	}
	if (alignType == Global) {
		for (t = 1; t <= k && t < tLen; t++) {
			scoreMat[rc2index(0, t + k , nCols)] = t * del;
			pathMat[rc2index(0, t + k , nCols)] = Left;
		}
	}
	if (alignType == QueryFit or alignType == Fit) {
		for (t = 1; t <= k && t < tLen; t++) {
			scoreMat[rc2index(0, t + k , nCols)] = 0;
			pathMat[rc2index(0, t + k , nCols)] = Left;
		}
	}
	if (alignType == TargetFit or alignType == Fit) {
		for (q = 1; q <= k && q < qLen; q++) {
			scoreMat[rc2index(q, 0, nCols)] = 0;
			pathMat[rc2index(q, 0, nCols)] = Up;
		}
	}

	//
	// Initialize the 0,0 position to be a match.
	//
	scoreMat[rc2index(0, k, nCols)] = 0;
	pathMat[rc2index(0, k, nCols)] = Diagonal;

	int matchScore, insScore, delScore;

	for (q = 1; q <= qLen; q++) {
		for (t = q - k; t < q + k + 1; t++) {
			if (t < 1)
				continue;
			if (t > tLen)
				continue;

			// On left boundary of k-band. 
			// do not allow deletions of t.
			if (t == q - k) {
				delScore = INF_INT;
			}
			else {
				// cur row = q
				// cur col = t - q 
				// prev col therefore t - q - 1
				// and offset from diagonal is k + t - q - 1
				delScore = scoreMat[rc2index(q, k + t - q - 1, nCols)] + scoreFn.Deletion(tSeq, (DNALength) t-1, qSeq, (DNALength)q-1);
			}

			// cur row = q
			// cur col = t - q
			
			// cur query index = q - 1
			// cur target index = t - 1
			// therefore match row (up) = q 
			//           match col (left, but since up shifted right) = t - q
			assert(rc2index(q - 1, k + t - q, nCols) < scoreMat.size());
			assert(t-1 >= 0);
			assert(q-1 >= 0);
			int  tmpMatchScore = scoreFn.Match(tSeq, t-1, qSeq, q-1);
			matchScore = scoreMat[rc2index(q - 1, k + t - q, nCols)] + tmpMatchScore;

			//
			//  Possibly on right boundary of k-band, in which
			//  case do not allow insertions from q.
			if (t == q + k ) {
				insScore = INF_INT;
			}
			else {
				// cur row = q
				// cur col = t - q
				// therefore insertion col = t - q + 1
				insScore = scoreMat[rc2index(q-1, k + t - q+1, nCols)] + scoreFn.Insertion(tSeq, (DNALength) t-1, qSeq, q-1);
			}
		
			int minScore = MIN(matchScore, MIN(insScore, delScore));
			int curIndex = rc2index(q, k + t - q, nCols);
			assert(curIndex < scoreMat.size());
			scoreMat[curIndex] = minScore;
      int nEqual = 0;
      (matchScore == minScore ? nEqual++ : nEqual );
      (insScore == minScore ? nEqual++ : nEqual );
      (delScore == minScore ? nEqual++ : nEqual );
      if (samplePaths == false or nEqual == 1) {
        if (minScore == matchScore) {
          pathMat[curIndex] = Diagonal;
        }
        else if (minScore == delScore) {
          pathMat[curIndex] = Left;
        }
        else {
          pathMat[curIndex] = Up;
        }
      }
      else {
        //
        // When there are paths of equal score reaching 
        //
        if (nEqual == 3) {
          int v = RandomInt(3);
          if (v == 0) { pathMat[curIndex] = Diagonal; }
          else if (v == 1) { pathMat[curIndex] = Left; }
          else if (v == 2) { pathMat[curIndex] = Up; }
        }
        else {
          assert(nEqual == 2);
          int v = RandomInt(2);
          if (matchScore == insScore) {
            if (v == 0) { pathMat[curIndex] = Diagonal; } else { pathMat[curIndex] = Up; }
          }
          else if (matchScore == delScore) {
            if (v == 0) { pathMat[curIndex] = Diagonal; } else { pathMat[curIndex] = Left; }
          }
          else if (delScore == insScore) {
            if (v == 0) { pathMat[curIndex] = Left; } else { pathMat[curIndex] = Up; }
          }
          else {
            std::cout << "ERROR, counted two values equal to the minimum but cannot find them." << std::endl;
            assert(0);
          }
        }
        alignment.nSampledPaths++;
      }
		}
	}

	//
	// Now create the alignment.
	//
	
	q = qLen ;
	t = k - (qLen - tLen);

	int globalMinScore;
	int minLastColScoreIndex=0, minLastRowScoreIndex=0;
	globalMinScore = scoreMat[rc2index(q,t,nCols)];
	int minLastColScore = globalMinScore, minLastRowScore = globalMinScore;

	if (alignType == QueryFit or alignType == Fit) {
		int q2,t2;
		q2 = qLen;
		t2 = k - (qLen - tLen);
		bool minScoreSet = false;
		for (t2 = q - k; t2 < q2 + k + 1; t2++) {
			if (t2 < 1)
				continue;
			if (t2 > tLen)
				continue;
      //      std::cout << t2 << " " << tLen << " " << " " << scoreMat[rc2index(q2, k+t2-q,nCols)] << " " << minLastRowScore <<std::endl;
			if (minScoreSet == false or scoreMat[rc2index(q2, k+t2-q,nCols)] < minLastRowScore){ 
				minScoreSet     = true;
				minLastRowScore = scoreMat[rc2index(q2,k+t2-q,nCols)];
				minLastRowScoreIndex   = t2;
			}
		}
		if (minScoreSet) {
      t2 = k - (q - minLastRowScoreIndex);
			t = t2; 
			q = q2;
		}
	}
	if (alignType == TargetFit or alignType == Fit) {
		//  Fit the target inside the query wh
		int q2,t2;
		q2 = qLen;
		t2 = k - (qLen - tLen);
		bool minScoreSet = false;
		for (q2 = qLen; q2 >= tLen - k and q2 > 0; q2--) {
			int index = rc2index(q2, k+tLen-q2, nCols);
			if (minScoreSet == false or scoreMat[index] < minLastColScore) {
				minLastColScore = scoreMat[index];
				minScoreSet = true;
				minLastColScoreIndex = q2;
			}
		}
		if (alignType == Fit) {
			if (minLastColScore < minLastRowScore) {
				t = t2;
				q = minLastColScoreIndex;
			}
		}
		else if (alignType == TargetFit) {
			t = t2;
			q = minLastColScoreIndex;
		}
	}

	std::vector<Arrow>  optAlignment;
	int optScore = scoreMat[rc2index(q, t, nCols)];
	Arrow arrow;
	/*
	PrintFlatMatrix(&pathMat[0], qLen + 1, nCols, debugOut);
	std::cout << std::endl;
  ofstream debugOut;
  stringstream debugOutName;
  debugOutName << "kband_" << kbandcounter << ".table";
  debugOut.open(debugOutName.str().c_str());
	PrintFlatMatrix(&scoreMat[0], qLen + 1, nCols, debugOut);
  kbandcounter++;
  */
  /*
	std::cout << std::endl;
	*/
	//
	// Use some logic to deal with unsigned types.  When t > k, t must
	// also be greater than 0, so it's not worth checking to see if it
	// hits a boundary.
	//
	if (alignType == Global or alignType == QueryFit) {
		while (q > 0 and (t < k ? (k - t != q) : true)) {
			arrow = pathMat[rc2index(q,t, nCols)];
			if (arrow == NoArrow) {
				break;
			}
			optAlignment.push_back(arrow);
			if (arrow == Diagonal) {
				q--;
			}
			else if (arrow == Up) {
				q--;
				t++;
			}
			else if (arrow == Left) {
				t--;
			}
		}
	}
	else if (alignType == Fit) {
		while (q > 0 and (t < k ? (k - t != q) : true)  and  ( q <= k ? k - q != t : true) ) {
			arrow = pathMat[rc2index(q,t, nCols)];
			if (arrow == NoArrow) {
				break;
			}
			optAlignment.push_back(arrow);
			if (arrow == Diagonal) {
				q--;
			}
			else if (arrow == Up) {
				q--;
				t++;
			}
			else if (arrow == Left) {
				t--;
			}
		}
	}
	else if (alignType == TargetFit) {
		while (q > 0 and ( q < k ? k - q != t : true) ) {
			arrow = pathMat[rc2index(q,t, nCols)];
			if (arrow == NoArrow) {
				break;
			}
			optAlignment.push_back(arrow);
			if (arrow == Diagonal) {
				q--;
			}
			else if (arrow == Up) {
				q--;
				t++;
			}
			else if (arrow == Left) {
				t--;
			}
		}
	}

	// remove the boundary condition.
	//optAlignment.pop_back();
	//	qSeq.Free();
	//	tSeq.Free();
	alignment.qPos = q;
	
	//
	// Use a little extra logic to deal with the unsignedness of indices.
	//
	if (t < k) {
		alignment.tPos = (k - t) - q;
	}
	else {
		alignment.tPos = (t - k) - q;
	}
	std::reverse(optAlignment.begin(), optAlignment.end());
	alignment.ArrowPathToAlignment(optAlignment);
	return optScore;
}

#endif // _BLASR_K_BAND_ALIGN_HPP_

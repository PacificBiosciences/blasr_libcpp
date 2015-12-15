#include "GuidedAlign.hpp"

int GuideRow::GetRowLength() {
    return tPost + tPre + 1;
}
		
int GetBufferIndexFunctor::operator()(Guide &guide, int seqRow, int seqCol, int &index) {
    //
    // Whenever a previous index was found, just use the one in the array to the right
    //
    int prevIndex = -1;
    if (index != -1) {
        index++;
        return 1;
    }

    int row, col;
    row = seqRow + 1 - seqRowOffset;
    col = seqCol; // use the variable here for standardized naming.
    if (row < 0 or row > guideSize ) {
        return 0;
    }
    if (col <= guide[row].t and guide[row].t - col <= guide[row].tPre) {
        index = guide[row].matrixOffset - (guide[row].t - col);
        assert(prevIndex == -1 or prevIndex == index);
        return 1;
    }
    else if (col > guide[row].t and col - guide[row].t  <= guide[row].tPost) {
        index = guide[row].matrixOffset + (col - guide[row].t);
        assert(prevIndex == -1 or prevIndex == index);
        return 1;
    }
    return 0;
}

int ComputeMatrixNElem(Guide &guide) {
    int totalSize = 0;
    for (size_t r = 0; r < guide.size(); r++) {
        totalSize += guide[r].GetRowLength();
        //    cout << r << " " << totalSize << endl;
        assert(guide[r].GetRowLength() >= 0);
    }
    return totalSize;
}

void StoreMatrixOffsets(Guide &guide) {
    int curMatrixSize = 0;
    for (size_t r = 0; r < guide.size(); r++) {
        guide[r].matrixOffset = guide[r].tPre + curMatrixSize;
        curMatrixSize += guide[r].GetRowLength();
    }
}

float QVToLogPScale(char qv) {
    return qv/-10.0;
}

void QVToLogPScale(QualityValueVector<QualityValue> &qualVect, size_t phredVectLength, std::vector<float> &lnVect) {
    if (phredVectLength > lnVect.size()) {
        lnVect.resize(phredVectLength);
    }
    for (size_t i = 0; i < phredVectLength; i++) {
        lnVect[i] = qualVect[i]/-10.0; //log(qualVect.ToProbability(i));
    }
}


int AlignmentToGuide(Alignment &alignment, Guide &guide, int bandSize)  {
	guide.clear();
	if (alignment.size() == 0) {
		// no blocks to make guide, exit.
		return 0;
	}
	
	int tStart, qStart, qEnd;
	int firstBlock   = 0;
	int lastBlock    = alignment.size() - 1;
	
	tStart = alignment.blocks[firstBlock].tPos;
	//int tEnd   = alignment.blocks[lastBlock].TEnd();
	qStart = alignment.blocks[firstBlock].qPos;
	qEnd   = alignment.blocks[lastBlock].QEnd();

	int qAlignLength = qEnd - qStart; 
	
	// Add one extra block for boundary conditions.
	guide.resize(qAlignLength+1);

	// Initilize the first (boundary condition) row.
	guide[0].t     = tStart - 1;
	guide[0].q     = qStart - 1;
	int drift = abs(tStart - qStart);	
	if (drift > bandSize) {
		guide[0].tPost = drift;
	}
	else {
		guide[0].tPost = bandSize;
	}
	guide[0].tPre  = 0;

	// The first row of the guide matches 
	int guideIndex = 1;

	for (size_t b = 0; b < alignment.blocks.size(); b++) {
		//
		// First add the match stored in block b, each block is a
		// diagonal, so that makes life easy.  
		//
		for (DNALength bp = 0; bp < alignment.blocks[b].length; bp++) {
			guide[guideIndex].t     = alignment.blocks[b].tPos + bp;
			guide[guideIndex].q     = alignment.blocks[b].qPos + bp;
			// 
			// This complicated logic is to determine how far back the band
			// should stretch.  The problem is that if the band stretches
			// back further than the previous row, it's possible for the
			// path matrix to go backwards into cells that should not be
			// touched.  
			//
			//int tDiff = guide[guideIndex].t - guide[guideIndex-1].t;
			if (bp == 0) {
				guide[guideIndex].tPre  = guide[guideIndex].t - 
					(guide[guideIndex-1].t - guide[guideIndex-1].tPre); 
				guide[guideIndex].tPost = bandSize + abs(drift);
			}
			else {
				//
				// Within aligned blocks, align around the band size.
				//
				int fullLengthTPre = (guide[guideIndex].t - 
									  (guide[guideIndex-1].t - 
									   guide[guideIndex-1].tPre));
				guide[guideIndex].tPre  = min(bandSize, fullLengthTPre);
				guide[guideIndex].tPost = min(MAX_BAND_SIZE, bandSize);
                //if (guide[guideIndex].tPre > 500 or guide[guideIndex].tPost > 500) {
                //cout << guideIndex << " " << guide[guideIndex].tPre << " " << guide[guideIndex].tPost << endl;
                //}
                //assert(guide[guideIndex].tPre >= 0);
                //assert(guide[guideIndex].tPost >= 0);
			}
			guideIndex++;
		}

		//
		// Now, widen k around regions where there is drift from the diagonal.
		//


		int diagonalLength;
		int qGap, tGap;

		if (b < alignment.blocks.size()-1) {
			qGap = alignment.blocks[b+1].qPos - alignment.blocks[b].QEnd();
			tGap = alignment.blocks[b+1].tPos - alignment.blocks[b].TEnd();
			// 
			// Drift is how far from the diagonal the next block starts at.
			//
			drift           = ComputeDrift(alignment.blocks[b], alignment.blocks[b+1]);
            //drift = min(drift, 100);
			diagonalLength  = std::min(qGap, tGap);
			
			int diagPos;
			int qPos, tPos;
			int qEnd;
			
			qPos = alignment.blocks[b].QEnd();
			tPos = alignment.blocks[b].TEnd();

			qEnd = alignment.blocks[b+1].qPos;
			//int tEnd = alignment.blocks[b+1].tPos;

			for (diagPos = 0; diagPos < diagonalLength; diagPos++, tPos++, qPos++) {
				guide[guideIndex].t     = tPos;
				guide[guideIndex].q     = qPos;
				guide[guideIndex].tPre  = min(MAX_BAND_SIZE, (guide[guideIndex].t - 
                                   (guide[guideIndex-1].t - guide[guideIndex-1].tPre)));
				guide[guideIndex].tPost = min(MAX_BAND_SIZE, bandSize + abs(drift));
/*        if (guide[guideIndex].tPre > 500 or guide[guideIndex].tPost > 500) {
          cout << guideIndex << " " << guide[guideIndex].tPre << " " << guide[guideIndex].tPost << endl;
        } */
				++guideIndex;
			}
			
			//
			// If the query gap is shorter than target (there is a deletion
			// of the target),  the guide must be extended down the side of
			// the gap.  See the figure below.
			//
			//  *****
			//  **  
			//  * *
			//  *  *
			//  *   * // extend down from here.
		//  *   *
			//  *   * 
			//

			while (qPos < qEnd) {
				guide[guideIndex].t = tPos;
				guide[guideIndex].q = qPos;
				// move q down.
				qPos++;
				// keep tPos fixed, the guide is straight down here.
				guide[guideIndex].tPre = min(MAX_BAND_SIZE, guide[guideIndex].t - 
					(guide[guideIndex-1].t - guide[guideIndex-1].tPre)); //bandSize + abs(drift);

				guide[guideIndex].tPost = min(MAX_BAND_SIZE, bandSize + abs(drift));
				guideIndex++;
			}
		}
	}
  //int i;
  //for (i = 0; i < guide.size(); i++) {
  //  guide[i].tPre = min(guide[i].tPre, 200);
  //  guide[i].tPost = min(guide[i].tPost, 200);
  // }
  return 1; // signal ok.

}

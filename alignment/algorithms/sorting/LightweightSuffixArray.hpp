#ifndef ALGORITHMS_SORTING_LIGHTWEIGHT_SUFFIX_ARRAY_H_
#define ALGORITHMS_SORTING_LIGHTWEIGHT_SUFFIX_ARRAY_H_

#include <algorithm>
#include "qsufsort.hpp"
#include "MultikeyQuicksort.hpp"
#include "DifferenceCovers.hpp"
#include "Types.h"

/*
 * a - b potentially may not fit into a signed type.  Use some logic
 * to get around that.
 */

UInt DiffMod(UInt a, UInt b, UInt d); 

void BuildDiffCoverReverseLookup(UInt diffCover[], UInt diffCoverLength,
        UInt reverseDiffCover[] // of size diffCoverSize
        );

UInt DiffCoverFindH(UInt diffCover[], UInt diffCoverLength, UInt diffCoverSize, UInt textSize); 

class DiffCoverMu {
public:
    UInt *diffCoverReverseLookup;
    UInt diffCoverLength;
    UInt diffCoverSize;
    UInt textSize;
    UInt h;
    UInt *diffCover;

    UInt compute(UInt i, UInt j); 

    UInt operator()(const UInt k); 

    DiffCoverMu(); 

    ~DiffCoverMu(); 

    void Initialize(UInt diffCoverP[], UInt diffCoverLengthP, UInt diffCoverSizeP, UInt textSizeP); 

};

void BuildDiffCoverLookup(UInt diffCover[], UInt diffCoverLength, UInt v, UInt diffCoverLookup[]); 

class DiffCoverDelta {
public:
    UInt *diffCoverLookup;
    UInt diffCoverSize;

    void Initialize(UInt diffCoverP[], UInt diffCoverLengthP, UInt diffCoverSizeP); 

    UInt operator()(UInt i, UInt j); 

    ~DiffCoverDelta();
};



UInt NCompareSuffices(unsigned char text[], UInt a, UInt b, UInt n); 

UInt ComputeDSetSize(UInt diffCover, UInt diffCoverLength, UInt diffCoverSize, UInt textSize); 

void ComputeSufVNaming(UInt diffCover[], UInt diffCoverLength, UInt diffCoverN, UInt textSize, UInt lexNaming[], 
        DiffCoverMu &mu,
        UInt sufVNaming[]); 

UInt IndexToDiffCoverIndex(UInt index, UInt diffCoverlookup[], UInt diffCoverSize, UInt diffCoverLength );

void DiffCoverComputeLOrder(UInt sufVNaming[], UInt sufVNamingLength, UInt maxVNaming, UInt textLength, DiffCoverMu &mu, UInt lOrder[]); 


/*
 * Build the lex naming of the v-ordered suffices.  
 *
 * Input: textVOrder - the v-ordering of a subset of the text.
 *        textSize   - the size of the v-order set.
 *        diffCoverLength - length of diffCover
 *        diffCoverSize - the size of the diff cover.
 * Output: lexNaming: the lex-naming of the v-order suffices.  The
 *        names are implemented as unsigned integers. 
 * Returns: the largest value of the lex-ordering.
 */
UInt DiffCoverBuildLexNaming( unsigned char text[], UInt textSize,
        UInt textVOrder[],
        UInt dSetSize, UInt diffCoverLength, UInt diffCoverSize, 
        UInt diffCoverLookup[],
        UInt lexNaming[]); 

class DiffCoverCompareSuffices {
public:
    UInt *lOrder;
    DiffCoverDelta *delta;
    UInt diffCoverSize;
    UInt diffCoverLength;
    UInt *diffCoverReverseLookup;
    int operator()(UInt a, UInt b); 
};

bool LightweightSuffixSort(unsigned char text[], UInt textLength, UInt *index, int diffCoverSize); 

#endif

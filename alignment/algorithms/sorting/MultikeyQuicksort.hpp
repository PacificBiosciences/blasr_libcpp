#ifndef _BLASR_MULTIKEY_QUICKSORT_HPP_
#define _BLASR_MULTIKEY_QUICKSORT_HPP_

/*
 * This is an implementation of MultiKey Quicksort, or ssort1 from
 * Bentley and Sedgewick, Fast Algorithms for Sorting and Searching
 * Strings, Proc. 8th Annual ACM-SIAM Symposium on Discrete Algorithms
 * (SODA), pages 360-369, January 1997. 
 *
 * The implementation here sorts lists of indices into an array of
 * substrings rather than the substrings themselves (for use in suffix
 * array generation).  Furthermore, it is made to be sse-alizeable to
 * speed up sorts on modern architectures. 
 */

#include <limits.h>
#include <vector>
#include <algorithm>
#include "NucConversion.hpp"
#include "FASTASequence.hpp"
#include "Types.h"

typedef unsigned int UInt;

void UIntSwap(unsigned int &a, unsigned int &b); 

void VecSwap( UInt i, UInt j, UInt n, UInt index[]); 

unsigned char ComputeMedianValue(unsigned char text[], UInt index[], int length, 
        UInt low, UInt high, int offset, UInt maxPossible, UInt *freq ); 

UInt FindFirstOf(unsigned char text[], UInt index[],
        UInt low, UInt high, int offset, Nucleotide nuc);

void SwapIndices(UInt &a, UInt &b);

void TransformSequenceForSorting(unsigned char text[], UInt textLength, int bound); 

void TransformBackSequence(Nucleotide text[], UInt textLength); 

/*
 * depth: the current depth of how much is sorted.
 * bound: how far to sort.
 */
void MediankeyBoundedQuicksort(unsigned char text[], UInt index[], UInt length,
        UInt low, UInt high, int depth, int bound, UInt maxChar= 0, UInt *freq=NULL); 

#endif // _BLASR_MULTIKEY_QUICKSORT_HPP_

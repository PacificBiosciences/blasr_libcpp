#ifndef _BLASR_TYPES_H_
#define _BLASR_TYPES_H_

#include <stdint.h>

typedef unsigned char Nucleotide;
// DNA sequence related definitions.
typedef uint32_t DNALength;
// Genome length, first step towards > 4G genome length, not there yet.
typedef uint64_t GenomeLength;

//
// Add definitions to handle 64/32 bit computing environments
//
// H5 Dataset length (e.g., total number of pulses in a pls.h5 file.)
typedef uint64_t DSLength; 
// used in alignment/tuple/BaseTuple class
typedef uint64_t TupleData;

typedef uint32_t VectorIndex;
typedef uint32_t UInt;
typedef uint8_t  Byte;
typedef uint8_t  UChar;
typedef uint16_t HalfWord;
typedef float MatchWeight;
#endif // _BLASR_TYPES_H_

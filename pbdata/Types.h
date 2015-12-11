#ifndef _BLASR_TYPES_H_
#define _BLASR_TYPES_H_

#include <stdint.h>

typedef unsigned long ULong;

// Move DNA sequence related definitions.
typedef uint32_t DNALength;
typedef unsigned char Nucleotide;

//
// Add definitions to handle 64/32 bit computing environments
//
// H5 Dataset length (e.g., total number of pulses in a pls.h5 file.)
typedef uint64_t DSLength; 
typedef uint32_t VectorIndex;
typedef uint32_t UInt;
typedef uint8_t  Byte;
typedef uint8_t  UChar;
typedef uint16_t HalfWord;
typedef float MatchWeight;
#endif // _BLASR_TYPES_H_

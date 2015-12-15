#ifndef _BLASR_FORMAT_SAMPRINTER_HPP_
#define _BLASR_FORMAT_SAMPRINTER_HPP_ 

#include <sstream>
#include <stdint.h>
#include "SMRTSequence.hpp"
#include "datastructures/alignment/AlignmentCandidate.hpp"
#include "datastructures/alignment/AlignmentContext.hpp"
#include "datastructures/alignment/Alignment.hpp"
#include "datastructures/alignment/Alignment.hpp"
#include "datastructures/alignmentset/SAMSupplementalQVList.hpp"


#define MULTI_SEGMENTS 0x1
#define ALL_SEGMENTS_ALIGNED 0x2
#define SEGMENT_UNMAPPED 0x4
#define NEXT_SEGMENT_UNMAPPED 0x8
#define SEQ_REVERSED 0x10
#define SEQ_NEXT_REVERSED 0x20
#define FIRST_SEGMENT 0x40
#define LAST_SEGMENT 0x80
#define SECONDARY_ALIGNMENT 0x100
#define NO_PASS_QUALITY 0x200
#define PCR_OR_OPTICAL_DUPLICATE 0x400

namespace SAMOutput {

enum Clipping {hard, soft, subread, none};

void BuildFlag(T_AlignmentCandidate &alignment, AlignmentContext &context, uint16_t &flag); 

//
// Trimming is used for both hard non-clipping
// so it is called trim instead of clip.
//
void CreateDNAString(DNASequence &seq, DNASequence &clippedSeq, 
    DNALength trimFront=0, DNALength trimEnd=0); 

void AddGaps(T_AlignmentCandidate &alignment, int gapIndex,
        std::vector<int> &opSize, std::vector<char> &opChar);

// Add sequence match/mismatch CIGAR string Ops for block b.
void AddMatchBlockCigarOps(DNASequence & qSeq, DNASequence & tSeq, blasr::Block & b,
        DNALength & qSeqPos, DNALength & tSeqPos,
        std::vector<int> & opSize, std::vector<char> & opChar);

// Merge adjacent indels and mismatches.
void MergeAdjacentIndels(std::vector<int> &opSize, std::vector<char> &opChar,
                         const char mismatchChar);

// If cigarUseSeqMatch is true, cigar string uses '=' and 'X' 
// instead of 'M' to represent sequence match and mismatch;
void CreateNoClippingCigarOps(T_AlignmentCandidate &alignment, 
        std::vector<int> &opSize, std::vector<char> &opChar,
        bool cigarUseSeqMatch = false, 
        const bool allowAdjacentIndels = true); 
//
// 
// The aligned sequence is either the sequence from the first
// aligned base to the last (hard and no clipping), or first high
// quality base to the last high quality base (soft clipping).
//
template<typename T_Sequence>
void SetAlignedSequence(T_AlignmentCandidate &alignment, T_Sequence &read,
    T_Sequence &alignedSeq, Clipping clipping = none); 

template<typename T_Sequence>
void SetSoftClip(T_AlignmentCandidate &alignment, T_Sequence &read, 
    DNALength hardClipPrefix, DNALength hardClipSuffix,
    DNALength &softClipPrefix, DNALength &softClipSuffix); 

template<typename T_Sequence>
void SetHardClip(T_AlignmentCandidate &alignment, T_Sequence &read, 
    DNALength &prefixClip, DNALength &suffixClip); 

void CigarOpsToString(std::vector<int> &opSize, std::vector<char> &opChar, 
        std::string &cigarString);

//
// Straight forward: create the cigar string allowing some clipping
// The read is provided to give length and hq information.
//
template<typename T_Sequence>
void CreateCIGARString(T_AlignmentCandidate &alignment, T_Sequence &read,
        std::string &cigarString, Clipping clipping, 
        DNALength &prefixSoftClip, DNALength &suffixSoftClip,
        DNALength &prefixHardClip, DNALength &suffixHardClip,
        bool cigarUseSeqMatch = false,
        const bool allowAdjacentIndels = true); 

template<typename T_Sequence>
void PrintAlignment(T_AlignmentCandidate &alignment, T_Sequence &read,
        std::ostream &samFile, AlignmentContext &context, 
        SupplementalQVList & qvList, Clipping clipping = none,
        bool cigarUseSeqMatch = false,
        const bool allowAdjacentIndels = true); 
}

#include "SAMPrinterImpl.hpp"
#endif // _BLASR_FORMAT_SAMPRINTER_HPP_

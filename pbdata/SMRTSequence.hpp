// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Mark Chaisson

#ifndef _BLASR_SMRT_SEQUENCE_HPP_
#define _BLASR_SMRT_SEQUENCE_HPP_

#include <cassert>
#include <iostream>
#include <sstream>

#include "Types.h"
#include "Enumerations.h"
#include "NucConversion.hpp"
#include "FASTQSequence.hpp"
#include "reads/RegionTable.hpp"
#include "reads/ZMWGroupEntry.hpp"


class SMRTSequence : public FASTQSequence {
private:
    enum SnrIndex4Base {A=0, C=1, G=2, T=3};
    float hqRegionSnr_[4]; // Always saved as 'ACGT'
 
public:
    int16_t xy[2];
    UInt holeNumber;

    float readScore;
    ZMWGroupEntry zmwData;
    PlatformId platform;

    HalfWord *preBaseFrames;
    HalfWord *widthInFrames;
    //
    // The following are fields that are read in from the pulse file.
    // Because they are not standard in bas.h5 files, these fields
    // should not be preallocated when resizing a SMRTSequence, and
    // memory should be managed separately.  For now, these fields all
    // have the same length as the number of bases, but this could
    // change so that all pulse values are stored in a SMRTSequence.
    //
    HalfWord *meanSignal, *maxSignal, *midSignal;
    float *classifierQV;
    unsigned int *startFrame;
    int *pulseIndex;
    DNALength lowQualityPrefix, lowQualitySuffix;
    int highQualityRegionScore; // High quality region score in region table.

public:
    DNALength subreadStart;
    DNALength subreadEnd;

protected:
    // read group id associated with each SMRTSequence
    std::string readGroupId; 

public:
    // Whether or not this is originally copied from a BamRecord.
    bool copiedFromBam;

    void SetNull(); 

    SMRTSequence();
    inline ~SMRTSequence();

    // Access to HQRegion SNRs must be done via public API.
    inline float HQRegionSnr(const char base) const;

    inline SMRTSequence & HQRegionSnr(const char base, float v);

    void Allocate(DNALength length); 

    void SetSubreadTitle(SMRTSequence &subread, DNALength subreadStart, 
        DNALength subreadEnd); 

    void SetSubreadBoundaries(SMRTSequence &subread, DNALength subreadStart, 
        DNALength subreadEnd); 

    void MakeSubreadAsMasked(SMRTSequence &subread, DNALength subreadStart = 0, 
        int subreadEnd = -1); 

    void MakeSubreadAsReference(SMRTSequence &subread, DNALength subreadStart = 0, 
        int subreadEnd = -1); 

    void Copy(const SMRTSequence &rhs); 

    void Copy(const SMRTSequence &rhs, int rhsPos, int rhsLength); 

    void Print(std::ostream &out); 

    SMRTSequence& operator=(const SMRTSequence &rhs); 

    void Free(); 

    bool StoreXY(int16_t xyP[]); 

    bool StorePlatformId(PlatformId pid); 

    bool StoreHoleNumber(UInt holeNumberP);

    bool StoreHoleStatus(unsigned char s); 

    bool StoreZMWData(ZMWGroupEntry &data); 

    bool GetXY(int xyP[]); 

    bool GetHoleNumber(UInt & holeNumberP);   

    inline UInt HoleNumber(void) const;

    // Get read group id for this sequence.
    std::string GetReadGroupId();

    // Set readGroup Id for this sequence.
    void SetReadGroupId(const std::string & rid);
    
#ifdef USE_PBBAM
public:
    // Copy read sequence, title, holeNumber, readGroupId, and QVs
    // (iq, dq, sq, mq, st, dt) from BamRecord to this SMRTSequence.
    // If copyAllQVs is false, also copy all QVs.
    void Copy(const PacBio::BAM::BamRecord & record, 
              bool copyAllQVs = false);

    // Keep track of BamRecord from which this SMRTSequence is 
    // originally copied. However, one should NOT assume
    // that this SMRTSequence has the same sequence, title, QVs as 
    // the BamRecord, because this SMRTSequence may be created by
    // MakeSubreadAsMasked(...) or MakeRC(...).
    PacBio::BAM::BamRecord bamRecord;
#endif 
};

inline SMRTSequence::~SMRTSequence(){
    SMRTSequence::Free();
}

inline UInt SMRTSequence::HoleNumber(void) const {
    return holeNumber;
}

inline float SMRTSequence::HQRegionSnr(const char base) const {
    if (::toupper(base) == 'A')      return hqRegionSnr_[SMRTSequence::SnrIndex4Base::A];
    else if (::toupper(base) == 'C') return hqRegionSnr_[SMRTSequence::SnrIndex4Base::C];
    else if (::toupper(base) == 'G') return hqRegionSnr_[SMRTSequence::SnrIndex4Base::G];
    else if (::toupper(base) == 'T') return hqRegionSnr_[SMRTSequence::SnrIndex4Base::T];
    else assert("Base must be in A, C, G, T" == 0);
}

inline SMRTSequence & SMRTSequence::HQRegionSnr(const char base, float v) {
    if (::toupper(base) == 'A')      hqRegionSnr_[SMRTSequence::SnrIndex4Base::A] = v;
    else if (::toupper(base) == 'C') hqRegionSnr_[SMRTSequence::SnrIndex4Base::C] = v;
    else if (::toupper(base) == 'G') hqRegionSnr_[SMRTSequence::SnrIndex4Base::G] = v;
    else if (::toupper(base) == 'T') hqRegionSnr_[SMRTSequence::SnrIndex4Base::T] = v;
    else assert("Base must be in A, C, G, T" == 0);
    return *this;
}
#endif  // _BLASR_SMRT_SEQUENCE_HPP_

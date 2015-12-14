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
friend class HDFZMWReader;
friend class HDFZMWWriter;
friend class HDFZMWMetricsWriter;

private:
    enum SnrIndex4Base {A=0, C=1, G=2, T=3};
    float hqRegionSnr_[4]; // Always saved as 'ACGT'

    DNALength subreadStart_;
    DNALength subreadEnd_;

    // read group id associated with each SMRTSequence
    std::string readGroupId_;

public:
    ZMWGroupEntry zmwData;

    DNALength lowQualityPrefix, lowQualitySuffix;
    int highQualityRegionScore; // High quality region score in region table.
    float readScore;

    // Whether or not this is originally copied from a BamRecord.
    bool copiedFromBam;

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

public:
    SMRTSequence();

    inline ~SMRTSequence();

    /// \name Sets and gets attributes.
    /// \{
    /// Set HoleNumber.
    /// \returns this SMRTSequence
    SMRTSequence & HoleNumber(UInt holeNumber);

    /// \reutrns HoleNumber
    UInt HoleNumber(void) const;

    /// Set HoleXY
    SMRTSequence & HoleXY(const int x, const int y);

    /// \returns HoleX
    UInt HoleX(void) const;

    /// \returns HoleY
    UInt HoleY(void) const;

    /// Set HoleStatus
    SMRTSequence & HoleStatus(const unsigned char);

    /// \returns HoleStatus
    unsigned char HoleStatus(void) const;

    /// \returns movie name parsed from sequence title
    std::string MovieName(void) const;

    /// \returns start pos of this sequence in coordinate of zmw polymerase sequence
    DNALength SubreadStart(void) const;

    /// Sets subreadStart.
    SMRTSequence & SubreadStart(const DNALength start);

    /// \returns subread end pos of this sequence in coordinate of zmw polymerase sequence
    DNALength SubreadEnd(void) const;

    /// Set subread end pos in coordinate of polymerase sequence.
    SMRTSequence & SubreadEnd(const DNALength end);

    /// A SMRTSequence's this->seq may point to sequence of a whole
    /// polymerase read, but only represents a subread [subreadStart_,  subreadEnd_).
    /// \returns subread length (SubreadEnd() - SubreadStart())
    DNALength SubreadLength(void) const;

    /// \returns read group id for this sequence.
    std::string ReadGroupId(void) const;

    /// Set readGroup Id for this sequence.
    SMRTSequence & ReadGroupId(const std::string & rid);

    /// Access to HQRegion SNRs must be done via public API.
    float HQRegionSnr(const char base) const;

    /// Set HQRegion SNR of base as v.
    SMRTSequence & HQRegionSnr(const char base, float v);

    /// \}

public:
    /// \name Clip subread
    /// \{
    SMRTSequence & Clip(const DNALength subreadStart, const DNALength subreadEnd);
    /// \}

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

    void Copy(const SMRTSequence &rhs, DNALength rhsPos, DNALength rhsLength); 

    void Print(std::ostream &out) const;

    SMRTSequence& operator=(const SMRTSequence &rhs); 

    void Free(); 
    
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

#endif  // _BLASR_SMRT_SEQUENCE_HPP_

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

#include <stdlib.h> 
#include "SMRTSequence.hpp"

using namespace std;

void SMRTSequence::SetNull() {
    pulseIndex    = NULL;
    preBaseFrames = NULL;
    widthInFrames = NULL;
    xy[0] = 0; xy[1] = 0;
    // These are not allocted by default.
    meanSignal = maxSignal = midSignal = NULL;
    classifierQV = NULL;
    startFrame   = NULL;
    platform     = NoPlatform;
    // By default, allow the entire read.
    lowQualityPrefix = lowQualitySuffix = 0;
    highQualityRegionScore = 0;
    // ZMWMetrics
    for (size_t i = 0; i < 4; i++) {
        hqRegionSnr_[i] = -1;
    }
    readScore = -1;
    holeNumber = static_cast<UInt>(-1);
    readGroupId = "";
    copiedFromBam = false;
#ifdef USE_PBBAM
    bamRecord = PacBio::BAM::BamRecord();
#endif
}

SMRTSequence::SMRTSequence() : FASTQSequence() {
    SetNull();
}

void SMRTSequence::Allocate(DNALength length) {
    // Assert *this has no allocated space.
    if (not (seq == NULL && preBaseFrames == NULL &&
             widthInFrames == NULL and pulseIndex == NULL)) {
        cout << "ERROR, trying to double-allocate memory for a SMRTSequence." << endl;
        exit(1);
    }

    FASTQSequence::AllocateRichQualityValues(length);
    seq           = new Nucleotide[length];
    this->length  = length;
    qual.Allocate(length);
    preBaseFrames = new HalfWord[length];
    widthInFrames = new HalfWord[length];
    pulseIndex    = new int[length];
    subreadEnd    = length;
    deleteOnExit  = true;
}

void SMRTSequence::SetSubreadTitle(SMRTSequence &subread, DNALength subreadStart, DNALength subreadEnd) {
    stringstream titleStream;
    titleStream << title << "/"<< subreadStart << "_" << subreadEnd;
    subread.CopyTitle(titleStream.str());
}    

void SMRTSequence::SetSubreadBoundaries(SMRTSequence &subread, DNALength subreadStart, DNALength subreadEnd) {
    if (subreadEnd == -1) {
        subreadEnd = length;
    }
    assert(subreadEnd - subreadStart <= length);
    subread.subreadStart= subreadStart;
    subread.subreadEnd  = subreadEnd;
    SetSubreadTitle(subread, subreadStart, subreadEnd);
}

void SMRTSequence::MakeSubreadAsMasked(SMRTSequence &subread, 
    DNALength subreadStart, int subreadEnd) {
    subread.Free();
    //
    // This creates the entire subread, but masks out the portions
    // that do not correspond to this insert.
    //
    static_cast<SMRTSequence*>(&subread)->Copy(*this);
    SetSubreadBoundaries(subread, subreadStart, subreadEnd);
    DNALength pos;
    for (pos = 0; pos < subreadStart; pos++) { subread.seq[pos] = 'N'; }
    for (pos = subreadEnd; pos < length; pos++) { subread.seq[pos] = 'N'; }
    // This is newly allocated memory, free it on exit.
    assert(subread.deleteOnExit);
}

void SMRTSequence::MakeSubreadAsReference(SMRTSequence &subread, 
    DNALength subreadStart, int subreadEnd) {
    subread.Free();
    //
    // Just create a reference to a substring of this read.  
    //
    static_cast<FASTQSequence*>(&subread)->ReferenceSubstring(*this, subreadStart, subreadEnd - subreadStart);
    SetSubreadBoundaries(subread, subreadStart, subreadEnd);
    // The subread references this read, protect the memory.
    assert(not subread.deleteOnExit);
}

void SMRTSequence::Copy(const SMRTSequence &rhs) {
    SMRTSequence::Copy(rhs, 0, rhs.length);
}

void SMRTSequence::Copy(const SMRTSequence &rhs, int rhsPos, int rhsLength) {
    // Sanity check
    CheckBeforeCopyOrReference(rhs, "SMRTSequence");
    
    // Free this SMRTSequence before copying anything from rhs.
    SMRTSequence::Free();

    FASTQSequence subseq; 
    // subseq.seq is referenced, while seq.title is not, we need to call 
    // subseq.Free() to prevent memory leak.
    static_cast<FASTQSequence*>(&subseq)->ReferenceSubstring(rhs, rhsPos, rhsLength);
    static_cast<FASTQSequence*>(&subseq)->CopyTitle(rhs.title, rhs.titleLength); 

    if (rhs.length == 0) {
        static_cast<FASTQSequence*>(this)->Copy(subseq);
        //
        // Make sure that no values of length 0 are allocated by returning here.
        //
    }
    else {
        assert(rhs.seq != seq);
        assert(rhsLength <= rhs.length);
        assert(rhsPos < rhs.length);

        // Copy seq, title and FASTQ QVs from subseq
        static_cast<FASTQSequence*>(this)->Copy(subseq); 

        // Copy SMRT QVs
        if (rhs.preBaseFrames != NULL) {
            preBaseFrames = new HalfWord[length];
            memcpy(preBaseFrames, rhs.preBaseFrames, length*sizeof(HalfWord));
        }
        if (rhs.widthInFrames != NULL) {
            widthInFrames = new HalfWord[length];
            memcpy(widthInFrames, rhs.widthInFrames, length*sizeof(HalfWord));
        }
        if (rhs.pulseIndex != NULL) {
            pulseIndex = new int[length];
            memcpy(pulseIndex, rhs.pulseIndex, sizeof(int) * length);
        }
    }

    // Copy other member variables from rhs
    subreadStart = rhs.subreadStart;
    subreadEnd   = rhs.subreadEnd;
    lowQualityPrefix = rhs.lowQualityPrefix;
    lowQualitySuffix = rhs.lowQualitySuffix;
    highQualityRegionScore = rhs.highQualityRegionScore;
    zmwData = rhs.zmwData;

    assert(deleteOnExit); // should have control over seq and all QVs

    subseq.Free();
    copiedFromBam = rhs.copiedFromBam;
#ifdef USE_PBBAM
    bamRecord = rhs.bamRecord;
#endif
}

void SMRTSequence::Print(ostream &out) {
    out << "SMRTSequence for zmw " << zmwData.holeNumber
        << ", [" << subreadStart << ", " << subreadEnd << ")" << endl;
    DNASequence::Print(out);
}

SMRTSequence& SMRTSequence::operator=(const SMRTSequence &rhs) {
    SMRTSequence::Copy(rhs);
    return *this;
}

void SMRTSequence::Free() {
    if (deleteOnExit == true) {
        if (preBaseFrames)  {
            delete[] preBaseFrames;
        }
        if (widthInFrames) {
            delete[] widthInFrames;
        }
        if (pulseIndex) {
            delete[] pulseIndex;
        }
        if (startFrame) {
            delete[] startFrame;
        }
        // meanSignal, maxSignal, midSignal and classifierQV
        // need to be handled separatedly.
    }

    // Reset SMRT QV pointers anyway
    preBaseFrames = NULL;
    widthInFrames = NULL;
    pulseIndex = NULL;
    startFrame = NULL;

    // Reset member variables
    xy[0] = 0; xy[1] = 0;
    lowQualityPrefix = lowQualitySuffix = 0;
    highQualityRegionScore = 0;
    holeNumber = static_cast<UInt>(-1);
    readGroupId = "";
    copiedFromBam = false;
#ifdef USE_PBBAM
    bamRecord = PacBio::BAM::BamRecord();
#endif 

    // Free seq, title and FASTQ QVs, also reset deleteOnExit.
    // Don't call FASTQSequence::Free() before freeing SMRT QVs.
    FASTQSequence::Free();
}

bool SMRTSequence::StoreXY(int16_t xyP[]) {
    xy[0] = xyP[0];
    xy[1] = xyP[1];
    return true;
}

bool SMRTSequence::StorePlatformId(PlatformId pid) {
    platform = pid;
    return true;
}

bool SMRTSequence::StoreHoleNumber(UInt holeNumberP){ 
    zmwData.holeNumber = holeNumber = holeNumberP;
    return true;
}

bool SMRTSequence::StoreHoleStatus(unsigned char s) {
    zmwData.holeStatus = s;
    return true;
}

bool SMRTSequence::StoreZMWData(ZMWGroupEntry &data) {
    zmwData = data;
    return true;
}

bool SMRTSequence::GetXY(int xyP[]) {
    xyP[0] = xy[0];
    xyP[1] = xy[1];
    return true;
}

bool SMRTSequence::GetHoleNumber(UInt & holeNumberP) {
    holeNumberP = holeNumber;
    return true;
}

std::string SMRTSequence::GetReadGroupId() {
    return readGroupId;
}

void SMRTSequence::SetReadGroupId(const std::string & rid) {
    readGroupId = rid;
}

#ifdef USE_PBBAM
void SMRTSequence::Copy(const PacBio::BAM::BamRecord & record,
                        bool copyAllQVs) {
    Free();

    copiedFromBam = true;

    bamRecord = PacBio::BAM::BamRecord(record);
    
    // Only copy insertionQV, deletionQV, substitutionQV, mergeQV, 
    // deletionTag and substitutionTag from BamRecord to SMRTSequence.
    // Do NOT copy other SMRTQVs such as startFrame, meanSignal...
    (static_cast<FASTQSequence*>(this))->Copy(record);

    if (copyAllQVs) {
        if (record.HasPreBaseFrames()) {
            std::vector<uint16_t> qvs = record.PreBaseFrames().DataRaw();
            assert(preBaseFrames == nullptr);
            preBaseFrames = new HalfWord[qvs.size()];
            std::memcpy(preBaseFrames, &qvs[0], qvs.size() * sizeof(HalfWord));
        }
    }

    // preBaseQVs are not included in BamRecord, and will not be copied.
    // Copy read group id from BamRecord.
    SetReadGroupId(record.ReadGroupId());

    // PacBio bam for secondary analysis does NOT carry zmw
    // info other than holeNumber, including holeStatus, holeX,
    // holeY, numEvents. 
    zmwData.holeNumber = static_cast<UInt> (record.HoleNumber()); 

    // Set hq region read score
    if (record.Impl().HasTag("rq")) {
        highQualityRegionScore = record.Impl().TagValue("rq").ToInt32();
        readScore = highQualityRegionScore / 1000.0;
    }

    // Set HQRegionSNR if record has the 'sn' tag 
    if (record.HasSignalToNoise()) {
        // Signal to noise ratio of ACGT (in that particular ORDER) over 
        // HQRegion from BAM: record.SignalToNoise() 
        std::vector<float> snrs = record.SignalToNoise();
        this->HQRegionSnr('A', snrs[0]) 
             .HQRegionSnr('C', snrs[1])
             .HQRegionSnr('G', snrs[2]) 
             .HQRegionSnr('T', snrs[3]);
    }
}

#endif

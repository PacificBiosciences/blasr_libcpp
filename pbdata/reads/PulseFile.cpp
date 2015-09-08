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

#include "PulseFile.hpp"

void PulseFile::CopySignal(HalfWord *signalData, // either a vector or matrix
                            int signalNDims,
                            int pulseStartPos,    // 0 if baseToPulseIndex maps to abs position
                            int *baseToPulseIndex,
                            Nucleotide *readSeq,
                            int readLength,
                            HalfWord *readData) {
    // if baseToPulseIndex maps bases to absolute pulse postions
    // pulseStartPos must be 0; 
    // otherwise, pulseStartPos is pulseStartPositions[holeIndex]

    std::map<char, size_t> baseMap = GetBaseMap();
    int i;
    if (signalNDims == 1) {
        for (i = 0; i < readLength; i++) {
            readData[i] = signalData[pulseStartPos + baseToPulseIndex[i] ];
        }
    }
    else {
        for (i = 0; i < readLength; i++) {
            readData[i] = signalData[baseToPulseIndex[i]*4 + baseMap[readSeq[i]]];
        }
    }
}

void PulseFile::CopyReadAt(uint32_t plsReadIndex, int *baseToPulseIndex, SMRTSequence &read) {
    int pulseStartPos = pulseStartPositions[plsReadIndex];
    bool allocResult;
    if (midSignal.size() > 0) {
    assert(midSignal.size() > pulseStartPos);
    allocResult = Realloc(read.midSignal, read.length);
    CopySignal(&midSignal[0], 
                midSignalNDims, 
                pulseStartPos,
                baseToPulseIndex,
                read.seq, read.length,
                read.midSignal);
    }

    if (maxSignal.size() > 0) {
        assert(maxSignal.size() > pulseStartPos); 
        Realloc(read.maxSignal, read.length);
        CopySignal(&maxSignal[0], 
                    maxSignalNDims, 
                    pulseStartPos,
                    baseToPulseIndex,
                    read.seq, read.length,
                    read.maxSignal);
    }

    if (meanSignal.size() > 0) {
        assert(meanSignal.size() > pulseStartPos); 
        Realloc(read.meanSignal, read.length);
        CopySignal(&meanSignal[0], 
                    meanSignalNDims, 
                    pulseStartPos,
                    baseToPulseIndex,
                    read.seq, read.length,
                    read.meanSignal);
    }
    if (plsWidthInFrames.size() > 0) {
        Realloc(read.widthInFrames, read.length);
        StoreField(plsWidthInFrames, baseToPulseIndex, read.widthInFrames, read.length);
    }
    if (classifierQV.size() > 0) {
        Realloc(read.classifierQV, read.length);
        StoreField(classifierQV, baseToPulseIndex, read.classifierQV, read.length);
    }
    if (startFrame.size() > 0) {
        Realloc(read.startFrame, read.length);
        StoreField(startFrame, baseToPulseIndex, read.startFrame, read.length);
    }
}

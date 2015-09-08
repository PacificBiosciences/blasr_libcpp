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

#ifndef DATASTRUCTURES_READS_SCAN_DATA_H_
#define DATASTRUCTURES_READS_SCAN_DATA_H_

#include <string>
#include <map>
#include "Enumerations.h"
#include "PacBioDefs.h"

class HDFScanDataReader;
class HDFScanDataWriter;

class ScanData {
friend class HDFScanDataReader;
friend class HDFScanDataWriter;
public:
    // Convert base map from a map<char, size_t> to a string. 
    // e.g., {{'A', 2}, {'C', 1}, {'T', 0}, {'G', 3}} --> TCAG
    static std::string BaseMapToStr(const std::map<char, size_t> & baseMap);

    // Convert base map from a string to a map<char, size_t>.
    // e.g., TCAG --> {{'A', 2}, {'C', 1}, {'T', 0}, {'G', 3}} 
    static std::map<char, size_t> StrToBaseMap(const std::string & baseMapStr);

    // A baseMap must contain maps from bases (e.g., ACGT) to indices (e.g., 0, 1, 2, 3). 
    static bool IsValidBaseMap(const std::map<char, size_t> & baseMap);
    
public:
    PlatformId platformId;
    float frameRate;
    unsigned int numFrames;
    std::string movieName, runCode;
    std::string whenStarted;
    std::map<char, size_t> baseMap;

    ScanData();
    std::string GetMovieName(); 

    ScanData & PlatformID(const PlatformId & id);
    ScanData & FrameRate(const float & rate);
    ScanData & NumFrames(const unsigned int & num);
    ScanData & MovieName(const std::string & name);
    ScanData & RunCode(const std::string & code);
    ScanData & WhenStarted(const std::string & when);
    ScanData & BaseMap(const std::map<char, size_t> & bmp);
    ScanData & BaseMap(const std::string & baseMapStr); 
    ScanData & SequencingKit(const std::string sequencingKit);
    ScanData & BindingKit(const std::string bindingKit);

    PlatformId PlatformID(void) const;
    float FrameRate(void) const;
    unsigned int NumFrames(void) const;
    std::string MovieName(void) const;
    std::string RunCode(void) const;
    std::string WhenStarted(void) const;
    std::map<char, size_t> BaseMap(void) const;
    std::string BaseMapStr(void) const;

    std::string SequencingKit(void) const;
    std::string BindingKit(void) const;
 

private:
    std::string sequencingKit_;
    std::string bindingKit_;
};

#endif

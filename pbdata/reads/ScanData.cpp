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

#include "ScanData.hpp"
#include <iostream>

std::string ScanData::BaseMapToStr(const std::map<char, size_t> & baseMap) {
    std::string baseMapStr = ""; //4 dye channels.
    if (not baseMap.empty()) {
        baseMapStr = "    ";
        for (auto it = baseMap.begin(); it != baseMap.end(); ++it){
            if (it->second != 4) {
                std::cout << "ERROR, must have exactly four dye channels."
                          << std::endl;
                exit(1);
            }
            baseMapStr[it->second]= it->first;
        }
    }
    return baseMapStr;
}

std::map<char, size_t> ScanData::StrToBaseMap(const std::string & baseMapStr) {
    std::map<char, size_t> ret;
    for (size_t i = 0; i < baseMapStr.size(); i++) {
        ret[baseMapStr[i]] = i;
    }
    return ret;
}

bool ScanData::IsValidBaseMap(const std::map<char, size_t> & baseMap) {
    const char X = 'x';
    std::string v(4, X);

    for(const char base : {'A', 'T', 'G', 'C'}) {
        size_t index = baseMap.find(base)->second;
        if (not (baseMap.find(base) != baseMap.end() and index <= 3))
            return false;
        else
            v[index] = 'o';
    }
    if (v.find(X) != std::string::npos) return false;
    else return true;
}

ScanData::ScanData(const AcqParams & acqParams) 
    : acqParams_(acqParams) 
{
    platformId = NoPlatform;
    frameRate = 0.0;
    numFrames = 0;
    movieName = runCode = whenStarted = "";
    baseMap.clear();
}

ScanData & ScanData::SetAcqParams(const AcqParams & acqParams) {
    acqParams_ = acqParams;
    return *this;
}

AcqParams ScanData::GetAcqParams(void) const {
    return acqParams_;
}

std::string ScanData::GetMovieName() {
    return movieName;
}

ScanData & ScanData::PlatformID(const PlatformId & id) {
    platformId = id;
    return *this;
}
ScanData & ScanData::FrameRate(const float & rate) {
    frameRate = rate;
    return *this;
}
ScanData & ScanData::NumFrames(const unsigned int & num) {
    numFrames = num;
    return *this;
}
ScanData & ScanData::MovieName(const std::string & name) {
    movieName = name;
    return *this;
}
ScanData & ScanData::RunCode(const std::string & code) {
    runCode = code;
    return *this;
}
ScanData & ScanData::WhenStarted(const std::string & when) {
    whenStarted = when;
    return *this;
}
ScanData & ScanData::BaseMap(const std::map<char, size_t> & bmp) {
    baseMap.clear();
    baseMap.insert(bmp.begin(), bmp.end());
    return *this;
}
ScanData & ScanData::BaseMap(const std::string & baseMapStr) {
   return this->BaseMap(ScanData::StrToBaseMap(baseMapStr));
}
ScanData & ScanData::SequencingKit(const std::string sequencingKit) {
    sequencingKit_ = sequencingKit;
    return *this;
}
ScanData & ScanData::BindingKit(const std::string bindingKit) {
    bindingKit_ = bindingKit;
    return *this;
}

PlatformId ScanData::PlatformID(void) const {
    return platformId;
}
float ScanData::FrameRate(void) const {
    return frameRate;
}
unsigned int ScanData::NumFrames(void) const {
    return numFrames;
}
std::string ScanData::MovieName(void) const {
    return movieName;
}
std::string ScanData::RunCode(void) const {
    return runCode;
}
std::string ScanData::WhenStarted(void) const {
    return whenStarted;
}
std::map<char, size_t> ScanData::BaseMap(void) const {
    return baseMap;
}

std::string ScanData::BaseMapStr(void) const {
    return ScanData::BaseMapToStr(baseMap);
}

std::string ScanData::SequencingKit(void) const {
    return sequencingKit_;
}
std::string ScanData::BindingKit(void) const {
    return bindingKit_;
}

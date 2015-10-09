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

// Author: Yuan Li

#include "HDFZMWWriter.hpp"

HDFZMWWriter::HDFZMWWriter(const std::string & filename, 
        HDFGroup & parentGroup, 
        bool hasHoleXY)
    : HDFWriterBase(filename)
    , parentGroup_(parentGroup)
    , hasHoleXY_(hasHoleXY)
{
    if (not parentGroup.groupIsInitialized)
        PARENT_GROUP_NOT_INITIALIZED_ERROR(PacBio::GroupNames::zmw);
    else {
        parentGroup_.AddGroup(PacBio::GroupNames::zmw); 

        if (zmwGroup_.Initialize(parentGroup_, PacBio::GroupNames::zmw) == 0)
            FAILED_TO_CREATE_GROUP_ERROR(PacBio::GroupNames::zmw);

        this->InitializeChildHDFGroups();
    }
}

HDFZMWWriter::~HDFZMWWriter() {
    this->_WriteAttributes();
    this->Close();
}

bool HDFZMWWriter::WriteOneZmw(const SMRTSequence & read) {
    int length_ = static_cast<int> (read.length);
    numEventArray_.Write(&length_, 1);

    UInt hn_ = read.HoleNumber(); 
    holeNumberArray_.Write(&hn_, 1);

    unsigned char hs_ = read.HoleStatus();
    holeStatusArray_.Write(&hs_, 1);

    if (HasHoleXY()) {
        int16_t xy[2] = {static_cast<int16_t>(read.HoleX()), 
                         static_cast<int16_t>(read.HoleY())};
        holeXYArray_.WriteRow(xy, 2);
    }
    return true;
}

void HDFZMWWriter::Flush(void) {
    numEventArray_.Flush();
    holeNumberArray_.Flush();
    holeStatusArray_.Flush();
    if (HasHoleXY())
        holeXYArray_.Flush();
}

void HDFZMWWriter::Close(void) {
    this->Flush();

    numEventArray_.Close();
    holeNumberArray_.Close();
    holeStatusArray_.Close();
    if (HasHoleXY())
        holeXYArray_.Close();
    zmwGroup_.Close();
}

bool HDFZMWWriter::InitializeChildHDFGroups(void) {
    bool OK = true;

    if (numEventArray_.Initialize(zmwGroup_, PacBio::GroupNames::numevent) == 0) { 
        FAILED_TO_CREATE_GROUP_ERROR(PacBio::GroupNames::numevent);
        OK = false;
    }

    if (holeNumberArray_.Initialize(zmwGroup_, PacBio::GroupNames::holenumber) == 0) {
        FAILED_TO_CREATE_GROUP_ERROR(PacBio::GroupNames::holenumber);
        OK = false;
    }

    if (holeStatusArray_.Initialize(zmwGroup_, PacBio::GroupNames::holestatus) == 0) {
        FAILED_TO_CREATE_GROUP_ERROR(PacBio::GroupNames::holestatus);
        OK = false;
    }

    if (HasHoleXY()) {
        if (holeXYArray_.Initialize(zmwGroup_, PacBio::GroupNames::holexy, 2) == 0) {
            FAILED_TO_CREATE_GROUP_ERROR(PacBio::GroupNames::holexy);
            OK = false;
        }
    }

    return OK;
}

void HDFZMWWriter::_WriteAttributes(void)
{
    if (holeNumberArray_.IsInitialized() and holeNumberArray_.size() > 0) {
        AddAttribute(holeNumberArray_, PacBio::AttributeNames::Common::description, PacBio::AttributeValues::ZMW::HoleNumber::description);
    }

    if (holeStatusArray_.IsInitialized() and holeStatusArray_.size() > 0) {
        AddAttribute(holeStatusArray_, PacBio::AttributeNames::Common::description, PacBio::AttributeValues::ZMW::HoleStatus::description);
        AddAttribute(holeStatusArray_, PacBio::AttributeNames::ZMW::HoleStatus::lookuptable, PacBio::AttributeValues::ZMW::HoleStatus::lookuptable);
    }

    if (holeXYArray_.IsInitialized()) {
        AddAttribute(holeXYArray_, PacBio::AttributeNames::Common::description, PacBio::AttributeValues::ZMW::HoleXY::description);
    }
}


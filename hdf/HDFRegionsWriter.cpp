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

#include "HDFRegionsWriter.hpp"

HDFRegionsWriter::HDFRegionsWriter(const std::string & filename,
                                   HDFGroup & parentGroup,
                                   const std::vector<std::string> & regionTypes) 
    : HDFWriterBase(filename)
    , parentGroup_(parentGroup)
    , regionTypes_(regionTypes)
    , curRow_(0)
{
    // Initialize the 'regions' group.
    regionsArray_.Initialize(parentGroup_, PacBio::GroupNames::regions, RegionAnnotation::NCOLS);
}

HDFRegionsWriter::~HDFRegionsWriter(void)
{
    WriteAttributes();
    Close();
}

bool HDFRegionsWriter::WriteAttributes(void) 
{
    if (curRow_ > 0) {
        AddAttribute(regionsArray_, PacBio::AttributeNames::Regions::columnnames, PacBio::AttributeValues::Regions::columnnames);
        AddAttribute(regionsArray_, PacBio::AttributeNames::Regions::regiontypes, regionTypes_);
        AddAttribute(regionsArray_, PacBio::AttributeNames::Regions::regiondescriptions, PacBio::AttributeValues::Regions::regiondescriptions);
        AddAttribute(regionsArray_, PacBio::AttributeNames::Regions::regionsources, PacBio::AttributeValues::Regions::regionsources);
        return true;
    } else {
        AddErrorMessage("Could not write attributes when Regions group is empty.");
        return false;
    }
}

bool HDFRegionsWriter::Write(const std::vector<RegionAnnotation> &annotations) {
    for (auto annotation: annotations)
        if (not Write(annotation))
            return false;
    return true;
}

bool HDFRegionsWriter::Write(const RegionAnnotation &annotation) {
    try {
        regionsArray_.WriteRow(annotation.row, HDFRegionsWriter::NCOLS);
    }
    catch (H5::Exception &e) {
        AddErrorMessage("Failed to write region annotation " + 
                annotation.GetHoleNumber());
        return false;
    }
    ++curRow_;
    return true;
}	

void HDFRegionsWriter::Flush(void) {
    regionsArray_.Flush();
}

void HDFRegionsWriter::Close(void) {
    Flush();
    regionsArray_.Close();
}

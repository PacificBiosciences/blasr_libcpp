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

#include "HDFBaxWriter.hpp"

HDFBaxWriter::HDFBaxWriter(const std::string & filename,
                           const ScanData & sd,
                           const std::string & basecallerVersion,
                           const std::vector<std::string> & qvsToWrite,
                           const std::vector<std::string> & regionTypes,
                           const H5::FileAccPropList & fileAccPropList)
    : HDFWriterBase(filename)
    , fileaccproplist_(fileAccPropList)
    , scandataWriter_(nullptr)
    , basecallsWriter_(nullptr) 
    , regionsWriter_(nullptr)
{
    // sanity check chemistry meta data. 
    SanityCheckChemistry(sd.BindingKit(),
                         sd.SequencingKit(), 
                         basecallerVersion);

    // open file 
    outfile_.Open(filename_, H5F_ACC_TRUNC, fileaccproplist_);

    // Add PulseData group to the root group '/'
    AddChildGroup(outfile_.rootGroup, pulseDataGroup_, PacBio::GroupNames::pulsedata);

    // Create a ScanData writer.
    scandataWriter_.reset(new HDFScanDataWriter(outfile_.rootGroup)); 
    scandataWriter_->Write(sd);

    // Create a BaseCaller writer.
    basecallsWriter_.reset(new HDFBaseCallsWriter(filename_, pulseDataGroup_, sd.BaseMap(), qvsToWrite));
    basecallsWriter_->WriteBaseCallerVersion(basecallerVersion);

    // Create a Regions writer.
    regionsWriter_.reset(new HDFRegionsWriter(filename_, pulseDataGroup_, regionTypes));
}

HDFBaxWriter::~HDFBaxWriter(void) {
    this->Close();
}

void HDFBaxWriter::Flush(void) {
    basecallsWriter_->Flush();
    regionsWriter_->Flush();
}

void HDFBaxWriter::Close(void) {
    basecallsWriter_->Close();
    scandataWriter_->Close();
    regionsWriter_->Close();
    outfile_.Close();
}

bool HDFBaxWriter::SanityCheckChemistry(const std::string & bindingKit,
                                        const std::string & sequencingKit,
                                        const std::string & basecallerVersion)
{
    bool OK = true;
    if (bindingKit.empty()) {
        OK = false;
        AddErrorMessage("Binding kit must be specified.");
    }
    if (sequencingKit.empty()) {
        OK = false;
        AddErrorMessage("Sequencing kit must be specified.");
    }
    if (basecallerVersion.empty()) {
        OK = false;
        AddErrorMessage("Base caller version must be specified.");
    }
    return OK;
}

bool HDFBaxWriter::WriteOneZmw(const SMRTSequence & seq) {
    return basecallsWriter_->WriteOneZmw(seq);
}

bool HDFBaxWriter::WriteOneZmw(const SMRTSequence & seq, 
                               const std::vector<RegionAnnotation> & regions) {
    if (not this->WriteOneZmw(seq)) {
        return false;
    }
    if (regions.size() == 0) {
        std::vector<RegionAnnotation> fake = {RegionAnnotation(seq.HoleNumber(), HQRegion, 0, 0, 0)};
        return regionsWriter_->Write(fake);
    } else {
        return regionsWriter_->Write(regions);
    }
}

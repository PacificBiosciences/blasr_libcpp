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

#include "HDFBaseCallsWriter.hpp"

HDFBaseCallsWriter::HDFBaseCallsWriter(const std::string & filename,
                                       HDFGroup & parentGroup,
                                       const std::map<char, size_t> & baseMap,
                                       const std::vector<std::string> & qvsToWrite)
    : HDFWriterBase(filename)
    , parentGroup_(parentGroup)
    , baseMap_(baseMap)
    , qvsToWrite_({}) // Input qvsToWrite must be checked.
    , zmwWriter_(nullptr)
    , zmwMetricsWriter_(nullptr)
{
    // Add BaseCalls as a child group to the parent group.
    AddChildGroup(parentGroup_, basecallsGroup_, PacBio::GroupNames::basecalls);

    // Initialize the 'basecall' group.
    basecallArray_.Initialize(basecallsGroup_, PacBio::GroupNames::basecall);

    // Sanity check QVs to write.
    if (SanityCheckQVs(qvsToWrite)) {
        // Initialize QV groups
        if (not InitializeQVGroups()) {
            AddErrorMessage("Failed to initialize QV Groups.");
        }
    }

    // Create a zmwWriter.
    zmwWriter_.reset(new HDFZMWWriter(Filename(), basecallsGroup_, true));

    // Create a zmwMetricsWriter.
    zmwMetricsWriter_.reset(new HDFZMWMetricsWriter(Filename(), basecallsGroup_, baseMap_));
}

std::vector<std::string> HDFBaseCallsWriter::Errors(void) const {
    std::vector<std::string> retErrors = this->errors_;
    std::vector<std::string> zmwErrors = zmwWriter_->Errors();
    std::vector<std::string> zmwMetricsErrors = zmwMetricsWriter_->Errors();

    retErrors.insert(retErrors.end(), zmwErrors.begin(), zmwErrors.end());
    retErrors.insert(retErrors.end(), zmwMetricsErrors.begin(), zmwMetricsErrors.end());
    return retErrors;
}

HDFBaseCallsWriter::~HDFBaseCallsWriter(void) {
    this->Close();
}

const std::vector<std::string> & HDFBaseCallsWriter::QVNamesToWrite(void) const {
    return qvsToWrite_;
}

const std::vector<std::string> & HDFBaseCallsWriter::ValidQVNames(void) const {
    return PacBio::GroupNames::BaxQVNames;
}

bool HDFBaseCallsWriter::InitializeQVGroups(void) {
    int ret = 1;
    if (_HasQV(PacBio::GroupNames::deletionqv)) 
        ret *= deletionQVArray_.Initialize(basecallsGroup_,      PacBio::GroupNames::deletionqv);
    if (_HasQV(PacBio::GroupNames::deletiontag))
        ret *= deletionTagArray_.Initialize(basecallsGroup_,     PacBio::GroupNames::deletiontag);
    if (_HasQV(PacBio::GroupNames::insertionqv))
        ret *= insertionQVArray_.Initialize(basecallsGroup_,     PacBio::GroupNames::insertionqv);
    if (_HasQV(PacBio::GroupNames::mergeqv))
        ret *= mergeQVArray_.Initialize(basecallsGroup_,         PacBio::GroupNames::mergeqv);
    if (_HasQV(PacBio::GroupNames::substitutionqv))
        ret *= substitutionQVArray_.Initialize(basecallsGroup_,  PacBio::GroupNames::substitutionqv);
    if (_HasQV(PacBio::GroupNames::substitutiontag))
        ret *= substitutionTagArray_.Initialize(basecallsGroup_, PacBio::GroupNames::substitutiontag);
    if (_HasQV(PacBio::GroupNames::prebaseframes))
        ret *= preBaseFramesArray_.Initialize(basecallsGroup_,   PacBio::GroupNames::prebaseframes);
    if (_HasQV(PacBio::GroupNames::widthinframes))
        ret *= widthInFramesArray_.Initialize(basecallsGroup_,   PacBio::GroupNames::widthinframes);
    return (ret != 0);
}

bool HDFBaseCallsWriter::SanityCheckQVs(const std::vector<std::string> & qvsToWrite) {
    bool allQVsToAddInSpec = true; 
    qvsToWrite_.clear();
    // Filter qvs which are not in format specification.
    const std::vector<std::string> & qvsInSpec = ValidQVNames();
    for(auto qv : qvsToWrite) {
        if (std::find(qvsInSpec.begin(), qvsInSpec.end(), qv) != qvsInSpec.end()) {
            if (std::find(qvsToWrite_.begin(), qvsToWrite_.end(), qv) == qvsToWrite_.end()) 
                qvsToWrite_.push_back(qv);
            // else redundant
        } else {
            allQVsToAddInSpec = false;
            AddErrorMessage(std::string("Unsupported quality value ") + qv);
        }
    }
    return allQVsToAddInSpec;
}

bool HDFBaseCallsWriter::WriteBaseCallerVersion(const std::string & basecallerVersion) {
    changeListIDAtom_.Create(basecallsGroup_.group, 
                             PacBio::AttributeNames::Common::changelistid,
                             basecallerVersion);
    return true;
}

bool HDFBaseCallsWriter::WriteOneZmw(const SMRTSequence & read) {
    bool OK = zmwWriter_->WriteOneZmw(read);
    OK = OK and zmwMetricsWriter_->WriteOneZmw(read);
    OK = OK and _WriteBasecall(read);
    OK = OK and _WriteDeletionQV(read);
    OK = OK and _WriteDeletionTag(read);
    OK = OK and _WriteInsertionQV(read);
    OK = OK and _WriteMergeQV(read);
    OK = OK and _WriteSubstitutionTag(read);
    OK = OK and _WriteSubstitutionQV(read);
    OK = OK and _WritePreBaseFrames(read);
    OK = OK and _WriteWidthInFrames(read);
    return OK;
}

bool HDFBaseCallsWriter::_WriteBasecall(const SMRTSequence & read) {
	basecallArray_.Write((const unsigned char*) read.seq, read.length);
    return true;
}

bool HDFBaseCallsWriter::_WriteDeletionQV(const SMRTSequence & read) {
    if (HasDeletionQV()) {
        if (read.deletionQV.Empty()) {
            AddErrorMessage(std::string(PacBio::GroupNames::deletionqv) + " absent in read " + read.GetTitle());
            return false;
        } else {
            deletionQVArray_.Write(read.deletionQV.data, read.length);
            return true;
        }
    }
    return true;
}

bool HDFBaseCallsWriter::_WriteDeletionTag(const SMRTSequence & read) {
    if (HasDeletionTag()) {
        if (read.deletionTag == nullptr) {
            AddErrorMessage(std::string(PacBio::GroupNames::deletiontag) + " absent in read " + read.GetTitle());
            return false;
        } else {
            deletionTagArray_.Write(read.deletionTag, read.length);
            return true;
        }
    }
    return true;
}

bool HDFBaseCallsWriter::_WriteInsertionQV(const SMRTSequence & read) {
    if (HasInsertionQV()) {
        if (read.insertionQV.Empty()) {
            AddErrorMessage(std::string(PacBio::GroupNames::insertionqv) + " absent in read " + read.GetTitle());
            return false;
        } else {
			insertionQVArray_.Write(read.insertionQV.data, read.length);
            return true;
        }
    }
    return true;
}

bool HDFBaseCallsWriter::_WriteSubstitutionTag(const SMRTSequence & read) {
    if (HasSubstitutionTag()) {
        if (read.substitutionTag == nullptr) {
            AddErrorMessage(std::string(PacBio::GroupNames::substitutiontag) + " absent in read " + read.GetTitle());
            return false;
        } else {
			substitutionTagArray_.Write(read.substitutionTag, read.length);
            return true;
        }
    }
    return true;
}

bool HDFBaseCallsWriter::_WriteSubstitutionQV(const SMRTSequence & read) {
    if (HasSubstitutionQV()) {
        if (read.substitutionQV.Empty()) {
            AddErrorMessage(std::string(PacBio::GroupNames::substitutionqv) + " absent in read " + read.GetTitle());
            return false;
        } else {
			substitutionQVArray_.Write(read.substitutionQV.data, read.length);
            return true;
        }
    }
    return true;
}

bool HDFBaseCallsWriter::_WriteMergeQV(const SMRTSequence & read) {
    if (HasMergeQV()) {
        if (read.mergeQV.Empty()) {
            AddErrorMessage(std::string(PacBio::GroupNames::mergeqv) + " absent in read " + read.GetTitle());
            return false;
        } else {
			mergeQVArray_.Write(read.mergeQV.data, read.length);
            return true;
        }
    }
    return true;
}

bool HDFBaseCallsWriter::_WritePreBaseFrames(const SMRTSequence & read) {
     if (HasPreBaseFrames()) {
        if (read.preBaseFrames == nullptr) {
            AddErrorMessage(std::string(PacBio::GroupNames::prebaseframes) + " absent in read " + read.GetTitle());
            return false;
        } else {
            preBaseFramesArray_.Write(read.preBaseFrames, read.length);
            return true;
        }
    }
    return true;
}

bool HDFBaseCallsWriter::_WriteWidthInFrames(const SMRTSequence & read) {
    if (HasWidthInFrames()) {
        if (read.widthInFrames == nullptr) {
            AddErrorMessage(std::string(PacBio::GroupNames::widthinframes) + " absent in read " + read.GetTitle());
            return false;
        } else {
            widthInFramesArray_.Write(read.widthInFrames, read.length);
            return true;
        }
    }
    return true;
}
 

void HDFBaseCallsWriter::Flush(void) {
    basecallArray_.Flush();

    if (HasDeletionQV())      deletionQVArray_.Flush();
    if (HasDeletionTag())     deletionTagArray_.Flush();
    if (HasInsertionQV())     insertionQVArray_.Flush();
    if (HasMergeQV())         mergeQVArray_.Flush();
    if (HasSubstitutionQV())  substitutionQVArray_.Flush();
    if (HasSubstitutionTag()) substitutionTagArray_.Flush();
    if (HasPreBaseFrames())   preBaseFramesArray_.Flush();
    if (HasWidthInFrames())   widthInFramesArray_.Flush();

    zmwWriter_->Flush();
    zmwMetricsWriter_->Flush();
}

void HDFBaseCallsWriter::Close(void) {
    this->Flush();

    basecallArray_.Close();

    if (HasDeletionQV())      deletionQVArray_.Close();
    if (HasDeletionTag())     deletionTagArray_.Close();
    if (HasInsertionQV())     insertionQVArray_.Close();
    if (HasMergeQV())         mergeQVArray_.Close();
    if (HasSubstitutionQV())  substitutionQVArray_.Close();
    if (HasSubstitutionTag()) substitutionTagArray_.Close();
    if (HasPreBaseFrames())   preBaseFramesArray_.Close();
    if (HasWidthInFrames())   widthInFramesArray_.Close();
}

// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.  //
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

#include "HDFWriterBase.hpp"

std::vector<std::string> HDFWriterBase::Errors(void) const {
    return errors_;
}

bool HDFWriterBase::AddChildGroup(HDFGroup & parentGroup, 
                                  HDFGroup & childGroup,
                                  const std::string & childGroupName) {
    parentGroup.AddGroup(childGroupName);
    if (childGroup.Initialize(parentGroup, childGroupName) == 0) {
        FAILED_TO_CREATE_GROUP_ERROR(childGroupName);
        return false;
    }
    return true;
}

bool HDFWriterBase::AddAttribute(HDFData & group, 
                                 const std::string & attributeName, 
                                 const std::vector<std::string> & attributeValues)
{
    try {
        HDFAtom<std::vector<std::string> > attributeAtom;
        attributeAtom.Create(group.dataset, std::string(attributeName), attributeValues);
        attributeAtom.Close();
    }
    catch (H5::Exception &e) {
        FAILED_TO_CREATE_ATTRIBUTE_ERROR(attributeName);
        return false;
    }
    return true;
}

void HDFWriterBase::AddErrorMessage(const std::string & errmsg) {
    errors_.push_back(errmsg);
}

void HDFWriterBase::FAILED_TO_CREATE_GROUP_ERROR(const std::string & groupName) {
    std::stringstream ss;
    ss << "Failed to create group " << groupName << " in " << filename_;
    AddErrorMessage(ss.str());
}

void HDFWriterBase::FAILED_TO_CREATE_ATTRIBUTE_ERROR(const std::string & attributeName) {
    std::stringstream ss;
    ss << "Failed to create attribute " << attributeName << " in " << filename_;
    AddErrorMessage(ss.str());
}

void HDFWriterBase::PARENT_GROUP_NOT_INITIALIZED_ERROR(const std::string & groupName) {
    std::stringstream ss;
    ss << "Parent hdf group of " << groupName << " in file " << filename_
       << " is not initialized.";
    AddErrorMessage(ss.str());
}

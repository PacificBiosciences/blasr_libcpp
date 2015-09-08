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

#ifndef _BLASR_HDFWRITERBASE_HPP_
#define _BLASR_HDFWRITERBASE_HPP_
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include "HDFGroup.hpp"
#include "HDFAtom.hpp"

class HDFWriterBase {
public:
    HDFWriterBase(const std::string & filename)
    : filename_(filename)
    {}

    ~HDFWriterBase() {}

public:
    /// \returns Target H5 filename.
    std::string Filename(void) {return filename_;}

    std::vector<std::string> Errors(void) const;

protected:
    std::string filename_;
    std::vector<std::string> errors_; 

    bool AddChildGroup(HDFGroup & parentGroup, 
                       HDFGroup & childGroup,
                       const std::string & childGroupName);

    bool AddAttribute(HDFData & group, 
                      const std::string & attributeName, 
                      const std::vector<std::string> & attributeValues);
    
    void AddErrorMessage(const std::string & errmsg);

    void FAILED_TO_CREATE_GROUP_ERROR(const std::string & groupName);
   
    void FAILED_TO_CREATE_ATTRIBUTE_ERROR(const std::string & attributeName);

    void PARENT_GROUP_NOT_INITIALIZED_ERROR(const std::string & groupName);

    virtual void Close(void) = 0;
};

#endif

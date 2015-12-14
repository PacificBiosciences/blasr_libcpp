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
#include "BufferedHDFArray.hpp"
#include "BufferedHDF2DArray.hpp"

class HDFWriterBase {
public:
    HDFWriterBase(const std::string & filename)
    : filename_(filename)
    {}

    virtual ~HDFWriterBase(void) = 0;

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

    bool AddAttribute(HDFGroup & group, 
                      const std::string & attributeName, 
                      const std::vector<std::string> & attributeValues);

    template<typename T>
    bool AddAttribute(HDFData & group, 
                      const std::string & attributeName, 
                      const T & attributeValue);

    template<typename T>
    bool AddAttribute(HDFGroup & group, 
                      const std::string & attributeName, 
                      const T & attributeValue);

    void AddErrorMessage(const std::string & errmsg);

    void FAILED_TO_CREATE_GROUP_ERROR(const std::string & groupName);
   
    void FAILED_TO_CREATE_ATTRIBUTE_ERROR(const std::string & attributeName);

    void PARENT_GROUP_NOT_INITIALIZED_ERROR(const std::string & groupName);

    virtual void Close(void) = 0;
};


template<typename T>
bool HDFWriterBase::AddAttribute(HDFData & group, 
                                 const std::string & attributeName, 
                                 const T & attributeValue) {
    try {
        HDFAtom<T> attributeAtom;
        attributeAtom.Create(group.dataset, std::string(attributeName));
        attributeAtom.Write(attributeValue);
        attributeAtom.Close();
    }
    catch (H5::Exception &e) {
        this->FAILED_TO_CREATE_ATTRIBUTE_ERROR(attributeName);
        return false;
    }
    return true;
}

template<typename T>
bool HDFWriterBase::AddAttribute(HDFGroup & group, 
                    const std::string & attributeName, 
                    const T & attributeValue) {
    try {
        HDFAtom<T> attributeAtom;
        attributeAtom.Create(group.group, std::string(attributeName));
        attributeAtom.Write(attributeValue);
        attributeAtom.Close();
    }
    catch (H5::Exception &e) {
        FAILED_TO_CREATE_ATTRIBUTE_ERROR(attributeName);
        return false;
    }
    return true;
}

/// \brief Write a dataset of name 'dsName' under h5 group 'dsGroup' with
///        length 'dsLength'. Just fill this dataset with buffer repeatedly.
template<typename T>
bool __WriteFakeDataSet(HDFGroup & dsGroup, const std::string & dsName,
                        const uint32_t dsLength, std::vector<T> & buffer) {
    BufferedHDFArray<T> dsArray_;
    if (dsArray_.Initialize(dsGroup, dsName) == 0) return false;
    uint32_t totalLength = 0; 
    while(totalLength < dsLength) {
        uint32_t thisLength = buffer.size();
        if (totalLength + thisLength <= dsLength) {
            totalLength = totalLength + thisLength;
        } else {
            thisLength = dsLength - totalLength;
            totalLength = dsLength;
        }
        dsArray_.Write(&buffer[0], thisLength);
        dsArray_.Flush();
    }
    dsArray_.Close();
    return true;
}

/// \brief Write a 2D dataset of name 'dSName' under h5 group 'dsGroup' with
///        row number 'rowNum', column number 'colNum'. Just fill this dataset
///        with value 'fillData'
template<typename T>
bool __WriteFake2DDataSet(HDFGroup & dsGroup, const std::string & dsName,
                           const uint32_t rowNum, const int colNum, const T & fillData) {
    BufferedHDF2DArray<T> dsArray;
    if (dsArray.Initialize(dsGroup, dsName, colNum) == 0) return false;
    T * data = new T[colNum];
    for(size_t i = 0; i < colNum; i++) {data[i] = fillData; }
    for(uint32_t i = 0; i< rowNum; i++) {
        dsArray.WriteRow(data, colNum);
    }
    dsArray.Close();
    delete []data;
    return true;
}

#endif

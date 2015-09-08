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


#ifndef _BLASR_HDF_HDFZMWWriter_HPP_
#define _BLASR_HDF_HDFZMWWriter_HPP_

#include "HDFWriterBase.hpp"
#include "BufferedHDFArray.hpp"
#include "BufferedHDF2DArray.hpp"
#include "SMRTSequence.hpp"

class HDFBaseCallerWriter;

class HDFZMWWriter: public HDFWriterBase {

friend class HDFBaseCallerWriter;

private:
    /// ZMW/NumEvent
   	BufferedHDFArray<int> numEventArray_;

    // ZMW/HoleNumber
	BufferedHDFArray<unsigned int> holeNumberArray_; 

    // ZMW/HoleStatus
	BufferedHDFArray<unsigned char> holeStatusArray_;

    // ZMW/HoleXY
	BufferedHDF2DArray<int16_t> holeXYArray_;

private:
	HDFGroup zmwGroup_;
    HDFGroup & parentGroup_;
    bool hasHoleXY_;

public:
    /// \name Constructors and Destructors
    /// \{
    HDFZMWWriter(const std::string & filename, 
                 HDFGroup & parentGroup, 
                 bool hasHoleXY = true);

    ~HDFZMWWriter() ;
    /// \}

    /// \name Public Methods
    /// \{
    
    /// \note Write info of a SMRTSequence to ZMW,
    ///       (1) add length (UInt) of the sequence to NumEvent,
    ///       (2) add zmw hole number (UInt) of the sequence as a UInt to HoleNumber,
    ///       (3) add hole status (unsigned char) to HoleStatus,
    ///       (4) add hole coordinate xy as (int16_t, int16_t) to HoleXY
    bool WriteOneZmw(const SMRTSequence & read);

    /// \returns Whether or not ZMW contains the HoleXY dataset.
    inline bool HasHoleXY(void) const;

    /// \note Flushes all data from cache to disc.
    void Flush(void);

    /// \note Closes this zmw group as well as child hdf groups.
    void Close(void);
 
    /// \}

private:
    /// \name Private Methods
    /// \{
    
    /// \note Initialize child hdf groups under ZMW, including
    ///       NumEvent, HoleNumber, HoleStatus, HoleXY
    /// \reutrns bool, whether or not child hdf groups successfully initialized.
    bool InitializeChildHDFGroups(void);

    /// \}
};

inline bool HDFZMWWriter::HasHoleXY(void) const 
{return hasHoleXY_;}

#endif

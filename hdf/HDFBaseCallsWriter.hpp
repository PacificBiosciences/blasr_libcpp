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

#ifndef _BLASR_HDF_BASECALLS_WRITER_HPP_
#define _BLASR_HDF_BASECALLS_WRITER_HPP_

//#include <boost/scoped_ptr.hpp>
#include <algorithm>
#include "HDFAtom.hpp"
#include "HDFWriterBase.hpp"
#include "HDFZMWWriter.hpp"
#include "HDFZMWMetricsWriter.hpp"

class HDFBaseCallsWriter: public HDFWriterBase {
    /// \name \{
public:
    HDFBaseCallsWriter(const std::string & filename,
                       HDFGroup & parentGroup,
                       const std::map<char, size_t> & baseMap,
                       const std::vector<std::string> & qvsToWrite = {});

    ~HDFBaseCallsWriter(void);

    /// \brief Write base caller version (changeListId)
    bool WriteBaseCallerVersion(const std::string & basecallerVersion);

    /// \brief Write a zmw read.
    bool WriteOneZmw(const SMRTSequence & read);

    /// \brief return a vector of QV name strings specified in file format specification.
    const std::vector<std::string> & ValidQVNames(void) const;

    /// \brief return a vector of QV name strings to write.
    const std::vector<std::string> & QVNamesToWrite(void) const;

    void Flush(void);

    void Close(void);

public:
    /// \brief Sanity check QVs to add. Remove QVs which are 
    ///        not included in file format specification, and
    ///        remove redundant QVs.
    /// \returns Whether or not a QV is not included in sepcification.
    bool SanityCheckQVs(const std::vector<std::string> & qvsToWrite);

    inline bool HasDeletionQV(void) const;
    inline bool HasDeletionTag(void) const;
    inline bool HasInsertionQV(void) const;
    inline bool HasSubstitutionTag(void) const;
    inline bool HasSubstitutionQV(void) const;
    inline bool HasMergeQV(void) const;
    inline bool HasPreBaseFrames(void) const;
    inline bool HasIPD(void) const;
    inline bool HasWidthInFrames(void) const;
    inline bool HasPulseWidth(void) const;

    std::vector<std::string> Errors(void) const;

private:
    inline bool _HasQV(const std::string & qvToQuery) const;

    bool _WriteBasecall(const SMRTSequence & read);
    bool _WriteDeletionQV(const SMRTSequence & read);
    bool _WriteDeletionTag(const SMRTSequence & read);
    bool _WriteInsertionQV(const SMRTSequence & read);
    bool _WriteSubstitutionTag(const SMRTSequence & read);
    bool _WriteSubstitutionQV(const SMRTSequence & read);
    bool _WriteMergeQV(const SMRTSequence & read);
    bool _WritePreBaseFrames(const SMRTSequence & read);
    bool _WriteWidthInFrames(const SMRTSequence & read);

private:
    /// \brief Create and initialize QV groups.
    /// \returns Whether or not QV groups initialized successfully.
    bool InitializeQVGroups(void);

private:
    HDFGroup & parentGroup_;
    std::map<char, size_t> baseMap_;
    std::vector<string> qvsToWrite_;
    std::unique_ptr<HDFZMWWriter> zmwWriter_;
    std::unique_ptr<HDFZMWMetricsWriter> zmwMetricsWriter_;
    //boost::scoped_ptr<HDFZMWWriter> zmwWriter_;
    //boost::scoped_ptr<HDFZMWMetricsWriter> zmwMetricsWriter_;
	HDFGroup basecallsGroup_;

private:
	HDFAtom<string> changeListIDAtom_;

    /// BaseCalls/Basecall group
	BufferedHDFArray<unsigned char> basecallArray_;

	/// \brief Define arrays for rich quality values.
    ///        DeletionQV         dq --> BaseCalls/DeletionQV
    ///        DeletionTag        dt --> BaseCalls/DeletionTag
    ///        InsertionQV        iq --> BaseCalls/InsertionQV
    ///        MergeQV            mq --> BaseCalls/MergeQV
    ///        SubstitutionQV     sq --> BaseCalls/SubstitutionQV
    ///        SubstitutionTag    st --> BaseCalls/SubstitutionTag
    ///        Ipd:Frames         ip --> BaseCalls/PreBaseFrames
    ///        PulseWidth:Frames  pw --> BaseCalls/WidthInFrames
	BufferedHDFArray<unsigned char> deletionQVArray_;
	BufferedHDFArray<unsigned char> deletionTagArray_;
	BufferedHDFArray<unsigned char> insertionQVArray_;
	BufferedHDFArray<unsigned char> mergeQVArray_;
	BufferedHDFArray<unsigned char> substitutionQVArray_;
	BufferedHDFArray<unsigned char> substitutionTagArray_;
	BufferedHDFArray<HalfWord> preBaseFramesArray_;
	BufferedHDFArray<HalfWord> widthInFramesArray_;
    /// \}

};

inline
bool HDFBaseCallsWriter::_HasQV(const std::string & qvToQuery) const {
    return (std::find(qvsToWrite_.begin(), qvsToWrite_.end(), qvToQuery) != qvsToWrite_.end());
}

inline
bool HDFBaseCallsWriter::HasDeletionQV(void) const
{return (_HasQV(PacBio::GroupNames::deletionqv) and 
        deletionQVArray_.IsInitialized());}

inline
bool HDFBaseCallsWriter::HasDeletionTag(void) const
{return (_HasQV(PacBio::GroupNames::deletiontag) and
        deletionTagArray_.IsInitialized());}

inline
bool HDFBaseCallsWriter::HasInsertionQV(void) const
{return (_HasQV(PacBio::GroupNames::insertionqv) and
        insertionQVArray_.IsInitialized());}

inline
bool HDFBaseCallsWriter::HasSubstitutionTag(void) const
{return (_HasQV(PacBio::GroupNames::substitutiontag) and
        substitutionTagArray_.IsInitialized());}

inline
bool HDFBaseCallsWriter::HasSubstitutionQV(void) const
{return (_HasQV(PacBio::GroupNames::substitutionqv) and 
        substitutionQVArray_.IsInitialized());}

inline
bool HDFBaseCallsWriter::HasMergeQV(void) const
{return (_HasQV(PacBio::GroupNames::mergeqv) and 
        mergeQVArray_.IsInitialized());}

inline
bool HDFBaseCallsWriter::HasPreBaseFrames(void) const
{return (_HasQV(PacBio::GroupNames::prebaseframes) and
        preBaseFramesArray_.IsInitialized());}

inline
bool HDFBaseCallsWriter::HasIPD(void) const
{return HasPreBaseFrames();}

inline
bool HDFBaseCallsWriter::HasWidthInFrames(void) const
{return (_HasQV(PacBio::GroupNames::widthinframes) and
        widthInFramesArray_.IsInitialized());}

inline
bool HDFBaseCallsWriter::HasPulseWidth(void) const
{return this->HasWidthInFrames();}
#endif

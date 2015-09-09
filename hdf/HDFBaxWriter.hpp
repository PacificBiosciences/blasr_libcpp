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


#ifndef _BLASR_HDF_BAX_WRITER_HPP_
#define _BLASR_HDF_BAX_WRITER_HPP_

#include <sstream>
//#include <boost/scoped_ptr.hpp>
#include "Enumerations.h"
#include "SMRTSequence.hpp"
#include "HDFFile.hpp"
#include "HDFWriterBase.hpp"
#include "HDFScanDataWriter.hpp"
#include "HDFBaseCallsWriter.hpp"
#include "HDFRegionsWriter.hpp"

using namespace H5;
using namespace std;

class HDFBaxWriter : public HDFWriterBase {
public:
    /// \name Constructor & Related Methods
    /// \{
    /// \brief Sets output h5 file name, scan data, base caller version
    ///        QVs to write, and h5 file access property list.
    /// \param[in] filename output h5 file name.
    /// \param[in] ScanData meta data string, must contain bindingKit and sequencingKit.
    /// \param[in] basecallerVersion meta data string
    /// \param[in] qvsToWrite Quality values to include in output h5 file. 
    /// \param[in] regionTypes, regionTypes as /Regions/RegionTypes
    /// \param[in] fileAccPropList H5 file access property list
    HDFBaxWriter(const std::string & filename,
                 const ScanData & sd,
                 const std::string & basecallerVersion,
                 const std::vector<std::string> & qvsToWrite,
                 const std::vector<std::string> & regionTypes = PacBio::AttributeValues::Regions::regiontypes,
                 const H5::FileAccPropList & fileAccPropList = H5::FileAccPropList::DEFAULT);

	~HDFBaxWriter(void);

    /// \brief Write one zmw sequence to output h5 file. 
    /// \param[in] seq, the SMRTSequence to write
    bool WriteOneZmw(const SMRTSequence & seq);

    /// \brief Write one zmw sequence and its region table to output h5 file. 
    /// \param[in] seq, the SMRTSequence to write
    /// \param[in] regions, region annotations of this zmw.
    bool WriteOneZmw(const SMRTSequence & seq, 
                     const std::vector<RegionAnnotation> & regions);

    /// \brief Flushes buffered data.
    void Flush(void);

    /// \}

private:
    /// \name Private Variables
    /// \{
	HDFFile outfile_;  ///< HDFFile file handler

    H5::FileAccPropList fileaccproplist_; ///< H5 file access property list

	HDFGroup pulseDataGroup_; ///< /PulseData group

private:
    /// Points to scan data writer.
    std::unique_ptr<HDFScanDataWriter>  scandataWriter_;
    //boost::scoped_ptr<HDFScanDataWriter>  scandataWriter_;
    /// Points to base caller writer.
    std::unique_ptr<HDFBaseCallsWriter> basecallsWriter_;
    //boost::scoped_ptr<HDFBaseCallsWriter> basecallsWriter_;
    /// Points to region table writer.
    std::unique_ptr<HDFRegionsWriter>   regionsWriter_;
    //boost::scoped_ptr<HDFRegionsWriter>   regionsWriter_;
    /// \}

public:
    /// \name Which QV will be written.
    /// \{
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
    /// \}
 
private:
    /// \name Private Methods.
    /// \{
    /// \brief Checks whether chemistry triple, including
    ///        binding kit, sequencing kit and base caller version
    ///        are set. 
    ///        If not, add error messages.
    bool SanityCheckChemistry(const std::string & bindingKit,
                              const std::string & sequencingKit,
                              const std::string & basecallerVersion);

    /// \brief Closes HDFBaxWriter.
    void Close(void);
    /// \}
};

inline bool HDFBaxWriter::HasDeletionQV(void) const 
{return basecallsWriter_->HasDeletionQV();}

inline bool HDFBaxWriter::HasDeletionTag(void) const 
{return basecallsWriter_->HasDeletionTag();}

inline bool HDFBaxWriter::HasInsertionQV(void) const
{return basecallsWriter_->HasInsertionQV();}

inline bool HDFBaxWriter::HasSubstitutionTag(void) const
{return basecallsWriter_->HasSubstitutionTag();}

inline bool HDFBaxWriter::HasSubstitutionQV(void) const
{return basecallsWriter_->HasSubstitutionQV();}

inline bool HDFBaxWriter::HasMergeQV(void) const
{return basecallsWriter_->HasMergeQV();}

inline bool HDFBaxWriter::HasPreBaseFrames(void) const
{return basecallsWriter_->HasPreBaseFrames();}

inline bool HDFBaxWriter::HasIPD(void) const 
{return this->HasPreBaseFrames();}

inline bool HDFBaxWriter::HasWidthInFrames(void) const
{return basecallsWriter_->HasWidthInFrames();}

inline bool HDFBaxWriter::HasPulseWidth(void) const 
{return this->HasWidthInFrames();}
#endif

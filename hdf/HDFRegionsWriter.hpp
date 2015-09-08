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

#ifndef _HDF_REGIONS_WRITER_HPP_
#define _HDF_REGIONS_WRITER_HPP_

#include <string>
#include "Enumerations.h"
#include "reads/RegionTable.hpp"
#include "HDFFile.hpp"
#include "HDFArray.hpp"
#include "HDF2DArray.hpp"
#include "HDFAtom.hpp"
#include "HDFWriterBase.hpp"

using namespace H5;
using namespace std;

class HDFRegionsWriter: public HDFWriterBase {
public:
    /// \name Constructor and destructor
    /// \{
    /// \param[in] filename, hdf file name
    /// \param[in] parentGroup, parent hdf group in hirarchy
    HDFRegionsWriter(const std::string & filename, 
                     HDFGroup & parentGroup,
                     const std::vector<std::string> & regionTypes = PacBio::AttributeValues::Regions::regiontypes);
    ~HDFRegionsWriter(void);
    /// \}

private:
    /// \name Private variables for hdf IO. 
    /// \{
    HDFGroup & parentGroup_; //< parent hdf group

    /// A vector of strings of region types for RegionTypeIndex to look up. Order matters!
    std::vector<std::string> regionTypes_; 

	HDF2DArray<int> regionsArray_; //< HDF2DArray for writing regions to hdf

	int curRow_; //< which row to write

    static const int NCOLS = 5; //< number of columns in Regions table.

    /// \brief Write attributes of the 'regions' group
    bool WriteAttributes(void);
    /// \}
       
public:
    /// \name Method to write region annotations. 
    /// \{
    /// \brief Append a vector of region annotations to 'regions'
    /// \param[in] annotations - region annotations to append. 
    /// \returns true if succeeded.
    bool Write(const std::vector<RegionAnnotation> &annotations);

    /// \brief Append a region annotation to 'regions' 
    /// \param[in] annotation - region annotation to append
    /// \returns true if succeeded.
    bool Write(const RegionAnnotation &annotation);

    void Flush(void);

    void Close(void);
};

#endif

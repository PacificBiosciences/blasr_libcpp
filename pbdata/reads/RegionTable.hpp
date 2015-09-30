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

// Author: Mark Chaisson

#ifndef _BLASR_REGION_TABLE_HPP_
#define _BLASR_REGION_TABLE_HPP_

#include <cassert>
#include <cstring>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <ostream>
#include "Types.h"
#include "Enumerations.h"
#include "PacBioDefs.h"
#include "RegionAnnotation.hpp"
#include "RegionAnnotations.hpp"


class RegionTable {
private:
    /// RegionTable reading from h5 file 'Regions' dataset.
    /// \name member variables
    /// \{
    /// Map zmw hole number to zmw RegionAnnotations.
    std::map<UInt, RegionAnnotations> map_;
    /// \}
    
    /// \name Region table attributes.
    std::vector<std::string> columnNames;
    std::vector<std::string> regionTypes;
    std::vector<std::string> regionDescriptions;
    std::vector<std::string> regionSources;
    std::vector<RegionType>  regionTypeEnums;
    /// \}

public:
    /// \name Constructor & destructor & reset
    /// \{
    RegionTable() {}

    ~RegionTable() {}

    /// Clears member variables in region table.
    /// \returns *this
    RegionTable& Reset();
    /// \}
    
    // Different region tables have different ways of encoding regions.
    // This maps from the way they are encoded in the rgn table to a
    // standard encoding.
    //
    /// \name Accessor functions to region table attributes.
    /// \{

    /// \returns *default PacBio* region types (order matters).
    static std::vector<RegionType> DefaultRegionTypes(void);

    /// \returns RegionType enums (order matters).
    std::vector<RegionType> RegionTypeEnums(void) const;

    /// \returns RegionType strings in order
    std::vector<std::string> RegionTypes(void) const;

    /// \returns column names.
    std::vector<std::string> ColumnNames(void) const;

    /// \returns region descriptions.
    std::vector<std::string> RegionDescriptions(void) const;

    /// \returns region sources.
    std::vector<std::string> RegionSources(void) const;

    /// Construct map_ (holeNumber --> RegionAnnotations) from table.
    /// \params[in] region table containing region annotations of all zmws
    /// \params[in] ordered region type strings, which maps region types
    ///             to region type indice.
    RegionTable & ConstructTable(std::vector<RegionAnnotation> & table,
                                 const std::vector<std::string> & regionTypeStrs);

    /// Note that the ORDER of region types does matter.
    /// Set region types (order matters).
    RegionTable & RegionTypes(const std::vector<std::string> & in);

    /// Set column names, e.g., 
    /// {"HoleNumber", "TypeIndex", "Start", "End", "Score"}
    RegionTable & ColumnNames(const std::vector<std::string> & in);

    /// Set region descriptions. e.g., 
    /// {"desc of holenumber", "desc of index", "desc of start", "desc of end", "desc of score"}
    RegionTable & RegionDescriptions(const std::vector<std::string> & in);

    /// Set region sources, e.g., 
    /// {"source of holenumber", "source of index", "source of start", "source of end", "source of score"}
    RegionTable & RegionSources(const std::vector<std::string> & in);
    /// \}

    /// \name Assessor functions to zmw region annotations.
    /// \{
    /// \returns Whether or not this region table has regions of a zmw.
    bool HasHoleNumber(const UInt holeNumber) const;

    /// Get zmw region annotaions given its hole number.
    /// Note that HasHoleNumber must be called first.
    /// \returns RegionAnnotations of a zmw.
    RegionAnnotations operator [] (const UInt holeNumber) const;
    /// \}
};

#endif // _BLASR_REGION_TABLE_HPP_

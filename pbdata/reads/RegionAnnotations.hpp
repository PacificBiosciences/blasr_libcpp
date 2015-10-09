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

#ifndef _PBDATA_READS_REGION_ANNOTATIONS_HPP_
#define _PBDATA_READS_REGION_ANNOTATIONS_HPP_

#include <algorithm>
#include <string>
#include <vector>
#include "ReadInterval.hpp"
#include "RegionAnnotation.hpp"

class RegionAnnotations {

    /// \name Region Annotations of a single ZMW
    /// \{
private:
    /// \name region table of a zmw
    std::vector<RegionAnnotation> table_;

    /// \name hole number of a zmw
    UInt holeNumber_;

    /// \name region types in order.
    std::vector<RegionType> types_;

public:
    RegionAnnotations(const UInt holeNumber,
                      const std::vector<RegionAnnotation> & annotations,
                      const std::vector<RegionType> & types);

    RegionAnnotations(const RegionAnnotations & rhs);

    ~RegionAnnotations() {}

    /// \returns zmw holeNumber.
    UInt HoleNumber(void) const;

    /// \returns sorted adapters of this zmw
    std::vector<RegionAnnotation> Adapters() const;

    /// \returns whether or not has HQ region specified in table.
    bool HasHQRegion() const;

    /// \returns exactly one HQ region of this zmw.
    /// \note If no HQ region exists, return a RegionAnnotation of length  0.
    ///       If more than one HQ region is found for this zmw, raise an assertion error.
    RegionAnnotation TheHQRegion() const;

    /// \returns HQ start position of this zmw.
    DNALength HQStart() const;

    /// \returns HQ end position of this zmw.
    DNALength HQEnd() const;

    /// \returns HQ score of this zmw.
    int HQScore() const;

    /// \returns sorted insert regions of this zmw.
    std::vector<RegionAnnotation> Inserts() const;

    /// \returns a vector of all adapters
    std::vector<ReadInterval> AdapterIntervals() const;

    /// \returns a vector of all subreads
    /// \param[in] wholeLength     Length of unrolled sequence of this zmw. Note that
    ///                            this piece of info does not exist in region table.
    /// \param[in] byAdapter       false: return inserts in region table directly.
    ///                            true : infer inserts according to adapters.
    /// \param[in] byHQRegion      false: inserts may contain both HQ and LQ regions
    ///                            true : inserts in HQ regions only.
    std::vector<ReadInterval> 
    SubreadIntervals(const DNALength wholeLength,
                     const bool byAdapter = true,
                     const bool byHQRegion = true) const;

private:
    /// \returns sorted vector of region annotations of a RegionType.
    std::vector<RegionAnnotation> 
    RegionAnnotationsOfType(RegionType type) const;

    /// \returns HQ regions of this zmw.
    std::vector<RegionAnnotation> HQRegions() const;

    /// \}
};

#endif

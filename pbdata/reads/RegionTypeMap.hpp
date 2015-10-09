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

#ifndef _BLASR_REGION_TYPE_MAP_HPP_
#define _BLASR_REGION_TYPE_MAP_HPP_

#include <cassert>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include "Types.h"
#include "Enumerations.h"


class RegionTypeMap {
public:
    /// \name Map region type to/from string and index
    /// \{
    static std::string ToString(RegionType rt);

    static RegionType ToRegionType(const std::string & str);

    /// \params[in] typeStr - query region type as a string
    /// \params[in] typeStrs - a vector region type strings in order
    /// \returns index of a region type as string in a vector of region type strings
    static int ToIndex(const std::string & typeStr,
                       const std::vector<std::string> & typeStrs);

    /// \params[in] rt - query region type
    /// \params[in] typeStrs - a vector region type strings in order
    /// \returns index of the query region type in a vector of region type strings
    static int ToIndex(RegionType rt,
                       const std::vector<std::string> & typeStrs);

    /// \params[in] rt - query region type
    /// \params[in] regionTypes - a vector region type strings in order
    /// \returns index of the query region type in a vector of region type enums 
    static int ToIndex(RegionType rt,
                       const std::vector<RegionType> & regionTypes);
private:
    // Map region type to string
    static const std::map<RegionType, std::string> RegionTypeToString;

    // Map string to region type 
    static const std::map<std::string, RegionType> StringToRegionType;
    /// \}
};

#endif // _BLASR_REGION_TYPE_MAP_HPP_

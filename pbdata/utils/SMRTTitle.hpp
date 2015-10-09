#ifndef _BLASR_SMRT_TITLE_HPP_
#define _BLASR_SMRT_TITLE_HPP_
#include <string>
#include <vector>
#include <sstream>
#include "Types.h"
#include "StringUtils.hpp"

class SMRTTitle {
public:
    std::string movieName;
    UInt holeNumber;
    DNALength start;
    DNALength end;
    // if input name is a smrt title.
    bool isSMRTTitle;

    SMRTTitle(const std::string & name);

    /// \returns smrt title movie/zmw/s_e, if input read is a smrt title;
    /// otherwise, return an empty string.
    std::string ToString();

public:
    inline std::string MovieName(void) const;
    inline UInt HoleNumber(void) const;
    inline DNALength Start(void) const;
    inline DNALength End(void) const;
    inline operator bool(void) const;
};

inline std::string SMRTTitle::MovieName(void) const {return movieName;}
inline UInt SMRTTitle::HoleNumber(void) const {return holeNumber;}
inline DNALength SMRTTitle::Start(void) const {return start;}
inline DNALength SMRTTitle::End(void) const {return end;}
inline SMRTTitle::operator bool(void) const {return isSMRTTitle;}

#endif

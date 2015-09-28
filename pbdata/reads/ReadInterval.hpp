#ifndef _BLASR_READ_INTERVAL_HPP_
#define _BLASR_READ_INTERVAL_HPP_

#include "RegionAnnotation.hpp"

class RegionAnnotation;

class ReadInterval {
public:
    int start;
    int end;
    int score;

    ReadInterval(int s, int e, int sc=0) : start(s), end(e), score(sc) {};

    ReadInterval& operator=(const ReadInterval &rhs) {
        start = rhs.start;
        end   = rhs.end;
        score = rhs.score;
        return *this;
    }

    ReadInterval(const RegionAnnotation & ra)
    : start(ra.GetStart())
    , end(ra.GetEnd())
    , score(ra.GetScore()) {}
};

#endif

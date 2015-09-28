#ifndef _BLASR_READ_INTERVAL_HPP_
#define _BLASR_READ_INTERVAL_HPP_

#include "RegionAnnotation.hpp"

class RegionAnnotation;

class ReadInterval {
public:
    int start;
    int end;
    int score;

    ReadInterval(int s=0, int e=0, int sc=0) : start(s), end(e), score(sc) {};

    ReadInterval(const RegionAnnotation & ra)
    : start(ra.GetStart())
    , end(ra.GetEnd())
    , score(ra.GetScore()) {}

    ReadInterval& operator=(const ReadInterval &rhs) {
        start = rhs.start;
        end   = rhs.end;
        score = rhs.score;
        return *this;
    }

    bool operator==(const ReadInterval &rhs) const {
        return (start == rhs.start and
                end   == rhs.end   and
                score == rhs.score);
    }
};

#endif

#ifndef _BLASR_WEIGHTED_INTERVAL_HPP_
#define _BLASR_WEIGHTED_INTERVAL_HPP_

#include <vector>
#include <set>
#include "DNASequence.hpp"
#include "datastructures/anchoring/MatchPos.hpp" 

class WeightedInterval {
public:
    DNALength size; // not necessarily end - start + 1
    DNALength start;
    DNALength end;
    DNALength qStart, qEnd;
    int readIndex;
    float pValue;
    std::vector<int> positions;
    std::vector<ChainedMatchPos> matches;
    float pValueVariance, pValueNStdDev;
    float sizeVariance, sizeNStdDev;
    int nAnchors;
    int totalAnchorSize;
    bool isOverlapping;

    WeightedInterval();

    WeightedInterval(int _size, int _start, int _end, 
        int _readIndex, float _pValue =0.0);

    WeightedInterval(int _size, int _start, int _end, 
        int _readIndex, float _pValue, int _qStart, int _qEnd);

    WeightedInterval(int _size, unsigned int _nAnchors, 
        unsigned int _totalAnchorSize, int _start, int _end, 
        int _readIndex, float _pValue, int _qStart, int _qEnd, 
        std::vector<ChainedMatchPos> &_matches);

    void Init(int _size, int _start, int _end, 
        int _readIndex, float _pValue);

    int GetStrandIndex() const;
    void SetPValueVariance(float v);
    void SetPValueNStdDev(float v);
    void SetSizeVariance(float v);
    void SetSizeNStdDev(float v);
    int operator<(const WeightedInterval &intv) const;
    int operator==(const WeightedInterval &intv) const;

    friend std::ostream & operator << (std::ostream & out, WeightedInterval & wi) {
        out << wi.size << " " << wi.start << " " << wi.end
            << " " << wi.qStart << " " << wi.qEnd << " " << wi.readIndex
            << " " << wi.pValue ;
        return out;
    }

    float PValue() const;
    int Size() const;
};

class CompareWeightedIntervalByPValue {
public:
    int operator()(const WeightedInterval& a, const WeightedInterval& b) const;
};

typedef std::vector<WeightedInterval> WeightedIntervalVector;

typedef std::multiset<WeightedInterval, CompareWeightedIntervalByPValue> 
T_WeightedIntervalMultiSet;

class WeightedIntervalSet: public T_WeightedIntervalMultiSet {
public:
    size_t maxSize;
    WeightedIntervalSet();
    WeightedIntervalSet(const size_t maxSizeP);
    bool insert(WeightedInterval &intv);
    friend std::ostream & operator << (std::ostream & out, WeightedIntervalSet & wis) {
        WeightedIntervalSet::iterator it;
        for (it = wis.begin(); it != wis.end(); it++) {
            out << *((WeightedInterval*)&(*it)) << std::endl;
        }
        return out;
    }
};

#endif

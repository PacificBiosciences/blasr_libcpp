#include "datastructures/anchoring/WeightedInterval.hpp"

WeightedInterval::WeightedInterval(){}

void WeightedInterval::Init(int _size, int _start, int _end, 
    int _readIndex, float _pValue) {
    size      = _size; 
    start     = _start; 
    end       = _end; 
    readIndex = _readIndex;
    pValue    = _pValue;
    qStart = 0;
    qEnd   = 0;
    nAnchors = 0;
    totalAnchorSize = 0;
    pValueVariance = 0;
    pValueNStdDev = 0;
    sizeVariance = 0;
    sizeNStdDev = 0;
}

WeightedInterval::WeightedInterval(int _size, int _start, int _end, 
    int _readIndex, float _pValue) {
    Init(_size, _start, _end, _readIndex, _pValue);
}

WeightedInterval::WeightedInterval(int _size, int _start, int _end, 
    int _readIndex, float _pValue, int _qStart, int _qEnd){
    Init(_size, _start, _end, _readIndex, _pValue);
    qStart    = _qStart;
    qEnd      = _qEnd;
}

WeightedInterval::WeightedInterval(int _size, unsigned int _nAnchors, 
    unsigned int _totalAnchorSize, int _start, int _end, int _readIndex, 
    float _pValue, int _qStart, int _qEnd, 
    std::vector<ChainedMatchPos> &_matches) {

    Init(_size, _start, _end, _readIndex, _pValue);
    qStart    = _qStart;
    qEnd      = _qEnd;
    matches   = _matches;
    nAnchors  = _nAnchors;
    totalAnchorSize = _totalAnchorSize;
}

float WeightedInterval::PValue() const {
    return pValue;
}

int WeightedInterval::Size() const {
    return size;
}

int WeightedInterval::GetStrandIndex() const {
    return readIndex;
}

void WeightedInterval::SetPValueVariance(float v) {
    pValueVariance = v;
}

void WeightedInterval::SetPValueNStdDev(float v) {
    pValueNStdDev = v;
}

void WeightedInterval::SetSizeVariance(float v) {
    sizeVariance = v;
}

void WeightedInterval::SetSizeNStdDev(float v) {
    sizeNStdDev = v;
}

int WeightedInterval::operator<(const WeightedInterval &intv) const {
    if (size == intv.size) {
        return start > intv.start;
    }
    else {
        return size < intv.size;
    }
}

int WeightedInterval::operator==(const WeightedInterval &intv) const {
    return size == intv.size;
}

// Functions of class CompareWeightedIntervalByPValue
int CompareWeightedIntervalByPValue::operator()(const WeightedInterval &a,
    const WeightedInterval &b) const {
    if (a.PValue() != b.PValue()) {
        return a.PValue() < b.PValue();
    }
    else {
        return a.start < b.start;
    }
}

// Functions of class WeightedIntervalSet
WeightedIntervalSet::WeightedIntervalSet()
    : maxSize (0) 
{ }

WeightedIntervalSet::WeightedIntervalSet(const size_t maxSizeP)
    : maxSize(maxSizeP) 
{ } 

bool WeightedIntervalSet::insert(WeightedInterval &intv) {
    intv.isOverlapping = false;

    //
    // Make sure this interval is not contained inside any other
    // weighted intervals.  
    //
    WeightedIntervalSet::iterator it = (*this).begin();
    WeightedIntervalSet::iterator endit = (*this).end();
    bool isContained = false;
    while (it != endit and isContained == false) {
        if (intv.qStart >= (*it).qStart and intv.qEnd <= (*it).qEnd and
            intv.start  >= (*it).start and intv.end <= (*it).end and 
            intv.readIndex == (*it).readIndex and
            intv.pValue >= (*it).pValue) {
            //
            // This already overlaps an existing interval, don't bother
            // trying to add it.
            //
            isContained = true;
            intv.isOverlapping = true;
        } else if((*it).start >= intv.start and 
                  (*it).end <= intv.end and
                  (*it).qStart >= intv.qStart and 
                  (*it).qEnd <= intv.qEnd and
                  (*it).readIndex == intv.readIndex and
                  (*it).pValue >= intv.pValue) {
            WeightedIntervalSet::iterator next = it;
            ++next;
            this->erase(it);
            it = next;
        } else {
            ++it;
        }
    }

    //
    // Take a peek to see if this interval is too low of a score to
    // bother attempting to add at all. 
    //
    if (size() >= maxSize and maxSize > 0) {
        WeightedIntervalSet::iterator last = (*this).end();
        last--;

        if (last->pValue < intv.pValue) {
            return false;
        }
    }

    if (isContained == false) {
        bool addInsert = false;
        if (size() == 0) {
            addInsert = true;
        }
        else {
            it = end();
            --it;
            if (size() < maxSize or (*it).pValue > intv.pValue) {
                addInsert = true;
                //
                // Keep the size of the stack the same if it is at the limit.
                //
                if (maxSize != 0 and size() >= maxSize and size() > 0) { 
                    erase(it);
                }
            }
        }
        if (addInsert) {
            ((T_WeightedIntervalMultiSet*)this)->insert(intv);
        }
        return true;
    }
    return false;
}


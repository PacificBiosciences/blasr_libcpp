#ifndef _BLASR_BASE_FILE_HPP_
#define _BLASR_BASE_FILE_HPP_

#include <stdint.h>
#include <vector>
#include "Enumerations.h"
#include "SMRTSequence.hpp"

#include "HoleXY.hpp"
#include "PulseBaseCommon.hpp"
#include "ScanData.hpp"

class CompareHoleXY {
public:
    bool operator()(const HoleXY &lhs, const HoleXY &rhs) const; 
};

class BaseFile : public PulseBaseCommon {
public:
    std::vector<unsigned char> baseCalls;
    std::vector<uint8_t> holeStatus;
    std::vector<HoleXY>   holeXY;
    std::vector<uint16_t> basWidthInFrames;
    std::vector<uint16_t> preBaseFrames;
    std::vector<int> pulseIndex;
    std::vector<unsigned char> qualityValues;
    std::vector<unsigned char> deletionQV;
    std::vector<unsigned char> deletionTag;
    std::vector<unsigned char> insertionQV;
    std::vector<unsigned char> substitutionQV;
    std::vector<unsigned char> substitutionTag;
    std::vector<unsigned char> mergeQV;
    std::vector<uint32_t>  readLengths;
    std::vector<int>  readStartPositions;
    int nReads;
    int nBases;

    bool LookupReadIndexByXY(uint16_t x, uint16_t y, int &index); 

    void CopyReadAt(uint32_t readIndex, SMRTSequence &read); 

    template<typename T>
    void CopyArray(std::vector<T> &fullArray, int pos, int length, T*dest); 
    
};

#include "BaseFileImpl.hpp"

#endif // _BLASR_BASE_FILE_HPP_

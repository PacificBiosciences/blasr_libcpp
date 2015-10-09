#include <algorithm>
#include "BaseFile.hpp"

bool CompareHoleXY::operator()(const HoleXY &lhs, const HoleXY &rhs) const {
    if (lhs.xy[0] == rhs.xy[0]) {
        return lhs.xy[1] < rhs.xy[1];
    }
    else {
        return lhs.xy[0] < rhs.xy[0];
    }
}

bool BaseFile::LookupReadIndexByXY(uint16_t x, uint16_t y, int &index) {
    int16_t xy[2];
    xy[0] = x; xy[1] = y;
    std::vector<HoleXY>::iterator holeIt;
    holeIt = lower_bound(holeXY.begin(), holeXY.end(), xy);
    if ((*holeIt).xy[0] == xy[0] and (*holeIt).xy[1] == xy[1]) {
        index = holeIt - holeXY.begin();
        return true;
    }
    else {
        return false;
    }
}


void BaseFile::CopyReadAt(uint32_t readIndex, SMRTSequence &read) {
    assert(holeNumbers.size() > readIndex);
    read.HoleNumber(holeNumbers[readIndex]);
    if (holeXY.size() > 0) {
        assert(holeXY.size() > readIndex);
        read.HoleXY(holeXY[readIndex].xy[0],
                    holeXY[readIndex].xy[1]);
    }

    int startPos = readStartPositions[readIndex];
    int readLength = readLengths[readIndex];
    read.length = readLength;
    read.Allocate(readLength);
    if (baseCalls.size() > 0) {
        assert(baseCalls.size() >= readLength + startPos);
        CopyArray(baseCalls, startPos, readLength, read.seq);
    }
    if (qualityValues.size() > 0) {
        assert(qualityValues.size() >= readLength + startPos);
        CopyArray(qualityValues, startPos, readLength, read.qual.data);
    }
    if (basWidthInFrames.size() > 0) {
        assert(basWidthInFrames.size() >= readLength + startPos);
        CopyArray(basWidthInFrames, startPos, readLength, read.widthInFrames);
    }
    if (deletionQV.size() > 0) {
        assert(deletionQV.size() >= readLength + startPos);
        CopyArray(deletionQV, startPos, readLength, read.deletionQV.data);
    }
    if (deletionTag.size() > 0) {
        assert(deletionTag.size() >= readLength + startPos);
        CopyArray(deletionTag, startPos, readLength, read.deletionTag);
    }
    if (insertionQV.size() > 0) {
        assert(insertionQV.size()  >= readLength + startPos);
        CopyArray(insertionQV, startPos, readLength, read.insertionQV.data);
    }
    if (substitutionQV.size() > 0) {
        assert(substitutionQV.size() >= readLength + startPos);
        CopyArray(substitutionQV, startPos, readLength, read.substitutionQV.data);
    }
    if (mergeQV.size() > 0) {
        assert(mergeQV.size() >= readLength + startPos);
        CopyArray(mergeQV, startPos, readLength, read.mergeQV.data);
    }
    if (substitutionTag.size() > 0) {
        assert(substitutionTag.size() >= readLength + startPos);
        CopyArray(substitutionTag, startPos, readLength, read.substitutionTag);
    }
    if (preBaseFrames.size() > 0) {
        assert(preBaseFrames.size() >= readLength + startPos);
        CopyArray(preBaseFrames, startPos, readLength, read.preBaseFrames);
    }

}

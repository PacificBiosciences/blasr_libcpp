#ifndef _BLASR_REGION_UTILS_IMPL_HPP 
#define _BLASR_REGION_UTILS_IMPL_HPP 


//FIXME: move all functions to class SMRTSequence
template<typename T_Sequence>
bool MaskRead(T_Sequence &fastaRead,
              ZMWGroupEntry &zmwData,
              RegionTable &regionTable)
{
    if (not regionTable.HasHoleNumber(zmwData.holeNumber)) {
        return false;
    } else {
        RegionAnnotations regions = regionTable[zmwData.holeNumber];

        // Mask off the low quality portion of this read.
        DNALength readPos;
        for (readPos = 0; readPos < std::min(regions.HQStart(), fastaRead.length); readPos++) {
            fastaRead.seq[readPos] = 'N';
        }
        for (readPos = regions.HQEnd(); readPos < fastaRead.length; readPos++) {
            fastaRead.seq[readPos] = 'N';
        }
        return regions.HasHQRegion();
    }
}

/// \params[in]  - fastaRead, zmwData, regionTable
/// \params[out] - readStart
/// \params[out] - readEnd
/// \params[out] - score
/// \returns Whether or not read coordinate trimmed according to HQRegion
template<typename T_Sequence>
bool GetReadTrimCoordinates(T_Sequence &fastaRead,
	ZMWGroupEntry &zmwData,	RegionTable &regionTable,
	DNALength &readStart ,DNALength &readEnd, int &score) {

    if (regionTable.HasHoleNumber(zmwData.holeNumber)) {
        RegionAnnotations regions = regionTable[zmwData.holeNumber];
        if (regions.HasHQRegion()) {
            readStart = regions.HQStart();
            readEnd   = regions.HQEnd();
            return true;
        }
    }

    readStart = 0;
    readEnd = fastaRead.length;
    return false;
}

#endif

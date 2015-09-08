#ifndef _BLASR_REGION_UTILS_IMPL_HPP 
#define _BLASR_REGION_UTILS_IMPL_HPP 


template<typename T_Sequence>
bool MaskRead(T_Sequence &fastaRead, ZMWGroupEntry &zmwData,
    RegionTable &regionTable) {
    int regionIndex;						 
    int regionLowIndex, regionHighIndex;
    regionLowIndex = regionHighIndex = 0;

    regionTable.LookupRegionsByHoleNumber(zmwData.holeNumber, 
        regionLowIndex, regionHighIndex);

    bool readHasGoodRegion = true;

    DNALength readPos;

    regionIndex = regionLowIndex;
    int lastHQRegionIndex;

    int hqRegionStart=0, hqRegionEnd=0, hqRegionScore = 0;
    readHasGoodRegion = LookupHQRegion(zmwData.holeNumber, 
        regionTable, hqRegionStart, hqRegionEnd, hqRegionScore);

    //
    // Mask off the low quality portion of this read.
    //
    for (readPos = 0; (readPos < hqRegionStart and
         readPos < fastaRead.length); readPos++) {
        fastaRead.seq[readPos] = 'N';
    }

    for (readPos = hqRegionEnd; readPos < fastaRead.length; readPos++) {
        fastaRead.seq[readPos] = 'N';
    }

    //
    // Look to see if there is region information provided, but the
    // entire read is bad.
    //
    if (hqRegionEnd == hqRegionStart) {
        //
        // This read is entirely bad, flag that.
        //
        readHasGoodRegion = false;
    }

    return readHasGoodRegion;
}


template<typename T_Sequence>
bool GetReadTrimCoordinates(T_Sequence &fastaRead,
	ZMWGroupEntry &zmwData,	RegionTable &regionTable,
	DNALength &readStart ,DNALength &readEnd, int &score) {

	int regionIndex;						 
	int regionLowIndex, regionHighIndex;
	regionLowIndex = regionHighIndex = 0;

	regionTable.LookupRegionsByHoleNumber(zmwData.holeNumber,
        regionLowIndex, regionHighIndex);

	bool readHasGoodRegion = true;

	DNALength readPos;

	regionIndex = regionLowIndex;
	int lastHQRegionIndex;
	
	while (regionIndex < regionHighIndex and 
		   regionTable.GetType(regionIndex) != HQRegion) {
	    regionIndex++;
	}
	
	if (regionIndex < regionHighIndex ) {
		readStart = regionTable.GetStart(regionIndex);
		readEnd   = regionTable.GetEnd(regionIndex);
        score     = regionTable.GetScore(regionIndex);
		return true;
	}
	else {
		readStart = 0;
		readEnd   = fastaRead.length;
		return false;
	}
}


template<typename T_Sequence>
bool TrimRead(T_Sequence &fastaRead, ZMWGroupEntry &zmwData,
    RegionTable &regionTable, T_Sequence &trimmedRead) {

	DNALength readStart, readEnd;

	GetReadTrimCoordinates(fastaRead, zmwData, regionTable,
        readStart, readEnd);

	if (readEnd - readStart > 0) {
		trimmedRead.CopySubsequence((FASTQSequence&)fastaRead, 
			readStart, readEnd);
		// signal that the read has a good region.
		return true;
	}
	else {
		//
		// There is no information for this read. Make it skipped.
		//
		trimmedRead.seq = NULL;
		trimmedRead.CopyTitle(fastaRead.title);
		// signal this read has no good region.
		return false;
	}
}

template<typename T_Sequence>
void CollectSubreadIntervals(T_Sequence &read, RegionTable *regionTablePtr, 
    std::vector<ReadInterval> &subreadIntervals, bool byAdapter) {

	int regionIndex;						 
	int regionLowIndex, regionHighIndex;
	regionLowIndex = regionHighIndex = 0;

	regionTablePtr->LookupRegionsByHoleNumber(read.zmwData.holeNumber,
        regionLowIndex, regionHighIndex);

	if (byAdapter == false) { 
        // read subreads (insert) directly from region table. 
		for (regionIndex = regionLowIndex; 
             regionIndex < regionHighIndex; regionIndex++) {
			if (regionTablePtr->GetType(regionIndex) ==  Insert) {

                RegionAnnotation & ra = regionTablePtr->table[regionIndex];
				subreadIntervals.push_back(ReadInterval(
                            ra.GetStart(),
                            ra.GetEnd(),
                            ra.GetScore()));
			}
		}
	}
	else { // Determine subreads according to adapters only.
        std::vector<int> adapterIntervalIndices;
		for (regionIndex = regionLowIndex; 
             regionIndex < regionHighIndex; regionIndex++) {
			if (regionTablePtr->GetType(regionIndex) == Adapter) {
				adapterIntervalIndices.push_back(regionIndex);
			}
		}
		// Sort indices so that the intervals appear in order on the read.
		SortRegionIndicesByStart(*regionTablePtr, adapterIntervalIndices);
		int curIntervalStart = 0;
		int i;
		if (adapterIntervalIndices.size() == 0) {
            // no adapter, this zmw has only one subread (pass)
			subreadIntervals.push_back(ReadInterval(0, read.length));
		}
		else {
            // The first subread covers [0, RegionStart of first adapter)
			subreadIntervals.push_back(ReadInterval(0, 
                regionTablePtr->table[adapterIntervalIndices[0]].
                    GetStart()));

            // The subread[i] covers (RegionEnd of i-1-th adapter, RegionStart of i-th adapter)
			for (i = 0; i + 1 < adapterIntervalIndices.size(); i++) {
				subreadIntervals.push_back(ReadInterval(
                    regionTablePtr->table[adapterIntervalIndices[i]].GetEnd(),
				    regionTablePtr->table[adapterIntervalIndices[i+1]].GetStart()));
			}
            // The last subread covers (RegionEnd of last adapter, end of read)
			subreadIntervals.push_back(
                ReadInterval(regionTablePtr->table[
                    adapterIntervalIndices[adapterIntervalIndices.size()-1]].GetEnd(),
				read.length));
		}
	}
	sort(subreadIntervals.begin(), subreadIntervals.end(), 
        OrderRegionsByReadStart());
}


// Get all adapter intervals of a ZMW.
// Input:
//   read - read.zmwData.holeNumber specifies the zmw.
//   regionTablePtr - a pointer to a region table.
// Output:
//   adapterIntervals - where to assign all adapter intervals of the zmw
template<typename T_Sequence>
void CollectAdapterIntervals(T_Sequence &read, 
    RegionTable *regionTablePtr, 
    std::vector<ReadInterval> &adapterIntervals) {

    assert(regionTablePtr != NULL);
    int regionLowIndex = 0, regionHighIndex = 0;

    regionTablePtr->LookupRegionsByHoleNumber(read.zmwData.holeNumber,
        regionLowIndex, regionHighIndex);

    for (int regionIndex = regionLowIndex; 
         regionIndex < regionHighIndex; regionIndex++) {

        if (regionTablePtr->GetType(regionIndex) ==  Adapter) {
            RegionAnnotation & ra = regionTablePtr->table[regionIndex];
            adapterIntervals.push_back(ReadInterval(
                        ra.GetStart(),
                        ra.GetEnd(), 
                        ra.GetScore()));
        }
    }
}

#endif

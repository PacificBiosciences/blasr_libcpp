#include "files/FragmentCCSIterator.hpp"


void FragmentCCSIterator::
Initialize(CCSSequence *_seqPtr, RegionTable *_regionTablePtr) {
    seqPtr         = _seqPtr;
    regionTablePtr = _regionTablePtr;
    curPass = 0;
    numPasses = 0;
    subreadIntervals.clear();
    readIntervalDirection.clear();

    int hqRegionStart, hqRegionEnd, hqRegionScore;
    hqRegionStart = hqRegionEnd = hqRegionScore = 0;

    bool hasHQRegion = LookupHQRegion(seqPtr->zmwData.holeNumber, 
        *regionTablePtr, hqRegionStart, hqRegionEnd, hqRegionScore);

    if (not hasHQRegion) {
        return; // Don't bother if there is no HQ region.
    }

    //
    // Since this iterator covers all passes, and not just those
    // included in the ccs, the the regions need to be loaded.
    subreadIntervals = (*regionTablePtr)[seqPtr->HoleNumber()].SubreadIntervals(seqPtr->unrolledRead.length, true);
    if (subreadIntervals.size() == 0) { return;}

    readIntervalDirection.resize(subreadIntervals.size());
    fill(readIntervalDirection.begin(), readIntervalDirection.end(), 2);

    //
    // Assign the read interval directions based on the pass direction
    // for the pass that has a similar start position.  This allows
    // some wiggle although in practice they coordinates of the pass
    // start base and the template should always match up. 
    //
    int i, j;
    for (i = 0; i < subreadIntervals.size(); i++) {
        for (j = 0; j < seqPtr->passStartBase.size(); j++) {
            if (abs( ((int)subreadIntervals[i].start)  - 
                     ((int)seqPtr->passStartBase[j]) ) < 10) {
                readIntervalDirection[i] = seqPtr->passDirection[j];
                break;
            }
        }
    }

    int firstAssignedSubread = 0;
    while (firstAssignedSubread < subreadIntervals.size() and 
           readIntervalDirection[firstAssignedSubread] == 2) { 
        firstAssignedSubread++; 
    }

    if (firstAssignedSubread == subreadIntervals.size()) {
        // None of the subread has been assigned a direction, guess.
        firstAssignedSubread = 0;
        readIntervalDirection[0] = 0;
    }

    // Assign directions to intervals to the left of the first assigned.
    if (firstAssignedSubread < subreadIntervals.size() and 
        subreadIntervals.size() > 0) {
        int curSubreadDir = readIntervalDirection[firstAssignedSubread];
        assert(curSubreadDir == 0 or curSubreadDir == 1);
        for (i = firstAssignedSubread - 1; i >= 0; i--) {
            curSubreadDir = (curSubreadDir==0)?1:0;
            readIntervalDirection[i] = curSubreadDir;
        }
    }

    // Assign directions to intervals which are to the right of the first 
    // assigned and whose direction is unknown.
    for (i = firstAssignedSubread + 1; i < subreadIntervals.size(); i++) {
        int & di = readIntervalDirection[i];
        int   dp = readIntervalDirection[i-1]; 
        if (di != 0 and di != 1) {
            di = (dp==0)?1:0; 
        }
    }

    //
    // So far, subreadIntervals have been sorted and each assigned 
    // a passDirection. But since all or part of a subreadInterval 
    // may not be in the HQ region, we need to trim low quality regions 
    // from subreads, remove subreads which do not have any high quality 
    // regions from subreadIntervals and their corresponding pass directions 
    // from readIntervalDirection. 
    //
    GetHighQualitySubreadsIntervals(subreadIntervals, 
            readIntervalDirection,
            hqRegionStart, hqRegionEnd);
    // Update number of passes. 
    numPasses = subreadIntervals.size();
}

int FragmentCCSIterator::
GetNext(int &direction, int &startBase, int &numBases) {
    if (curPass >= subreadIntervals.size()) {
        return 0;
    }
    direction = int(readIntervalDirection[curPass]);
    startBase = int(subreadIntervals[curPass].start);
    numBases  = int(subreadIntervals[curPass].end - subreadIntervals[curPass].start);
    ++curPass;
    return 1;
}


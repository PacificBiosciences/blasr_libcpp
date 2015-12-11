#include "Types.h"
#include "HDFZMWReader.hpp"

HDFZMWReader::HDFZMWReader() {
    closeFileOnExit = false;
    readHoleNumber  = false;
    readHoleXY      = false;
    readNumEvent    = false;
    readHoleStatus  = false;
    nZMWEntries = curZMW = 0;
    parentGroupPtr = NULL;
}

int HDFZMWReader::Initialize(HDFGroup *parentGroupP) {
    parentGroupPtr = parentGroupP;
    closeFileOnExit = false;
    return Initialize();
}

int HDFZMWReader::Initialize() {

    //
    // Make sure we can open the component containing the zmw information.
    //
    if (parentGroupPtr->ContainsObject("ZMW") == 0 or
            zmwGroup.Initialize(parentGroupPtr->group, "ZMW") == 0) {
        return 0;
    }

    //
    // Now open all the important datasets in the zmw group.  Some of
    // these are optional, so flags must be set if they do not exist.
    //
    if (zmwGroup.ContainsObject("HoleNumber")) {
        if (holeNumberArray.Initialize(zmwGroup, "HoleNumber") == 0) {
            return 0;
        }
        readHoleNumber = true;
    }
    else {
        readHoleNumber = false;
    }

    if (zmwGroup.ContainsObject("HoleStatus")) {
        if (holeStatusArray.Initialize(zmwGroup, "HoleStatus") == 0) {
            return 0;
        }
        readHoleStatus = true;
    }
    else {
        readHoleStatus = false;
    }


    if (zmwGroup.ContainsObject("HoleXY")) {
        if (xyArray.Initialize(zmwGroup, "HoleXY") == 0) {
            return 0;
        }
        readHoleXY = true;
    }
    else {
        readHoleXY = false;
    }
    if (numEventArray.Initialize(zmwGroup, "NumEvent") == 0) {
        return 0;
    }
    nZMWEntries = numEventArray.arrayLength;
    readNumEvent = true;
    curZMW      = 0;
    return 1;
}

int HDFZMWReader::Advance(UInt nSteps) {
    if (curZMW >= nZMWEntries) {
        return 0;
    }
    else {
        curZMW += nSteps;
        return 1;
    }
}

bool HDFZMWReader::GetNext(ZMWGroupEntry &groupEntry) {
    if (curZMW == nZMWEntries) {
        return false;
    }
    if (readHoleNumber) {
        holeNumberArray.Read(curZMW, curZMW+1, &groupEntry.holeNumber);
    }
    if (readHoleStatus) {
        holeStatusArray.Read(curZMW, curZMW+1, &groupEntry.holeStatus);
    }
    if (readHoleXY){ 
        int16_t holeXY[2];
        xyArray.Read(curZMW, curZMW+1, holeXY);
        groupEntry.x = holeXY[0];
        groupEntry.y = holeXY[1];
    }
    numEventArray.Read(curZMW, curZMW+1, &groupEntry.numEvents);
    curZMW++;
    return true;
}

void HDFZMWReader::Close() {
    if (readHoleNumber) {
        holeNumberArray.Close();
    }
    if (readHoleStatus) {
        holeStatusArray.Close();
    }
    if (readHoleXY) {
        xyArray.Close();
    }
    if (readNumEvent) {
        numEventArray.Close();
    }

    if (closeFileOnExit == true) {
        //
        // This instance is owner of it's reader.  Close the reader file.
        //
        hdfPlsFile.close();
    }
    zmwGroup.Close();

}

bool HDFZMWReader::GetHoleNumberAt(UInt index, UInt &holeNumber) {
    if (index >= nZMWEntries) { return false; }
    holeNumberArray.Read(index, index + 1, &holeNumber); 
    return true;
}

HDFZMWReader::~HDFZMWReader() {
    Close();
}

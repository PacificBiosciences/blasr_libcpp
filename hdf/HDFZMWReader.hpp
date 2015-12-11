#ifndef _BLASR_HDF_ZMW_READER_HPP_
#define _BLASR_HDF_ZMW_READER_HPP_

#include <cstdint>
#include "H5Cpp.h"
#include "reads/ZMWGroupEntry.hpp"
#include "HDFArray.hpp"
#include "HDF2DArray.hpp"
#include "HDFGroup.hpp"

class HDFZMWReader {
public:
    HDFGroup *parentGroupPtr;
    HDFGroup zmwGroup;
    HDFArray<UInt> holeNumberArray;
    HDFArray<unsigned char> holeStatusArray;
    HDF2DArray<int16_t> xyArray;
    HDFArray<DNALength> numEventArray;
    bool readHoleNumber, readHoleStatus;
    bool readHoleXY;
    bool readNumEvent;
    UInt curZMW;
    UInt nZMWEntries;
    bool  closeFileOnExit;
    H5::H5File hdfPlsFile;

    HDFZMWReader(); 

    int Initialize(HDFGroup *parentGroupP); 

    int Initialize(); 

    int Advance(UInt nSteps); 

    bool GetNext(ZMWGroupEntry &groupEntry); 

    void Close(); 
    
    // Return true if get hole number at ZMW/HoleNumber[index].
    bool GetHoleNumberAt(UInt index, UInt & holeNumber);

    ~HDFZMWReader(); 
};

#endif

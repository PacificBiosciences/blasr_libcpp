#ifndef _BLASR_HDF_PULSE_DATA_FILE_HPP_
#define _BLASR_HDF_PULSE_DATA_FILE_HPP_

#include <string>
#include <vector>
#include <H5Cpp.h>
#include "HDFGroup.hpp"
#include "HDFZMWReader.hpp"
#include "HDFScanDataReader.hpp"

class HDFPulseDataFile {
public:
    H5::H5File hdfBasFile;
    HDFGroup pulseDataGroup;
    HDFGroup rootGroup;
    HDFGroup *rootGroupPtr;
    string pulseDataGroupName;
    HDFScanDataReader scanDataReader;
    bool useScanData;
    bool closeFileOnExit;

    // Max number of rows allowed to read into memory at a time.
    DSLength  maxAllocNElements;

    HDFZMWReader zmwReader;
    std::vector<DSLength> eventOffset;
    UInt nReads;
    bool preparedForRandomAccess;

    DSLength GetAllReadLengths(std::vector<DNALength> &readLengths);

    void CheckMemoryAllocation(long allocSize, long allocLimit,
        const char *fieldName = NULL); 

    HDFPulseDataFile();

    void PrepareForRandomAccess();

    int OpenHDFFile(std::string fileName,
        const H5::FileAccPropList & fileAccPropList=H5::FileAccPropList::DEFAULT);

    //
    // All pulse data files contain the "PulseData" group name.
    // 
    //
    int InitializePulseDataFile(std::string fileName,
        const H5::FileAccPropList & fileAccPropList=H5::FileAccPropList::DEFAULT);

    int Initialize(std::string fileName, 
        const H5::FileAccPropList & fileAccPropList=H5::FileAccPropList::DEFAULT);

    //
    // Initialize inside another open group.
    //
    int Initialize(HDFGroup *rootGroupP);

    //
    // Initialize all fields
    //
    int Initialize();

    int InitializePulseGroup();

    size_t GetAllHoleNumbers(std::vector<unsigned int> &holeNumbers);

    void Close();

};

#endif

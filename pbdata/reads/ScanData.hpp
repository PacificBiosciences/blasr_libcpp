#ifndef DATASTRUCTURES_READS_SCAN_DATA_H_
#define DATASTRUCTURES_READS_SCAN_DATA_H_

#include <string>
#include <map>
#include "Enumerations.h"

class HDFScanDataReader;
class HDFScanDataWriter;

class ScanData {
    friend class HDFScanDataReader;
    friend class HDFScanDataWriter;

public:
    PlatformId platformId;
    float frameRate;
    unsigned int numFrames;
    std::string movieName, runCode;
    std::string whenStarted;
    std::map<char, int> baseMap;
    ScanData();
    std::string GetMovieName(); 

    ScanData & PlatformID(const PlatformId & id);
    ScanData & FrameRate(const float & rate);
    ScanData & NumFrames(const unsigned int & num);
    ScanData & MovieName(const std::string & name);
    ScanData & RunCode(const std::string & code);
    ScanData & WhenStarted(const std::string & when);
    ScanData & BaseMap(const std::map<char, int> bmp);
    ScanData & SequencingKit(const std::string sequencingKit);
    ScanData & BindingKit(const std::string bindingKit);

    PlatformId PlatformID(void) const;
    float FrameRate(void) const;
    unsigned int NumFrames(void) const;
    std::string MovieName(void) const;
    std::string RunCode(void) const;
    std::string WhenStarted(void) const;
    std::map<char, int> BaseMap(void) const;

    std::string SequencingKit(void) const;
    std::string BindingKit(void) const;

private:
    std::string sequencingKit_;
    std::string bindingKit_;
};

#endif

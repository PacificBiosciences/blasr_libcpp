#ifndef DATA_HDF_HDF_SCAN_DATA_WRITER_H_
#define DATA_HDF_HDF_SCAN_DATA_WRITER_H_

#include <string>
#include <iostream>
#include "HDFFile.hpp"
#include "HDFGroup.hpp"
#include "HDFAtom.hpp"
#include "Enumerations.h"
#include "reads/ScanData.hpp"

class HDFScanDataWriter {
private:
    HDFGroup * rootGroupPtr;
	HDFGroup scanDataGroup;
	HDFGroup acqParamsGroup;
    HDFGroup dyeSetGroup;
	HDFGroup runInfoGroup;

	HDFAtom<std::string> whenStartedAtom;
	HDFAtom<float> frameRateAtom;
	HDFAtom<unsigned int> numFramesAtom;

    HDFAtom<std::string> baseMapAtom;
    HDFAtom<unsigned int> numAnalogAtom;

	HDFAtom<std::string> movieNameAtom;
	HDFAtom<std::string> runCodeAtom;

	HDFAtom<std::string> bindingKitAtom;
	HDFAtom<std::string> sequencingKitAtom;

	HDFAtom<unsigned int> platformIdAtom;
	HDFAtom<std::string> platformNameAtom;

    void CreateAcqParamsGroup();

    void CreateDyeSetGroup();

    void CreateRunInfoGroup();

public:
	HDFScanDataWriter(HDFFile & _outFile);

    HDFScanDataWriter(HDFGroup & _rootGroup);

    ~HDFScanDataWriter();
    
    int Initialize(HDFGroup & _rootGroup);
      
    void Write(ScanData & scanData);
   
	void WriteFrameRate(float frameRate);

    void WriteNumFrames(unsigned int numFrames);

    void WriteWhenStarted(const std::string whenStarted);

    std::string BaseMapToStr(std::map<char, int> & baseMap);
   
    void WriteBaseMap(const std::string baseMapStr);
   
    void WriteNumAnalog(const unsigned int numAnalog);

    void WritePlatformId(const PlatformId id);
   
    void WriteMovieName(const std::string movieName);

    void WriteRunCode(const std::string runCode);

	void Close();

private:
    void WriteBindingKit(const std::string & bindingKit);

    void WriteSequencingKit(const std::string & sequencingKit);
};

#endif

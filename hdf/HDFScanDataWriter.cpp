#include "HDFScanDataWriter.hpp"

void HDFScanDataWriter::CreateAcqParamsGroup() {
    if (acqParamsGroup.Initialize(scanDataGroup, "AcqParams") == 0) {
        std::cout << "ERROR could not create /ScanData/AcqParams." << std::endl;
        exit(1);
    }
    frameRateAtom.Create(acqParamsGroup.group, "FrameRate");
    numFramesAtom.Create(acqParamsGroup.group, "NumFrames");
    whenStartedAtom.Create(acqParamsGroup.group, "WhenStarted");
}

void HDFScanDataWriter::CreateDyeSetGroup(){
    if (dyeSetGroup.Initialize(scanDataGroup, "DyeSet") == 0) {
        std::cout << "ERROR could not create /ScanData/DyeSet." << std::endl;
        exit(1);
    }
    baseMapAtom.Create(dyeSetGroup.group, "BaseMap");
    numAnalogAtom.Create(dyeSetGroup.group, "NumAnalog");
}

void HDFScanDataWriter::CreateRunInfoGroup(){
    if (runInfoGroup.Initialize(scanDataGroup, "RunInfo") == 0) {
        std::cout << "ERROR, could not create /ScanDta/RunInfo." << std::endl;
        exit(1);
    }
    movieNameAtom.Create(runInfoGroup.group, "MovieName");
    platformIdAtom.Create(runInfoGroup.group, "PlatformId");
    platformNameAtom.Create(runInfoGroup.group, "PlatformName");
    runCodeAtom.Create(runInfoGroup.group, "RunCode");
    bindingKitAtom.Create(runInfoGroup.group, "BindingKit");
    sequencingKitAtom.Create(runInfoGroup.group, "SequencingKit");
}

HDFScanDataWriter::HDFScanDataWriter(HDFFile & _outFile) {
    Initialize(_outFile.rootGroup);
}

HDFScanDataWriter::HDFScanDataWriter(HDFGroup & _rootGroup) {
    Initialize(_rootGroup);
}

HDFScanDataWriter::~HDFScanDataWriter() { 
    this->Close();
}

int HDFScanDataWriter::Initialize(HDFGroup & _rootGroup) {
    rootGroupPtr = &(_rootGroup);
    rootGroupPtr->AddGroup("ScanData"); 
    if (scanDataGroup.Initialize(*(rootGroupPtr), "ScanData") == 0) {
        std::cout << "ERROR, could not create /ScanData group." << std::endl;
        exit(1);
    }
    scanDataGroup.AddGroup("AcqParams");
    scanDataGroup.AddGroup("DyeSet");
    scanDataGroup.AddGroup("RunInfo");

    CreateAcqParamsGroup();
    CreateDyeSetGroup();
    CreateRunInfoGroup();

    return 1;
}

void HDFScanDataWriter::Write(ScanData & scanData) {
    WriteFrameRate((scanData.frameRate==0)?
            (75):(scanData.frameRate));
    WriteNumFrames((scanData.numFrames==0)?
            (1000000):(scanData.numFrames));
    WriteWhenStarted((scanData.whenStarted.empty())?
            ("2013-01-01T01:01:01"):(scanData.whenStarted));
    std::string baseMapStr = BaseMapToStr(scanData.baseMap);
    WriteBaseMap((baseMapStr == "")?("TGAC"):baseMapStr);
    WriteNumAnalog(4);

    WriteMovieName((scanData.movieName.empty()?
                ("simulated_movie"):scanData.movieName));
    WriteRunCode((scanData.runCode.empty())?
            "simulated_runcode":(scanData.runCode));
    WritePlatformId((scanData.platformId==NoPlatform)?
            (Springfield):(scanData.platformId));
    WriteBindingKit(scanData.BindingKit());
    WriteSequencingKit(scanData.SequencingKit());
}

void HDFScanDataWriter::WriteFrameRate(float frameRate) {
    // Write /ScanData/AcqParams/FrameRate attribute.
    frameRateAtom.Write(frameRate);
}

void HDFScanDataWriter::WriteNumFrames(unsigned int numFrames) {
    // Write /ScanData/AcqParams/NumFrames attribute.
    numFramesAtom.Write(numFrames);
}

void HDFScanDataWriter::WriteWhenStarted(const std::string whenStarted) {
    // Write /ScanData/AcqParams/WhenStarted attribute.
    whenStartedAtom.Write(whenStarted);
}

std::string HDFScanDataWriter::BaseMapToStr(std::map<char, int> & baseMap) {
    std::string baseMapStr = ""; //4 dye channels.
    if (not baseMap.empty()) {
        baseMapStr = "    ";
        map<char, int>::iterator it;
        for (it = baseMap.begin(); it != baseMap.end(); ++it){
            if (it->second > 4 or it->second < 0) {
                std::cout << "ERROR, there are more than four dye channels."
                          << std::endl;
                exit(1);
            }
            baseMapStr[it->second]= it->first;
        }
    }
    return baseMapStr;
}

void HDFScanDataWriter::WriteBaseMap(const std::string baseMapStr) {
    //Write /ScanData/DyeSet/BaseMap attribute.
    baseMapAtom.Write(baseMapStr);
}

void HDFScanDataWriter::WriteNumAnalog(const unsigned int numAnalog) {
    //Write /ScanData/DyeSet/NumAnalog attribute.
    numAnalogAtom.Write(numAnalog);
}

void HDFScanDataWriter::WritePlatformId(const PlatformId id) {
    //Write /ScanData/RunInfo/Platform attribute.
    std::string name = (id == Springfield)?"Springfield":"Astro";
    platformIdAtom.Write(id);
    platformNameAtom.Write(name);
}

void HDFScanDataWriter::WriteMovieName(const std::string movieName) {
    //Write /ScanData/RunInfo/MovieName attribute.
    movieNameAtom.Write(movieName);
}

void HDFScanDataWriter::WriteRunCode(const std::string runCode) {
    //Write /ScanData/RunInfo/MovieName attribute.
    runCodeAtom.Write(runCode);
}

void HDFScanDataWriter::WriteBindingKit(const std::string & bindingKit) {
    bindingKitAtom.Write(bindingKit);
}

void HDFScanDataWriter::WriteSequencingKit(const std::string & sequencingKit) {
    sequencingKitAtom.Write(sequencingKit);
}

void HDFScanDataWriter::Close() {
    // Close /ScanData/AcqParams attributes.
    whenStartedAtom.Close();
    frameRateAtom.Close();
    numFramesAtom.Close();

    // Close /ScanData/DyeSet attributes.
    baseMapAtom.Close();

    // Close /ScanData/RunInfo attributes.
    movieNameAtom.Close();
    runCodeAtom.Close();
    platformIdAtom.Close();
    platformNameAtom.Close();
    bindingKitAtom.Close();
    sequencingKitAtom.Close();

    // Close /ScanData/AcqParams|DyeSet|RunInfo.
    acqParamsGroup.Close();
    dyeSetGroup.Close();
    runInfoGroup.Close();

    // Close /ScanData
    scanDataGroup.Close();
}

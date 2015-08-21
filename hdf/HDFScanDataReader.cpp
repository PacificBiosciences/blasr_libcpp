#include "HDFScanDataReader.hpp"

using namespace std;

HDFScanDataReader::HDFScanDataReader() {
    //
    // Assume the file is written without a movie name.  This is
    // flipped when a movie name is found.
    //
    Reset();
}

void HDFScanDataReader::Reset() {
    useMovieName    = false;
    useRunCode      = false;
    useWhenStarted  = false;
    fileHasScanData = false;
    movieName = "";
    runCode   = "";
    platformId      = NoPlatform;
    initializedAcqParamsGroup = initializedRunInfoGroup = false;
}

int HDFScanDataReader::InitializeAcqParamsAtoms() {
    if (frameRateAtom.Initialize(acqParamsGroup.group, "FrameRate") == 0) { return 0; }
    if (numFramesAtom.Initialize(acqParamsGroup.group, "NumFrames") == 0) { return 0; }
    if (acqParamsGroup.ContainsAttribute("WhenStarted")) {
        if (whenStartedAtom.Initialize(acqParamsGroup.group, "WhenStarted") == 0) { return 0; }
        useWhenStarted = true;
    }
    return 1;
}

//
// This is created on top of a file that is already opened, so
// instead of initializing by opening a file, it is initialized by
// passing the root group that should contain the ScanData group.  
// When shutting down, no file handle needs to be closed.
//
int HDFScanDataReader::Initialize(HDFGroup *pulseDataGroup) {

    //
    // Initiailze groups for reading data.
    //

    initializedAcqParamsGroup = false;
    initializedRunInfoGroup   = false;
    if (pulseDataGroup->ContainsObject("ScanData") == 0 or 
            scanDataGroup.Initialize(pulseDataGroup->group, "ScanData") == 0) {
        return 0;
    }
    fileHasScanData = true;

    if (scanDataGroup.ContainsObject("DyeSet") == 0 or
            dyeSetGroup.Initialize(scanDataGroup.group, "DyeSet") == 0) {
        return 0;
    }

    if (scanDataGroup.ContainsObject("AcqParams") == 0 or
            acqParamsGroup.Initialize(scanDataGroup.group, "AcqParams") == 0) {
        return 0;
    }
    initializedAcqParamsGroup = true;

    if (scanDataGroup.ContainsObject("RunInfo") == 0 or
            runInfoGroup.Initialize(scanDataGroup.group, "RunInfo") == 0) {
        return 0;
    }
    initializedRunInfoGroup = true;
    if (InitializeAcqParamsAtoms() == 0) {
        return 0;
    }

    //
    // Read in the data that will be used later on either per read or
    // when the entire bas/pls file is read.
    //

    if (ReadPlatformId(platformId) == 0) {
        return 0;
    }

    if (runInfoGroup.ContainsAttribute("RunCode") and
            runCodeAtom.Initialize(runInfoGroup, "RunCode")) {
        useRunCode = true;
    }

    //
    // Load baseMap which maps bases (ATGC) to channel orders.
    // This should always be present.
    //
    if (LoadBaseMap(baseMap) == 0)
        return 0;

    //
    // Attempt to load the movie name.  This is not always present.
    //
    LoadMovieName(movieName);

    return 1;
}

string HDFScanDataReader::GetMovieName() {
    // If this object is correctly initialized, movieName
    // is guaranteed to be loaded if it exists, no need to reload.
    return movieName;
}

string HDFScanDataReader::GetRunCode() {
    return runCode;
}

int HDFScanDataReader::Read(ScanData &scanData) {
    // All parameters below are required.
    if (ReadPlatformId(scanData.platformId) == 0) return 0;
    LoadMovieName(scanData.movieName);
    LoadBaseMap(scanData.baseMap);

    if (useRunCode) {
        runCodeAtom.Read(scanData.runCode);
    }
    frameRateAtom.Read(scanData.frameRate);
    numFramesAtom.Read(scanData.numFrames);

    if (useWhenStarted) {
        whenStartedAtom.Read(scanData.whenStarted);
    }

    ReadSequencingKit(scanData.sequencingKit_);

    ReadBindingKit(scanData.bindingKit_);

    return 1;
}

void HDFScanDataReader::ReadWhenStarted(string &whenStarted) {
    whenStartedAtom.Read(whenStarted);
}

PlatformId HDFScanDataReader::GetPlatformId() {
    return platformId;
}

int HDFScanDataReader::ReadPlatformId(PlatformId &pid) {
    if (runInfoGroup.ContainsAttribute("PlatformId")) {
        if (platformIdAtom.Initialize(runInfoGroup, "PlatformId") == 0) {
            return 0;
        }
        platformIdAtom.Read((unsigned int&)pid);
    }
    else {
        pid = Astro;
    }
    return 1;
}

int HDFScanDataReader::ReadStringAttribute(std::string & attributeValue,
        const std::string & attributeName, HDFGroup & group,
        HDFAtom<std::string> & atom) {

    if (group.ContainsAttribute(attributeName) and
        (atom.isInitialized or atom.Initialize(group, attributeName))) {
        atom.Read(attributeValue);
        return 1;
    } else {
        return 0;
    }
}

int HDFScanDataReader::ReadBindingKit(std::string &bindingKit)
{
    return ReadStringAttribute(bindingKit, "BindingKit", runInfoGroup, bindingKitAtom);
}

int HDFScanDataReader::ReadSequencingKit(std::string &sequencingKit)
{
    return ReadStringAttribute(sequencingKit, "SequencingKit", runInfoGroup, sequencingKitAtom);
}

int HDFScanDataReader::LoadMovieName(string &movieNameP) {
    // Groups for building read names
    if (ReadStringAttribute(movieNameP, "MovieName", runInfoGroup, movieNameAtom) == 0) {
        // Internal analysis may manually edit the movie name and set STRSIZE to a value
        // which != movie name length. Handle this case. 
        movieNameP = string(movieNameP.c_str()); 
        return 0;
    } else {
        useMovieName = true;
        int e = movieNameP.size() - 1;
        while (e > 0 and movieNameP[e] == ' ') e--;
        movieNameP = movieNameP.substr(0, e+1);
        movieNameP = string(movieNameP.c_str());
        return 1;
    }
}

int HDFScanDataReader::LoadBaseMap(map<char, int> & baseMap) {
    // Map bases to channel order in hdf pls file.
    if (dyeSetGroup.ContainsAttribute("BaseMap") and
            baseMapAtom.Initialize(dyeSetGroup, "BaseMap")) {
        string baseMapStr;
        baseMapAtom.Read(baseMapStr);
        if (baseMapStr.size() != 4) {
            cout << "ERROR, there are more than four types of bases "
                << "according to /ScanData/DyeSet/BaseMap." << endl;
            exit(1);
        }
        baseMap.clear();
        for(size_t i = 0; i < baseMapStr.size(); i++) {
            baseMap[toupper(baseMapStr[i])] = i;
            baseMap[tolower(baseMapStr[i])] = i;
        }
        return 1;
    }
    return 0;
}

void HDFScanDataReader::Close() {
    if (useMovieName) {
        movieNameAtom.Close();
    }
    if (useRunCode) {
        runCodeAtom.Close();
    }
    if (useWhenStarted) {
        whenStartedAtom.Close();
    }

    baseMapAtom.Close();
    platformIdAtom.Close();
    frameRateAtom.Close();
    numFramesAtom.Close();
    sequencingKitAtom.Close();
    bindingKitAtom.Close();

    scanDataGroup.Close();
    dyeSetGroup.Close();
    acqParamsGroup.Close();
    runInfoGroup.Close();
    Reset();
}


std::string HDFScanDataReader::GetMovieName_and_Close(std::string & fileName) {
    HDFFile file;
    file.Open(fileName, H5F_ACC_RDONLY);

    fileHasScanData = false;
    if (file.rootGroup.ContainsObject("ScanData") == 0 or 
        scanDataGroup.Initialize(file.rootGroup, "ScanData") == 0) {
        return ""; 
    }
    fileHasScanData = true;

    initializedRunInfoGroup = false;
    if (scanDataGroup.ContainsObject("RunInfo") == 0 or
        runInfoGroup.Initialize(scanDataGroup.group, "RunInfo") == 0) {
        return "";
    }
    initializedRunInfoGroup = true;
    
    string movieName;
    LoadMovieName(movieName);
    Close();
    file.Close();
    return movieName;
}

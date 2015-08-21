#include "ScanData.hpp"

ScanData::ScanData() {
    platformId = NoPlatform;
    frameRate = 0.0;
    numFrames = 0;
    movieName = runCode = whenStarted = "";
    baseMap.clear();
}

std::string ScanData::GetMovieName() {
    return movieName;
}

ScanData & ScanData::PlatformID(const PlatformId & id) {
    platformId = id;
    return *this;
}
ScanData & ScanData::FrameRate(const float & rate) {
    frameRate = rate;
    return *this;
}
ScanData & ScanData::NumFrames(const unsigned int & num) {
    numFrames = num;
    return *this;
}
ScanData & ScanData::MovieName(const std::string & name) {
    movieName = name;
    return *this;
}
ScanData & ScanData::RunCode(const std::string & code) {
    runCode = code;
    return *this;
}
ScanData & ScanData::WhenStarted(const std::string & when) {
    whenStarted = when;
    return *this;
}
ScanData & ScanData::BaseMap(const std::map<char, int> bmp) {
    baseMap = bmp;
    return *this;
}
ScanData & ScanData::SequencingKit(const std::string sequencingKit) {
    sequencingKit_ = sequencingKit;
    return *this;
}
ScanData & ScanData::BindingKit(const std::string bindingKit) {
    bindingKit_ = bindingKit;
    return *this;
}

PlatformId ScanData::PlatformID(void) const {
    return platformId;
}
float ScanData::FrameRate(void) const {
    return frameRate;
}
unsigned int ScanData::NumFrames(void) const {
    return numFrames;
}
std::string ScanData::MovieName(void) const {
    return movieName;
}
std::string ScanData::RunCode(void) const {
    return runCode;
}
std::string ScanData::WhenStarted(void) const {
    return whenStarted;
}
std::map<char, int> ScanData::BaseMap(void) const {
    return baseMap;
}
std::string ScanData::SequencingKit(void) const {
    return sequencingKit_;
}
std::string ScanData::BindingKit(void) const {
    return bindingKit_;
}

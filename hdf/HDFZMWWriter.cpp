#include "libconfig.h"
#ifdef USE_PBBAM

#include "HDFZMWWriter.hpp"

HDFZMWWriter::HDFZMWWriter(const std::string & filename, 
                           HDFGroup & parentGroup,
                           const bool inPulseCalls,
                           const std::map<char, size_t> & baseMap)
    : HDFWriterBase(filename)
    , parentGroup_(parentGroup)
    , baseMap_(baseMap)
    , inPulseCalls_(inPulseCalls)
    , arrayLength_(0)
{ 
    if (not parentGroup.groupIsInitialized)
        PARENT_GROUP_NOT_INITIALIZED_ERROR(PacBio::GroupNames::zmw);
    else {
        parentGroup_.AddGroup(PacBio::GroupNames::zmw);

        if (zmwGroup_.Initialize(parentGroup_, PacBio::GroupNames::zmw) == 0)
            FAILED_TO_CREATE_GROUP_ERROR(PacBio::GroupNames::zmw);

        this->InitializeChildHDFGroups();
    }
}

HDFZMWWriter::HDFZMWWriter(const std::string & filename, 
                           HDFGroup & parentGroup) 
    : HDFZMWWriter(filename, parentGroup, false, {})
{ }

HDFZMWWriter::~HDFZMWWriter(void) {
    this->_WriteAttributes();
    this->Close();
}

bool HDFZMWWriter::InPulseCalls(void) const {
    return inPulseCalls_;
}

bool HDFZMWWriter::HasBaseLineSigma(void) const {
    return false; // base line sigma won't be provided util 3.0.1 bug 29486 
    //return InPulseCalls() and baseLineSigmaArray_.IsInitialized();
}
    
bool HDFZMWWriter::WriteOneZmw(const SMRTSequence & read) {
    _WriteNumEvent(read.length);
    _WriteHoleNumber(read.HoleNumber());
    _WriteHoleXY(static_cast<int16_t>(read.HoleX()),
                 static_cast<int16_t>(read.HoleY()));
    _WriteHoleStatus(read.HoleStatus());
    arrayLength_++;
    return Errors().empty();
}

bool HDFZMWWriter::WriteOneZmw(const PacBio::BAM::BamRecord & read) {
    if (InPulseCalls()) {
        if (not read.HasPulseCall()) {
            AddErrorMessage(std::string("PulseCall absent in read ") + read.FullName());
        } else {
            _WriteNumEvent(read.PulseCall().size());
        }
    } else {
        _WriteNumEvent(read.Sequence().size());
    }
    arrayLength_++;

    uint32_t hn_ = read.HoleNumber();
    _WriteHoleNumber(hn_);
    _WriteHoleXY(static_cast<int16_t>(hn_ & 0x0000FFFF), 
                 static_cast<int16_t>(hn_ >> 16));
    _WriteHoleStatus(PacBio::AttributeValues::ZMW::HoleStatus::sequencingzmw);
    _WriteBaseLineSigma(read);
    return Errors().empty();
}

bool HDFZMWWriter::WriteFakeDataSets(void) {
    // Fake BaselineLevel, BaselineSigma, SignalLevel, SignalSigma
    //std::vector<float> buffer(arrayLength_);
    //std::fill(buffer.begin(), buffer.end(), 0.0);
    return __WriteFake2DDataSet<float>(zmwGroup_, PacBio::GroupNames::baselinelevel, arrayLength_, 4, 0) and
           __WriteFake2DDataSet<float>(zmwGroup_, PacBio::GroupNames::baselinesigma, arrayLength_, 4, 0) and
           __WriteFake2DDataSet<float>(zmwGroup_, PacBio::GroupNames::signallevel, arrayLength_, 4, 0) and
           __WriteFake2DDataSet<float>(zmwGroup_, PacBio::GroupNames::signalsigma, arrayLength_, 4, 0);
}

bool HDFZMWWriter::_WriteNumEvent(const uint32_t numEvent) {
    numEventArray_.Write(&numEvent, 1);
    return true;
}

bool HDFZMWWriter::_WriteHoleNumber(const uint32_t holeNumber) {
    holeNumberArray_.Write(&holeNumber, 1);
    return true;
}

bool HDFZMWWriter::_WriteHoleXY(const int16_t holeX, const int16_t holeY) {
    int16_t xy[2] = {holeX, holeY};
    holeXYArray_.WriteRow(xy, 2);
    return true;
}

bool HDFZMWWriter::_WriteHoleStatus(const unsigned char holeStatus) {
    // NOTE that: We assume that all zmws in BAM are SEQUENCING zmws.
    holeStatusArray_.Write(&holeStatus, 1);
    return true;
}

bool HDFZMWWriter::_WriteBaseLineSigma(const PacBio::BAM::BamRecord & read) {
    // FIXME: pbbam should provide APIs: HasBaseLineSigma and BaseLineSigma().
    if (HasBaseLineSigma()) {
        if (read.Impl().HasTag("bs")) {
            const PacBio::BAM::Tag & tag = read.Impl().TagValue("bs");
            std::vector<float> data = tag.ToFloatArray();
            if (data.size() != 4) {
                AddErrorMessage(std::string("Tag BaseLineSigma must have 4 values per each record in read " + read.FullName()));
            }
            float bls[4]; 
            bls[baseMap_['A']] = data[0]; bls[baseMap_['C']] = data[1];
            bls[baseMap_['G']] = data[2]; bls[baseMap_['T']] = data[3];
            baseLineSigmaArray_.WriteRow(bls, 4);
        } else {
            AddErrorMessage(std::string("Tag BaseLineSigma is absent in read ") + read.FullName());
        }
    } else {
        //BaseLineSigma won't be available uill v 3.0.1.
        //AddErrorMessage("BaseLineSigma array is not initialized.");
    }
    return Errors().empty();
}

void HDFZMWWriter::Flush(void) {
    // Mandatory metrics
    numEventArray_.Flush();
    holeNumberArray_.Flush();
    holeStatusArray_.Flush();
    holeXYArray_.Flush();

    // Optional metrics
    if (HasBaseLineSigma()) baseLineSigmaArray_.Flush();
}

void HDFZMWWriter::Close(void) {
    this->Flush();

    // Mandatory metrics
    numEventArray_.Close();
    holeNumberArray_.Close();
    holeStatusArray_.Close();
    holeXYArray_.Close();

    // Optional metrics
    if (HasBaseLineSigma()) baseLineSigmaArray_.Close();

    zmwGroup_.Close();
}

bool HDFZMWWriter::InitializeChildHDFGroups(void) {
    // Mandatory metrics
    if (numEventArray_.Initialize(zmwGroup_, PacBio::GroupNames::numevent) == 0) { 
        FAILED_TO_CREATE_GROUP_ERROR(PacBio::GroupNames::numevent);
    }

    if (holeNumberArray_.Initialize(zmwGroup_, PacBio::GroupNames::holenumber) == 0) {
        FAILED_TO_CREATE_GROUP_ERROR(PacBio::GroupNames::holenumber);
    }

    if (holeStatusArray_.Initialize(zmwGroup_, PacBio::GroupNames::holestatus) == 0) {
        FAILED_TO_CREATE_GROUP_ERROR(PacBio::GroupNames::holestatus);
    }

    if (holeXYArray_.Initialize(zmwGroup_, PacBio::GroupNames::holexy, 2) == 0) {
        FAILED_TO_CREATE_GROUP_ERROR(PacBio::GroupNames::holexy);
    }

    /* BaseLineSigma not available until 3.0.1
    if (InPulseCalls()) {
        if (baseLineSigmaArray_.Initialize(zmwGroup_, PacBio::GroupNames::baselinesigma, 4) == 0) {
            FAILED_TO_CREATE_GROUP_ERROR(PacBio::GroupNames::baselinesigma);
        }
    }
    */
    return Errors().empty();
}

void HDFZMWWriter::_WriteAttributes(void)
{
    // NumEvent has no attributes.

    if (holeNumberArray_.IsInitialized() and holeNumberArray_.size() > 0) {
        AddAttribute<std::string>(holeNumberArray_, PacBio::AttributeNames::Common::description, PacBio::AttributeValues::ZMW::HoleNumber::description);
    }

    if (holeStatusArray_.IsInitialized() and holeStatusArray_.size() > 0) {
        AddAttribute<std::string>(holeStatusArray_, PacBio::AttributeNames::Common::description, PacBio::AttributeValues::ZMW::HoleStatus::description);
        AddAttribute<std::vector<std::string>>(holeStatusArray_, PacBio::AttributeNames::ZMW::HoleStatus::lookuptable, PacBio::AttributeValues::ZMW::HoleStatus::lookuptable);
    }

    if (holeXYArray_.IsInitialized() and holeXYArray_.GetNRows() > 0) {
        AddAttribute<std::string>(holeXYArray_, PacBio::AttributeNames::Common::description, PacBio::AttributeValues::ZMW::HoleXY::description);
    }

    if (HasBaseLineSigma() and baseLineSigmaArray_.GetNRows() > 0) {
        AddAttribute<std::string>(baseLineSigmaArray_, PacBio::AttributeNames::Common::description, PacBio::AttributeValues::ZMW::BaseLineSigma::description);
    }
}
#endif

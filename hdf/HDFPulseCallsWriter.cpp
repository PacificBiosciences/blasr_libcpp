#include "../pbdata/libconfig.h"
#ifdef USE_PBBAM
#include <ctype.h>
#include <sstream>
#include <algorithm>
#include "HDFPulseCallsWriter.hpp"
#include "../pbdata/utils/TimeUtils.hpp"

const std::vector<PacBio::BAM::BaseFeature>  HDFPulseCallsWriter::ValidQVEnums = 
{     PacBio::BAM::BaseFeature::PULSE_CALL
    , PacBio::BAM::BaseFeature::LABEL_QV
    , PacBio::BAM::BaseFeature::PKMID
    , PacBio::BAM::BaseFeature::PKMEAN
    , PacBio::BAM::BaseFeature::PULSE_MERGE_QV
    , PacBio::BAM::BaseFeature::ALT_LABEL
    , PacBio::BAM::BaseFeature::ALT_LABEL_QV
    , PacBio::BAM::BaseFeature::START_FRAME
    , PacBio::BAM::BaseFeature::PULSE_CALL_WIDTH 
};

std::vector<PacBio::BAM::BaseFeature> 
HDFPulseCallsWriter::WritableQVs(const std::vector<PacBio::BAM::BaseFeature> & qvsToWrite) {
    std::vector<PacBio::BAM::BaseFeature> ret;
    for(auto qv : qvsToWrite) {
        // Filter qvs which are not in format specification.
        if (std::find(ValidQVEnums.begin(), ValidQVEnums.end(), qv) != ValidQVEnums.end()) {
            if (std::find(ret.begin(), ret.end(), qv) == ret.end()) ret.push_back(qv); 
        }
    }
    return ret;
}

HDFPulseCallsWriter::HDFPulseCallsWriter(const std::string & filename,
                                         HDFGroup & parentGroup,
                                         const std::map<char, size_t> & baseMap,
                                         const std::string & basecallerVersion,
                                         const std::vector<PacBio::BAM::BaseFeature> & qvsToWrite)
    : HDFWriterBase(filename)
    , parentGroup_(parentGroup)
    , baseMap_(baseMap)
    , qvsToWrite_({}) // Input qvsToWrite must be checked.
    , zmwWriter_(nullptr)
    , arrayLength_(0)
    , basecallerVersion_(basecallerVersion)
{
    // Add PulseCalls as a child group to the parent group.
    AddChildGroup(parentGroup_, pulsecallsGroup_, PacBio::GroupNames::pulsecalls);

    if (basecallerVersion_.empty()) {
        AddErrorMessage("BaseCallerVersion must not be empty!");
        return;
    }
 
    this->qvsToWrite_ = WritableQVs(qvsToWrite);

    // Any QVs to write?
    if (qvsToWrite_.size() == 0) {
        AddErrorMessage("No QVs to write.");
        return;
    }

    // Reject if it does not have StartFrame
    if (not _HasQV(PacBio::BAM::BaseFeature::START_FRAME)) {
        AddErrorMessage("Bam input must contain tag 'StartFrame' using PPA 2.0.0 or later.");
        return;
    }

    // Initialize QV groups
    if (not InitializeQVGroups()) {
        AddErrorMessage("Failed to initialize QV Groups.");
        return;
    }

    // Create a zmwWriter.
    zmwWriter_.reset(new HDFZMWWriter(Filename(), pulsecallsGroup_, true, baseMap));

    inverseGain_ = 1.0f;

    // Note: ignore /PulseCalls/ZMWMetrics none of its metrics exist in BAM.
}

uint32_t HDFPulseCallsWriter::NumZMWs(void) const {
    if (zmwWriter_) return zmwWriter_->NumZMWs();
    else return 0;
}

void HDFPulseCallsWriter::Content(std::vector<std::string> & names,
                                  std::vector<std::string> & types) const {
    // Print order matters;
    const bool fakeChi2 =true, fakeMaxSignal = true, fakeMidStdDev = true;

    const std::string uint8_t_str        = "uint8_t";
    const std::string uint16_t_str       = "uint16_t";
    const std::string uint32_t_str       = "uint32_t";

    const std::string channeltype        = uint8_t_str;
    const std::string ispulsetype        = uint8_t_str;
    const std::string labelqvtype        = uint8_t_str;
    const std::string mergeqvtype        = uint8_t_str;
    const std::string altlabeltype       = uint8_t_str;
    const std::string altlabelqvtype     = uint8_t_str;

    const std::string chi2type           = uint16_t_str;
    const std::string maxsignaltype      = uint16_t_str;
    const std::string meansignaltype     = uint16_t_str;
    const std::string midsignaltype      = uint16_t_str;
    const std::string midstddevtype      = uint16_t_str;
    const std::string widthinframestype  = uint16_t_str;

    const std::string startframetype     = uint32_t_str;

    names.clear(); types.clear();

    if (HasAltLabel()) {
        names.push_back(PacBio::GroupNames::altlabel);
        types.push_back(altlabeltype);
    }

    if (HasAltLabelQV()) {
        names.push_back(PacBio::GroupNames::altlabelqv);
        types.push_back(altlabelqvtype);
    }

    if (HasPulseCall()) {
        names.push_back(PacBio::GroupNames::channel);
        types.push_back(channeltype);
    }

    if (fakeChi2) { // fake chi2
        names.push_back(PacBio::GroupNames::chi2);
        types.push_back(chi2type);
    }

    if (HasIsPulse()) {
        names.push_back(PacBio::GroupNames::ispulse);
        types.push_back(ispulsetype);
    }

    if (HasLabelQV()) {
        names.push_back(PacBio::GroupNames::labelqv);
        types.push_back(labelqvtype);
    }

    if (fakeMaxSignal) { // fake maxsignal
        names.push_back(PacBio::GroupNames::maxsignal);
        types.push_back(maxsignaltype);
    }

    if (HasPkmean()) {
        names.push_back(PacBio::GroupNames::meansignal);
        types.push_back(meansignaltype);
    }

    if (HasPulseMergeQV()) {
        names.push_back(PacBio::GroupNames::mergeqv);
        types.push_back(mergeqvtype);
    }

    if (HasPkmid()) {
        names.push_back(PacBio::GroupNames::midsignal);
        types.push_back(midsignaltype);
    }

    if (fakeMidStdDev) { // fake MidStdDev
        names.push_back(PacBio::GroupNames::midstddev);
        types.push_back(midstddevtype);
    }

    if (HasStartFrame()) {
        names.push_back(PacBio::GroupNames::startframe);
        types.push_back(startframetype);
    }

    if (HasPulseCallWidth()) {
        names.push_back(PacBio::GroupNames::widthinframes);
        types.push_back(widthinframestype);
    }
}

void HDFPulseCallsWriter::SetInverseGain(const float igain) {
    inverseGain_ = igain;
}

bool HDFPulseCallsWriter::_WriteAttributes(void) {
    std::vector<std::string> content_names, content_types;
    Content(content_names, content_types);

    // ChangeListID
    bool OK = 
        AddAttribute<std::string>(pulsecallsGroup_, 
                     PacBio::AttributeNames::Common::changelistid,
                     basecallerVersion_) and

    // Content
        AddAttribute<std::vector<std::string>>(pulsecallsGroup_,
                     PacBio::AttributeNames::Common::content,
                     content_names) and

    // ContentStored
        AddAttribute<uint32_t>(pulsecallsGroup_,
                     PacBio::AttributeNames::Common::contentstored,
                     NumZMWs()) and

    // DateCreated
        AddAttribute<std::string>(pulsecallsGroup_,
                     PacBio::AttributeNames::Common::datacreated,
                     GetTimestamp()) and

    // SchemaRevision
        AddAttribute<std::string>(pulsecallsGroup_,
                     PacBio::AttributeNames::Common::schemarevision,
                     PacBio::AttributeValues::Common::schemarevision);

    return OK;
}


std::vector<std::string> HDFPulseCallsWriter::Errors(void) const {
    std::vector<std::string> retErrors = this->errors_;
    if (zmwWriter_) {
        const std::vector<std::string> & zmwErrors = zmwWriter_->Errors();
        retErrors.insert(retErrors.end(), zmwErrors.begin(), zmwErrors.end());
    }
    return retErrors;
}

HDFPulseCallsWriter::~HDFPulseCallsWriter(void) { 
    this->Close();
}

bool HDFPulseCallsWriter::InitializeQVGroups(void) {
    int ret = 1;
    if (_HasQV(PacBio::BAM::BaseFeature::PULSE_CALL)) 
        ret *= pulseCallArray_.Initialize(pulsecallsGroup_,      PacBio::GroupNames::channel);
        ret *= isPulseArray_.Initialize(pulsecallsGroup_,        PacBio::GroupNames::ispulse);
    if (_HasQV(PacBio::BAM::BaseFeature::LABEL_QV))
        ret *= labelQVArray_.Initialize(pulsecallsGroup_,        PacBio::GroupNames::labelqv);
    if (_HasQV(PacBio::BAM::BaseFeature::PKMEAN))
        ret *= pkmeanArray_.Initialize(pulsecallsGroup_,         PacBio::GroupNames::meansignal, 4);
    if (_HasQV(PacBio::BAM::BaseFeature::PULSE_MERGE_QV))
        ret *= pulseMergeQVArray_.Initialize(pulsecallsGroup_,   PacBio::GroupNames::mergeqv);
    if (_HasQV(PacBio::BAM::BaseFeature::PKMID))
        ret *= pkmidArray_.Initialize(pulsecallsGroup_,          PacBio::GroupNames::midsignal);
    if (_HasQV(PacBio::BAM::BaseFeature::START_FRAME))
        ret *= startFrameArray_.Initialize(pulsecallsGroup_,     PacBio::GroupNames::startframe);
    if (_HasQV(PacBio::BAM::BaseFeature::PULSE_CALL_WIDTH))
        ret *= pulseCallWidthArray_.Initialize(pulsecallsGroup_, PacBio::GroupNames::widthinframes);
    if (_HasQV(PacBio::BAM::BaseFeature::ALT_LABEL))
        ret *= altLabelArray_.Initialize(pulsecallsGroup_,       PacBio::GroupNames::altlabel);
    if (_HasQV(PacBio::BAM::BaseFeature::ALT_LABEL_QV))
        ret *= altLabelQVArray_.Initialize(pulsecallsGroup_,     PacBio::GroupNames::altlabelqv);
    return (ret != 0);
}

bool HDFPulseCallsWriter::WriteOneZmw(const SMRTSequence & read) {
    const PacBio::BAM::BamRecord & record = read.bamRecord;

    if (zmwWriter_) 
        zmwWriter_->WriteOneZmw(record);

    _WritePulseCall(record); // Write PulseCall and IsPulse
    _WriteLabelQV(record);
    _WritePkmean(record);
    _WritePulseMergeQV(record);
    _WritePkmid(record);
    _WriteStartFrame(record);
    _WritePulseCallWidth(record);
    _WriteAltLabel(record);
    _WriteAltLabelQV(record); 

    arrayLength_ +=  record.PulseCall().size();

    return Errors().empty();
}

bool HDFPulseCallsWriter::WriteFakeDataSets() {
    uint32_t block_sz = 65536; // This is a data buffer.
    std::vector<uint16_t> buffer_uint16_5M_0(block_sz, 0);
    // Write 2D Array: Chi2, 1D Arrays: MaxSignal, MidStdDev 
    bool OK = __WriteFake2DDataSet<uint16_t>(pulsecallsGroup_, PacBio::GroupNames::chi2, arrayLength_, 4, 0) and 
              __WriteFakeDataSet<uint16_t>(pulsecallsGroup_, PacBio::GroupNames::maxsignal, arrayLength_, buffer_uint16_5M_0) and 
              __WriteFakeDataSet<uint16_t>(pulsecallsGroup_, PacBio::GroupNames::midstddev, arrayLength_, buffer_uint16_5M_0);

    if (zmwWriter_)
        return OK and zmwWriter_->WriteFakeDataSets();

    return OK;
}

bool HDFPulseCallsWriter::_CheckRead(const PacBio::BAM::BamRecord & read, 
                                     const uint32_t qvLength,
                                     const std::string & qvName) {
    // Tag 'pb' will not be available until 3.0.1, bug 29486
    (void)(read); (void)(qvLength); (void)(qvName);
    return true;
    /*
    // FIXME: pbbam should provide HasPulseBlockSize() and PulseBlockSize().
    if (read.Impl().HasTag("pb")) {
        if (qvLength != read.Impl().TagValue("pb").ToUInt32()) {
            AddErrorMessage(std::string(qvName) + "'s size does not match PulseBlockSize in read " + read.FullName());
            return false;
        } else {
            return true;
        }
    } else {
        AddErrorMessage(std::string("Tag PulseBlockSize is absent in read ") + read.FullName());
        return false;
    }
    */
}

bool HDFPulseCallsWriter::_WritePulseCall(const PacBio::BAM::BamRecord & read) {
    if (HasPulseCall()) {
        if (read.HasPulseCall()) {
            // Write both PulseCall (i.e., Channel) and IsPulse
            const std::string & pulsecall = read.PulseCall();
            const unsigned int length = pulsecall.size();
            _CheckRead(read, length, "PulseCall");
            std::vector<uint8_t> channel; channel.resize(length);
            uint32_t num_bases = 0;
            for (size_t i = 0; i < length; i++) {
                const char base = pulsecall[i];
                if (base == 'A' or base == 'C' or base == 'G' or base == 'T') {
                    channel[i] = static_cast<uint8_t>(baseMap_[std::toupper(base)]);
                    num_bases += 1;
                } else if (base == 'a' or base == 'c' or base == 'g' or base == 't') {
                    channel[i] = static_cast<uint8_t>(baseMap_[std::toupper(base)]);
                } else {
                    AddErrorMessage(std::string("Unrecognizable base in PulseCall of " + read.FullName()));
                }
            }
            if (num_bases != read.Sequence().size()) {
                AddErrorMessage(std::string("Number of bases in PulseCall does not match BaseCall in read ") + read.FullName());
            }
            pulseCallArray_.Write(&channel[0], length);

            std::vector<uint8_t> ispulse(length, 1); // all pulses in BAM are real pulses
            isPulseArray_.Write(&ispulse[0], length);
        } else {
            AddErrorMessage(std::string("PulseCall is absent in read " + read.FullName()));
        }
    }
    return Errors().empty();
}

bool HDFPulseCallsWriter::_WriteLabelQV(const PacBio::BAM::BamRecord & read) {
    if (HasLabelQV()) {
        if (read.HasLabelQV()) {
            const PacBio::BAM::QualityValues & qvs = read.LabelQV();
            std::vector<uint8_t> data(qvs.cbegin(), qvs.cend());
            _CheckRead(read, data.size(), "LabelQV");
            labelQVArray_.Write(&data[0], data.size());
        } else {
            AddErrorMessage(std::string("LabelQV is absent in read " + read.FullName()));
        }
    }
    return Errors().empty();
}

bool HDFPulseCallsWriter::_WritePkmean(const PacBio::BAM::BamRecord & read) {
    if (HasPkmean()) {
        if (read.HasPkmean()) {
            std::vector<float> pkmids = read.Pkmean();

            // convert from photoE to counts
            std::for_each(pkmids.begin(), pkmids.end(), [&](float& x){x*=inverseGain_;});

            // down-convert to ushorts (required by pulse file specification)
            std::vector<uint16_t> data(pkmids.begin(), pkmids.end());

            const std::string & pulsecall = read.PulseCall();
            for (size_t i = 0; i < pulsecall.size(); i++) {
                uint16_t pkm[4] = {0, 0, 0, 0};
                pkm[baseMap_[std::toupper(pulsecall[i])]] = data[i];
                pkmeanArray_.WriteRow(pkm, 4);
            }
        } else {
            AddErrorMessage(std::string("Pkmean is absent in read " + read.FullName()));
        }
    }
    return Errors().empty();
}

bool HDFPulseCallsWriter::_WritePulseMergeQV(const PacBio::BAM::BamRecord & read) {
    if (HasPulseMergeQV()) {
        if (read.HasPulseMergeQV()) {
            const PacBio::BAM::QualityValues & qvs = read.PulseMergeQV();
            std::vector<uint8_t> data(qvs.cbegin(), qvs.cend());
            _CheckRead(read, data.size(), "PulseMergeQV");
            pulseMergeQVArray_.Write(&data[0], data.size());
        } else {
            AddErrorMessage(std::string("PulseMergeQV is absent in read " + read.FullName()));
        }
    }
    return Errors().empty();
}

bool HDFPulseCallsWriter::_WritePkmid(const PacBio::BAM::BamRecord & read) {
    if (HasPkmid()) {
        if (read.HasPkmid()) {
            std::vector<float> pkmids = read.Pkmid();

            // convert from photoE to counts
            std::for_each(pkmids.begin(), pkmids.end(), [&](float& x){x*=inverseGain_;});

            // down-convert to ushorts (required by pulse file specification)
            std::vector<uint16_t> data(pkmids.begin(), pkmids.end());

            pkmidArray_.Write(&data[0], data.size());
        } else {
            AddErrorMessage(std::string("Pkmid is absent in read " + read.FullName()));
        }
    }
    return Errors().empty();
}

bool HDFPulseCallsWriter::_WriteStartFrame(const PacBio::BAM::BamRecord & read) {
    if (HasStartFrame()) {
        if (read.HasStartFrame()) {
            const std::vector<uint32_t> data = read.StartFrame();
            _CheckRead(read, data.size(), "StartFrame");
            startFrameArray_.Write(&data[0], data.size());
        } else {
            AddErrorMessage(std::string("StartFrame is absent in read " + read.FullName()));
        }
    }
    return Errors().empty();
}

bool HDFPulseCallsWriter::_WritePulseCallWidth(const PacBio::BAM::BamRecord & read) {
    if (HasPulseCallWidth()) {
        if (read.HasPulseCallWidth()) {
            //FIXME: pbbam PulseCallWidth().Data() returns incorrect vector (e.g., size)
            const PacBio::BAM::Tag & tag = read.Impl().TagValue("px");
            std::vector<uint16_t> data = tag.ToUInt16Array();
            _CheckRead(read, data.size(), "PulseCallWidth");
            pulseCallWidthArray_.Write(&data[0], data.size());
        } else {
            AddErrorMessage(std::string("PulseCallWidth is absent in read " + read.FullName()));
        }
    }
    return Errors().empty();
}

bool HDFPulseCallsWriter::_WriteAltLabel(const PacBio::BAM::BamRecord & read) {
    if (HasAltLabel()) {
        if (read.HasAltLabelTag()) {
            const std::string & tags = read.AltLabelTag();
            std::vector<unsigned char> data(tags.begin(), tags.end());
            _CheckRead(read, data.size(), "AltLabel");
            altLabelArray_.Write(&data[0], tags.size());
        } else {
            AddErrorMessage(std::string("AltLabel is absent in read " + read.FullName()));
        }
    }
    return Errors().empty();
}

bool HDFPulseCallsWriter::_WriteAltLabelQV(const PacBio::BAM::BamRecord & read) {
    if (HasAltLabelQV()) {
        if (read.HasAltLabelQV()) {
            const PacBio::BAM::QualityValues & qvs = read.AltLabelQV();
            std::vector<uint8_t> data(qvs.begin(), qvs.end());
            _CheckRead(read, data.size(), "AltLabelQV");
            altLabelQVArray_.Write(&data[0], data.size());
        } else {
            AddErrorMessage(std::string("AltLabelQV is absent in read " + read.FullName()));
        }
    }
    return Errors().empty();
}

void HDFPulseCallsWriter::Flush(void) {
    if (HasPulseCall())      pulseCallArray_.Flush();
    if (HasIsPulse())        isPulseArray_.Flush();
    if (HasLabelQV())        labelQVArray_.Flush();
    if (HasPkmean())         pkmeanArray_.Flush();
    if (HasPulseMergeQV())   pulseMergeQVArray_.Flush();
    if (HasPkmid())          pkmidArray_.Flush();
    if (HasStartFrame())     startFrameArray_.Flush();
    if (HasPulseCallWidth()) pulseCallWidthArray_.Flush();
    if (HasAltLabel())       altLabelArray_.Flush();
    if (HasAltLabelQV())     altLabelQVArray_.Flush();

    if (zmwWriter_)          zmwWriter_->Flush();
}

void HDFPulseCallsWriter::Close(void) {
    this->Flush();

    // Write attributes to pulsecallsGroup
    try { _WriteAttributes(); }
    catch (H5::Exception e) {
        AddErrorMessage("Failed to write attributes to " +
                        PacBio::GroupNames::pulsecalls);
    }

    if (HasPulseCall())      pulseCallArray_.Close();
    if (HasIsPulse())        isPulseArray_.Close();
    if (HasLabelQV())        labelQVArray_.Close();
    if (HasPkmean())         pkmeanArray_.Close();
    if (HasPulseMergeQV())   pulseMergeQVArray_.Close();
    if (HasPkmid())          pkmidArray_.Close();
    if (HasStartFrame())     startFrameArray_.Close();
    if (HasPulseCallWidth()) pulseCallWidthArray_.Close();
    if (HasAltLabel())       altLabelArray_.Close();
    if (HasAltLabelQV())     altLabelQVArray_.Close();
}

#endif

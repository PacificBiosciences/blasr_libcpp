#include  "libconfig.h"
#ifdef USE_PBBAM
#include "HDFBaseCallsWriter.hpp"

const std::vector<PacBio::BAM::BaseFeature> HDFBaseCallsWriter::ValidQVEnums = 
{     PacBio::BAM::BaseFeature::DELETION_QV
    , PacBio::BAM::BaseFeature::DELETION_TAG
    , PacBio::BAM::BaseFeature::INSERTION_QV
    , PacBio::BAM::BaseFeature::MERGE_QV
    , PacBio::BAM::BaseFeature::SUBSTITUTION_QV
    , PacBio::BAM::BaseFeature::SUBSTITUTION_TAG
    , PacBio::BAM::BaseFeature::IPD
    , PacBio::BAM::BaseFeature::PULSE_WIDTH
    , PacBio::BAM::BaseFeature::PULSE_CALL
};

std::vector<PacBio::BAM::BaseFeature> 
HDFBaseCallsWriter::WritableQVs(const std::vector<PacBio::BAM::BaseFeature> & qvsToWrite) {
    std::vector<PacBio::BAM::BaseFeature> ret;
    for(auto qv : qvsToWrite) {
        // Filter qvs which are not in format specification.
        if (std::find(ValidQVEnums.begin(), ValidQVEnums.end(), qv) != ValidQVEnums.end()) {
            if (std::find(ret.begin(), ret.end(), qv) == ret.end()) ret.push_back(qv); 
        }
    }
    return ret;
}

HDFBaseCallsWriter::HDFBaseCallsWriter(const std::string & filename,
                                       HDFGroup & parentGroup,
                                       const std::map<char, size_t> & baseMap,
                                       const std::string & basecallerVersion,
                                       const std::vector<PacBio::BAM::BaseFeature> & qvsToWrite)
    : HDFWriterBase(filename)
    , parentGroup_(parentGroup)
    , baseMap_(baseMap)
    , qvsToWrite_({})
    , basecallerVersion_(basecallerVersion)
    , arrayLength_(0)
    , zmwWriter_(nullptr)
    , zmwMetricsWriter_(nullptr)
{
    // Add BaseCalls as a child group to the parent group.
    AddChildGroup(parentGroup_, basecallsGroup_, PacBio::GroupNames::basecalls);

    if (basecallerVersion.empty()) {
        AddErrorMessage("BaseCallerVersion must not be empty!");
        return;
    }

    // Initialize the 'basecall' group.
    basecallArray_.Initialize(basecallsGroup_, PacBio::GroupNames::basecall);

    qvsToWrite_ = HDFBaseCallsWriter::WritableQVs(qvsToWrite);

    // It's ok to only write basecall.
    if (qvsToWrite_.size() > 0 and not InitializeQVGroups()) {
        AddErrorMessage("Failed to initialize QV Groups.");
        return;
    }

    // Create a zmwWriter.
    zmwWriter_.reset(new HDFZMWWriter(filename, basecallsGroup_));

    // Create a zmwMetricsWriter.
    zmwMetricsWriter_.reset(new HDFZMWMetricsWriter(filename, basecallsGroup_, baseMap));
}

std::vector<std::string> HDFBaseCallsWriter::Errors(void) const {
    std::vector<std::string> retErrors = this->errors_;
    if (zmwWriter_) {
        const std::vector<std::string> & zmwErrors = zmwWriter_->Errors();
        retErrors.insert(retErrors.end(), zmwErrors.begin(), zmwErrors.end());
    }
    if (zmwMetricsWriter_) {
        const std::vector<std::string> & zmwMetricsErrors = zmwMetricsWriter_->Errors();
        retErrors.insert(retErrors.end(), zmwMetricsErrors.begin(), zmwMetricsErrors.end());
    }
    return retErrors;
}

HDFBaseCallsWriter::~HDFBaseCallsWriter(void) {
    this->Close();
}

bool HDFBaseCallsWriter::InitializeQVGroups(void) {
    int ret = 1;

    // normal datasets
    if (_HasQV(PacBio::BAM::BaseFeature::DELETION_QV)) 
        ret *= deletionQVArray_.Initialize(basecallsGroup_,      PacBio::GroupNames::deletionqv);
    if (_HasQV(PacBio::BAM::BaseFeature::DELETION_TAG))
        ret *= deletionTagArray_.Initialize(basecallsGroup_,     PacBio::GroupNames::deletiontag);
    if (_HasQV(PacBio::BAM::BaseFeature::INSERTION_QV))
        ret *= insertionQVArray_.Initialize(basecallsGroup_,     PacBio::GroupNames::insertionqv);
    if (_HasQV(PacBio::BAM::BaseFeature::MERGE_QV))
        ret *= mergeQVArray_.Initialize(basecallsGroup_,         PacBio::GroupNames::mergeqv);
    if (_HasQV(PacBio::BAM::BaseFeature::SUBSTITUTION_QV))
        ret *= substitutionQVArray_.Initialize(basecallsGroup_,  PacBio::GroupNames::substitutionqv);
    if (_HasQV(PacBio::BAM::BaseFeature::SUBSTITUTION_TAG))
        ret *= substitutionTagArray_.Initialize(basecallsGroup_, PacBio::GroupNames::substitutiontag);
    if (_HasQV(PacBio::BAM::BaseFeature::IPD))
        ret *= ipdArray_.Initialize(basecallsGroup_,             PacBio::GroupNames::prebaseframes);
    if (_HasQV(PacBio::BAM::BaseFeature::PULSE_WIDTH))
        ret *= pulseWidthArray_.Initialize(basecallsGroup_,      PacBio::GroupNames::widthinframes);
    if (_HasQV(PacBio::BAM::BaseFeature::PULSE_CALL))
        ret *= pulseIndexArray_.Initialize(basecallsGroup_,      PacBio::GroupNames::pulseindex);
    return (ret != 0);
}

bool HDFBaseCallsWriter::_WriteAttributes(void) {
    // SchemaRevision
    bool OK = 
        AddAttribute<std::string>(basecallsGroup_,
                     PacBio::AttributeNames::Common::schemarevision,
                     PacBio::AttributeValues::Common::schemarevision) and
    // ChangeListID
        AddAttribute<std::string>(basecallsGroup_,
                     PacBio::AttributeNames::Common::changelistid,
                     basecallerVersion_);
    return OK;
}

bool HDFBaseCallsWriter::WriteOneZmw(const SMRTSequence & read) {
    bool OK = true;
    if (zmwWriter_)
        OK = OK and zmwWriter_->WriteOneZmw(read);
    if (zmwMetricsWriter_)
        OK = OK and zmwMetricsWriter_->WriteOneZmw(read);

    OK = OK and _WriteBasecall(read);

    OK = OK and _WriteDeletionQV(read);
    OK = OK and _WriteDeletionTag(read);
    OK = OK and _WriteInsertionQV(read);
    OK = OK and _WriteMergeQV(read);
    OK = OK and _WriteSubstitutionTag(read);
    OK = OK and _WriteSubstitutionQV(read);
    OK = OK and _WriteIPD(read);
    OK = OK and _WritePulseWidth(read);
    OK = OK and _WritePulseIndex(read);

    arrayLength_ +=  read.length;
    return OK;
}

bool HDFBaseCallsWriter::_WriteBasecall(const SMRTSequence & read) {
	basecallArray_.Write((const unsigned char*) read.seq, read.length);
    return true;
}

/*
bool HDFBaseCallsWriter::_WriteQualityValue(const SMRTSequence & read) {
    if (FakeQualityValue()) {
        if (read.length <= 0) {
            AddErrorMessage(read.GetTitle() + std::string(" is empty."));
            return false;
        }
        if (not read.deletionQV.Empty()) {
            // Use deletionQV to fake QualityValue if possible.
            qualityValueArray_.Write(read.deletionQV.data, read.length);
        } else { // otherwise, fill with 255.
            QualityValueVector<QualityValue> fakedata;
            fakedata.Allocate(read.length);
            memset(fakedata.data, MAX_QUALITY_VALUE, read.length * sizeof(QualityValue));
            qualityValueArray_.Write(fakedata.data, read.length);
            fakedata.Free();
        }
    }
    return true;
}
*/

bool HDFBaseCallsWriter::_WriteDeletionQV(const SMRTSequence & read) {
    if (HasDeletionQV()) {
        if (read.deletionQV.Empty()) {
            AddErrorMessage(std::string(PacBio::GroupNames::deletionqv) + " absent in read " + read.GetTitle());
            return false;
        } else {
            deletionQVArray_.Write(read.deletionQV.data, read.length);
            return true;
        }
    }
    return true;
}

bool HDFBaseCallsWriter::_WriteDeletionTag(const SMRTSequence & read) {
    if (HasDeletionTag()) {
        if (read.deletionTag == nullptr) {
            AddErrorMessage(std::string(PacBio::GroupNames::deletiontag) + " absent in read " + read.GetTitle());
            return false;
        } else {
            deletionTagArray_.Write(read.deletionTag, read.length);
            return true;
        }
    }
    return true;
}

bool HDFBaseCallsWriter::_WriteInsertionQV(const SMRTSequence & read) {
    if (HasInsertionQV()) {
        if (read.insertionQV.Empty()) {
            AddErrorMessage(std::string(PacBio::GroupNames::insertionqv) + " absent in read " + read.GetTitle());
            return false;
        } else {
			insertionQVArray_.Write(read.insertionQV.data, read.length);
            return true;
        }
    }
    return true;
}

bool HDFBaseCallsWriter::_WriteSubstitutionTag(const SMRTSequence & read) {
    if (HasSubstitutionTag()) {
        if (read.substitutionTag == nullptr) {
            AddErrorMessage(std::string(PacBio::GroupNames::substitutiontag) + " absent in read " + read.GetTitle());
            return false;
        } else {
			substitutionTagArray_.Write(read.substitutionTag, read.length);
            return true;
        }
    }
    return true;
}

bool HDFBaseCallsWriter::_WriteSubstitutionQV(const SMRTSequence & read) {
    if (HasSubstitutionQV()) {
        if (read.substitutionQV.Empty()) {
            AddErrorMessage(std::string(PacBio::GroupNames::substitutionqv) + " absent in read " + read.GetTitle());
            return false;
        } else {
			substitutionQVArray_.Write(read.substitutionQV.data, read.length);
            return true;
        }
    }
    return true;
}

bool HDFBaseCallsWriter::_WriteMergeQV(const SMRTSequence & read) {
    if (HasMergeQV()) {
        if (read.mergeQV.Empty()) {
            AddErrorMessage(std::string(PacBio::GroupNames::mergeqv) + " absent in read " + read.GetTitle());
            return false;
        } else {
			mergeQVArray_.Write(read.mergeQV.data, read.length);
            return true;
        }
    }
    return true;
}

bool HDFBaseCallsWriter::_WriteIPD(const SMRTSequence & read) {
     if (HasIPD()) {
        if (read.preBaseFrames == nullptr) {
            AddErrorMessage(std::string(PacBio::GroupNames::prebaseframes) + " absent in read " + read.GetTitle());
            return false;
        } else {
            ipdArray_.Write(read.preBaseFrames, read.length);
            return true;
        }
    }
    return true;
}

bool HDFBaseCallsWriter::_WritePulseWidth(const SMRTSequence & read) {
    if (HasPulseWidth()) {
        if (read.widthInFrames == nullptr) {
            AddErrorMessage(std::string(PacBio::GroupNames::widthinframes) + " absent in read " + read.GetTitle());
            return false;
        } else {
            pulseWidthArray_.Write(read.widthInFrames, read.length);
            return true;
        }
    }
    return true;
}

bool HDFBaseCallsWriter::_WritePulseIndex(const SMRTSequence & read) {
    if (HasPulseIndex()) {
        if (read.copiedFromBam) {
            const PacBio::BAM::BamRecord & record = read.bamRecord;
            if (record.HasPulseCall()) {
                const std::string & pulsecall = record.PulseCall();
                std::vector<uint16_t> data(read.length, 0);
                size_t indx = 0;
                for(size_t i = 0; i < pulsecall.size(); i++) {
                    const char base = pulsecall[i];
                    if (base == 'A' or base == 'C' or base == 'G' or base == 'T') {
                        assert(indx < read.length);
                        data[indx] = static_cast<uint16_t>(i);
                        indx ++;
                    } else {
                        assert(base == 'a' or base == 'c' or base == 'g' or base == 't');
                    }
                }
                assert(indx == read.length);
                pulseIndexArray_.Write(&data[0], read.length);
                return true;
            }
        }
        AddErrorMessage(std::string(PacBio::GroupNames::pulseindex) + " absent in read " + read.GetTitle());
        return false;
    }
    return true;
}

bool HDFBaseCallsWriter::WriteFakeDataSets() {   
    // Fake QualityValue with 255
    uint32_t block_sz = 65536; // This is a data buffer.
    std::vector<uint8_t> buffer_uint16_5M_0(block_sz);
    std::fill(buffer_uint16_5M_0.begin(), buffer_uint16_5M_0.end(), 255);
    bool OK = __WriteFakeDataSet<uint8_t>(basecallsGroup_, PacBio::GroupNames::qualityvalue, arrayLength_, buffer_uint16_5M_0);
    return OK;
}
 
void HDFBaseCallsWriter::Flush(void) {
    basecallArray_.Flush();

    if (HasDeletionQV())      deletionQVArray_.Flush();
    if (HasDeletionTag())     deletionTagArray_.Flush();
    if (HasInsertionQV())     insertionQVArray_.Flush();
    if (HasMergeQV())         mergeQVArray_.Flush();
    if (HasSubstitutionQV())  substitutionQVArray_.Flush();
    if (HasSubstitutionTag()) substitutionTagArray_.Flush();
    if (HasIPD())             ipdArray_.Flush();
    if (HasPulseWidth())      pulseWidthArray_.Flush();
    if (HasPulseIndex())      pulseIndexArray_.Flush();

    if (zmwWriter_)           zmwWriter_->Flush();
    if (zmwMetricsWriter_)    zmwMetricsWriter_->Flush();
}

void HDFBaseCallsWriter::Close(void) {
    this->Flush();

    try { _WriteAttributes(); }
    catch (H5::Exception e) {
        AddErrorMessage("Failed to write attributes to " +
                        PacBio::GroupNames::basecalls);
    }

    basecallArray_.Close();

    if (HasDeletionQV())      deletionQVArray_.Close();
    if (HasDeletionTag())     deletionTagArray_.Close();
    if (HasInsertionQV())     insertionQVArray_.Close();
    if (HasMergeQV())         mergeQVArray_.Close();
    if (HasSubstitutionQV())  substitutionQVArray_.Close();
    if (HasSubstitutionTag()) substitutionTagArray_.Close();
    if (HasIPD())             ipdArray_.Close();
    if (HasPulseWidth())      pulseWidthArray_.Close();
    if (HasPulseIndex())      pulseIndexArray_.Close();
}
#endif

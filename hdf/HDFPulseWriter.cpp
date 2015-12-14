#include "libconfig.h"
#ifdef USE_PBBAM
#include "HDFPulseWriter.hpp"

HDFPulseWriter::HDFPulseWriter(const std::string & filename,
                           const ScanData & sd,
                           const std::string & basecallerVersion,
                           const std::vector<PacBio::BAM::BaseFeature> & qvsToWrite,
                           const H5::FileAccPropList & fileAccPropList)
    : HDFWriterBase(filename)
    , fileaccproplist_(fileAccPropList)
    , scandataWriter_(nullptr)
    , basecallsWriter_(nullptr) 
    , pulsecallsWriter_(nullptr)
    , regionsWriter_(nullptr)
{
    // sanity check chemistry meta data. 
    SanityCheckChemistry(sd.BindingKit(),
                         sd.SequencingKit(), 
                         basecallerVersion);

    // open file 
    outfile_.Open(filename_, H5F_ACC_TRUNC, fileaccproplist_);

    // Add PulseData group to the root group '/'
    AddChildGroup(outfile_.rootGroup, pulseDataGroup_, PacBio::GroupNames::pulsedata);

    // Create a ScanData writer.
    scandataWriter_.reset(new HDFScanDataWriter(outfile_.rootGroup)); 
    scandataWriter_->Write(sd);

    // Create a BaseCaller writer.
    basecallsWriter_.reset(new HDFBaseCallsWriter(filename_, pulseDataGroup_, sd.BaseMap(), basecallerVersion, qvsToWrite));

    // Create a PulseCalls writer
    pulsecallsWriter_.reset(new HDFPulseCallsWriter(filename_, pulseDataGroup_, sd.BaseMap(), basecallerVersion, qvsToWrite));
}

HDFPulseWriter::HDFPulseWriter(const std::string & filename,
                           const ScanData & sd,
                           const std::string & basecallerVersion,
                           const std::vector<PacBio::BAM::BaseFeature> & qvsToWrite,
                           const std::vector<std::string> & regionTypes,
                           const H5::FileAccPropList & fileAccPropList)
    :HDFPulseWriter(filename, sd, basecallerVersion, qvsToWrite, fileAccPropList)
{
    // Create a Regions writer.
    regionsWriter_.reset(new HDFRegionsWriter(filename_, pulseDataGroup_, regionTypes));
}
 

HDFPulseWriter::~HDFPulseWriter(void) {
    this->Close();
}

void HDFPulseWriter::Flush(void) {
    basecallsWriter_->Flush();
    pulsecallsWriter_->Flush();
    if (HasRegions()) regionsWriter_->Flush();
}

std::vector<std::string> HDFPulseWriter::Errors(void) {
    std::vector<std::string> errors = errors_;

    //for (auto error: scandataWriter_->Errors())
    //    errors.emplace_back(error);

    for (auto error: basecallsWriter_->Errors())
        errors.emplace_back(error);

    for (auto error: pulsecallsWriter_->Errors())
        errors.emplace_back(error);

    if (HasRegions()) {
        for (auto error: regionsWriter_->Errors())
            errors.emplace_back(error);
    }

    return errors;
}

void HDFPulseWriter::Close(void) {
    if (basecallsWriter_) basecallsWriter_.reset();
    if (pulsecallsWriter_) pulsecallsWriter_.reset();
    if (scandataWriter_) scandataWriter_.reset();
    if (HasRegions() and regionsWriter_) regionsWriter_.reset();
    outfile_.Close();
}

bool HDFPulseWriter::SanityCheckChemistry(const std::string & bindingKit,
                                        const std::string & sequencingKit,
                                        const std::string & basecallerVersion)
{
    bool OK = true;
    if (bindingKit.empty()) {
        OK = false;
        AddErrorMessage("Binding kit must be specified.");
    }
    if (sequencingKit.empty()) {
        OK = false;
        AddErrorMessage("Sequencing kit must be specified.");
    }
    if (basecallerVersion.empty()) {
        OK = false;
        AddErrorMessage("Base caller version must be specified.");
    }
    return OK;
}

bool HDFPulseWriter::WriteOneZmw(const SMRTSequence & seq) {
    bool OK = basecallsWriter_->WriteOneZmw(seq);
    OK = OK and pulsecallsWriter_->WriteOneZmw(seq);
    return OK;
}

bool HDFPulseWriter::WriteOneZmw(const SMRTSequence & seq, 
                                 const std::vector<RegionAnnotation> & regions) {
    if (not this->WriteOneZmw(seq)) {
        return false;
    }
    if (HasRegions()) {
        if (regions.size() == 0) {
            std::vector<RegionAnnotation> fake = {RegionAnnotation(seq.HoleNumber(), HQRegion, 0, 0, 0)};
            return regionsWriter_->Write(fake);
        } else {
            return regionsWriter_->Write(regions);
        }
    }
    return true;
} 

bool HDFPulseWriter::WriteFakeDataSets() {
    return pulsecallsWriter_->WriteFakeDataSets();
}
#endif // ifndef USE_PBBAM

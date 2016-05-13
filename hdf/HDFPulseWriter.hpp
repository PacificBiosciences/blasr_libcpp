#ifndef _BLASR_HDF_PULSE_WRITER_HPP_
#define _BLASR_HDF_PULSE_WRITER_HPP_
#include "../pbdata/libconfig.h"
#ifdef USE_PBBAM

#include <sstream>
#include <memory>
// pbdata/
#include "../pbdata/Enumerations.h"
#include "../pbdata/SMRTSequence.hpp"

#include "HDFWriterBase.hpp"
#include "HDFScanDataWriter.hpp"
#include "HDFBaseCallsWriter.hpp"
#include "HDFPulseCallsWriter.hpp"
#include "HDFRegionsWriter.hpp"

using namespace H5;
using namespace std;

class HDFPulseWriter : public HDFWriterBase {
public:

    /// \name Constructor & Related Methods
    /// \{
    /// \brief Sets output h5 file name, scan data, base caller version
    ///        QVs to write, and h5 file access property list.
    /// \note Set /PulseData/Regions
    /// \param[in] filename output h5 file name.
    /// \param[in] basecallerVersion meta data string
    /// \param[in] qvsToWrite Quality values to include in output h5 file. 
    /// \param[in] regionTypes, regionTypes as /Regions/RegionTypes
    /// \param[in] fileAccPropList H5 file access property list
    HDFPulseWriter(const std::string & filename,
                   const std::string & basecallerVersion,
                   const std::map<char, size_t>& baseMap,
                   const std::vector<PacBio::BAM::BaseFeature> & qvsToWrite,
                   const std::vector<std::string> & regionTypes,
                   const H5::FileAccPropList & fileAccPropList = H5::FileAccPropList::DEFAULT);

    /// \note No /PulseData/Regions
    HDFPulseWriter(const std::string & filename,
                   const std::string & basecallerVersion,
                   const std::map<char, size_t>& baseMap,
                   const std::vector<PacBio::BAM::BaseFeature> & qvsToWrite,
                   const H5::FileAccPropList & fileAccPropList = H5::FileAccPropList::DEFAULT);

	~HDFPulseWriter(void);

    /// \brief Write one zmw sequence to output h5 file. 
    /// \param[in] seq, the SMRTSequence to write
    bool WriteOneZmw(const SMRTSequence & seq);

    /// \brief Write one zmw sequence and its region table to output h5 file. 
    /// \param[in] seq, the SMRTSequence to write
    /// \param[in] regions, region annotations of this zmw.
    bool WriteOneZmw(const SMRTSequence & seq, 
                     const std::vector<RegionAnnotation> & regions);

    /// \brief Write fake datasets to create POC/PulseRecognizer compatible plx file.
    bool WriteFakeDataSets();

    inline bool HasRegions(void) const;

    /// \brief Flushes buffered data.
    void Flush(void);

    /// \returns all errors from all writers.
    std::vector<std::string> Errors(void);

    /// Pass-through method for setting the inverse gain to PulseCallsWriter
    void SetInverseGain(float igain);

    /// \}

private:
    H5::FileAccPropList fileaccproplist_; ///< H5 file access property list
	HDFGroup pulseDataGroup_; ///< /PulseData group

private:
    /// Points to base caller writer.
    std::unique_ptr<HDFBaseCallsWriter> basecallsWriter_;
    /// Points to pulse calls writer.
    std::unique_ptr<HDFPulseCallsWriter> pulsecallsWriter_;
    /// Points to region table writer.
    std::unique_ptr<HDFRegionsWriter>   regionsWriter_;
    /// \}

private:
    /// \name Private Methods.
    /// \{

    /// \brief Closes HDFPulseWriter.
    void Close(void);
    /// \}
};

inline bool HDFPulseWriter::HasRegions(void) const
{ return bool(regionsWriter_); }

#endif  // USE_PBBAM
#endif

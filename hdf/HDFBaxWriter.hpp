#ifndef _BLASR_HDF_BAX_WRITER_HPP_
#define _BLASR_HDF_BAX_WRITER_HPP_
#include "../pbdata/libconfig.h"
#ifdef USE_PBBAM

#include <sstream>
#include <memory>
// pbdata/
#include "../pbdata/Enumerations.h"
#include "../pbdata/SMRTSequence.hpp"

#include "HDFFile.hpp"
#include "HDFWriterBase.hpp"
#include "HDFScanDataWriter.hpp"
#include "HDFBaseCallsWriter.hpp"
#include "HDFRegionsWriter.hpp"

using namespace H5;
using namespace std;

class HDFBaxWriter : public HDFWriterBase {
public:
    /// \name Constructor & Related Methods
    /// \{
    /// \brief Sets output h5 file name, scan data, base caller version
    ///        QVs to write, and h5 file access property list.
    /// \note Should not create /PulseData/Regions.
    /// \param[in] filename output h5 file name.
    /// \param[in] basecallerVersion meta data string
    /// \param[in] qvsToWrite Quality values to include in output h5 file. 
    /// \param[in] fileAccPropList H5 file access property list
    HDFBaxWriter(const std::string & filename,
                 const std::string & basecallerVersion,
                 const std::map<char, size_t>& baseMap,
                 const std::vector<PacBio::BAM::BaseFeature> & qvsToWrite,
                 const H5::FileAccPropList & fileAccPropList = H5::FileAccPropList::DEFAULT);

    /// \brief Sets output h5 file name, scan data, base caller version
    ///        QVs to write, regions types and h5 file access property list.
    /// \note /PulseData/Regions dataset should be created.
    /// \param[in] regionTypes, regionTypes as /Regions/RegionTypes
    HDFBaxWriter(const std::string & filename,
                 const std::string & basecallerVersion,
                 const std::map<char, size_t>& baseMap,
                 const std::vector<PacBio::BAM::BaseFeature> & qvsToWrite,
                 const std::vector<std::string> & regionTypes,
                 const H5::FileAccPropList & fileAccPropList = H5::FileAccPropList::DEFAULT);
    
	~HDFBaxWriter(void);

    /// \returns whether or not write region table.
    inline bool HasRegions(void) const;

    /// \brief Write one zmw sequence to output h5 file. 
    /// \param[in] seq, the SMRTSequence to write
    bool WriteOneZmw(const SMRTSequence & seq);

    /// \brief Write one zmw sequence and its region table to output h5 file. 
    /// \param[in] seq, the SMRTSequence to write
    /// \param[in] regions, region annotations of this zmw.
    bool WriteOneZmw(const SMRTSequence & seq, 
                     const std::vector<RegionAnnotation> & regions);

    /// \brief Write fake datasets for POC compatible bax file.
    bool WriteFakeDataSets();

    /// \brief Flushes buffered data.
    void Flush(void);

    /// \returns all errors from all writers.
    std::vector<std::string> Errors(void);

    /// \}

private:
    H5::FileAccPropList fileaccproplist_; ///< H5 file access property list
	HDFGroup pulseDataGroup_; ///< /PulseData group

private:
    /// Points to scan data writer.
    std::unique_ptr<HDFScanDataWriter>  scandataWriter_;
    /// Points to base caller writer.
    std::unique_ptr<HDFBaseCallsWriter> basecallsWriter_;
    /// Points to region table writer.
    std::unique_ptr<HDFRegionsWriter>   regionsWriter_;
    /// \}

private:
    /// \name Private Methods.
    /// \{
    /// \brief Closes HDFBaxWriter.
    void Close(void);
    /// \}
};

inline bool HDFBaxWriter::HasRegions(void) const 
{ return bool(regionsWriter_); }

#endif // end of #ifdef USE_PBBAM

#endif // end of #ifdef _BLASR_HDF_BAX_WRITER_HPP_

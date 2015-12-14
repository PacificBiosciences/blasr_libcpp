#ifndef _BLASR_HDF_HDFZMWWriter_HPP_
#define _BLASR_HDF_HDFZMWWriter_HPP_

#include "libconfig.h"
#ifdef USE_PBBAM

#include "HDFWriterBase.hpp"
#include "BufferedHDFArray.hpp"
#include "BufferedHDF2DArray.hpp"
#include "SMRTSequence.hpp"
#include <pbbam/BamRecord.h>

class HDFBaseCallerWriter;
class HDFPulseCallerWriter;

class HDFZMWWriter: public HDFWriterBase {

friend class HDFBaseCallerWriter;
friend class HDFPulseCallerWriter;

public:
    /// \name Constructors and Destructors
    /// \{
    HDFZMWWriter(const std::string & filename, 
                 HDFGroup & parentGroup);

    /// \params[in] filename
    /// \params[in] parentGroup
    /// \params[in] inPulseCalls, true if this ZMW is within PulseCalls. 
    /// \params[in] baseMap, base to channel index in H5.
    HDFZMWWriter(const std::string & filename, 
                 HDFGroup & parentGroup,
                 const bool inPulseCalls,
                 const std::map<char, size_t> & baseMap);

    ~HDFZMWWriter(void);
    /// \}

    /// \name Public Methods
    /// \{
    /// \returns true if this ZMW is within PulseCalls.
    bool InPulseCalls(void) const;

    /// \returns Whether or not to write BaseLineSigma
    bool HasBaseLineSigma(void) const;
    
    /// \note Write info of a SMRTSequence to ZMW,
    ///       (1) add number of pulses to NumEvent if InPulseCalls();
    ///           otherwise, add number of bases
    ///       (2) add zmw hole number (UInt) of the sequence as a UInt to HoleNumber,
    ///       (3) add hole status (unsigned char) to HoleStatus,
    ///       (4) add hole coordinate xy as (int16_t, int16_t) to HoleXY
    bool WriteOneZmw(const PacBio::BAM::BamRecord & read);

    bool WriteOneZmw(const SMRTSequence & read);

    /// \brief Write fake datasets under /PulseCalls/ZMW
    bool WriteFakeDataSets(void);

    uint32_t NumZMWs(void) const {return arrayLength_;}

    /// \note Flushes all data from cache to disc.
    void Flush(void);

    /// \note Closes this zmw group as well as child hdf groups.
    void Close(void);
 
    /// \}

private:
    /// \name Private Data 
    /// ZMW/NumEvent
    BufferedHDFArray<uint32_t>      numEventArray_;
    /// ZMW/HoleNumber
    BufferedHDFArray<unsigned int>  holeNumberArray_;
    /// ZMW/HoleStatus
    BufferedHDFArray<unsigned char> holeStatusArray_;
    /// ZMW/HoleXY
    BufferedHDF2DArray<int16_t>     holeXYArray_;
    /// ZMW/BaseLineSigma
    BufferedHDF2DArray<float>       baseLineSigmaArray_;

    /// ZMW group
    HDFGroup zmwGroup_;
    /// Parent group PulseCalls or BaseCalls
    HDFGroup & parentGroup_;

    /// Map 'ACGT' to channel indices, defined in /ScanData/RunInfo/BaseMap
    std::map<char, size_t> baseMap_;
    /// true if parent gropu is PulseCalls.
    bool inPulseCalls_;
    /// dataset (array) size.
    uint32_t arrayLength_;
    /// \}

private:
    /// \name Private Methods
    /// \{
    
    /// \note Initialize child hdf groups under ZMW, including
    ///       NumEvent, HoleNumber, HoleStatus, HoleXY, and BaseLineSigma
    /// \reutrns bool, whether or not child hdf groups successfully initialized.
    bool InitializeChildHDFGroups(void);

    /// \note Write number of bases to ZMW/NumEvent if not InPulseCalls();
    ///       Otherwise, write number of pulses.
    bool _WriteNumEvent(const uint32_t numEvent);

    /// \note Write HoleNumber
    bool _WriteHoleNumber(const uint32_t holeNumber);

    /// \note Write HoleXY
    bool _WriteHoleXY(const int16_t holeX, const int16_t holeY);

    /// \note Always write 'SEQUENCINGZMW' to HoleStatus
    bool _WriteHoleStatus(const unsigned char holeStatus);

    /// \note Write BaseLineSigma if it is required to write and exists in read
    bool _WriteBaseLineSigma(const PacBio::BAM::BamRecord & read);

    /// \note Add attributes to HoleNumber, HoleXY, HoleStatus, HoleXY and BaseLineSigma
    void _WriteAttributes(void);

    /// \}
};

#endif // end of #ifdef USE_PBBAM

#endif // end of #ifndef _BLASR_HDF_HDFZMWWriter_HPP_

#ifndef _BLASR_PULSE_FILE_HPP_
#define _BLASR_PULSE_FILE_HPP_

#include <map>
#include <vector>
#include "Types.h"
#include "Enumerations.h"
#include "DNASequence.hpp"
#include "SMRTSequence.hpp"
#include "PulseBaseCommon.hpp"
#include "ScanData.hpp"

class PulseFile : public PulseBaseCommon {
 public:
    unsigned int numFrames;
    PlatformId platformId;
    std::vector<unsigned int> startFrame;
    std::vector<uint16_t>  plsWidthInFrames;
    int midSignalNDims, maxSignalNDims, meanSignalNDims;
    std::vector<uint16_t>  midSignal;
    std::vector<uint16_t>  maxSignal;
    std::vector<uint16_t>  meanSignal; 
    std::vector<DNALength> numEvent;
    std::vector<DSLength>  pulseStartPositions;
    std::vector<float>     classifierQV;

    PulseFile(){numFrames = 0; platformId = Springfield;}

    // Copy signals (QVs) from signalData (which has dimension signalNDims)
    // starting from position pulseStartPos with length readLength to readData.
    // If signal diamension is 2 (2D Array e.g., Chi2), then bases in readSeq at a 
    // particular location is used to determine which channel to copy the signal from.
    // e.g., if basemap is 'GCTA', and base at position 3 is T, each base corresponds 
    // to 4 channels and then we need to copy signal from channel 2 of position 3.
    void CopySignal(HalfWord *signalData,   // either a vector or matrix
                    int signalNDims,        // diamension of signal, either 1 or 2.
                    DSLength pulseStartPos, // 0 if baseToPulseIndex maps to abs position
                    int *baseToPulseIndex,  // Map base to pulse index.
                    Nucleotide *readSeq,
                    DNALength readLength,   // read length (num of signals to copy)
                    HalfWord *readData);
    
    // plsReadIndex: index of this hole number in /PulseCalls/ZMW/HoleNumber.
    // baseToPulseIndex: index from pulse to base from the beginning of the read.
    // read: a SMRTSequence.
    void CopyReadAt(UInt plsReadIndex, int *baseToPulseIndex, SMRTSequence &read);

    template<typename T_FieldType>
    void StoreField(std::vector<T_FieldType> &source, int *basToPlsIndex, T_FieldType *dest, const DSLength destLength);
    
    template <typename T>
    bool Realloc(T *&ptr, const DNALength length);
};

#include "PulseFileImpl.hpp"

#endif

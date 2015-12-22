#include "MovieAlnIndexLookupTable.hpp"


MovieAlnIndexLookupTable::MovieAlnIndexLookupTable() {
    skip = true;
    //skip this alignment unless we can find all the information
}

void MovieAlnIndexLookupTable::SetValue(
                const bool & skipP,            
                const size_t  & movieAlignmentIndexP,
                const UInt & alignmentIndexP,  
                const size_t  & refGroupIndexP,   
                const size_t  & readGroupIndexP,  
                const UInt & holeNumberP,
                const UInt & offsetBeginP,     
                const UInt & offsetEndP,
                const UInt & queryStartP,      
                const UInt & queryEndP,
                const size_t  & readIndexP,
                const size_t  & readStartP,       
                const int  & readLengthP,
                const size_t  & plsReadIndexP) {
    skip = skipP; 
    movieAlignmentIndex = movieAlignmentIndexP;
    alignmentIndex      = alignmentIndexP;
    refGroupIndex       = refGroupIndexP;
    readGroupIndex      = readGroupIndexP;
    holeNumber          = holeNumberP;
    offsetBegin         = offsetBeginP;
    offsetEnd           = offsetEndP;
    queryStart          = queryStartP;
    queryEnd            = queryEndP;
    readIndex           = readIndexP;
    readStart           = readStartP;
    readLength          = readLengthP;
    plsReadIndex        = plsReadIndexP;
}


void MovieAlnIndexLookupTable::print() {
    // Print this lookup table for debugging . 
    if (skip) 
        std::cout << "skip = True, ";
    else 
        std::cout << "skip = False, ";
    std::cout << "movieAlnIndex    = " << alignmentIndex << ", refGroupIndex = " << refGroupIndex
            << ", readGroupIndex = " << readGroupIndex << ", holeNumber    = " << holeNumber
            << ", offsetBegin    = " << offsetBegin    << ", offsetEnd     = " << offsetEnd
            << ", queryStart     = " << queryStart     << ", queryEnd      = " << queryEnd
            << ", readIndex      = " << readIndex      << ", readStart     = " << readStart
            << ", readLength     = " << readLength     << ", plsReadIndex  = " << plsReadIndex
            << std::endl;
}

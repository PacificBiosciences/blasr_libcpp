/*
 * =====================================================================================
 *
 *       Filename:  MovieAlnIndexLookupTable.h
 *
 *    Description:  For cpp/pbihdfutils/LoadPulses.cpp use only.
 *                  Cache all the information required for loading metric from pls/bas
 *                  to cmp file.
 *                  mapping movieAlignmentIndex in movieIndexSets[movieIndex] to:
 *                  alignmentIndex: index of this alignment in cmpFile.alnInfo.alignments
 *                  refGroupIndex : index of this alignment in cmpReader.refAlignGroups
 *                  readGroupIndex: index of this alignment in cmpReader.refAlignGroups\
 *                                  [refGroupIndex]->readGroupds[readGroupIndex]
 *                  offsetBegin   : offset begin for this alignment in cmpFile 
 *                                  dataset /ref/movie/AlnArray
 *                  offsetEnd     : offset enda
 *                  queryStart    : start position of query of this alignment 
 *                  queryEnd      : end position of query of this alignment 
 *                  readIndex     : index of this alignment in baseFile.readStartPositions
 *                  readStart     : start position of this alignment in baseFile
 *                  readLength    : read length of this alignment in baseFile
 *                  plsReadIndex  : index of this alignment in pulseFile.pulseStartPositions
 *                  alignedSequence
 *                                : aligned sequence of this alignment
 *                  
 *        Version:  1.0
 *        Created:  12/19/2012 03:50:21 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#ifndef _LOADPULSES_MovieAlnIndexLookupTable_HPP_
#define _LOADPULSES_MovieAlnIndexLookupTable_HPP_
#include <iostream>
#include <string>

#include "../Types.h"

class MovieAlnIndexLookupTable {
public: 
    bool  skip;
    // as movies may split into multiple parts, skip=true if 
    // this alignment is not found in the movie
     
    size_t   movieAlignmentIndex;
    // movieIndexSets[movieIndex][toFrom[movieAlignmentIndex]] 
    
    UInt  alignmentIndex;
    // cmpFile.alnInfo.alignments[alignmentIndex]
    
    UInt  holeNumber;
    // holeNumber corresponding to this alignment in baseFile
    
    size_t   refGroupIndex;
    // cmpReader.refAlignGroups[refGroupIndex]
    
    size_t   readGroupIndex;
    // cmpReader.refAlignGroups[refGroupIndex]->readGroups[readGroupIndex]
    
    UInt  offsetBegin, offsetEnd; 
    // offset begin and end for this alignment in /ref/movie/AlnArray
    // = cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetBegin/End()
    
    UInt   queryStart,  queryEnd;
    // start/end for this read = cmpFile.alnInfo.
    // alignments[alignmentIndex].GetQueryStart/End()

    size_t   readIndex;
    // index of this alignment in baseFile.readStartPositions
    // = index of this hole number in BaseCalls/ZMW/HoleNumber
    // baseFile.LookupReadIndexByHoleNumber(holeNumber, out=readIndex)
    
    size_t   readStart;
    // start pos of this alignment in baseFile
    // = baseFile.readStartPositions[readIndex]

    int   readLength;
    // read length of this alignment in baseFile
    // = baseFile.readStartPositions[readIndex+1] - readStart 

    size_t   plsReadIndex;
    // index of this alignment in pulseFile.pulseStartPositions
    // = index of this hole number in PulseCalls/ZMW/HoleNumbers
    // = pulseFile.LookupReadIndexByHoleNumber(holeNumber, out=plsReadIndex)

    // vector<int> baseToAlignmentMap; 
    // keep all the baseToAlignmentMap in memory for now
    // Note that baseToAlignmentMap is not initialized when 
    // BuildLookupTable is called.
   
    std::string alignedSequence;

    MovieAlnIndexLookupTable();

    void SetValue(const bool & skipP,
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
                  const size_t  & plsReadIndexP);

    void print();
};
#endif

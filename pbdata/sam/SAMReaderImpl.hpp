#ifndef _BLASR_SAM_READER_IMPL_HPP_
#define _BLASR_SAM_READER_IMPL_HPP_
#include <iostream>
#include "SAMReader.hpp"

template<typename T_ReferenceSequence, typename T_ReadGroup, typename T_SAMAlignment>
bool SAMReader<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment>::Initialize(std::string samFileName) {
  if(samFileName != "stdin") {
    CrucialOpen(samFileName, samFile, std::ios::in);
    samFilePtr = &samFile;
  } else {
    samFilePtr = &std::cin;
  }
  return true;
}

template<typename T_ReferenceSequence, typename T_ReadGroup, typename T_SAMAlignment>
void SAMReader<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment>::Close() {
  if(samFile.is_open()) {
    samFile.close();
  }
}

template<typename T_ReferenceSequence, typename T_ReadGroup, typename T_SAMAlignment>
void SAMReader<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment>::GetLine(std::istream &in, std::string &line) {
  getline(in, line);
}

template<typename T_ReferenceSequence, typename T_ReadGroup, typename T_SAMAlignment>
bool SAMReader<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment>::LineTypeIsHeader(LineType lineType) {
  return (lineType == HSHeader or
          lineType == HSSequence or
          lineType == HSReadGroup or 
          lineType == HSProgram or 
          lineType == HSComment);
}

template<typename T_ReferenceSequence, typename T_ReadGroup, typename T_SAMAlignment>
bool SAMReader<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment>::PeekLineIsHeader(std::istream &in) {
  if (in and in.peek() == '@') {
    return true;
  }
  else {
    return false;
  }
}

template<typename T_ReferenceSequence, typename T_ReadGroup, typename T_SAMAlignment>
typename SAMReader<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment>::LineType
SAMReader<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment>::GetLineType(std::string &line) {
  if (line.length() == 0) {
    return Blank;
  }
  else if (line[0] == '@') {
    std::stringstream strm(line);
    std::string tag;
    strm >> tag;
    if (tag == "@HD") { return HSHeader; }
    else if (tag == "@SQ") { return HSSequence; }
    else if (tag == "@RG") { return HSReadGroup; }
    else if (tag == "@PG") { return HSProgram; }
    else if (tag == "@CO") { return HSComment; }
    else { return Error; }
  }
  else {
    return Alignment;
  }
}

template<typename T_ReferenceSequence, typename T_ReadGroup, typename T_SAMAlignment>
void SAMReader<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment>::StoreKVPairs(std::string line,
                             std::vector<SAMKeywordValuePair> &kvPairs) {
  //
  // Split on tab delineated line.
  //
  std::vector<std::string> kvPairStrings;
  Splice(line, "\t", kvPairStrings);
  KeywordValueStringsToPairs(kvPairStrings, kvPairs);
}

template<typename T_ReferenceSequence, typename T_ReadGroup, typename T_SAMAlignment>
void SAMReader<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment>::StoreHeader(std::vector<SAMKeywordValuePair> &kvPairs,
                           AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments) {
  alignments.header.StoreValues(kvPairs, lineNumber);
}

template<typename T_ReferenceSequence, typename T_ReadGroup, typename T_SAMAlignment>
void SAMReader<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment>::StoreReferenceSequence(std::vector<SAMKeywordValuePair> &kvPairs,
                                       AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments) {
  alignments.references.push_back(T_ReferenceSequence());
  int lastRefIndex = alignments.references.size() - 1;
  alignments.references[lastRefIndex].StoreValues(kvPairs, lineNumber);
}

template<typename T_ReferenceSequence, typename T_ReadGroup, typename T_SAMAlignment>
void SAMReader<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment>::StoreReadGroup(std::vector<SAMKeywordValuePair> &kvPairs,
                               AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments) {
  alignments.readGroups.push_back(T_ReadGroup());
  int lastReadGroupIndex = alignments.readGroups.size() - 1;
  alignments.readGroups[lastReadGroupIndex].StoreValues(kvPairs, lineNumber);
}

template<typename T_ReferenceSequence, typename T_ReadGroup, typename T_SAMAlignment>
void SAMReader<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment>::StoreAlignment(std::string & line,
                               AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments) {
  alignments.alignments.push_back(T_SAMAlignment());   
  int lastAlignmentIndex = alignments.alignments.size() - 1;
  alignments.alignments[lastAlignmentIndex].StoreValues(line, lineNumber);
}

template<typename T_ReferenceSequence, typename T_ReadGroup, typename T_SAMAlignment>
void SAMReader<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment>::Read(std::string samFileName, AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments) {
  Initialize(samFileName);
  Read(alignments);
}

template<typename T_ReferenceSequence, typename T_ReadGroup, typename T_SAMAlignment>
std::vector<std::string> SAMReader<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment>::ReadHeader(AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments) {
  std::vector<std::string> allHeaders;
  std::string line;
  LineType lineType;
  lineNumber = 0;
  while (*samFilePtr and PeekLineIsHeader(*samFilePtr)) {
    getline(*samFilePtr, line);
    lineType = GetLineType(line);
    if (LineTypeIsHeader(lineType)) {
      allHeaders.push_back(line);
      std::stringstream strm(line);
      std::string tag;
      strm >> tag;
      std::string remainder;
      getline(strm, remainder);
      std::vector<SAMKeywordValuePair> kvPairs;
      StoreKVPairs(remainder, kvPairs);
      if (lineType == HSHeader) {
        StoreHeader(kvPairs, alignments);
      }
      else if (lineType == HSSequence) {
        StoreReferenceSequence(kvPairs, alignments);
      }
      else if (lineType == HSReadGroup) {
        StoreReadGroup(kvPairs, alignments);
      }
      else if (lineType == HSProgram) {
        StoreProgram(kvPairs, alignments);
      }
      else if (lineType == HSComment) {
        // do nothing with comments for now.
      }
    }
    ++lineNumber;
  }
  return allHeaders;
}

template<typename T_ReferenceSequence, typename T_ReadGroup, typename T_SAMAlignment>
void SAMReader<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment>::Read(AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments) {
  std::string line;
  LineType lineType;
  lineNumber = 0;
  ReadHeader(alignments);
  while (getline(*samFilePtr, line)) {
    lineType = GetLineType(line);
    if (LineTypeIsHeader(lineType)) {
      std::cout << "ERROR! Header line found outside of the header at " << lineNumber << std::endl;
      exit(1);
    }
    else if (lineType == Alignment) {
      StoreAlignment(line, alignments);
    }
    else {
      std::cout << "Error, line type unknown at " << lineNumber << std::endl;
      std::cout << line << std::endl;
      exit(1);
    }
    ++lineNumber;
  }
}

template<typename T_ReferenceSequence, typename T_ReadGroup, typename T_SAMAlignment>
bool SAMReader<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment>::GetNextAlignment(SAMAlignment &alignment) {
  if (*samFilePtr) {
    std::string line;
    if (getline(*samFilePtr, line)) {
      alignment.StoreValues(line, lineNumber);
      ++lineNumber;
      return true;
    }
    else {
      return false;
    }
  }
  else {
    return false;
  }
}

#endif

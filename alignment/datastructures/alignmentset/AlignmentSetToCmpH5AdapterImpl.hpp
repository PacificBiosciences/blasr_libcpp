#ifndef _BLASR_ALIGNMENT_SET_TO_CMPH5_ADAPTER_IMPL_HPP_
#define _BLASR_ALIGNMENT_SET_TO_CMPH5_ADAPTER_IMPL_HPP_

#include "AlignmentSetToCmpH5AdapterImpl.hpp"

template<typename T_CmpFile>
unsigned int AlignmentSetToCmpH5Adapter<T_CmpFile>::StoreMovieInfo(
        std::string movieName, T_CmpFile &cmpFile) {

  std::map<std::string, unsigned int>::iterator mapIt;
  mapIt = knownMovies.find(movieName);
  if (mapIt != knownMovies.end()) {
      return mapIt->second;
  }
  else {
      unsigned int id = cmpFile.movieInfoGroup.AddMovie(movieName);
      knownMovies[movieName] = id;
      return id;
  }
}

template<typename T_CmpFile>
template<typename T_Reference>
void AlignmentSetToCmpH5Adapter<T_CmpFile>::StoreReferenceInfo(
    std::vector<T_Reference> &references, T_CmpFile &cmpFile) {
  for (int r = 0; r < references.size(); r++) {
    std::string sequenceName, md5;
    sequenceName = references[r].GetSequenceName();
    md5 = references[r].GetMD5();
    unsigned int length = references[r].GetLength();

    // Add this reference to /RefInfo. 
    // Don't create /ref0000x and register it in /RefGroup at this point, 
    // because this reference may not map to any alignments at all. 
    cmpFile.AddRefInfo(sequenceName, length, md5);

    // Update refNameToRefInfoIndex
    if (refNameToRefInfoIndex.find(sequenceName) != refNameToRefInfoIndex.end()) {
      cout << "ERROR. Reference name " << sequenceName 
            << " is not unique." << endl;
      exit(1);
    }
    refNameToRefInfoIndex[sequenceName] = r;
  }
}

template<typename T_CmpFile>
unsigned int AlignmentSetToCmpH5Adapter<T_CmpFile>::StoreRefGroup(
        std::string refName, T_CmpFile & cmpFile) {
  // Find out whether there is a refGroup associated with refName.
  std::map<std::string, RefGroupNameId>::iterator mapIt;
  mapIt = refNameToRefGroupNameandId.find(refName);
  if (mapIt != refNameToRefGroupNameandId.end()) {
    // An existing refGroup is associated with this refName. 
    return mapIt->second.id;
  } else {
    // No refGroup is associated with refName, create one.
    int refInfoIndex = refNameToRefInfoIndex[refName];
    unsigned int refInfoId = refInfoIndex + 1;

    std::string refGroupName;
    unsigned int refGroupId = cmpFile.AddRefGroup(refName, refInfoId, refGroupName);  

    // Update refNameToRefGroupNameandId.
    refNameToRefGroupNameandId[refName] = RefGroupNameId(refGroupName, refGroupId);
    return refGroupId;
  }
}

template<typename T_CmpFile>
unsigned int AlignmentSetToCmpH5Adapter<T_CmpFile>::StorePath(
      std::string & path, T_CmpFile &cmpFile) {
  if (knownPaths.find(path) != knownPaths.end()) {
    return knownPaths[path];
  }
  else {
    unsigned int id = cmpFile.alnGroupGroup.AddPath(path);
    knownPaths[path] = id;
    return id;
  }
}

template<typename T_CmpFile>
void AlignmentSetToCmpH5Adapter<T_CmpFile>::RemoveGapsAtEndOfAlignment(AlignmentCandidate<> &alignment) {
  int numEndDel = 0, numEndIns = 0;
  if (alignment.gaps.size() > 0) {
    int lastGap = alignment.gaps.size() - 1;
    int g;
    for (g = 0; g < alignment.gaps[lastGap].size(); g++) {
      if (alignment.gaps[lastGap][g].seq == Gap::Target) {
        numEndIns += alignment.gaps[lastGap][g].length;
      }
      else if (alignment.gaps[lastGap][g].seq == Gap::Query) {
        numEndDel += alignment.gaps[lastGap][g].length;
      }
    }
  }
  alignment.qAlignedSeqLength -= numEndIns;
  alignment.tAlignedSeqLength -= numEndDel;
}

template<typename T_CmpFile>
void AlignmentSetToCmpH5Adapter<T_CmpFile>::StoreAlignmentCandidate(
    AlignmentCandidate<> &alignment, 
    int alnSegment,
    T_CmpFile &cmpFile,
    int moleculeNumber,
    bool copyQVs) {
  //
  // Find out where the movie is going to get stored.
  //
  std::string movieName;
  int holeNumber = 0;
  bool nameParsedProperly;
  
  nameParsedProperly = ParsePBIReadName(alignment.qName, movieName, holeNumber);
  if (!nameParsedProperly) {
    std::cout <<"ERROR. Attempting to store a read with name " 
          << alignment.qName << " that does not " << std::endl
          << "appear to be a PacBio read." << std::endl;
    exit(1);
  }

  unsigned int movieId = StoreMovieInfo(movieName, cmpFile);

  // Check whether the reference is in /RefInfo.
  std::map<std::string, int>::iterator mapIt;
  mapIt = refNameToRefInfoIndex.find(alignment.tName);
  if (mapIt == refNameToRefInfoIndex.end()) {
    std::cout << "ERROR. The reference name " << alignment.tName 
          << " was not found in the list of references." << std::endl;
    std::cout << "Perhaps a different reference file was aligned to than " << std::endl
          << "what was provided for SAM conversion. " << std::endl;
    exit(1);
  } 

  // Store refGroup
  unsigned int refGroupId = StoreRefGroup(alignment.tName, cmpFile);
  std::string refGroupName = refNameToRefGroupNameandId[alignment.tName].name; 
  assert(refGroupId  == refNameToRefGroupNameandId[alignment.tName].id);

  if (cmpFile.refGroupIdToArrayIndex.find(refGroupId) == cmpFile.refGroupIdToArrayIndex.end()) {
    std::cout << "ERROR. The reference ID is not indexed. " 
          << "This is an internal inconsistency." << std::endl;
    exit(1);
  }

  int    refGroupIndex= cmpFile.refGroupIdToArrayIndex[refGroupId];
  assert(refGroupIndex + 1 == refGroupId);

  std::string path = "/" + refGroupName + "/" + movieName;
  unsigned int pathId = StorePath(path, cmpFile);

  vector<unsigned int> alnIndex;
  alnIndex.resize(22);

  RemoveGapsAtEndOfAlignment(alignment);

  /*
    * Store the alignment string
    */
  vector<unsigned char> byteAlignment;
  AlignmentToByteAlignment(alignment, 
                            alignment.qAlignedSeq, alignment.tAlignedSeq,
                            byteAlignment);

  unsigned int offsetBegin, offsetEnd;
  cmpFile.StoreAlnArray(byteAlignment, alignment.tName, movieName, offsetBegin, offsetEnd);
  // Copy QVs into cmp.h5
  if (copyQVs) {
    std::vector<std::string> optionalQVs;
    alignment.CopyQVs(&optionalQVs);
    for (int qv_i=0; qv_i<optionalQVs.size(); qv_i++) {
      std::string *qvName = &alignment.optionalQVNames[qv_i];
      std::string *qvString = &optionalQVs[qv_i];
      
      // If the qvString is empty, then the alignment is missing the quality
      // value
      if (qvString->size() == 0) {
        continue;
      }

      unsigned int qvOffsetBegin, qvOffsetEnd;
      if (qvName->compare(qvName->size() - 3, 3, "Tag") == 0) {
        std::vector<char> qvVector;
        QVsToCmpH5QVs(*qvString, byteAlignment, true, &qvVector);
        cmpFile.StoreTags(qvVector, alignment.tName, *qvName,
                          movieName, &qvOffsetBegin, &qvOffsetEnd);
      } else {
        std::vector<UChar> qvVector;
        QVsToCmpH5QVs(*qvString, byteAlignment, false, &qvVector);
        cmpFile.StoreQVs(qvVector, alignment.tName, *qvName,
                         movieName, &qvOffsetBegin, &qvOffsetEnd);
      }
      assert(qvOffsetBegin == offsetBegin);
      assert(qvOffsetEnd == offsetEnd);
    }
  }

  numAlignments++;

  DistanceMatrixScoreFunction<DNASequence, DNASequence> distScoreFn;
  //distScoreFn does not matter since the score is not stored.
  ComputeAlignmentStats(alignment, alignment.qAlignedSeq.seq, alignment.tAlignedSeq.seq, distScoreFn);

  /*
    The current AlnIndex column names:
    (0): "AlnID", "AlnGroupID", "MovieID", "RefGroupID", "tStart",
    (5): "tEnd", "RCRefStrand", "HoleNumber", "SetNumber",
    (9): "StrobeNumber", "MoleculeID", "rStart", "rEnd", "MapQV", "nM",
    (15): "nMM", "nIns", "nDel", "Offset_begin", "Offset_end",
    (20): "nBackRead", "nReadOverlap"
  */
  if (moleculeNumber == -1) {
    moleculeNumber =  numZMWsPerMovieSpringField * (movieId - 1) + holeNumber;
  }
  alnIndex[0]  = numAlignments;  // AlnId
  alnIndex[1]  = pathId;        // AlnGroupID
  alnIndex[2]  = movieId;    // MovieID
  alnIndex[3]  = refGroupId; // RefGroupID
  alnIndex[4]  = alignment.tAlignedSeqPos; // tStart
  alnIndex[5]  = alignment.tAlignedSeqPos +  alignment.tAlignedSeqLength; // tEnd
  alnIndex[6]  = alignment.tStrand; // RCRefStrand
  alnIndex[7]  = holeNumber;
  alnIndex[8]  = 0; // SET NUMBER -- parse later!!!!
  alnIndex[9]  = alnSegment; // strobenumber
  alnIndex[10] = moleculeNumber;
  alnIndex[11] = alignment.qAlignedSeqPos; 
  alnIndex[12] = alignment.qAlignedSeqPos + alignment.qAlignedSeqLength;
  alnIndex[13] = alignment.mapQV;
  alnIndex[14] = alignment.nMatch;
  alnIndex[15] = alignment.nMismatch;
  alnIndex[16] = alignment.nIns;
  alnIndex[17] = alignment.nDel;
  alnIndex[18] = offsetBegin;
  alnIndex[19] = offsetEnd;
  alnIndex[20] = 0;
  alnIndex[21] = 0;
  cmpFile.alnInfoGroup.WriteAlnIndex(alnIndex);
}

template<typename T_CmpFile>
void AlignmentSetToCmpH5Adapter<T_CmpFile>::StoreAlignmentCandidateList(
    std::vector<AlignmentCandidate<> > &alignments,
    T_CmpFile &cmpFile,
    int moleculeNumber,
    bool copyQVs) {
  
  int a;
  for (a = 0; a < alignments.size(); a++) {
    StoreAlignmentCandidate(alignments[a], a, cmpFile, moleculeNumber, copyQVs);
  }
}

#endif

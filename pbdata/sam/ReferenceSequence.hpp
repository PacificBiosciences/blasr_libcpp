#ifndef _BLASR_SAM_REFERENCE_SEQUENCE_HPP_
#define _BLASR_SAM_REFERENCE_SEQUENCE_HPP_

#include <string>

#include "sam/SAMKeywordValuePair.hpp"

class SAMReferenceSequence {
 public:
  std::string sequenceName;
  uint64_t length;
  std::string GetSequenceName() {
    return sequenceName;
  }
  //
  // By definining accessor functions here but no data, we can
  // economize on the amount of space used for each element.   This is
  // no big deal for references, but for pairwise alignments, it is
  // big.
  //
  std::string GetMD5() {
    return "";
  }
  uint64_t GetLength() {
    return length;
  }
  std::string GetGenomeAssembly() {
    return "";
  }
  std::string GetSpecies() {
    return "";
  }
  std::string GetURI() {
    return "";
  }

  enum SAMReferenceSequenceRequiredFields { SQ_SN, SQ_LN };

  static const char* SAMReferenceSequenceFieldNames[];

  void StoreValues(std::vector<SAMKeywordValuePair> &kvPairs,
                   uint64_t lineNumber=0);
};

class SAMFullReferenceSequence : public SAMReferenceSequence {
 public:
  std::string md5;
  std::string species;
  std::string uri;
  std::string genomeAssembly;
  std::string GetMD5() {
    return md5;
  }
  std::string GetSpecies() {
    return species;
  }
  std::string GetURI() {
    return uri;
  }
  std::string GetGenomeAssembly() {
    return genomeAssembly;
  }

  enum FullReferenceSequenceRequiredFields { SQ_AS, SQ_M5, SQ_SP, SQ_UR };

  static const char* SAMFullReferenceSequenceFieldNames[];

  void StoreValues(std::vector<SAMKeywordValuePair> &kvPairs, uint64_t lineNumber=0);
};


#endif

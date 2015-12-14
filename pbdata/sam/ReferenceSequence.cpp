#include "ReferenceSequence.hpp"

const char* SAMReferenceSequence::SAMReferenceSequenceFieldNames[] = {"SN", "LN"};

void SAMReferenceSequence::StoreValues(std::vector<SAMKeywordValuePair> &kvPairs, uint64_t lineNumber) {
  (void)(lineNumber);
  std::vector<bool> usedFields;
  usedFields.resize(SQ_LN);
  std::fill(usedFields.begin(), usedFields.end(), false);
  for (size_t i = 0; i < kvPairs.size(); i++) {
    if (kvPairs[i].key == "SN") {
      sequenceName = kvPairs[i].value;
      usedFields[SQ_SN] = true;
    }
    else if (kvPairs[i].key == "LN") {
      StoreValue(kvPairs[i].value, length);
      usedFields[SQ_SN] = true;
    }
  }
  for (size_t i = 0; i < usedFields.size(); i++) {
    if (usedFields[i] == false) {
      std::cout << "SQ specifier missing " << SAMReferenceSequenceFieldNames[i] << std::endl;
   }
  }
}

const char* SAMFullReferenceSequence::SAMFullReferenceSequenceFieldNames[] = {"AS", "M5", "SP", "UR"};

void SAMFullReferenceSequence::StoreValues(std::vector<SAMKeywordValuePair> &kvPairs,
                                           uint64_t lineNumber) {
  SAMReferenceSequence::StoreValues(kvPairs, lineNumber);
  for (size_t i = 0; i < kvPairs.size(); i++ ){
    if (kvPairs[i].key == "AS") {
      genomeAssembly = kvPairs[i].value;
    }
    else if (kvPairs[i].key == "M5") {
      md5 = kvPairs[i].value;
    }
    else if (kvPairs[i].key == "SP") {
      species = kvPairs[i].value;
    }
    else if (kvPairs[i].key == "UR") {
      uri = kvPairs[i].value;
    }
  }
}

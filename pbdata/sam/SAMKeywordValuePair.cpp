#include "SAMKeywordValuePair.hpp"


bool SplitSAMKeyValuePair(std::string &kvPair, std::string &key, std::string &value) {
  size_t sepIndex = kvPair.find_first_of(":");
  if (sepIndex == kvPair.npos) {
    return false;
  }
  else {
    key = kvPair.substr(0, sepIndex);
    value = kvPair.substr(sepIndex+1);
    return true;
  }
}

bool SplitSAMTypedKeyValuePair(std::string kvPair, std::string &key,
                               std::string &kvType, std::string &value) {
  std::vector<std::string> strValues;
  ParseSeparatedList(kvPair, strValues, ':', 3);
  if (strValues.size() != 3) {
    return false;
  }
  else {
    key = strValues[0];
    kvType = strValues[1];
    value  = strValues[2];
    return true;
  }
}

bool TypedKeywordValuePair::Separate(std::string &kvPair, std::string &kvKey,
                                     std::string &kvType, std::string &kvValue) {
  if (SplitSAMTypedKeyValuePair(kvPair, kvKey, kvType, kvValue) == false) {
    return false;
  }
  return true;
}

void KeywordValueStringsToPairs(std::vector<std::string> &kvStrings, std::vector<SAMKeywordValuePair> &kvPairs) {

  kvPairs.resize(kvStrings.size());

  if (kvStrings.size() == 0) {
    return;
  }

  for (size_t i = 0; i < kvStrings.size(); i++ ) {
    SplitSAMKeyValuePair(kvStrings[i], kvPairs[i].key, kvPairs[i].value);
  }
}

#ifndef _BLASR_SAM_READ_GROUP_HPP_
#define _BLASR_SAM_READ_GROUP_HPP_

#include <string>

/*
 * Minimal components of read group.  Define only required items. Use
 * if memory is scarce.  
 */
#include "SAMKeywordValuePair.hpp"

class SAMReadGroup {
 public:
  std::string id;
  void StoreValues(std::vector<SAMKeywordValuePair> &kvPairs,
                   uint64_t lineNumber = 0);
};

/*
 * Full read group. Use when all data are required.
 */
class SAMFullReadGroup : public SAMReadGroup {
 public:
  std::string centerName;
  std::string description;
  std::string date;
  std::string flowOrder;
  std::string keySequence;
  std::string library;
  std::string processingProgram;
  uint64_t    insertSize;
  std::string platform;
  std::string platformUnit;
  std::string sample;

  void StoreValues(std::vector<SAMKeywordValuePair> &kvPairs,
                   uint64_t lineNumber = 0);
};

#endif

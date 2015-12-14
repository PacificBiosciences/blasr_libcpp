#ifndef _BLASR_SAM_ALIGNMENT_HPP_
#define _BLASR_SAM_ALIGNMENT_HPP_

#include <string>
#include <sstream>

#include "sam/SAMKeywordValuePair.hpp"
#include "sam/CigarString.h"

class SAMAlignment {
 public:
  enum SAMAlignmentRequiredFields {
    S_QNAME,
    S_FLAG,
    S_RNAME,
    S_POS,
    S_MAPQV,
    S_CIGAR,
    S_RNEXT,
    S_PNEXT,
    S_TLEN,
    S_SEQ,
    S_QUAL};

  static const char* SAMAlignmentRequiredFieldNames[];

  std::string qName;
  unsigned int flag;
  std::string rName;
  unsigned int pos;
  short mapQV;
  CigarString cigar;
  std::string rNext;
  int pNext;
  int tLen;
  std::string seq;
  std::string qual;

  // Optional tags defined in blasr:
  // "RG" read group Id
  // "AS" alignment score
  // "XS" read alignment start position without counting previous soft clips (1 based) 
  // "XE" read alignment end position without counting previous soft clips (1 based) 
  // "XL" aligned read length 
  // "XQ" query read length
  // "XT" # of continues reads, always 1 for blasr 
  // "NM" # of subreads 
  // "FI" read alignment start position (1 based) 
  //
  float score;
  int xs, xe;
  int xl;
  int xq;

  std::string rg;
  int as;
  int xt;
  int nm;
  int fi;
  std::string optTagStr; 
  //
  // Quality values.
  //
  std::string iq, dq, sq, mq, st, dt;	
  static const char* optionalQVTags[];
  static const char* optionalQVNames[];

  //
  // Initialize all optional fields.  Required fields will be
  // assigned a value later.
  //
  SAMAlignment();

  void PrintSAMAlignment(std::ostream &out = std::cout);

  // Find position of the nth character in a string.
  size_t FindPosOfNthChar(std::string str, int n, char c);

  // Trim the end '\n\r' characters from a string.
  std::string TrimStringEnd(std::string str);

  bool StoreValues(std::string &line,  int lineNumber=0);
  
  // CopyQVs writes the strings from the optional QV tags to a vector. The
  // order of QVs in the vector is given by optionalQVNames[]
  void CopyQVs(std::vector<std::string> *optionalQVs);
};


class SAMPosAlignment : public SAMAlignment {
 public:
  unsigned int qStart, qEnd;
  unsigned int tStart, tEnd;
  int qStrand, tStrand;
};

#endif

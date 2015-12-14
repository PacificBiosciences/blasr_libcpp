#include "SAMAlignment.hpp"

SAMAlignment::SAMAlignment() {
  score = xs = xe = as = xt = xq = nm = fi = xl = 0;
  rg = optTagStr = "";
}

void SAMAlignment::PrintSAMAlignment(std::ostream & out) {
    out << qName << "\t" << flag  << "\t" << rName << "\t"
        << pos   << "\t" << mapQV << "\t" << cigar << "\t"
        << rNext << "\t" << pNext << "\t" << tLen  << "\t"
        << seq   << "\t" << qual  << "\t" << optTagStr 
        << std::endl;
}

size_t SAMAlignment::FindPosOfNthChar(std::string str, int n, char c) {
  if (n < 1) {
    std::cout << "Nth should be a positive integer." << std::endl;
    exit(0);
  }
  int count = 1;
  size_t pos = str.find(c, 0);
  // pos is the position of the first character c;
  while(count < n and pos != std::string::npos) {
      pos = str.find(c, pos+1);
      count ++;
  }
  return pos;
}

std::string SAMAlignment::TrimStringEnd(std::string str) {
  std::string newStr = str;
  while(newStr[newStr.size()-1] == '\r' or
        newStr[newStr.size()-1] == '\n') {
      newStr.erase(newStr.size()-1, 1);
  }
  return newStr;
}

bool SAMAlignment::StoreValues(std::string &line,  int lineNumber) {
  std::stringstream strm(line);
  std::vector<bool> usedFields;
  usedFields.resize(S_QUAL);
  fill(usedFields.begin(), usedFields.end(), false);
  std::string kvPair;
  bool parseError = false;
  SAMAlignmentRequiredFields field;
  //
  // Define a temporary mapqv value that gets over a GMAP bug that prints a mapqv < 0.
  //
  int tmpMapQV;
  if (!(strm >> qName)) {
    parseError = true;
    field = S_QNAME;
  }
  else if (! (strm >> flag) ){ 
    parseError = true;
    field = S_FLAG;
  }
  else if (! (strm >> rName) ) {
    parseError = true;
    field = S_RNAME;
  }
  else if (! (strm >> pos) ) {
    parseError = true;
    field = S_POS;
  }
  else if (! (strm >> tmpMapQV)) {
    parseError = true; field = S_MAPQV;
  }
  else if (! (strm >> cigar)) {
    parseError = true; field = S_CIGAR;
  }

  else if (! (strm >> rNext)) {
    parseError = true; field = S_RNEXT;
  }

  else if (! (strm >> pNext)) {
    parseError = true; field = S_PNEXT;
  }
  else if (! (strm >> tLen)) {
    parseError = true; field = S_TLEN;
  }
  else if (! (strm >> seq)) {
    parseError = true; field = S_SEQ;
  }
  else if (! (strm >> qual)) {
    parseError = true; field = S_QUAL;
  }

  mapQV = (unsigned char) tmpMapQV;

  // Find posisition of the 11th tab.
  size_t optTagsStartPos = FindPosOfNthChar(strm.str(), 11, '\t');
  // Save all optional tags.
  if (optTagsStartPos != std::string::npos) {
      optTagStr = strm.str().substr(optTagsStartPos+1);
      optTagStr = TrimStringEnd(optTagStr);
  } 

  //
  // If not aligned, stop trying to read in elements from the sam string.
  //
  if (rName == "*") {
    return true;
  }

  if (parseError) {
      std::cout << "Error parsing alignment line " << lineNumber << ". Missing or error in field " << SAMAlignmentRequiredFieldNames[field] << std::endl;
      exit(1);
  }
  
  //
  // Now parse optional data.
  //
  while (strm) {
    std::string kvName, kvType, kvValue;
    std::string typedKVPair;
    if (!(strm >> typedKVPair)) {
      break;
    }
    if (TypedKeywordValuePair::Separate(typedKVPair, kvName, kvType, kvValue)) {
      std::stringstream strm(kvValue);
      if (kvName == "RG") {
        rg = kvValue;
      }
      else if (kvName == "AS") {
        strm >> as;
      }
      else if (kvName == "XS") {
        strm >> xs;
      }
      else if (kvName == "XE") {
        strm >> xe;
      }
      else if (kvName == "XL") {
        strm >> xl;
      }
      else if (kvName == "XT") {
        strm >> xt;
      }
      else if (kvName == "NM") {
        strm >> nm;
      }
      else if (kvName == "FI") {
        strm >> fi;
      }
      else if (kvName == "XQ") {
        strm >> xq;
      } // Add quality values, including QualityValue?, 
        // InsertionQV, DeletionQV, SubstitutionQV, 
        // MergeQV and SubstitutionTag and DeletionTag
      else if (kvName == "iq") {
          strm >> iq;
      }
      else if (kvName == "dq") {
          strm >> dq;
      }
      else if (kvName == "sq") {
          strm >> sq;
      }
      else if (kvName == "mq") {
          strm >> mq;
      }
      else if (kvName == "st") {
          strm >> st;
      }
      else if (kvName == "dt") {
          strm >> dt;
      }
    }
    else {
      std::cout << "ERROR. Could not parse typed keyword value " << typedKVPair << std::endl;
      exit(0);
    }
  }
  return true;
}

void SAMAlignment::CopyQVs(std::vector<std::string> *optionalQVs) {  
  optionalQVs->clear();
  optionalQVs->push_back(iq);
  optionalQVs->push_back(dq);
  optionalQVs->push_back(sq);
  optionalQVs->push_back(mq);
  optionalQVs->push_back(st);
  optionalQVs->push_back(dt);
}

const char* SAMAlignment::SAMAlignmentRequiredFieldNames[] = { "QNAME", "FLAG", 
                                                               "RNAME", "POS", 
                                                               "MAPQV", "CIGAR", 
                                                               "RNEXT", "PNEXT", 
                                                               "TLEN", "SEQ", "QUAL"} ;
const char* optionalQVTags[] = {"iq", "dq", "sq", "mq", "st", "dt"};
const char* optionalQVNames[] = {"InsertionQV", "DeletionQV", "SubstitutionQV",
                                 "MergeQV", "SubstitutionTag", "DeletionTag"};

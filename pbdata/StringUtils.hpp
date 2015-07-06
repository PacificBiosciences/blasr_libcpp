#ifndef UTILS_STRING_UTILS_H_
#define UTILS_STRING_UTILS_H_

#include <cassert>
#include <string>
#include <vector>
#include <sstream>
#include "Types.h"
#include "MD5Utils.hpp"
#include "reads/ReadType.hpp"

using namespace std;

int ExactPatternMatch(string orig, string pattern); 

void MakeMD5(const char *data, unsigned int dataLength, string &md5Str, int nChars = 0); 

void MakeMD5(string &data, string &md5Str, int nChars=0); 

int IsWhitespace(char c); 

int IsSpace(char c); 

int ToWords(string &orig, vector<string> &words); 

int Splice(const string & orig, const string & pattern, vector<string> &tokens); 

void ParseSeparatedList(const string &csl, vector<string> &values, char delim=','); 

int AssignUntilFirstSpace(char *orig, int origLength, string &result); 

template<typename T_Value>
void ParseSeparatedList(const string &csl, vector<T_Value> &values, char delim, int maxVals) {
  //Parse up to 'maxVals' lines of lists.
  stringstream cslStrm(csl);
  T_Value val;
  string valString;
  string next;
  int valIndex = 0;
  do {
    if (maxVals == 0 or valIndex < maxVals - 1 ) {
      getline(cslStrm, valString, delim);
    }
    else {
      // If on last value, get the rest of the line.
      getline(cslStrm, valString);
    }
    if (cslStrm and valString.size() > 0) {
      stringstream valStrm(valString);
      if (! (valStrm  >> val) ) {
        cout << "Error, value " << valString << " is malformatted." << endl;
      }
      else {
        values.push_back(val);
      }
    }
    valIndex++;
  }
  while (cslStrm);
}

template<typename T_Value>
void ParseSeparatedList(const string &csl, vector<T_Value> &values, char delim=',') {
	ParseSeparatedList(csl, values, delim, 0);
}

string RStrip(string & fileName);

/// \returns Read Group Id, given movieName and readType.
string MakeReadGroupId(const string & movieName, const ReadType::ReadTypeEnum & readType);
#endif

#include <string>
#include <stdint.h>
#include <vector>
#include <sstream>
#include "Types.h"
#include "MD5Utils.hpp"
#include "StringUtils.hpp"

int ExactPatternMatch(string orig, string pattern) {
    string::size_type pos = orig.find(pattern);
    if (pos == orig.npos) {
        return 0;
    }
    else {
        return 1;
    }
}

void MakeMD5(const char *data, unsigned int dataLength, string &md5Str, int nChars) {

    MD5 md5engine;
    md5engine.update((unsigned char*) data, dataLength);
    md5engine.finalize();

    char *md5c_str = md5engine.hex_digest();
    assert(md5c_str != NULL);
    if (nChars == 0) {
        nChars = 32;
    }
    md5Str.assign(md5c_str, nChars);
    delete[] md5c_str;
}


void MakeMD5(string &data, string &md5Str, int nChars) {
    MakeMD5(data.c_str(), data.size(), md5Str, nChars);
}


int IsWhitespace(char c) {
    return (c == ' ' or c == '\t' or c == '\n' or c == '\r' or c == '\0');
}

int IsSpace(char c) {
    return (c == ' ' or c == '\t');
}

int ToWords(string &orig, vector<string> &words) {
    int curWordStart, curWordEnd;
    curWordStart = 0;
    while(curWordStart < orig.size()) {
        while (curWordStart < orig.size() and IsSpace(orig[curWordStart])) { curWordStart++; }
        curWordEnd = curWordStart;
        while (curWordEnd < orig.size() and !IsSpace(orig[curWordEnd])) { curWordEnd++; }
        string word;
        if (curWordEnd != curWordStart) {
            word.assign(orig, curWordStart, curWordEnd - curWordStart);
            words.push_back(word);
        }
        curWordStart = curWordEnd;
    }
    return words.size();
}

// Splice a string by pattern and save to a vector of token strings.
int Splice(const string & orig, const string & pattern, vector<string> & tokens) {
    assert(pattern.size() > 0);

    tokens.clear();
    size_t search_start = 0;
    size_t find_pos = orig.find(pattern, search_start);
    while(find_pos != string::npos) {
        string x = orig.substr(search_start, find_pos - search_start);
        tokens.push_back(x);
        search_start = find_pos + pattern.size();
        find_pos = orig.find(pattern, search_start);
    }
    tokens.push_back(orig.substr(search_start));
    return tokens.size();
}

void ParseSeparatedList(const string &csl, vector<string> &values, char delim) {
    stringstream cslStrm(csl);
    string valString;
    string next;
    do {
        if (getline(cslStrm, valString, delim)) {
            if (valString.size() > 0) {
                values.push_back(valString);
            }
        }
    }
    while (cslStrm);
}

int AssignUntilFirstSpace(char *orig, int origLength, string &result) {
    int i;
    for (i = 0; i < origLength; i++ ){ 
        if (orig[i] == ' ' or orig[i] == '\t' or orig[i] == '\n' or orig[i] == '\r' or orig[i] == '\0') {
            break;
        }
    }
    result.assign(orig, i);
    return i;
}

string RStrip(string & fileName) {
    // Remove right-ended spaces
    int i = fileName.size();
    if (i == 0) {
        return "";
    }
    while (i >= 1) {
        i--;
        if (not IsWhitespace(fileName[i])) {
            break;
        }
    }
    return fileName.substr(0, i + 1);
}

string MakeReadGroupId(const string & movieName, const ReadType::ReadTypeEnum & readType) {
    // PBBAM spec 3.0b5:
    // Read Group Id is computed as MD5(${movieName}//${readType})[0:8], where
    // movieName is PacBio platform unit id, e.g., (m140905_042...77_s1_X0),
    // readtype is SUBREAD, CCS or UNKNOWN, 
    // CCS reads for a movie named "movie32" would have 
    // RGID STRING = "f5b4ffb6"
    string seed = movieName + "//" + ReadType::ToString(readType);
    string readGroupId;
    MakeMD5(seed, readGroupId, 8);
    return readGroupId;
}

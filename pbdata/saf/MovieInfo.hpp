#ifndef DATASTRUCTURES_SAF_HDF_MOVIE_INFO_H_
#define DATASTRUCTURES_SAF_HDF_MOVIE_INFO_H_

#include "Types.h"
#include <string>
#include <vector>

using namespace std;

class MovieInfo {
 public:
	vector<string> name;
	vector<UInt>   run;
	vector<UInt>   experiment;
	vector<UInt>   id;
	int FindMovie(UInt idKey, string &nameVal);
};

#endif

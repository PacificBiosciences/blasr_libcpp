#include "MovieInfo.hpp"


int MovieInfo::FindMovie(UInt idKey, string &nameVal) {
    for (size_t i = 0; i < id.size(); i++) {
        if (id[i] == idKey) {
            nameVal = name[i];
            return 1;
        }
    }
    return 0;
}

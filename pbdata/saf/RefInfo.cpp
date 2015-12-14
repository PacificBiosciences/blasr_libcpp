#include "RefInfo.hpp"

OneRefInfo::OneRefInfo() {
    fullName = md5 = "";
    id = length = 0;
}

bool RefInfo::RefIdToIndex(uint32_t qid, int &index) {
/*
 * Translate from external id to position in array.
 */
    for (size_t i = 0; i < id.size(); i++) {
        if (id[i] == qid) {
        index = i;
        return true;
        }
    }
    return false;
}

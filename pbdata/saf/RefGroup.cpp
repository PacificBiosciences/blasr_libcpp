#include "RefGroup.hpp"


bool RefGroup::IdToIndex(uint32_t idKey, uint32_t &idIndex) {
    for (size_t i = 0; i < refInfoId.size(); i++) {
        if (refInfoId[i] == idKey) {
        idIndex = i; return true;
        }
    }
    return false;
}
	

int RefGroup::FindPath(uint32_t idKey, string &pathVal, string &groupNameVal) {
    for (size_t i = 0; i < id.size(); i++) {
        if (id[i] == idKey) {
            pathVal = path[i];
            groupNameVal = refGroupName[i];
            return 1;
        }
    }
    return 0;
}

#include "AlnGroup.hpp"

int AlnGroup::FindPath(unsigned int idKey, std::string &val) {
    for (size_t i = 0; i < id.size(); i++) {
        if (idKey == id[i]) {
            val = path[i];
            return 1;
        }
    }
    return 0;
}

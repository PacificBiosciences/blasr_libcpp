#include <cassert>
#include "HDFAttributable.hpp"

using namespace std;
using namespace H5;

void CallStoreAttributeName(H5Location &obj, string attrName, void *attrList){ 
    ((vector<string>*)attrList)->push_back(string(attrName));
}

void HDFAttributable::StoreAttributeNames(H5Location &thisobject, 
    const std::vector<std::string> &attributeNames) {
    int nAttr = thisobject.getNumAttrs();
    unsigned int bounds[2];
    bounds[0] = 0;
    bounds[1] = nAttr;
    attributeNameList.clear();
    thisobject.iterateAttrs(&CallStoreAttributeName, 
            bounds, (void*) &attributeNames);
}

H5Location* HDFAttributable::GetObject() {
    return NULL;
}

int HDFAttributable::ContainsAttribute(const string & attributeName) {
    size_t i;
    std::vector<std::string> tmpAttributeNames;
    try{
        H5Location *obj = GetObject();
        assert(obj != NULL);
        StoreAttributeNames(*obj, tmpAttributeNames);
        for (i = 0; i < tmpAttributeNames.size(); i++) {
            if (tmpAttributeNames[i] == attributeName) return true;
        }
    } catch (H5::Exception e) {
        //Failed to read attribute // e.printError();
    }
    return false;
}


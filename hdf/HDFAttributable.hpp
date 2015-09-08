#ifndef _BLASR_HDF_ATTRIBUTABLE_HPP_
#define _BLASR_HDF_ATTRIBUTABLE_HPP_

#include <string>
#include <vector>
#include "H5Cpp.h"

void CallStoreAttributeName(H5::H5Location &obj, std::string attrName, void *attrListPtr);

class HDFAttributable {
public:
    std::vector<std::string> attributeNameList;

    void StoreAttributeNames(H5::H5Location &thisobject, 
        const std::vector<std::string> &attributeNames); 

    virtual H5::H5Location* GetObject(); 

    int ContainsAttribute(const std::string & attributeName); 

};

#endif

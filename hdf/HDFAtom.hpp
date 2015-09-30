#ifndef _BLASR_HDF_ATOM_HPP_
#define _BLASR_HDF_ATOM_HPP_

#include <cassert>
#include <string>
#include <vector>
#include <cstdint>
#include "H5Cpp.h"

#include "HDFConfig.hpp"
#include "HDFGroup.hpp"
#include "HDFData.hpp"

template<typename T>
class HDFAtom : public HDFData {
public:
    H5::Attribute attribute;

    HDFAtom() {
        isInitialized = false;
    }
    ~HDFAtom() {
        if (IsInitialized()) {
            attribute.close();
        }
    }

    H5::H5Location *GetObject() {
        return NULL;
    }

    int Initialize(H5::H5Location &object, const std::string & attributeName, bool createIfMissing=false) {
        attribute = object.openAttribute(attributeName.c_str());
        isInitialized = true;
        return 1;
    }

    int Initialize(HDFGroup &group, const std::string & attributeName, bool createIfMissing=false) {
        return Initialize(group.group, attributeName);
    }

    int Initialize(HDFData &data, const std::string & attributeName, bool createIfMissing=false) {
        return Initialize(data.dataset, attributeName);
    }

    int Initialize(H5::Group &object, const std::string & attributeName, bool createIfMissing=false) {
        try {
            attribute = object.openAttribute(attributeName.c_str());
        }
        catch (H5::Exception e) {
            cout << "ERROR. Could not open attribute " << attributeName << endl;
            exit(1);
        }
        isInitialized  = true;
        return 1;
    }

    int Initialize(H5::H5File &hdfFile, const std::string & groupName, const std::string & attributeName) {
        HDFGroup group;
        group.Initialize(hdfFile, groupName);
        attribute = group.group.openAttribute(attributeName.c_str());
        isInitialized = true;
        return 1;
    }

    //
    // This handles creation of all non-std::string types.  A specialization
    // for std::strings is provided below.
    //
    void Create(H5::H5Location &object, const std::string & atomName) {
        hsize_t defaultDims[] = {1};
        H5::DataSpace defaultDataSpace(1, defaultDims);
        TypedCreate(object, atomName, defaultDataSpace);
    }


    void Create(H5::H5Location &object, const std::string & name, const std::string & value) {
        H5::StrType strType(0, value.size());
        attribute = object.createAttribute(name.c_str(), strType, H5::DataSpace(0,NULL));
        isInitialized = true;
        attribute.write(strType, value.c_str());
    }

    void Create(H5::H5Location &object, const std::string & name, std::vector<int> &vect) {
        hsize_t length = vect.size();
        H5::ArrayType arrayDataType(H5::PredType::NATIVE_INT, 1, &length);
        attribute = object.createAttribute(name.c_str(), H5::PredType::NATIVE_INT, H5::DataSpace(1, &length));
        attribute.write(H5::PredType::NATIVE_INT, &((vect)[0]));    
    }

    void Create(H5::H5Location &object, const std::string & name, const std::vector<std::string> &vect) {
        hsize_t length = vect.size();
        H5::StrType strType(0,H5T_VARIABLE);
        H5::ArrayType arrayDataType(strType, 1, &length);
        attribute = object.createAttribute(name.c_str(), strType, H5::DataSpace(1, &length));
        attribute.write(strType, &((vect)[0]));    
    }

    void TypedCreate(H5::H5Location &object, const std::string &atomName, H5::DataSpace &dataSpace) {
        assert("Calling HDFAtom<T>::typedCreate on an unsupported type" == 0);
    }

    void Write(T value) {
        assert("Calling HDFAtom<T>::Write on an unsupported type" == 0);
    }

    void Read(T& value) {
        assert("Calling read on an unsupported type!" == 0);
    }

};

//
// Special create for std::strings.  Since this uses a H5::StrType for the
// typename rather than specifying a H5::PredType, it mertis its own
// function.
//

template<>
void HDFAtom<std::string>::Create(H5::H5Location &object, const std::string & atomName);

#define DECLARE_TYPED_CREATE_ATOM(T, predtype) template<> \
	void HDFAtom<T>::TypedCreate(H5::H5Location &object, const std::string & atomName, H5::DataSpace &defaultDataSpace);

DECLARE_TYPED_CREATE_ATOM(int, H5::PredType::NATIVE_INT)
DECLARE_TYPED_CREATE_ATOM(unsigned int, H5::PredType::NATIVE_UINT)
DECLARE_TYPED_CREATE_ATOM(unsigned char, H5::PredType::NATIVE_UINT8)
DECLARE_TYPED_CREATE_ATOM(char, H5::PredType::NATIVE_INT8)
DECLARE_TYPED_CREATE_ATOM(float, H5::PredType::NATIVE_FLOAT)
DECLARE_TYPED_CREATE_ATOM(uint64_t, H5::PredType::STD_I64LE)


template<>
void HDFAtom<std::vector<int> >::Write(const std::vector<int> vect);

template<>
void HDFAtom<std::string>::Write(std::string value);

template<>
void HDFAtom<uint64_t>::Write(uint64_t value);

template<>
void HDFAtom<int>::Write(int value);

template<>
void HDFAtom<unsigned int>::Write(unsigned int value);

template<>
void HDFAtom<unsigned char>::Write(unsigned char value);

template<>
void HDFAtom<char>::Write(char value);

template<>
void HDFAtom<float>::Write(float value);

template<>
void HDFAtom<std::string>::Read(std::string &value);

template<>
void HDFAtom<int>::Read(int &value);

template<>
void HDFAtom<uint64_t>::Read(uint64_t &value);

template<>
void HDFAtom<unsigned int>::Read(unsigned int &value);

template<>
void HDFAtom<float>::Read(float &value);

template<>
void HDFAtom<std::vector<std::string> >::Read(std::vector<std::string> &values);

#endif

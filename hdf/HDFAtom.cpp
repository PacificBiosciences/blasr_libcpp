#include "HDFAtom.hpp"

template<>
void HDFAtom<std::string>::Create(H5::H5Location &object, const std::string & atomName) {
    H5::StrType strType(0, H5T_VARIABLE);
    hsize_t defaultDims[] = {1};
    H5::DataSpace defaultDataSpace(1, defaultDims);
    attribute = object.createAttribute(atomName.c_str(), strType, H5::DataSpace(H5S_SCALAR));
    initialized= true;
}


#define DEFINE_TYPED_CREATE_ATOM(T, predtype) template<> \
	void HDFAtom<T>::TypedCreate(H5::H5Location &object, const std::string &atomName, H5::DataSpace &defaultDataSpace) {				\
  attribute = object.createAttribute(atomName.c_str(), (predtype), defaultDataSpace );	\
}


DEFINE_TYPED_CREATE_ATOM(int, H5::PredType::NATIVE_INT)
DEFINE_TYPED_CREATE_ATOM(unsigned int, H5::PredType::NATIVE_UINT)
DEFINE_TYPED_CREATE_ATOM(unsigned char, H5::PredType::NATIVE_UINT8)
DEFINE_TYPED_CREATE_ATOM(char, H5::PredType::NATIVE_INT8)
DEFINE_TYPED_CREATE_ATOM(float, H5::PredType::NATIVE_FLOAT)
DEFINE_TYPED_CREATE_ATOM(uint64_t, H5::PredType::STD_I64LE)


template<>
void HDFAtom<std::vector<int> >::Write(const std::vector<int> vect) {
  hsize_t  length = vect.size();
  H5::DataType baseType = H5::PredType::NATIVE_INT;
  H5::ArrayType arrayDataType(baseType, 1, &length);
  attribute.write(arrayDataType, &((vect)[0]));
}


template<>
void HDFAtom<std::string>::Write(std::string value) {
  //	H5::StrType strType(0, value.size());
	H5::StrType strType(0, H5T_VARIABLE);
	attribute.write(strType, std::string(value.c_str()));
}

template<>
void HDFAtom<uint64_t>::Write(uint64_t value) {
	attribute.write( H5::PredType::STD_I64LE, &value);
}

template<>
void HDFAtom<int>::Write(int value) {
	attribute.write( H5::PredType::NATIVE_INT, &value);
}

template<>
void HDFAtom<unsigned int>::Write(unsigned int value) {
	attribute.write( H5::PredType::NATIVE_INT, &value);
}

template<>
void HDFAtom<unsigned char>::Write(unsigned char value) {
	attribute.write( H5::PredType::NATIVE_UINT8, &value);
}

template<>
void HDFAtom<char>::Write(char value) {
	attribute.write( H5::PredType::NATIVE_INT8, &value);
}

template<>
void HDFAtom<float>::Write(float value) {
	attribute.write( H5::PredType::NATIVE_FLOAT, &value);
}

template<>
void HDFAtom<std::string>::Read(std::string &value) {
	/*
	 * Read in a std::string that has been stored either as an array or a
	 * variable length std::string.  To decide which, query the
	 * isVariableStr() option.
	 */
	H5::StrType stringType = attribute.getStrType();
	bool stringIsVariableLength = stringType.isVariableStr();
	if (stringIsVariableLength) 
		attribute.read(stringType, value);
	else {
		hsize_t stsize = attribute.getStorageSize();
		value.resize(stsize);
		//		char *valueStr = new char[stsize+1];
		attribute.read(stringType, &value[0]);
		if (stsize > 0 and value[stsize-1] == '\0') {
			value.resize(stsize-1);
		}
	}
}

template<>
void HDFAtom<int>::Read(int &value) {
	H5::DataType intType(H5::PredType::NATIVE_INT);
	attribute.read(intType, &value);
}

template<>
void HDFAtom<uint64_t>::Read(uint64_t &value) {
	H5::DataType intType(H5::PredType::STD_I64LE);
	attribute.read(intType, &value);
}

template<>
void HDFAtom<unsigned int>::Read(unsigned int &value) {
	H5::DataType uintType(H5::PredType::NATIVE_UINT);
	attribute.read(uintType, &value);
}

template<>
void HDFAtom<float>::Read(float &value) {
	H5::DataType type(H5::PredType::NATIVE_FLOAT);
	attribute.read(type, &value);
}

template<>
void HDFAtom<std::vector<std::string> >::Read(std::vector<std::string> &values) {
	std::string value;
	/*
	 * This attribute is an array of std::strings. They are read in by
	 * storing pointers to std::strings in memory controlled by HDF.  To read
	 * the std::strings, read the pointers into a temporary array, then copy
	 * those std::strings to the values array. This way when the values array
	 * is destroyed, it will not try and get rid of space that is under
	 * HDF control.
	 */
    H5::DataSpace attributeSpace = attribute.getSpace();
	hsize_t nPoints;
	nPoints = attributeSpace.getSelectNpoints();
	H5::DataType attrType = attribute.getDataType(); // necessary for attr.read()

	// Declare and initialize std::vector of pointers to std::string attribute list.
	std::vector<char*> ptrsToHDFControlledMemory;
	ptrsToHDFControlledMemory.resize(nPoints);
	// Copy the pointers.
	attribute.read(attrType, &ptrsToHDFControlledMemory[0]);
	// Copy the std::strings into memory the main program has control over.
	unsigned int i;
	for (i = 0; i < ptrsToHDFControlledMemory.size(); i++ ){
		values.push_back(ptrsToHDFControlledMemory[i]);
		free(ptrsToHDFControlledMemory[i]);
	}
}

// Explict instantiate classes
template class HDFAtom<int>;
template class HDFAtom<unsigned int>;
template class HDFAtom<char>;
template class HDFAtom<unsigned char>;
template class HDFAtom<uint64_t>;
template class HDFAtom<float>;
template class HDFAtom< std::string >;
template class HDFAtom< std::vector <int> >;
template class HDFAtom< std::vector <unsigned int> >;
template class HDFAtom< std::vector <string> >;


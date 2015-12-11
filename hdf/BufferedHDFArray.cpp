#include "BufferedHDFArray.hpp"

/*
 *
 * Implementation of a 1-D array for IO from an HDF array.
 * This is templated, but specialized for a few data types, so that 
 * the HDF data types do not need to be specified by anybody.
 *
 *  Two examples of the usage of this class follow:
 *
 *	HDFArray<int> nElemArray;
 * 	nElemArray.Initialize(hdfFile, "PulseData/BaseCalls/ZMW/NumEvent");
 *  nElemArray.Read(i, i+1, &nElem);
 *
 * 	HDFArray<unsigned char> qualArray;
 *	qualArray.Initialize(hdfFile, "PulseData/BaseCalls/QualityValue");
 *  qualArray.Read(cur, cur + nElem, qual);
 *
 */

//
// Type specializations for some standard types. Use the macro for
// vanilla specializations (that only require the HDF type ID to be
// specified). 
//
#define DEFINE_TYPED_READ_ARRAY(T, Pred) template<>  \
   	void BufferedHDFArray<T>::Read(DSLength start, DSLength end, T* dest) { \
   	Read(start,end, Pred, dest); \
	}

DEFINE_TYPED_READ_ARRAY(int, H5::PredType::NATIVE_INT)
DEFINE_TYPED_READ_ARRAY(char, H5::PredType::NATIVE_INT8)
DEFINE_TYPED_READ_ARRAY(unsigned char, H5::PredType::NATIVE_UINT8)
DEFINE_TYPED_READ_ARRAY(unsigned int, H5::PredType::NATIVE_UINT)
DEFINE_TYPED_READ_ARRAY(uint16_t, H5::PredType::NATIVE_UINT16)
DEFINE_TYPED_READ_ARRAY(float, H5::PredType::NATIVE_FLOAT)
DEFINE_TYPED_READ_ARRAY(char*, H5::PredType::C_S1)


#define DEFINE_TYPED_READ_DATASET(T, Pred) template<> \
	void BufferedHDFArray<T>::ReadDataset(std::vector<T> &dest) { \
	dest.resize(arrayLength); \
  Read(0, arrayLength, Pred, &dest[0]); \
}

DEFINE_TYPED_READ_DATASET(int, H5::PredType::NATIVE_INT)
DEFINE_TYPED_READ_DATASET(char, H5::PredType::NATIVE_INT8)
DEFINE_TYPED_READ_DATASET(unsigned char, H5::PredType::NATIVE_UINT8)
DEFINE_TYPED_READ_DATASET(unsigned int, H5::PredType::NATIVE_UINT)
DEFINE_TYPED_READ_DATASET(uint16_t, H5::PredType::NATIVE_UINT16)
DEFINE_TYPED_READ_DATASET(float, H5::PredType::NATIVE_FLOAT)
DEFINE_TYPED_READ_DATASET(char*, H5::PredType::C_S1)


#define DEFINE_TYPED_CREATE_ARRAY(T, Pred) template<> \
	void BufferedHDFArray<T>::TypedCreate(H5::DataSpace &fileSpace, H5::DSetCreatPropList &cparms) { \
	T zero; zero = 0;\
	cparms.setFillValue(Pred, &zero);\
	dataset = container->createDataSet(datasetName.c_str(), Pred, fileSpace, cparms); \
}

DEFINE_TYPED_CREATE_ARRAY(int, H5::PredType::NATIVE_INT)
DEFINE_TYPED_CREATE_ARRAY(char, H5::PredType::NATIVE_INT8)
DEFINE_TYPED_CREATE_ARRAY(char*, H5::StrType(0,H5T_VARIABLE))
DEFINE_TYPED_CREATE_ARRAY(unsigned char, H5::PredType::NATIVE_UINT8)
DEFINE_TYPED_CREATE_ARRAY(unsigned int, H5::PredType::NATIVE_UINT)
DEFINE_TYPED_CREATE_ARRAY(float, H5::PredType::NATIVE_FLOAT)
DEFINE_TYPED_CREATE_ARRAY(uint16_t, H5::PredType::NATIVE_UINT16)

template<>
void BufferedHDFArray<string>::TypedCreate(H5::DataSpace &space, H5::DSetCreatPropList &cparms) {
    H5::StrType varStrType(0,H5T_VARIABLE);
    dataset = container->createDataSet(datasetName.c_str(), varStrType, space, cparms);
}


#define DEFINE_TYPED_WRITE_ARRAY(T, Pred) template<> \
void BufferedHDFArray<T>::TypedWrite(const T *data, \
    const H5::DataSpace &memorySpace, const H5::DataSpace &fileSpace) { \
    dataset.write(data, Pred, memorySpace, fileSpace); \
}

DEFINE_TYPED_WRITE_ARRAY(int, H5::PredType::NATIVE_INT)
DEFINE_TYPED_WRITE_ARRAY(unsigned int, H5::PredType::NATIVE_UINT)
DEFINE_TYPED_WRITE_ARRAY(unsigned char, H5::PredType::NATIVE_UINT8)
DEFINE_TYPED_WRITE_ARRAY(char, H5::PredType::NATIVE_INT8)
DEFINE_TYPED_WRITE_ARRAY(float, H5::PredType::NATIVE_FLOAT)
DEFINE_TYPED_WRITE_ARRAY(uint16_t, H5::PredType::NATIVE_UINT16)
DEFINE_TYPED_WRITE_ARRAY(char*, H5::StrType(0,H5T_VARIABLE))
DEFINE_TYPED_WRITE_ARRAY(std::string, H5::StrType(0,H5T_VARIABLE))



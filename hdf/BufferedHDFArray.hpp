#ifndef _BLASR_BUFFERED_HDF_ARRAY_HPP_
#define _BLASR_BUFFERED_HDF_ARRAY_HPP_

#include <cassert>
#include <cstdint>
#include <string>
#include <vector>

// HDF5 library includes
#include "hdf5.h"
#include "H5Cpp.h"

#include "Types.h"
#include "HDFConfig.hpp"
#include "HDFData.hpp"
#include "HDFGroup.hpp"
#include "HDFWriteBuffer.hpp"
#include "HDFFile.hpp"
#include "DNASequence.hpp"
#include "FASTQSequence.hpp"

/*
 *
 * Implementation of a 1-D array for IO from an HDF array.
 * This is templated, but specialized for a few data types, so that 
 * the HDF data types do not need to be specified by anybody.
 *
 *  Two examples of the usage of this class follow:
 *
 *  HDFArray<int> nElemArray;
 *  nElemArray.Initialize(hdfFile, "PulseData/BaseCalls/ZMW/NumEvent");
 *  nElemArray.Read(i, i+1, &nElem);
 *
 *  HDFArray<unsigned char> qualArray;
 *  qualArray.Initialize(hdfFile, "PulseData/BaseCalls/QualityValue");
 *  qualArray.Read(cur, cur + nElem, qual);
 *
 */

template<typename T>
class BufferedHDFArray : public HDFData, public HDFWriteBuffer<T> {
public:
    hsize_t   nDims;
    hsize_t   *dimSize;
    int       maxDims;
    DSLength  arrayLength;

    /*
     * Constructor meant to be used for data that will be written.  
     * This allocates the write buffer.
     */
    BufferedHDFArray(int pBufferSize=32768);

    BufferedHDFArray(H5::CommonFG* _container, std::string _datasetName);

    ~BufferedHDFArray(); 

    void SetBufferSize(int _bufferSize); 

    void Write(const T *data, DSLength dataLength, bool append=true, 
        DSLength writePos = 0);

    void Flush(bool append=true, DSLength writePos = 0); 

    void TypedWrite(const char **data, const H5::DataSpace &memorySpace, 
        const H5::DataSpace &extendedSpace); 

    void TypedWrite(const T*data, const H5::DataSpace &memorySpace, 
        const H5::DataSpace &extendedSpace); 

    void TypedCreate(H5::DataSpace &fileSpace, H5::DSetCreatPropList &cparms); 

    void Create(HDFGroup &parentGroup, std::string _datasetName); 

    void Create(H5::CommonFG* _container, std::string _datasetName); 

    /*
     * Initialize for reading.
     *
     * Open a dataset in an hdf file. Only call this on datasets that
     * exist, since this currently handles errors with opening datasets
     * by ungracefully exiting the program. 
     */

    int InitializeForReading(HDFGroup &parentGroup, 
        const std::string datasetName); 

    int Initialize(HDFGroup &parentGroup, const std::string &datasetName); 

    int Initialize(HDFGroup &parentGroup, const std::string &datasetName,
        bool createIfMissing, const DSLength newArrayLength=0); 

    int UpdateH5Dataspace(); 

    int Resize(const DSLength newArrayLength); 

    void Close(); 

    DSLength size(); 

    /*
     * Unspecialized form of read.
     * Read cannot be called on a type T* that does not have a
     * specialized template definition.  This is all determined at
     * compile time.  To ensure this, the following
     * default definition is provided that gives a nasty warning and
     * exits the code.
     */

    virtual void Read(DSLength start, DSLength end, T* dest); 

    /*
     * Read in type T from the opened dataset from the interval (start,
     * end].
     */
    void ReadDataset(std::vector<T> &dest); 

    void Read(DSLength start, DSLength end, H5::DataType typeID, T* dest); 

    void ReadCharArray(DSLength start, DSLength end, std::string* dest); 
};

/*
 * Type specializations for some standard types. Use the macro for
 * vanilla specializations (that only require the HDF type ID to be
 * specified). 
 */
#define DECLARE_TYPED_READ_ARRAY(T, Pred) template<>  \
void BufferedHDFArray<T>::Read(DSLength start, DSLength end, T* dest);

DECLARE_TYPED_READ_ARRAY(int, H5::PredType::NATIVE_INT)
DECLARE_TYPED_READ_ARRAY(char, H5::PredType::NATIVE_INT8)
DECLARE_TYPED_READ_ARRAY(unsigned char, H5::PredType::NATIVE_UINT8)
DECLARE_TYPED_READ_ARRAY(unsigned int, H5::PredType::NATIVE_UINT)
DECLARE_TYPED_READ_ARRAY(uint16_t, H5::PredType::NATIVE_UINT16)
DECLARE_TYPED_READ_ARRAY(float, H5::PredType::NATIVE_FLOAT)
DECLARE_TYPED_READ_ARRAY(char*, H5::PredType::C_S1)


#define DECLARE_TYPED_READ_DATASET(T, Pred) template<>  \
void BufferedHDFArray<T>::ReadDataset(vector<T>  &dest);

DECLARE_TYPED_READ_DATASET(int, H5::PredType::NATIVE_INT)
DECLARE_TYPED_READ_DATASET(char, H5::PredType::NATIVE_INT8)
DECLARE_TYPED_READ_DATASET(unsigned char, H5::PredType::NATIVE_UINT8)
DECLARE_TYPED_READ_DATASET(unsigned int, H5::PredType::NATIVE_UINT)
DECLARE_TYPED_READ_DATASET(uint16_t, H5::PredType::NATIVE_UINT16)
DECLARE_TYPED_READ_DATASET(float, H5::PredType::NATIVE_FLOAT)
DECLARE_TYPED_READ_DATASET(char*, H5::PredType::C_S1)


#define DECLARE_TYPED_CREATE_ARRAY(T, Pred) template <> \
void BufferedHDFArray<T>::TypedCreate(H5::DataSpace &fileSpace,  \
    H5::DSetCreatPropList &cparms);

DECLARE_TYPED_CREATE_ARRAY(int, H5::PredType::NATIVE_INT)
DECLARE_TYPED_CREATE_ARRAY(char, H5::PredType::NATIVE_INT8)
DECLARE_TYPED_CREATE_ARRAY(char*, H5::StrType(0,H5T_VARIABLE))
DECLARE_TYPED_CREATE_ARRAY(unsigned char, H5::PredType::NATIVE_UINT8)
DECLARE_TYPED_CREATE_ARRAY(unsigned int, H5::PredType::NATIVE_UINT)
DECLARE_TYPED_CREATE_ARRAY(float, H5::PredType::NATIVE_FLOAT)
DECLARE_TYPED_CREATE_ARRAY(uint16_t, H5::PredType::NATIVE_UINT16)

template<>
void BufferedHDFArray<std::string>::TypedCreate(H5::DataSpace &space, H5::DSetCreatPropList &cparms);

#define DECLARE_TYPED_WRITE_ARRAY(T, Pred) template <> \
void BufferedHDFArray<T>::TypedWrite(const T *data, \
    const H5::DataSpace &memorySpace, const H5::DataSpace &fileSpace);

DECLARE_TYPED_WRITE_ARRAY(int, H5::PredType::NATIVE_INT)
DECLARE_TYPED_WRITE_ARRAY(unsigned int, H5::PredType::NATIVE_UINT)
DECLARE_TYPED_WRITE_ARRAY(unsigned char, H5::PredType::NATIVE_UINT8)
DECLARE_TYPED_WRITE_ARRAY(char, H5::PredType::NATIVE_INT8)
DECLARE_TYPED_WRITE_ARRAY(float, H5::PredType::NATIVE_FLOAT)
DECLARE_TYPED_WRITE_ARRAY(uint16_t, H5::PredType::NATIVE_UINT16)
DECLARE_TYPED_WRITE_ARRAY(char*, H5::StrType(0,H5T_VARIABLE))
DECLARE_TYPED_WRITE_ARRAY(std::string, H5::StrType(0,H5T_VARIABLE))

//
// This is a nonstandard definition because it requires the creation
// of a special datatype for variable length string type.
//
//
// Use this as the base class for other lists.
//
typedef BufferedHDFArray<void> BaseHDFArray;

#define DEFINE_TYPED_CLASS(CLASS_NAME, TEMPLATE_TYPE) \
class CLASS_NAME : public BufferedHDFArray< TEMPLATE_TYPE  > { \
public: \
    void Read(DSLength start, DSLength end, TEMPLATE_TYPE *dest) { \
        BufferedHDFArray<TEMPLATE_TYPE>::Read(start, end, dest); \
    } \
}; 

DEFINE_TYPED_CLASS(HDFIntArray, int)
DEFINE_TYPED_CLASS(HDFUIntArray, unsigned int)
DEFINE_TYPED_CLASS(HDFUCharArray, unsigned char)
DEFINE_TYPED_CLASS(HDFCharArray, char)
DEFINE_TYPED_CLASS(HDFUShortArray, uint16_t )
//DEFINE_TYPED_CLASS(HDFStringArray, std::string)
DEFINE_TYPED_CLASS(HDFFloatArray, float)

#include "BufferedHDFArrayImpl.hpp"

#endif

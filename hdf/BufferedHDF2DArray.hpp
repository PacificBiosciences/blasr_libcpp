#ifndef _BLASR_HDF_BUFFERED_HDF_2D_ARRAY_HPP_
#define _BLASR_HDF_BUFFERED_HDF_2D_ARRAY_HPP_

#include <cstdlib>
#include <iostream>
#include <string>

#include "H5Cpp.h"

#include "Types.h"
#include "HDFConfig.hpp"
#include "HDFData.hpp"
#include "HDFGroup.hpp"
#include "HDFWriteBuffer.hpp"

/*
 *
 * Implementation of a 2-D array for IO from an HDF array.
 * This is templated, but specialized for a few data types, so that 
 * the HDF data types do not need to be specified by somebody when reading.
 *
 * Currently no support exists for reading non-contiguous blocks of data, and
 * the main intended use is to read in increments of rows.

 int main(int argc, char* argv[]) {
	if (argc < 1) {
		cout << "usage: testHDFReading hdfFile" << endl;
		exit(0);
	}

	std::string hdfFileName = argv[1];
	
	H5File hdfFile;
	hdfFile.openFile(hdfFileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	BufferedHDF2DArray<uint16_t> xyArray;
	xyArray.Initialize(hdfFile, "PulseData/BaseCalls/ZMW/HoleXY");
	int curX = 0;
	xyArray.Read(curX, curX + 1, 0, 2, holeXY);

	or, to read a row:
	xyArray.Read(curX, curX+1, holeXY);

 *
 */
template<typename T>
class BufferedHDF2DArray : public HDFData, public HDFWriteBuffer<T> {

private:
    hsize_t   nDims;
    hsize_t   *dimSize;
    DSLength  maxDims;
    DSLength  rowLength,  colLength;

public:

    BufferedHDF2DArray(H5::CommonFG *_container, std::string _datasetName); 

    BufferedHDF2DArray();

    DSLength GetNRows();

    DSLength GetNCols(); 

    void Close(); 

    ~BufferedHDF2DArray(); 

    int InitializeForReading(HDFGroup &group, std::string datasetName);

    /*
     * Initialize HDF2D for reading.  No write buffer initialization is
     * required.  The assumption is that the dataspace is in two
     * dimensions, and this exits without grace if it is not. 
     */
    int Initialize(HDFGroup &group, std::string datasetName, DSLength _rowLength=0, 
        int _bufferSize=0, bool createIfMissing=true); 

    DSLength size(); 

    /*
     * Read rows in the range (startX, endX] in to dest.
     */
    void Read(DSLength startX, DSLength endX, H5::DataType typeID, T*dest); 

    void Read(DSLength startX, DSLength endX, T*dest); 

    /*
     * This is the non-specialized definition.  Since this should only
     * operate on specialized types, report an error and bail.
     */
    void Read(DSLength startX, DSLength endX, DSLength startY, DSLength endY, T* dest); 

    void Read(DSLength startX, DSLength endX, DSLength startY, DSLength endY, H5::DataType typeID, 
        T *dest); 

    void Create(H5::CommonFG *_container, std::string _datasetName, DSLength _rowLength); 

    void TypedCreate(H5::DataSpace &fileSpace, H5::DSetCreatPropList &cparms); 
    
    // Append
    void TypedWriteRow(const T*, const H5::DataSpace &memoryDataSpace, 
        const H5::DataSpace &fileDataSpace); 

    /*
     * This code is copied directly form BufferedHDFArray.  I'm not sure
     * how to set up the objects nicely to share the code between the
     * two since the Flush() function is different.  There probably is a
     * design pattern or simply better way to engineer this, but for now
     * it's 15 lines of code.
     */
    void WriteRow(const T *data, DSLength dataLength, DSLength destRow=static_cast<DSLength>(-1)); 

    void Flush(DSLength destRow = static_cast<DSLength>(-1)); 
};

DSLength GetDatasetNDim(H5::CommonFG &parentGroup, std::string datasetName);

#define DECLARE_TYPED_WRITE_ROW(T, Pred) template<>\
void BufferedHDF2DArray<T>::TypedWriteRow(const T *data, \
    const H5::DataSpace &memorySpace, const H5::DataSpace &fileSpace) ;


DECLARE_TYPED_WRITE_ROW(int, H5::PredType::NATIVE_INT)
DECLARE_TYPED_WRITE_ROW(unsigned int, H5::PredType::NATIVE_UINT)
DECLARE_TYPED_WRITE_ROW(unsigned char, H5::PredType::NATIVE_UINT8)
DECLARE_TYPED_WRITE_ROW(uint16_t, H5::PredType::NATIVE_UINT16)
DECLARE_TYPED_WRITE_ROW(int16_t, H5::PredType::NATIVE_INT16)
DECLARE_TYPED_WRITE_ROW(float, H5::PredType::NATIVE_FLOAT)


#define DECLARE_TYPED_READ_ROW(T, Pred) template<>\
void BufferedHDF2DArray<T>::Read(DSLength startX, DSLength endX, DSLength startY, DSLength endY, T* dest);


DECLARE_TYPED_READ_ROW(int, H5::PredType::NATIVE_INT)
DECLARE_TYPED_READ_ROW(unsigned int, H5::PredType::NATIVE_UINT)
DECLARE_TYPED_READ_ROW(char, H5::PredType::NATIVE_INT8)
DECLARE_TYPED_READ_ROW(unsigned char, H5::PredType::NATIVE_UINT8)
DECLARE_TYPED_READ_ROW(uint16_t, H5::PredType::NATIVE_UINT16)
DECLARE_TYPED_READ_ROW(int16_t, H5::PredType::NATIVE_INT16)
DECLARE_TYPED_READ_ROW(float, H5::PredType::NATIVE_FLOAT)

#define DECLARE_TYPED_CREATE_ROW(T, Pred)template<>\
void BufferedHDF2DArray<T>::TypedCreate(H5::DataSpace &fileSpace, \
    H5::DSetCreatPropList &cparms) ;

DECLARE_TYPED_CREATE_ROW(int, H5::PredType::NATIVE_INT)
DECLARE_TYPED_CREATE_ROW(unsigned int, H5::PredType::NATIVE_UINT)
DECLARE_TYPED_CREATE_ROW(char, H5::PredType::NATIVE_INT8)
DECLARE_TYPED_CREATE_ROW(unsigned char, H5::PredType::NATIVE_UINT8)
DECLARE_TYPED_CREATE_ROW(uint16_t, H5::PredType::NATIVE_UINT16)
DECLARE_TYPED_CREATE_ROW(int16_t, H5::PredType::NATIVE_INT16)
DECLARE_TYPED_CREATE_ROW(float, H5::PredType::NATIVE_FLOAT)

#include "BufferedHDF2DArrayImpl.hpp"

#endif

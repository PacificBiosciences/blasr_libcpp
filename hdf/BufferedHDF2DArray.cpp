#include "BufferedHDF2DArray.hpp"
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

	string hdfFileName = argv[1];
	
	H5File hdfFile;
	hdfFile.openFile(hdfFileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	BufferedHDF2DArray<uint16_t> xyArray;
	xyArray.Initialize(hdfFile, "PulseData/BaseCalls/ZMW/HoleXY");
	DSLength curX = 0;
	xyArray.Read(curX, curX + 1, 0, 2, holeXY);

	or, to read a row:
	xyArray.Read(curX, curX+1, holeXY);

 */

DSLength GetDatasetNDim(H5::CommonFG &parentGroup, std::string datasetName) {
    HDFData tmpDataset;
    tmpDataset.InitializeDataset(parentGroup, datasetName);
    H5::DataSpace dataspace = tmpDataset.dataset.getSpace();
    DSLength nDims = dataspace.getSimpleExtentNdims();
    dataspace.close();
    tmpDataset.dataset.close();
    return nDims;
}

#define DEFINE_TYPED_WRITE_ROW(T, Pred) template<>\
void BufferedHDF2DArray<T>::TypedWriteRow(const T *data, const H5::DataSpace &memorySpace, const H5::DataSpace &fileSpace) {\
	dataset.write(data, Pred, memorySpace, fileSpace);\
}

DEFINE_TYPED_WRITE_ROW(int, H5::PredType::NATIVE_INT)
DEFINE_TYPED_WRITE_ROW(unsigned int, H5::PredType::NATIVE_UINT)
DEFINE_TYPED_WRITE_ROW(unsigned char, H5::PredType::NATIVE_UINT8)
DEFINE_TYPED_WRITE_ROW(uint16_t, H5::PredType::NATIVE_UINT16)
DEFINE_TYPED_WRITE_ROW(int16_t, H5::PredType::NATIVE_INT16)
DEFINE_TYPED_WRITE_ROW(float, H5::PredType::NATIVE_FLOAT)


#define DEFINE_TYPED_READ_ROW(T, Pred) template<>\
void BufferedHDF2DArray<T>::Read(DSLength startX, DSLength endX, DSLength startY, DSLength endY, T* dest) {\
	Read(startX, endX, startY, endY, Pred, dest);\
}

DEFINE_TYPED_READ_ROW(int, H5::PredType::NATIVE_INT)
DEFINE_TYPED_READ_ROW(unsigned int, H5::PredType::NATIVE_UINT)
DEFINE_TYPED_READ_ROW(char, H5::PredType::NATIVE_INT8)
DEFINE_TYPED_READ_ROW(unsigned char, H5::PredType::NATIVE_UINT8)
DEFINE_TYPED_READ_ROW(uint16_t, H5::PredType::NATIVE_UINT16)
DEFINE_TYPED_READ_ROW(int16_t, H5::PredType::NATIVE_INT16)
DEFINE_TYPED_READ_ROW(float, H5::PredType::NATIVE_FLOAT)

#define DEFINE_TYPED_CREATE_ROW(T, Pred) template<>\
void BufferedHDF2DArray<T>::TypedCreate(H5::DataSpace &fileSpace, H5::DSetCreatPropList &cparms) {\
	dataset = container->createDataSet(datasetName.c_str(), Pred, fileSpace, cparms);\
}

DEFINE_TYPED_CREATE_ROW(int, H5::PredType::NATIVE_INT)
DEFINE_TYPED_CREATE_ROW(unsigned int, H5::PredType::NATIVE_UINT)
DEFINE_TYPED_CREATE_ROW(char, H5::PredType::NATIVE_INT8)
DEFINE_TYPED_CREATE_ROW(unsigned char, H5::PredType::NATIVE_UINT8)
DEFINE_TYPED_CREATE_ROW(uint16_t, H5::PredType::NATIVE_UINT16)
DEFINE_TYPED_CREATE_ROW(int16_t, H5::PredType::NATIVE_INT16)
DEFINE_TYPED_CREATE_ROW(float, H5::PredType::NATIVE_FLOAT)


#ifndef _BLASR_HDF_BUFFERED_HDF_2D_ARRAY_IMPL_HPP_
#define _BLASR_HDF_BUFFERED_HDF_2D_ARRAY_IMPL_HPP_

#include <cstring>
#include <cassert>
#include "utils.hpp"

template<typename T>
BufferedHDF2DArray<T>::BufferedHDF2DArray(H5::CommonFG *_container, 
    std::string _datasetName) : HDFData(_container, _datasetName) {}

template<typename T>
BufferedHDF2DArray<T>::BufferedHDF2DArray() : HDFData() {
    maxDims = 0;
    nDims = 2;
    dimSize =NULL;
    rowLength = -1;
    colLength = -1;
}

template<typename T>
DSLength BufferedHDF2DArray<T>::GetNRows() {
    return rowLength;
}

template<typename T>
DSLength BufferedHDF2DArray<T>::GetNCols() {
    return colLength;
}

template<typename T>
void BufferedHDF2DArray<T>::Close() {

    //
    // Clean up the write buffer.
    //
    //		Flush();
    if (dimSize != NULL) {
        delete[] dimSize;
        dimSize = NULL;
    }
    this->HDFWriteBuffer<T>::Free();
}

template<typename T>
BufferedHDF2DArray<T>::~BufferedHDF2DArray() {
    Close();
}

template<typename T>
int BufferedHDF2DArray<T>::InitializeForReading(HDFGroup& group, std::string datasetName)
{
    return Initialize(group, datasetName, 0, 0, false);
}

/*
 * Initialize HDF2D for reading.  No write buffer initialization is
 * required.  The assumption is that the dataspace is in two
 * dimensions, and this exits without grace if it is not. 
 */
template<typename T>
int BufferedHDF2DArray<T>::Initialize(HDFGroup &group, std::string datasetName,
    DSLength _rowLength, int _bufferSize, bool createIfMissing) {

    bool groupContainsDataset = group.ContainsObject(datasetName);
    if (groupContainsDataset == false) {
        //
        // Do some error checking.
        //
        if (createIfMissing == false) {
            std::cout << "ERROR! Could not open dataset " << datasetName << std::endl;
            exit(1);
        }
        if (_rowLength == 0) {
            std::cout << "ERROR!  Improper usage of BufferedHDF2DArray::Initialize.  The 2D Array "<<std::endl
                << "is being created but is given a number of columns of 0." << std::endl;
            exit(1);
        }
        Create(&group.group, datasetName, _rowLength);
    }
    else {
        InitializeDataset(group.group, datasetName);
        try {
            dataspace = dataset.getSpace();
        }
        catch(H5::DataSetIException &e) { 
            std::cout << e.getDetailMsg() << std::endl;
            exit(1);
        }

        maxDims   = MAX_DIMS;
        try {
            nDims     = dataspace.getSimpleExtentNdims();
            /*
             * Prevent abuse of this class for multidimensional IO.
             */
            if (nDims != 2) {
                std::cout << "ERROR in HDF format: dataset: " << datasetName << " should be 1-D, but it is not." << std::endl;
                exit(1);
            }

            /*
             * Load in the size of this dataset, and make a map to the whole thing.
             */
            if (dimSize) {
                delete [] dimSize;
            }
            dimSize = ProtectedNew<hsize_t>(nDims);
            dataspace.getSimpleExtentDims(dimSize);
            rowLength = dimSize[0];
            colLength = dimSize[1];
            if (rowLength == 0) {
                dataspace.close();
                return 1;
            }
            fullSourceSpace = H5::DataSpace(2, dimSize);
            dataspace.close();
        }
        catch(H5::Exception &e) {
            std::cout << e.getDetailMsg() << std::endl;
            exit(1);
        }
    }
    return 1;
}

template<typename T>
DSLength BufferedHDF2DArray<T>::size() {
    dataspace.getSimpleExtentDims(dimSize);
    return dimSize[0];
}

/*
 * Read rows in the range (startX, endX] in to dest.
 */
template<typename T>
void BufferedHDF2DArray<T>::Read(DSLength startX, DSLength endX, H5::DataType typeID, T*dest) {
    Read(startX, endX, 0, dimSize[1], typeID, dest);
}

template<typename T>
void BufferedHDF2DArray<T>::Read(DSLength startX, DSLength endX, T*dest) {
    Read(startX, endX, 0, dimSize[1], dest);
}
/*
 * This is the non-specialized definition.  Since this should only
 * operate on specialized types, report an error and bail.
 */
template<typename T>
void BufferedHDF2DArray<T>::Read(DSLength startX, DSLength endX, DSLength startY, DSLength endY, T* dest) {
    assert("ERROR, calling Read with an unsupported type. Use Read(startx,endx, starty,endy,datatype, dest) instead." == 0);
    exit(1);
}

template<typename T>
void BufferedHDF2DArray<T>::Read(DSLength startX, DSLength endX, DSLength startY, DSLength endY, H5::DataType typeID, T *dest) {
    hsize_t memSpaceSize[2] = {0, 0};
    memSpaceSize[0] = endX - startX;
    memSpaceSize[1] = endY - startY;
    hsize_t sourceSpaceOffset[2] = {0, 0};
    sourceSpaceOffset[0] = startX;
    sourceSpaceOffset[1] = startY;

    H5::DataSpace destSpace(2, memSpaceSize);		
    fullSourceSpace.selectHyperslab(H5S_SELECT_SET, memSpaceSize, sourceSpaceOffset);
    dataset.read(dest, typeID, destSpace, fullSourceSpace);
    destSpace.close();
}

template<typename T>
void BufferedHDF2DArray<T>::Create(H5::CommonFG *_container, string _datasetName, DSLength _rowLength) {
    container   = _container;
    datasetName = _datasetName;
    rowLength   = _rowLength;
    //
    // Make life easy if the buffer is too small to fit a row --
    // resize it so that rows may be copied and written out in an
    // atomic unit.
    //
    if (this->bufferSize < rowLength) {
        // When the buffer size is greater than 0, the write buffer
        // should exist.
        if (this->bufferSize > 0) {
            assert(this->writeBuffer != NULL);
            delete[] this->writeBuffer;
        }
        this->writeBuffer = ProtectedNew<T>(rowLength);
        this->bufferSize = rowLength;
    }

    hsize_t dataSize[2]    = {0, hsize_t(rowLength)};
    hsize_t maxDataSize[2] = {H5S_UNLIMITED, hsize_t(rowLength)};
    H5::DataSpace fileSpace(2, dataSize, maxDataSize);
    H5::DSetCreatPropList cparms;

    /*
     * For some reason, chunking must be enabled when creating a dataset
     * that  has an unlimited dimension.  Of course, this is not
     * mentioned in the hdf5 c++ documentation, because that
     * docuemntation was written for people who enjoy learning how to
     * use an API by reading comments in source code.
     */
    hsize_t chunkDims[2] = {16384, hsize_t(rowLength)};
    cparms.setChunk( 2, chunkDims );
    TypedCreate(fileSpace, cparms);
    fileSpace.close();

    //
    // Set some flags that indicate this dataset is ready for writing.
    //
    fileDataSpaceInitialized = true;
    isInitialized = true;
}

template<typename T>
void BufferedHDF2DArray<T>::TypedCreate(H5::DataSpace &fileSpace, 
    H5::DSetCreatPropList &cparms) {

    assert("Error, calling HDF2DArray<T>::TypedCreate on an unsupported type.  A specialization must be written in HDF2DArray.h" == 0);
}

// Append
template<typename T>
void TypedWriteRow(const T*, const H5::DataSpace &memoryDataSpace, 
    const H5::DataSpace &fileDataSpace) {

    assert("Error, calling HDF2DArray<T>::TypedWriteRow on an unsupported type.  A specialization must be written in HDF2DArray.h" == 0);
}


/*
 * This code is copied directly form BufferedHDFArray.  I'm not sure
 * how to set up the objects nicely to share the code between the
 * two since the Flush() function is different.  There probably is a
 * design pattern or simply better way to engineer this, but for now
 * it's 15 lines of code.
 */

template<typename T>
void BufferedHDF2DArray<T>::WriteRow(const T *data, DSLength dataLength, DSLength destRow) {
    // Fill the buffer with data. When there is overflow, write
    // that out to disk.
    //
    DSLength dataIndex = 0;
    int bufferCapacity;
    int bufferFillSize = 0;
    bool flushBuffer;
    while(dataIndex < dataLength) {
        //
        // Compute the capacity of this buffer to fit an integral number
        // of rows into it.
        //
        bufferCapacity = (this->bufferSize / rowLength)*rowLength - this->bufferIndex;
        flushBuffer = false;
        if (bufferCapacity  > dataLength - dataIndex) {
            bufferFillSize = dataLength - dataIndex;
        }
        else {
            bufferFillSize = bufferCapacity;
            flushBuffer = true;
        }
        memcpy((void*) &this->writeBuffer[this->bufferIndex], (void*) &data[dataIndex], sizeof(T)*bufferFillSize);
        dataIndex   += bufferFillSize;
        this->bufferIndex += bufferFillSize;
        if (flushBuffer) {
            Flush(destRow);
        }
        //
        //  When not appending, increment the position of where the data
        //  is to be written.
        //
        if (destRow != -1) {
            destRow += this->bufferIndex / rowLength;
        }
    }
}

template<typename T>
void BufferedHDF2DArray<T>::Flush(DSLength destRow) {

    //
    // A default writeRow of -1 implies append
    //
    DSLength numDataRows;
    //
    // this->bufferIndex points after the end of the last data in the
    // buffer (full rows), so this->bufferIndex / rowLength is the
    // number of number of rows to create.
    //
    numDataRows = this->bufferIndex / rowLength;

    if (numDataRows > 0) {
        assert(fileDataSpaceInitialized);

        H5::DataSpace fileSpace;
        fileSpace = dataset.getSpace();

        //
        // Load the current size of the array on disk.
        //
        hsize_t fileArraySize[2], fileArrayMaxSize[2], blockStart[2];
        fileSpace.getSimpleExtentDims(fileArraySize, fileArrayMaxSize);

        // Save this for later to determine the offsets
        blockStart[0] = fileArraySize[0];
        blockStart[1] = fileArraySize[1];

        //
        // Calculate the number of rows to create.  This is dependent
        // on the current file size, the destination of where the data
        // will go, and how much to write.
        //

        if (destRow == -1) {
            fileArraySize[0] += numDataRows;
        }
        else {
            // If the data cannot fit in the current file size, extend
            // it,  otherwise, do not toch the file array size.
            if (destRow + numDataRows > fileArraySize[0]) {
                fileArraySize[0] = destRow + numDataRows;
            }
        }

        //
        // Make room in the file for the array.
        //
        dataset.extend(fileArraySize);

        H5::DataSpace extendedSpace = dataset.getSpace();
        //
        // Store the newly dimensioned dataspaces.
        //
        fileSpace.getSimpleExtentDims(fileArraySize, fileArrayMaxSize);			
        //
        // Configure the proper addressing to append to the array.
        //
        hsize_t dataSize[2];
        dataSize[0] = numDataRows;
        dataSize[1] = rowLength;
        hsize_t offset[2];
        //
        // Determine which row to write to.
        //
        if (destRow == -1) {
            offset[0] = blockStart[0];
        }
        else {
            offset[0] = destRow;
        }
        offset[1] = 0;
        extendedSpace.selectHyperslab(H5S_SELECT_SET, dataSize, offset);
        H5::DataSpace memorySpace(2, dataSize);

        //
        // Finally, write out the data.  
        // This uses a generic function which is specialized with
        // templates later on to t
        // memorySpace addresses the entire array in linear format
        // fileSpace addresses the last dataLength blocks of dataset.
        //
        TypedWriteRow(this->writeBuffer, memorySpace, extendedSpace);
        memorySpace.close();
        extendedSpace.close();
        fileSpace.close();
    }
    this->ResetWriteBuffer();
}

#endif // _BLASR_HDF_BUFFERED_HDF_2D_ARRAY_IMPL_HPP_

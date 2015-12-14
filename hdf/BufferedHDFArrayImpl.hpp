#ifndef _BLASR_HDF_BUFFERED_HDF_ARRAY_IMPL_HPP_
#define _BLASR_HDF_BUFFERED_HDF_ARRAY_IMPL_HPP_

#include <cstdlib>
#include <iostream>
#include <cstring>
#include "utils.hpp"
#include "BufferedHDFArray.hpp"

template<typename T>
BufferedHDFArray<T>::BufferedHDFArray(int pBufferSize) : HDFData() {
    nDims = 0;
    maxDims = 0;
    arrayLength = 0;
    dimSize = NULL;
    this->bufferIndex = 0;
    this->InitializeBuffer(pBufferSize);
}

template<typename T>
BufferedHDFArray<T>::BufferedHDFArray(H5::CommonFG* _container, 
    std::string _datasetName) : HDFData(_container, _datasetName) {
    // no-op
}

template<typename T>
BufferedHDFArray<T>::~BufferedHDFArray() {
    //
    // Clean up the write buffer.
    //
    if (dimSize != NULL) {
        delete[] dimSize;
        dimSize = NULL;
    }
    this->Free();
}

template<typename T>
void BufferedHDFArray<T>::SetBufferSize(int _bufferSize) {
    this->InitializeBuffer(_bufferSize);
}

template<typename T>
void BufferedHDFArray<T>::Write(const T *data, DSLength dataLength, bool append, 
    DSLength writePos) {

    // Fill the buffer with data. When there is overflow, write
    // that out to disk.
    //
    DSLength dataIndex = 0;
    int bufferCapacity;
    int bufferFillSize = 0;
    bool flushBuffer;
    while(dataIndex < dataLength) {
        bufferCapacity = this->bufferSize - this->bufferIndex;
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
            Flush(append, writePos);
        }
    }
}

template<typename T>
void BufferedHDFArray<T>::Flush(bool append, DSLength writePos) {
    //
    // Flush contents of current buffer to the file.
    //
    if (this->WriteBufferEmpty()) {
        // 
        // There is no data in the buffer, so nothing can be written.
        // HDF does not support empty arrays (as far as I can tell), so
        // don't even bother trying to create the dataspace.
        // 
        return;
    }

    // fetch the current size of the dataspace
    if (fileDataSpaceInitialized == false) {
        std::cout << "ERROR, trying to flush a dataset that has not been ";
        std::cout << "created or initialized" << std::endl;
        exit(1);
        fileDataSpaceInitialized = true;
    }

    H5::DataSpace fileSpace;
    fileSpace = dataset.getSpace();

    //
    // Load the current size of the array on disk.
    //
    hsize_t fileArraySize[1], blockStart;
    fileArraySize[0] = fileSpace.getSimpleExtentNpoints();
    if (append) {
        blockStart = fileSpace.getSimpleExtentNpoints();
        fileArraySize[0] += this->bufferIndex;
        //
        // Make room in the file for the array.
        //
        dataset.extend(fileArraySize);
    }
    else {
        blockStart = writePos;
        if (blockStart + this->bufferIndex > fileArraySize[0]) {
            fileArraySize[0] = blockStart + this->bufferIndex;
            dataset.extend(fileArraySize);
        }
    }

    H5::DataSpace extendedSpace = dataset.getSpace();
    //int extendedSize        = extendedSpace.getSimpleExtentNpoints(); // FIXME(yli): why isn't this used?
    //
    // Configure the proper addressing to append to the array.
    //
    hsize_t dataSize[1];
    hsize_t offset[1];
    dataSize[0] = this->bufferIndex;
    offset[0]   = blockStart;
    extendedSpace.selectHyperslab(H5S_SELECT_SET, dataSize, offset);
    H5::DataSpace memorySpace(1, dataSize);

    //
    // Finally, write out the data.  
    // This uses a generic function which is specialized with
    // templates later on to t
    // memorySpace addresses the entire array in linear format
    // fileSpace addresses the last dataLength blocks of dataset.
    //
    try {
        TypedWrite(this->writeBuffer, memorySpace, extendedSpace);
    }
    catch(H5::DataSetIException e) {
        std::cout <<"ERROR! Could not write HDF5 data." << std::endl;
        e.printError();
        exit(1);
    }
    memorySpace.close();
    extendedSpace.close();
    fileSpace.close();

    // Clear the buffer.
    this->ResetWriteBuffer();
}

template<typename T>
void BufferedHDFArray<T>::TypedWrite(const char **data, const H5::DataSpace &memorySpace, 
    const H5::DataSpace &extendedSpace) {

    H5::StrType varStrType(0,H5T_VARIABLE);
    dataset.write(data, varStrType, memorySpace, extendedSpace);
}

template<typename T>
void BufferedHDFArray<T>::TypedWrite(const T* data, const H5::DataSpace &memorySpace, 
    const H5::DataSpace &extendedSpace) {

    assert("Calling TypedWrite on an unsupported type" == 0);
}

template<typename T>
void BufferedHDFArray<T>::TypedCreate(H5::DataSpace &fileSpace, 
    H5::DSetCreatPropList &cparms) {

    std::cout << "DEFAULT typed create " << std::endl;
}


template<typename T>
void BufferedHDFArray<T>::Create(HDFGroup &parentGroup, std::string _datasetName) {
    return Create(&parentGroup.group, _datasetName);
}


template<typename T>
void BufferedHDFArray<T>::Create(H5::CommonFG* _container, std::string _datasetName) {
    //
    // Initialize where the dataset will go.
    container   = _container;
    datasetName = _datasetName;

    hsize_t dataSize[]    = {0};
    hsize_t maxDataSize[] = {H5S_UNLIMITED};
    H5::DataSpace fileSpace(1, dataSize, maxDataSize);
    H5::DSetCreatPropList cparms;

    /*
     * For some reason, chunking must be enabled when creating a dataset
     * that  has an unlimited dimension.  Of course, this is not
     * mentioned in the hdf5 c++ documentation, because that
     * docuemntation was written for people who enjoy learning how to
     * use an API by reading comments in source code.
     */
    hsize_t      chunk_dims[1] = { 16384 };
    cparms.setChunk( 1, chunk_dims );
    TypedCreate(fileSpace, cparms);

    //
    // Since TypedCreate created an assigned a dataset, this array is
    // now initialized.  Do the bookkeeping here.
    //

    isInitialized = true;
    fileDataSpaceInitialized = true;
    fileSpace.close();
}



/*
 * Initialize for reading.
 *
 * Open a dataset in an hdf file. Only call this on datasets that
 * exist, since this currently handles errors with opening datasets
 * by ungracefully exiting the program. 
 */

template<typename T>
int BufferedHDFArray<T>::InitializeForReading(HDFGroup &parentGroup, 
    const std::string datasetName) {

    return Initialize(parentGroup, datasetName, false);
}

template<typename T>
int BufferedHDFArray<T>::Initialize(HDFGroup &parentGroup, 
    const std::string &datasetName) {

    return Initialize(parentGroup, datasetName, true);
}

template<typename T>
int BufferedHDFArray<T>::Initialize(HDFGroup &parentGroup, 
    const std::string &datasetName, bool createIfMissing, DSLength newArrayLength) {
    //
    // For writing to this dataset, start at the first position in the
    // write buffer.
    //
    this->bufferIndex = 0;
    //
    // It's possible that the group may be asked to initialize this
    // dataset when the dataset does not exist.  Check that here.
    //
    bool parentHasObject = parentGroup.ContainsObject(datasetName);
    if ( parentHasObject and InitializeDataset(parentGroup, datasetName) == 0) {
        // 
        // The parent group already contains this dataset.  Try to
        // initialize this dataset and if it does not exist, flag fail.
        //
        return 0;
    }

    //
    // This is a hack to create in read/write mode.  If the parent
    // does not have the object, try and create it.  The problem with
    // trying to open a dataset in append mode is it will fail if the
    // dataset does not exist.
    //
    if (parentHasObject == false) {
        if (createIfMissing) {
            Create(parentGroup, datasetName);
        }
        else {
            //
            // Trying to open a dataset to read only, but it does not
            // exist.  Bail.
            //
            return 0;
        }
    }
    int ret = UpdateH5Dataspace();
    if (newArrayLength > 0) {
        ret *= Resize(newArrayLength);
    }
    return ret;
}

template<typename T>
int BufferedHDFArray<T>::UpdateH5Dataspace() {
    try {
        dataspace = dataset.getSpace();
    }
    catch(H5::DataSetIException &e) { 
        e.printError();
        return 0;
    }
    maxDims = MAX_DIMS;
    try {
        nDims = dataspace.getSimpleExtentNdims();
        /*
         * Prevent abuse of this class for multidimensional IO.
         */
        if (nDims != 1) {
            std::cout << "ERROR in HDF format: dataset: ";
            std::cout << datasetName << " should be 1-D, but it is not."; 
            std::cout << std::endl;
            exit(1);
        }

        /*
         * Load in the size of this dataset, and make a map to the whole thing.
         */
        if (dimSize != NULL) {
            delete [] dimSize;
            dimSize = NULL;
        }
        dimSize = ProtectedNew<hsize_t>(nDims);

        dataspace.getSimpleExtentDims(dimSize);
        arrayLength = dimSize[0];
        if (dimSize[0] == 0) {
            // DONT create a real dataspace if the size is 0
            // cout << "WARNING, trying to open a zero sized dataspace." << endl;
            dataspace.close();
            return 1;
        }
        fullSourceSpace = H5::DataSpace(1, dimSize);
        dataspace.close();
    }
    catch(H5::Exception& e) {
        e.printError();
        return 0;
    }
    return 1;
}

template<typename T>
int BufferedHDFArray<T>::Resize(const DSLength newArrayLength) {
    //
    // Resize this dataset. May or may not allocate space in file.
    // May or may not write fill value.
    //
    try{
        H5::DataSpace fileSpace;
        fileSpace = dataset.getSpace();

        hsize_t fileArraySize[1];
        fileArraySize[0] = newArrayLength;
        arrayLength = newArrayLength;
        dataset.extend(fileArraySize);
        fileSpace.close();
    } catch(H5::DataSetIException &e) { 
        e.printError();
        return 0;
    }
    return 1;
}


template<typename T>
void BufferedHDFArray<T>::Close() {
    if (dimSize != NULL) {
        delete[] dimSize;
        dimSize = NULL;
        HDFData::Close();
    }
}

template<typename T>
DSLength BufferedHDFArray<T>::size() {
    dataspace = dataset.getSpace();
    hsize_t dimSizeArray[1];
    dataspace.getSimpleExtentDims(dimSizeArray);
    dataspace.close();
    return dimSizeArray[0];
}

/*
 * Unspecialized form of read.
 * Read cannot be called on a type T* that does not have a
 * specialized template definition.  This is all determined at
 * compile time.  To ensure this, the following
 * default definition is provided that gives a nasty warning and
 * exits the code.
 */

template<typename T>
void BufferedHDFArray<T>::Read(DSLength start, DSLength end, T* dest) {
    (void)(start);
    (void)(end);
    (void)(dest);
    assert("ERROR, calling Read with an unsupported type. Use Read(start,end,datatype, dest) instead." == 0);
    exit(1); // this is in case the assert statement is removed.
}

/*
 * Read in type T from the opened dataset from the interval (start,
 * end].
 */

template<typename T>
void BufferedHDFArray<T>::ReadDataset(std::vector<T> &dest) {
    assert("ERROR, calling ReadDataset with an unsupported type.");
    exit(1); // this is in case the assert statement is removed.
}

template<typename T>
void BufferedHDFArray<T>::Read(DSLength start, DSLength end, H5::DataType typeID, T *dest) {
    if (end - start == 0) {
        return;
    }
    hsize_t memSpaceSize[] = {0};
    memSpaceSize[0] = end - start;
    hsize_t sourceSpaceOffset[] = {0};
    sourceSpaceOffset[0] = start;
    H5::DataSpace destSpace(1, memSpaceSize);
    fullSourceSpace.selectHyperslab(H5S_SELECT_SET, memSpaceSize, sourceSpaceOffset);
    dataset.read(dest, typeID, destSpace, fullSourceSpace);
    destSpace.close();
}

template<typename T>
void BufferedHDFArray<T>::ReadCharArray(DSLength start, DSLength end, std::string* dest) {
    hsize_t memSpaceSize[] = {0};
    memSpaceSize[0] = end - start;
    hsize_t sourceSpaceOffset[] = {0};
    sourceSpaceOffset[0] = start;
    H5::DataSpace destSpace(1, memSpaceSize);
    H5::StrType strType(0, H5T_VARIABLE);
    fullSourceSpace.selectHyperslab(H5S_SELECT_SET, memSpaceSize, sourceSpaceOffset);
    std::vector<char*> tmpStringArray;
    tmpStringArray.resize(end-start);
    dataset.read(&tmpStringArray[0], strType, destSpace, fullSourceSpace);
    for (size_t i = 0; i < tmpStringArray.size(); i++) {
        dest[i] = tmpStringArray[i];
    }
    destSpace.close();
}

#endif // _BLASR_HDF_BUFFERED_HDF_ARRAY_IMPL_HPP_

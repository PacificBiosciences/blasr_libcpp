#ifndef _BLASR_HDF_DATA_HPP_
#define _BLASR_HDF_DATA_HPP_

#include <string>
#include "H5Cpp.h"
#include "HDFConfig.hpp"
#include "HDFGroup.hpp"
#include "HDFAttributable.hpp"

class HDFData : public HDFAttributable {
public:
    H5::DataSet   dataset;
    H5::DataSpace dataspace;
    H5::DataSpace sourceSpace, destSpace;
    H5::DataSpace fullSourceSpace;
    bool      fileDataSpaceInitialized;
    H5::CommonFG  *container;
    std::string    datasetName;
    bool      isInitialized;

    H5::H5Location* GetObject(); 

    HDFData(H5::CommonFG* _container, const std::string & _datasetName); 

    HDFData(); 

    bool IsInitialized() const; 

    //
    // Allow derived classes to be initialized generically.
    //
    // FIXME(yli): This method and the Initialize methods in
    // subclasses cause compilation warnings under clang.  For
    // example, see the discussion here:
    //
    // http://stackoverflow.com/questions/18515183/c-overloaded-virtual-function-warning-by-clang
    //
    // To fix the problem here may require rethinking the initialize
    // method more carefully.
    //
    virtual int Initialize(HDFGroup &parentGroup, const std::string &datasetName);

    int BaseInitializeDataset(H5::CommonFG &hdfFile, const std::string & _datasetName); 

    int InitializeDataset(HDFGroup &group, const std::string & _datasetName); 

    int InitializeDataset(H5::CommonFG &hdfFile, const std::string & _datasetName); 

    void Close(); 
};

#endif

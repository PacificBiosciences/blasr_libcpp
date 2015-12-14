#include "HDFData.hpp"

using namespace std;
using namespace H5;

H5Location* HDFData::GetObject() {
    return &dataset;
}

HDFData::HDFData(CommonFG* _container, const string & _datasetName) {
    container   = _container;
    datasetName = _datasetName;
    fileDataSpaceInitialized = false;
    isInitialized = false;
}

HDFData::HDFData() {
    container = NULL;
    fileDataSpaceInitialized = false;
    isInitialized = false;
}

bool HDFData::IsInitialized() const {
    return isInitialized;
}

//
// Allow derived classes to be initialized generically.
//
int HDFData::Initialize(HDFGroup &parentGroup, const string &datasetName) { 
    (void)(parentGroup);
    (void)(datasetName);
    cout << "ERROR! Only a subclass should call this." << endl;
    exit(1);
}

int HDFData::BaseInitializeDataset(CommonFG &hdfFile, const string & _datasetName) {
    dataset   = hdfFile.openDataSet(_datasetName.c_str());
    isInitialized = true;
    fileDataSpaceInitialized = true;
    return 1;
}

int HDFData::InitializeDataset(HDFGroup &group, const string & _datasetName) {
    return InitializeDataset(group.group, _datasetName);
}

int HDFData::InitializeDataset(CommonFG &hdfFile, const string & _datasetName) {
    try {
        datasetName = _datasetName;
        dataset   = hdfFile.openDataSet(_datasetName.c_str());
        isInitialized = true;
        fileDataSpaceInitialized = true;
    }
    catch(FileIException &e) {
        cerr << e.getDetailMsg() <<endl;
        return 0;
    }
    catch(GroupIException &e) {
        cerr << e.getDetailMsg() << endl;
        return 0;
    }
    catch(H5::Exception &e) {
        cerr << e.getDetailMsg() << endl;
        return 0;
    }
    return 1;
}

void HDFData::Close() {
    if (isInitialized) {
        dataspace.close();
        dataset.close();
        isInitialized = false;
    }
}

#include "utils/FileOfFileNames.hpp"
#include "HDFNewBasReader.hpp"
#include <cstdlib>

void
FileOfFileNames::StoreFileOrFileList(std::string fileName, 
    std::vector<std::string> &fofnList) {

    std::vector<std::string> tmpList;
    if (IsFOFN(fileName)) {
        FOFNToList(fileName, tmpList);
    }
    else {
        tmpList.push_back(fileName);
    }
    for (int i = 0; i < int(tmpList.size()); i++) {
        if (FileOfFileNames::IsFOFN(tmpList[i])) {
            std::cout << "ERROR. Nested File of File Names are not allowed. "
                      << std::endl;
            exit(1);
        } else if (FileOfFileNames::IsBasH5(tmpList[i])) {
            std::vector<std::string> baxFNs = FileOfFileNames::Bas2Bax(tmpList[i]);
            fofnList.insert(fofnList.end(), baxFNs.begin(), baxFNs.end());
        } else {
            fofnList.push_back(tmpList[i]);
        }
    }
}

void 
FileOfFileNames::FOFNToList(std::string &fofnFileName, 
    std::vector<std::string> &fofnList) {
    std::ifstream fofnIn;
    CrucialOpen(fofnFileName, fofnIn);
    while(fofnIn) {
        std::string name;
        getline(fofnIn, name);
        if (name.size() > 0) {
            fofnList.push_back(name);
        }
    }
}

bool 
FileOfFileNames::IsFOFN(std::string &fileName) {
    std::string::size_type dotPos = fileName.rfind(".");
    if (dotPos != std::string::npos) {
        std::string extension;
        extension.assign(fileName, dotPos+1, fileName.size() - (dotPos+1));
        if (extension == "fofn") {
            return true;
        }
    }
    return false;
}

bool 
FileOfFileNames::IsBasH5(std::string & fileName) {
    // Return true if file ends with bas.h5
    if (fileName.size() > 6) {
        if (fileName.rfind("bas.h5") == fileName.size() - 6) {
            return true;
        }
    }
    return false;
}

std::vector<std::string> 
FileOfFileNames::Bas2Bax(std::string & basFN) {
    // There are two types of bas.h5 files.
    // Before SMRT 2.0, bas.h5 files contain all the /PulseData data,
    // in this case, return the bas.h5.
    // After SMRT 2.0, bas.h5 files have been changed to only contain
    // paths to bax.h5 files (in the /MultiPart/Parts group), while
    // all base calls and QVs are in bax.h5 files. In this case,
    // return path to the bax.h5 files. Assumption is that bax.h5
    // files are in the same directory as bas.h5 file.
    vector<string> baxFNs;
    HDFNewBasReader reader;
    if (reader.Initialize(basFN) != 0) {
        baxFNs = reader.GetBaxFileNames();
    } else {
        baxFNs.push_back(basFN);
    }
    reader.Close();
    return baxFNs;
}


int
FileOfFileNames::ExpandFileNameList(std::vector<std::string> &fileNames) {
    int rfn;
    std::vector<std::string> expandedFileNames;
    for (rfn = 0; rfn < static_cast<int>(fileNames.size()); rfn++) {
        std::vector<std::string> tmpList;
        FileOfFileNames::StoreFileOrFileList(fileNames[rfn], tmpList);
        expandedFileNames.insert(expandedFileNames.end(),
                                 tmpList.begin(), tmpList.end());
    }
    fileNames = expandedFileNames;
    return fileNames.size();
}

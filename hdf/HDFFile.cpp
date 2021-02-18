#include <hdf/HDFFile.hpp>

using namespace H5;

HDFFile::HDFFile() {}

//
//  Open a file.  By default, if the file already exists, open it in
//  read/write mode.  The only other flag that is allowed is
//  H5F_ACC_TRUNC, which will truncate the file to zero size.
//
void HDFFile::Open(std::string fileName, unsigned int flags, const FileAccPropList& fileAccPropList)
{

    assert(flags == H5F_ACC_RDWR || flags == H5F_ACC_TRUNC || flags == H5F_ACC_RDONLY);
    std::ifstream testIn(fileName.c_str());
    bool fileExists = static_cast<bool>(testIn);
    bool flagsIsNotTrunc = flags != H5F_ACC_TRUNC;

    if (fileExists and H5File::isHdf5(fileName.c_str()) and flagsIsNotTrunc) {
        try {
            hdfFile.openFile(fileName.c_str(), flags, fileAccPropList);
        } catch (const FileIException& e) {
            std::cout << "Error opening file " << fileName << std::endl;
            std::exit(EXIT_FAILURE);
        }
    } else {
        try {
            //
            // Open a new file with TRUNC permissions, always read/write.
            //
            FileCreatPropList filePropList;
            filePropList.setUserblock(512);
            hdfFile = H5File(fileName.c_str(), H5F_ACC_TRUNC, filePropList);
        } catch (const FileIException& fileException) {
            std::cout << "Error creating file " << fileName << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    if (rootGroup.Initialize(hdfFile, "/") != 1) {
        std::cout << "Error initializing the root group for file " << fileName << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

void HDFFile::Close() { hdfFile.close(); }

#include <cstdlib>
#include "BaseSequenceIO.hpp"

using namespace std;

void BaseSequenceIO::SetFiles(FileType &pFileType, std::string &pFileName) {
    fileType = pFileType;
    fileName = pFileName;
}

FileType BaseSequenceIO::GetFileType() {
    return fileType;
}

int BaseSequenceIO::DetermineFileTypeByExtension(string &fileName, 
    FileType &type, bool exitOnFailure) {

    string::size_type dotPos = fileName.rfind(".");
    if (dotPos != string::npos) {
        string extension;
        extension.assign(fileName, dotPos+1, fileName.size() - (dotPos+1));
        if (extension == "fasta" or
                extension == "fa" or
                extension == "fas" or
                extension == "fsta" or
                extension == "screen"
           ) {
            type = FileType::Fasta;
            return 1;
        }
        else if (extension == "h5") {
            dotPos = fileName.rfind(".", dotPos-1); 
            extension.assign(fileName, dotPos+1, fileName.size() - (dotPos + 1));
            if (extension == "pls.h5" or extension == "plx.h5" ) {
                type = FileType::HDFPulse;
                return 1;
            }
            else if (extension == "bas.h5" or extension == "bax.h5") {
                type = FileType::HDFBase;
                return 1;
            }
            else if (extension == "ccs.h5") {
                type = FileType::HDFCCSONLY;
                return 1;
            }
            else {
                type = FileType::None;
                return 0;
            }
        }
        else if (extension == "fastq" or
                extension == "fq") {
            type = FileType::Fastq;
            return 1;
        }
        else if (extension == "4bit" or
                extension == "fourbit") {
            type = FileType::Fourbit;
            assert("Four bit reading is not yet implemented for the reader agglomerate!" == 0);
            return 1;
        }
        else if (extension == "bam") {
            type = FileType::PBBAM;
            return 1;
        }
        else if (extension == "xml") {
            type = FileType::PBDATASET;
            return 1;
        }
        else {
            type = FileType::None;
            if (exitOnFailure) {
                cout << "ERROR, file type '." << extension << "' is not understood to be one of pls.h5, fasta, fastq, nor bam. " << endl;
                exit(1);
            }
            return 0;
        }
        return 0;
    }
    return 0;
}

// Author: Yuan Li

#ifndef _BLASR_HDFWRITERBASE_HPP_
#define _BLASR_HDFWRITERBASE_HPP_
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include "pbdata/SMRTSequence.hpp"
#include "pbdata/reads/RegionAnnotation.hpp"
#include "HDFFile.hpp"
#include "HDFGroup.hpp"
#include "HDFAtom.hpp"
#include "BufferedHDFArray.hpp"
#include "BufferedHDF2DArray.hpp"
#include "HDFScanDataWriter.hpp"

class HDFWriterBase {
public:
    HDFWriterBase(const std::string & filename)
    : filename_(filename)
    {}

    virtual ~HDFWriterBase(void) = 0;

public:
    /// \returns Target H5 filename.
    std::string Filename(void) {return filename_;}

    std::vector<std::string> Errors(void) const;

    /// Copy the given group and write to the output
    /// \returns Object that was copied
    void CopyObject(HDFFile& src, const char* path); 

    void WriteScanData(const ScanData& scanData);

    virtual bool WriteFakeDataSets();

    virtual void Flush();

    virtual bool WriteOneZmw(const SMRTSequence& seq,
                             const std::vector<RegionAnnotation>& regions);
   
protected:
    std::string filename_;
    std::vector<std::string> errors_; 

    bool AddChildGroup(HDFGroup & parentGroup, 
                       HDFGroup & childGroup,
                       const std::string & childGroupName);

    bool AddAttribute(HDFData & group, 
                      const std::string & attributeName, 
                      const std::vector<std::string> & attributeValues);

    bool AddAttribute(HDFGroup & group, 
                      const std::string & attributeName, 
                      const std::vector<std::string> & attributeValues);

    template<typename T>
    bool AddAttribute(HDFData & group, 
                      const std::string & attributeName, 
                      const T & attributeValue);

    template<typename T>
    bool AddAttribute(HDFGroup & group, 
                      const std::string & attributeName, 
                      const T & attributeValue);

    /// \brief Checks whether chemistry triple, including
    ///        binding kit, sequencing kit and base caller version
    ///        are set. 
    ///        If not, add error messages.
    void SanityCheckChemistry(const std::string & bindingKit,
                              const std::string & sequencingKit);

    void AddErrorMessage(const std::string & errmsg);

    void FAILED_TO_CREATE_GROUP_ERROR(const std::string & groupName);
   
    void FAILED_TO_CREATE_ATTRIBUTE_ERROR(const std::string & attributeName);

    void PARENT_GROUP_NOT_INITIALIZED_ERROR(const std::string & groupName);

    virtual void Close(void) = 0;

	HDFFile outfile_;  ///< HDFFile file handler

};


template<typename T>
bool HDFWriterBase::AddAttribute(HDFData & group, 
                                 const std::string & attributeName, 
                                 const T & attributeValue) {
    try {
        HDFAtom<T> attributeAtom;
        attributeAtom.Create(group.dataset, std::string(attributeName));
        attributeAtom.Write(attributeValue);
        attributeAtom.Close();
    }
    catch (H5::Exception &e) {
        this->FAILED_TO_CREATE_ATTRIBUTE_ERROR(attributeName);
        return false;
    }
    return true;
}

template<typename T>
bool HDFWriterBase::AddAttribute(HDFGroup & group, 
                    const std::string & attributeName, 
                    const T & attributeValue) {
    try {
        HDFAtom<T> attributeAtom;
        attributeAtom.Create(group.group, std::string(attributeName));
        attributeAtom.Write(attributeValue);
        attributeAtom.Close();
    }
    catch (H5::Exception &e) {
        FAILED_TO_CREATE_ATTRIBUTE_ERROR(attributeName);
        return false;
    }
    return true;
}

/// \brief Write a dataset of name 'dsName' under h5 group 'dsGroup' with
///        length 'dsLength'. Just fill this dataset with buffer repeatedly.
template<typename T>
bool __WriteFakeDataSet(HDFGroup & dsGroup, const std::string & dsName,
                        const uint32_t dsLength, std::vector<T> & buffer) {
    BufferedHDFArray<T> dsArray_;
    if (dsArray_.Initialize(dsGroup, dsName) == 0) return false;
    uint32_t totalLength = 0; 
    while(totalLength < dsLength) {
        uint32_t thisLength = buffer.size();
        if (totalLength + thisLength <= dsLength) {
            totalLength = totalLength + thisLength;
        } else {
            thisLength = dsLength - totalLength;
            totalLength = dsLength;
        }
        dsArray_.Write(&buffer[0], thisLength);
        dsArray_.Flush();
    }
    dsArray_.Close();
    return true;
}

/// \brief Write a 2D dataset of name 'dSName' under h5 group 'dsGroup' with
///        row number 'rowNum', column number 'colNum'. Just fill this dataset
///        with value 'fillData'
template<typename T>
bool __WriteFake2DDataSet(HDFGroup & dsGroup, const std::string & dsName,
                           const uint32_t rowNum, const size_t colNum, const T & fillData) {
    BufferedHDF2DArray<T> dsArray;
    if (dsArray.Initialize(dsGroup, dsName, colNum) == 0) return false;
    T * data = new T[colNum];
    for(size_t i = 0; i < colNum; i++) {data[i] = fillData; }
    for(uint32_t i = 0; i< rowNum; i++) {
        dsArray.WriteRow(data, colNum);
    }
    dsArray.Close();
    delete []data;
    return true;
}

#endif

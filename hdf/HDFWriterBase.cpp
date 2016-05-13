// Author: Yuan Li

#include "HDFWriterBase.hpp"


HDFWriterBase::~HDFWriterBase(void) 
{ }

std::vector<std::string> HDFWriterBase::Errors(void) const {
    return errors_;
}

void HDFWriterBase::CopyObject(HDFFile& src, const char* path) {
    H5Ocopy(src.hdfFile.getId(), path, outfile_.hdfFile.getId(), path, H5P_DEFAULT, H5P_DEFAULT);
}

void HDFWriterBase::WriteScanData(const ScanData& scanData) {
    // sanity check chemistry meta data. 
    SanityCheckChemistry(scanData.BindingKit(), scanData.SequencingKit());
    HDFScanDataWriter writer(outfile_.rootGroup);
    writer.Write(scanData);
}

bool HDFWriterBase::WriteFakeDataSets() { return true; }

void HDFWriterBase::Flush() { }

bool HDFWriterBase::WriteOneZmw(const SMRTSequence & seq,
                               const std::vector<RegionAnnotation> & regions) {
    // XXX: unused variables
    return true;
}

bool HDFWriterBase::AddChildGroup(HDFGroup & parentGroup, 
                                  HDFGroup & childGroup,
                                  const std::string & childGroupName) {
    parentGroup.AddGroup(childGroupName);
    if (childGroup.Initialize(parentGroup, childGroupName) == 0) {
        FAILED_TO_CREATE_GROUP_ERROR(childGroupName);
        return false;
    }
    return true;
}

template<>
bool HDFWriterBase::AddAttribute(HDFData & group, 
                                 const std::string & attributeName, 
                                 const std::vector<std::string> & attributeValues)
{
    try {
        HDFAtom<std::vector<std::string> > attributeAtom;
        attributeAtom.Create(group.dataset, std::string(attributeName), attributeValues);
        attributeAtom.Close();
    }
    catch (H5::Exception &e) {
        FAILED_TO_CREATE_ATTRIBUTE_ERROR(attributeName);
        return false;
    }
    return true;
}

template<>
bool HDFWriterBase::AddAttribute(HDFGroup & group, 
                                 const std::string & attributeName, 
                                 const std::vector<std::string> & attributeValues)
{
    try {
        HDFAtom<std::vector<std::string> > attributeAtom;
        attributeAtom.Create(group.group, std::string(attributeName), attributeValues);
        attributeAtom.Close();
    }
    catch (H5::Exception &e) {
        FAILED_TO_CREATE_ATTRIBUTE_ERROR(attributeName);
        return false;
    }
    return true;
}

void HDFWriterBase::SanityCheckChemistry(const std::string & bindingKit,
                                         const std::string & sequencingKit)
{
    if (bindingKit.empty()) {
        AddErrorMessage("Binding kit must be specified.");
    }
    if (sequencingKit.empty()) {
        AddErrorMessage("Sequencing kit must be specified.");
    }
}


void HDFWriterBase::AddErrorMessage(const std::string & errmsg) {
    errors_.push_back(errmsg);
}

void HDFWriterBase::FAILED_TO_CREATE_GROUP_ERROR(const std::string & groupName) {
    std::stringstream ss;
    ss << "Failed to create group " << groupName << " in " << filename_;
    AddErrorMessage(ss.str());
}

void HDFWriterBase::FAILED_TO_CREATE_ATTRIBUTE_ERROR(const std::string & attributeName) {
    std::stringstream ss;
    ss << "Failed to create attribute " << attributeName << " in " << filename_;
    AddErrorMessage(ss.str());
}

void HDFWriterBase::PARENT_GROUP_NOT_INITIALIZED_ERROR(const std::string & groupName) {
    std::stringstream ss;
    ss << "Parent hdf group of " << groupName << " in file " << filename_
       << " is not initialized.";
    AddErrorMessage(ss.str());
}

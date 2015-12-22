#ifndef HDF_FILE_LOG_GROUP_H_
#define HDF_FILE_LOG_GROUP_H_
#include "HDFGroup.hpp"
#include "HDFArray.hpp"

class HDFFileLogGroup {
 public:
  HDFGroup group;
  HDFStringArray commandLineArray;
  HDFStringArray versionArray;
  HDFStringArray timestampArray;
  HDFArray<unsigned int> idArray;
  HDFStringArray logArray;
  HDFStringArray programArray;

  int Initialize(HDFGroup &parentGroup) {
    if (group.Initialize(parentGroup.group, "FileLog") == 0) { return 0; }
    int ret = 1;
    ret *= commandLineArray.Initialize(group, "CommandLine");
    ret *= versionArray.Initialize(group, "Version");
    ret *= timestampArray.Initialize(group, "Timestamp");
    ret *= idArray.Initialize(group, "ID");
    ret *= logArray.Initialize(group, "Log");
    ret *= programArray.Initialize(group, "Program");
    return ret;
  }

  
  void AddEntry(string command, string log, string program, string timestamp, string version) {
    commandLineArray.Write(&command, 1);
    versionArray.Write(&version, 1);
    timestampArray.Write(&timestamp, 1);
    programArray.Write(&program, 1);
    logArray.Write(&log, 1);
    
    unsigned int id = idArray.size();
    id = id + 1;
    idArray.Write(&id, 1);
    
  }

  bool Create(HDFGroup &parent) {
    parent.AddGroup("FileLog");
    if (group.Initialize(parent.group, "FileLog") == 0) { return false; }
    commandLineArray.Create(group, "CommandLine");
    versionArray.Create(group, "Version");
    timestampArray.Create(group, "Timestamp");
    programArray.Create(group, "Program");
    logArray.Create(group, "Log");
    idArray.Create(group, "ID");
    return true;
  }    

};

#endif

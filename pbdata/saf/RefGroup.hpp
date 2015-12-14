#ifndef _BLASR_REF_GROUP_HPP_
#define _BLASR_REF_GROUP_HPP_

#include <string>
#include <vector>

#include <stdint.h>

using namespace std;
class RefGroup {
 public:
	vector<uint32_t> id;
	vector<string> path;
	vector<string> refGroupName;
	vector<uint32_t> refInfoId;

    bool IdToIndex(uint32_t idKey, uint32_t &idIndex);
	int FindPath(uint32_t idKey, string &pathVal, string &groupNameVal);
};

#endif

#ifndef _BLASR_ZMW_GROUP_ENTRY_HPP_
#define _BLASR_ZMW_GROUP_ENTRY_HPP_

#include "Types.h"

class ZMWGroupEntry {
 public:
	UInt holeNumber;
	UInt x;
	UInt y;
	DNALength numEvents;
    unsigned char holeStatus;
    ZMWGroupEntry(); 
};


#endif // _BLASR_ZMW_GROUP_ENTRY_HPP_

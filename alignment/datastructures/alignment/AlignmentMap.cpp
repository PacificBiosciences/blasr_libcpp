#include <vector>
#include <string> 
#include "AlignmentMap.hpp"

using namespace std;

void 
CreateSequenceToAlignmentMap(const string & alignedSequence,
        vector<int> & baseToAlignmentMap) {
    baseToAlignmentMap.resize(alignedSequence.size());
    size_t alignedPos, unalignedPos;
    for (alignedPos=0, unalignedPos=0; 
         alignedPos < alignedSequence.size();
         alignedPos++) {
        if (not (alignedSequence[alignedPos] == ' ' or 
            alignedSequence[alignedPos] == '-')) {
            baseToAlignmentMap[unalignedPos] = alignedPos;
            unalignedPos++;
        }
    }
    baseToAlignmentMap.resize(unalignedPos);
}

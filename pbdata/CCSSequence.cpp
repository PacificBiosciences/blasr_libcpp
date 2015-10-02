#include "CCSSequence.hpp"

void CCSSequence::Free() {
    numPasses = 0;
    numConsensusBases = 0;
    SMRTSequence::Free();
    unrolledRead.Free();
}

int CCSSequence::GetStorageSize() {
    return SMRTSequence::GetStorageSize() + unrolledRead.GetStorageSize();
}

UInt CCSSequence::HoleNumber(void) const {
    return SMRTSequence::HoleNumber();
}

CCSSequence & CCSSequence::HoleNumber(const UInt holeNumber) {
    SMRTSequence::HoleNumber(holeNumber);
    unrolledRead.HoleNumber(holeNumber);
    return *this;
}

//
// In the first iteration, Explode simply pulls the subreads out
// that are used in the ccs.   Eventually, it will pull out all
// high-quality subreads.
// 
void CCSSequence::Explode(std::vector<SMRTSequence> &subreads) {
    subreads.resize(numPasses);
    int subreadIndex;
    for (subreadIndex = 0; subreadIndex < numPasses; subreadIndex++) {
        subreads[subreadIndex].ReferenceSubstring(this->unrolledRead, passStartBase[subreadIndex], passNumBases[subreadIndex]);
        subreads[subreadIndex].zmwData = unrolledRead.zmwData;
    }
}

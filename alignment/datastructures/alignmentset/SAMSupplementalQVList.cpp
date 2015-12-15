#include "SAMSupplementalQVList.hpp"

const char* SupplementalQVList::qvNames[] = {"InsertionQV", "DeletionQV", "SubstitutionQV", "MergeQV", "SubstitutionTag", "DeletionTag"};
const char* SupplementalQVList::qvTags[] = {"iq", "dq", "sq", "mq", "st", "dt"};

// Only the first 4 tags are quality values.
int SupplementalQVList::nqvTags = 4;
int SupplementalQVList::nTags = 6;

int SupplementalQVList::UseQV(std::vector<std::string> &qvList) {
    useqv = 0;
    for (size_t i = 0; i < qvList.size(); i++) {
        int j;
        for (j = 0; j < nTags; j++) {
            if (qvList[i] == qvNames[j]) {
                useqv |= 1 << j;
                break;
            }
        }
        if (j == nTags) {
            return 1;
        }
    }
    return 0;
}

void SupplementalQVList::FormatQVOptionalFields(SMRTSequence &alignedSubsequence) {
    int i;
    for (i = 0; i < nqvTags; i++) {
        if (alignedSubsequence.GetQVPointerByIndex(i+1)->data == NULL) {
            // mask off this quality value since it does not exist
            useqv = useqv & ~(1 << i);
        }
    }
    for (i = 0; i < nqvTags; i++) {
        if (useqv & (1 << i)) {
            QualityVectorToPrintable(alignedSubsequence.GetQVPointerByIndex(i+1)->data, alignedSubsequence.length);
        }
    }
}
	
void SupplementalQVList::PrintQVOptionalFields(SMRTSequence &alignedSubsequence, std::ostream &out) {
    int i = 0;
    for (i = 0; i < nqvTags; i++) {
        if (alignedSubsequence.GetQVPointerByIndex(i+1)->data == NULL) {
            // mask off this quality value since it does not exist
            useqv = useqv & ~(1 << i);
        }
    }
    for (i = 0; i < nTags; i++) {
        if (alignedSubsequence.GetQVPointerByIndex(i+1) != NULL and (useqv & (1 << i)) ) {
            out << "\t" << qvTags[i] << ":Z:";
            alignedSubsequence.PrintAsciiRichQuality(out, i + 1, 0);
        }
    }
    if (alignedSubsequence.substitutionTag != NULL and (useqv & SubstitutionTag)) {
        out << "\t" << qvTags[I_SubstitutionTag-1] << ":Z:";
        alignedSubsequence.PrintAsciiRichQuality(out, I_SubstitutionTag, 0);
    }
    if (alignedSubsequence.deletionTag != NULL and (useqv & DeletionTag)) {
        out << "\t" << qvTags[I_DeletionTag-1] << ":Z:";
        alignedSubsequence.PrintAsciiRichQuality(out, I_DeletionTag, 0);
    }
}

void SupplementalQVList::clear() {
    for (int j = 0; j < nTags; j++) {
        useqv &= 0 << j;
    }
}

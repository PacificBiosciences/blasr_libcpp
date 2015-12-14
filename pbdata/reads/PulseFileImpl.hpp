#ifndef _BLASR_PULSE_FILE_IMPL_HPP_
#define _BLASR_PULSE_FILE_IMPL_HPP_

#include "utils.hpp"

template<typename T_FieldType>
void PulseFile::StoreField(std::vector<T_FieldType> &source, int *basToPlsIndex, T_FieldType *dest, const DSLength destLength) {
    for (DSLength i = 0 ; i < destLength; i++) {
        dest[i] = source[basToPlsIndex[i]];
    }
}

template <typename T>
bool PulseFile::Realloc(T *&ptr, const DNALength length) {
    if (ptr != NULL) {
        delete[] ptr;
    }
    ptr = ProtectedNew<T>(length);
    return ptr != NULL;
}

#endif

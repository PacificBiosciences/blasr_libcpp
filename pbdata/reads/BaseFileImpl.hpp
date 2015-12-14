#ifndef _BLASR_BASE_FILE_IMPL_HPP_
#define _BLASR_BASE_FILE_IMPL_HPP_

#include "BaseFile.hpp"

template<typename T>
void BaseFile::CopyArray(std::vector<T> &fullArray, DSLength pos, DNALength length, T*dest) {
    memcpy(dest, &fullArray[pos], sizeof(T) * length);
}

#endif

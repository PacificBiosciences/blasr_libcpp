#ifndef _BLASR_UTIL_IMPL_HPP_
#define _BLASR_UTIL_IMPL_HPP_
#include <stdlib.h>
#include <cstdlib>   // abort()
#include <new>       // bad_alloc
#include <iostream>  // cout/cerr


template<typename t_file>
void CrucialOpen(std::string &fileName, t_file &file, std::ios_base::openmode mode) {
	if (mode==0)
		file.open(fileName.c_str());
	else
		file.open(fileName.c_str(), mode);

	if (!file.good()) {
		std::cout << "Could not open " << fileName << std::endl;
		exit(1);
	}
}
template<typename T_Int>
T_Int CeilOfFraction(T_Int num, T_Int denom) {
	return num / denom + ((num % denom) && 1);
}

template<typename T>
inline T* ProtectedNew(uint64_t size) {
    T * ptr = nullptr;
    try {
        ptr = new T[size];
    } catch (std::bad_alloc & ba) {
        std::cout << "ERROR, allocating " << size * sizeof(T) << " bytes."
                  << ba.what() << std::endl;
        abort();
    }
    return ptr;
}

template<typename T>
inline T* ProtectedNew(void) {
    T * ptr = nullptr;
    try {
       ptr = new T;
    } catch (std::bad_alloc & ba) {
        std::cout << "ERROR, allocating " << sizeof(T) << " bytes."
                  << ba.what() << std::endl;
        abort();
    }
    return ptr;
}

#endif

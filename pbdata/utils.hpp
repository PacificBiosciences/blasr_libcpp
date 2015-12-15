#ifndef _BLASR_UTILS_HPP_
#define _BLASR_UTILS_HPP_
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>

template<typename t_file>
void CrucialOpen(std::string &fileName, t_file &file, 
    std::ios_base::openmode mode=(std::ios_base::openmode)0 ); 

template<typename T_Int>
T_Int CeilOfFraction(T_Int num, T_Int denom);

template<typename T>
inline T* ProtectedNew(uint64_t size); 

template<typename T>
inline T* ProtectedNew(void); 

#include "utilsImpl.hpp"

#endif // _BLASR_UTILS_HPP_

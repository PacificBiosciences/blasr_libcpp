#ifndef _BLASR_DISTANCE_MATRIX_SCORE_FUNCTION_IMPL_HPP_
#define _BLASR_DISTANCE_MATRIX_SCORE_FUNCTION_IMPL_HPP_
#include <string>
#include <iostream>
#include <stdint.h>
#include <string>
#include <cstring>
#include <ostream>
// pbdata
#include "../../../pbdata/Types.h"
#include "../../../pbdata/NucConversion.hpp"
#include "../../../pbdata/Enumerations.h"
#include "../../../pbdata/DNASequence.hpp"
#include "../../../pbdata/FASTASequence.hpp"
#include "../../../pbdata/FASTQSequence.hpp"

#include "DistanceMatrixScoreFunction.hpp"

template<typename T_RefSequence, typename T_QuerySequence>
void DistanceMatrixScoreFunction<T_RefSequence,T_QuerySequence>::InitializeScoreMatrix(int scoreMatrixP[5][5]) {
    int i, j;
    for (i = 0; i < 5; i++ ){ 
        for (j = 0; j < 5; j++ ){
            scoreMatrix[i][j] = scoreMatrixP[i][j];
        }
    }
}

template<typename T_RefSequence, typename T_QuerySequence>
DistanceMatrixScoreFunction<T_RefSequence,T_QuerySequence>::DistanceMatrixScoreFunction() : BaseScoreFunction() {
}

template<typename T_RefSequence, typename T_QuerySequence>
DistanceMatrixScoreFunction<T_RefSequence,T_QuerySequence>::DistanceMatrixScoreFunction(int scoreMatrixP[5][5], int insertionP, int deletionP) : BaseScoreFunction() {
    InitializeScoreMatrix(scoreMatrixP);
    ins = insertionP;
    del = deletionP;
}

template<typename T_RefSequence, typename T_QuerySequence>
int DistanceMatrixScoreFunction<T_RefSequence,T_QuerySequence>::Deletion(
    T_RefSequence &seq, DNALength pos, T_QuerySequence &querySeq, 
    DNALength queryPos) {
    (void)(seq); (void)(pos); (void)(querySeq); (void)(queryPos);
    return del;
}

template<typename T_RefSequence, typename T_QuerySequence>
int DistanceMatrixScoreFunction<T_RefSequence,T_QuerySequence>::Insertion(
    T_RefSequence &seq, DNALength pos, T_QuerySequence &querySeq, 
    DNALength queryPos) {
    (void)(seq); (void)(pos); (void)(querySeq); (void)(queryPos);
    return ins;
}

template<typename T_RefSequence, typename T_QuerySequence>
int DistanceMatrixScoreFunction<T_RefSequence,T_QuerySequence>::Deletion(
    T_RefSequence &seq, DNALength pos) {
    (void)(seq); (void)(pos);
    return del;
}


template<typename T_RefSequence, typename T_QuerySequence>
int DistanceMatrixScoreFunction<T_RefSequence,T_QuerySequence>::Match(
    T_RefSequence &ref, DNALength refPos, T_QuerySequence &query, 
    DNALength queryPos) {
    return scoreMatrix[ThreeBit[ref[refPos]]][ThreeBit[query[queryPos]]];		
}

//
// Define the score function on dereferenced pointers for speed.
//
template<typename T_RefSequence, typename T_QuerySequence>
int DistanceMatrixScoreFunction<T_RefSequence,T_QuerySequence>::Match(
    Nucleotide ref, Nucleotide query) {
    return scoreMatrix[ThreeBit[ref]][ThreeBit[query]];
}

template<typename T_RefSequence, typename T_QuerySequence>
int DistanceMatrixScoreFunction<T_RefSequence,T_QuerySequence>::Insertion(
    T_QuerySequence &seq, DNALength pos) {
    (void)(seq); (void)(pos);
    return ins;
}

template<typename T_RefSequence, typename T_QuerySequence>
float 
DistanceMatrixScoreFunction<T_RefSequence,T_QuerySequence>
::NormalizedMatch(T_RefSequence &ref, DNALength refPos, 
    T_QuerySequence &query, DNALength queryPos) {
    (void)(ref); (void)(refPos); (void)(query); (void)(queryPos);
    return 0;
}

template<typename T_RefSequence, typename T_QuerySequence>
float DistanceMatrixScoreFunction<T_RefSequence,T_QuerySequence>
::NormalizedInsertion(T_RefSequence &ref, DNALength refPos, 
    T_QuerySequence &query, DNALength queryPos) {
    (void)(ref); (void)(refPos); (void)(query); (void)(queryPos);
    return 0;
}

template<typename T_RefSequence, typename T_QuerySequence>
float DistanceMatrixScoreFunction<T_RefSequence,T_QuerySequence>::NormalizedDeletion(T_RefSequence &ref, DNALength refPos, T_QuerySequence &query, DNALength queryPos) {
    (void)(ref); (void)(refPos); (void)(query); (void)(queryPos);
    return 0;
}


#endif // _BLASR_DISTANCE_MATRIX_SCORE_FUNCTION_IMPL_HPP_

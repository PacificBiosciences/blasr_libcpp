#include "QualityValueScoreFunction.hpp"

int QualityScoreTypeNotSpecified(const std::string func) {
    std::cout << func << " must be specialized with FASTQ or "
              << "FASTA sequencecs." << std::endl;
    exit(1);
    return 0;
}

template<typename T_RefSequence, typename T_QuerySequence>
int QualityValueScoreFunction<T_RefSequence, T_QuerySequence>::Deletion(
    T_RefSequence &seq, DNALength refPos, T_QuerySequence &querySeq, 
    DNALength queryPos) {
    return QualityScoreTypeNotSpecified("Deletion");
}
		
template<typename T_RefSequence, typename T_QuerySequence>
int QualityValueScoreFunction<T_RefSequence, T_QuerySequence>::Deletion(
    T_RefSequence &seq, DNALength pos) {
    return QualityScoreTypeNotSpecified("Deletion");
}

template<typename T_RefSequence, typename T_QuerySequence>
int QualityValueScoreFunction<T_RefSequence, T_QuerySequence>::Match(
    T_RefSequence &ref, DNALength refPos, T_QuerySequence &query, 
    DNALength queryPos) {
    return QualityScoreTypeNotSpecified("Match");
}

template<typename T_RefSequence, typename T_QuerySequence>
int QualityValueScoreFunction<T_RefSequence, T_QuerySequence>::Insertion(
    T_RefSequence &ref, DNALength refPos, T_QuerySequence &seq, 
    DNALength pos) {
    return QualityScoreTypeNotSpecified("Insertion");
}

template<typename T_RefSequence, typename T_QuerySequence>
int QualityValueScoreFunction<T_RefSequence, T_QuerySequence>::Insertion(
    T_QuerySequence &seq, DNALength pos) {
    return QualityScoreTypeNotSpecified("Insertion");
}

/*
 * Define all specializations for a FASTA reference and FASTQSequence 
 * for the query, or FASTA sequence for query.
 */
template<>
int QualityValueScoreFunction<FASTASequence, FASTQSequence>::Deletion(
    FASTASequence &ref, DNALength pos) {
    (void)(ref); (void)(pos);
    return del; // For now there is no global deletion quality value.
}

template<>
int QualityValueScoreFunction<DNASequence, FASTQSequence>::Deletion(
    DNASequence &ref, DNALength pos) {
    (void)(ref); (void)(pos);
    return del; // For now there is no global deletion quality value.
}

template<>
int QualityValueScoreFunction<DNASequence, FASTQSequence>::Deletion(
    DNASequence &seq, DNALength refPos, FASTQSequence &querySeq,
    DNALength queryPos) {
    (void)(querySeq); (void)(queryPos);
    return Deletion(seq, refPos);
}

template<>
int QualityValueScoreFunction<DNASequence, FASTQSequence>::Insertion(
    FASTQSequence &query, DNALength pos) {
    // Positive value for quality value penalizes the alignment.
    return query.qual[pos];
}

template<>
int QualityValueScoreFunction<DNASequence, FASTQSequence>::Insertion(
    DNASequence &ref, DNALength refPos, FASTQSequence &query, 
    DNALength pos) {
    (void)(ref); (void)(refPos); (void)(query); (void)(pos);
    // Positive value for quality value penalizes the alignment.
    // return query.qual[pos];
    // return Insertion(query, pos);
    return ins;
}

template<>
int QualityValueScoreFunction<DNASequence, FASTQSequence>::Match(
    DNASequence &ref, DNALength refPos, FASTQSequence &query, 
    DNALength queryPos) {
    // Positive value for quality value penalizes the alignment.
    return (QVDistanceMatrix[ThreeBit[query.seq[queryPos]]][ThreeBit[ref.seq[refPos]]]
            * query.qual[queryPos]);
}

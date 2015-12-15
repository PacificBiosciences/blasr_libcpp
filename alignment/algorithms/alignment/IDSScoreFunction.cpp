#include "IDSScoreFunction.hpp"

float SumAsValidPhred(float v1, float v2, float v3) {
    float sum = 0;
    if (v1 > 0) {
        sum = std::pow(10,v1/-10.0);
    }
    if (v2 > 0) {
        sum += std::pow(10,v2/-10.0);
    }
    if (v3 > 0) {
        sum += std::pow(10,v3/-10.0);
    }
    return sum;
}

template<>
int IDSScoreFunction<DNASequence,FASTQSequence>::Deletion(DNASequence &ref, DNALength refPos,
        FASTQSequence &query, DNALength queryPos) {

    if (query.deletionQV.Empty() == false and query.deletionTag != NULL) {
        if (query.deletionTag[queryPos] == 'N') {
            return globalDeletionPrior; //query.deletionQV[queryPos] ;
        }
        else {
            if (query.deletionTag[queryPos] == ref.seq[refPos]) {
                return query.deletionQV[queryPos];
            }
            else {
                return globalDeletionPrior; 
            }
        }
    }
    else {
        return del;
    }
}

template<>
int IDSScoreFunction<DNASequence, FASTQSequence>::Deletion(FASTQSequence &query, DNALength queryPos) {
    return query.deletionQV[queryPos];
}


template<>
int IDSScoreFunction<DNASequence, DNASequence>::Deletion(DNASequence &query, DNALength pos) {
    (void)(query);
    (void)(pos);
    return del; // For now there is no global deletion quality value.
}

template<>
int IDSScoreFunction<DNASequence, FASTQSequence>::Insertion(DNASequence &refSeq, DNALength refPos, 
        FASTQSequence &query, DNALength pos) {
    (void)(refSeq); (void)(refPos);
    return query.insertionQV[pos] ;
}

template<>
int IDSScoreFunction<DNASequence, FASTQSequence>::Insertion(FASTQSequence &query, DNALength pos) {
    return query.insertionQV[pos];
}

template<>
int IDSScoreFunction<DNASequence, FASTQSequence>::Match(DNASequence &ref, DNALength refPos, 
        FASTQSequence &query, DNALength queryPos) {

    if (query.seq[queryPos] == ref.seq[refPos]) {
        return 0;
    } else if (query.substitutionTag != NULL) {
        if (query.substitutionTag[queryPos] == ref.seq[refPos]) {
            return query.substitutionQV[queryPos];
        }
    }
    return substitutionPrior;
}


template<>
float IDSScoreFunction<DNASequence, FASTQSequence>::NormalizedMatch(DNASequence &ref, DNALength refPos, 
        FASTQSequence &query, DNALength queryPos) {
    /*
     * Return the match probability normalized such that the probability
     * of transitioning from refPos, queryPos is 1.
     */

    float matchScore = Match(ref, refPos, query, queryPos);

    float delScore   = -1;
    if (refPos > 0) {
        delScore = Deletion(ref, refPos-1, query, queryPos);
    }

    float insScore = -1;
    if (queryPos > 0) {
        insScore = Insertion(ref, refPos, query, queryPos-1);
    }

    float sumScore = SumAsValidPhred(matchScore, delScore, insScore);
    if (sumScore > 0) {
        float numerator = std::pow(10, matchScore/-10.0);
        return -10*std::log10( numerator / sumScore);
    }
    else {
        return 0;
    }
}



template<>
float IDSScoreFunction<DNASequence, FASTQSequence>::NormalizedInsertion(DNASequence &ref, DNALength refPos, 
        FASTQSequence &query, DNALength queryPos) {

    float insScore = Insertion(ref, refPos, query, queryPos);

    float delScore = -1;
    float matchScore = -1;  
    if (refPos < ref.length - 1) {
        matchScore = Match(ref, refPos + 1, query, queryPos);
        if (queryPos > 0) {
            delScore = Deletion(ref, refPos + 1, query, queryPos - 1);
        }
    }

    float sum = SumAsValidPhred(insScore, delScore, matchScore);
    if (sum > 0) {
        float numerator = std::pow(10,insScore/-10.0);
        return -10*std::log10( numerator / sum);
    }
    else {
        return 0;
    }
}


template<>
float IDSScoreFunction<DNASequence, FASTQSequence>::NormalizedDeletion(DNASequence &ref, DNALength refPos, 
        FASTQSequence &query, DNALength queryPos) {

    float delScore = Deletion(ref, refPos, query, queryPos);

    float matchScore = -1;
    float insScore = -1;

    if (queryPos < query.length - 1) {
        matchScore = Match(ref, refPos, query, queryPos + 1);

        if (refPos > 0) {
            insScore = Insertion(ref, refPos - 1, query, queryPos + 1);
        }
    }
    float sum = SumAsValidPhred(delScore, matchScore, insScore);
    if (sum > 0) {
        float numerator= std::pow(10, delScore/-10.0);
        return -10*std::log10( numerator / sum);
    }
    else {
        return 0;
    }
}

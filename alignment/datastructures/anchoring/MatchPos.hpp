#ifndef _BLASR_MATCH_POS_HPP_
#define _BLASR_MATCH_POS_HPP_

#include <vector>
#include <algorithm>
#include <ostream>
#include "Types.h"
#include "DNASequence.hpp"

class MatchPos {
public:
    DNALength t, q;
    MatchWeight w;
    DNALength l;
    int m; // multiplicity

    inline MatchPos(DNALength pt, DNALength pq, DNALength pl, int pm = 0); 

    MatchPos (const MatchPos &rhs); 

    inline MatchWeight Size(); 

    inline MatchPos(); 

    inline MatchPos& operator=(const MatchPos &rhs); 

    inline DNALength GetLength() const; 
    
    int GetMultiplicity() const; 
    
    MatchWeight GetWeight() const; 

    inline DNALength GetX() const; 
    
    inline DNALength GetY() const; 
    
    UInt GetT(); 

    UInt GetQ(); 

    UInt GetW(); 

    friend std::ostream& operator<<(std::ostream & out, MatchPos &p); 
};

inline MatchPos::MatchPos(DNALength pt, DNALength pq, DNALength pl, int pm) :
    t(pt), q(pq), w(0), l(pl), m(pm){ }

inline MatchPos::MatchPos() {
    t = q = -1;
    w = 0;
    l = 0;
    m = 0;
}

inline MatchWeight MatchPos::Size() {
    return l;
}

inline MatchPos& MatchPos::operator=(const MatchPos &rhs) {
    t = rhs.t; q = rhs.q; w = rhs.w;
    l = rhs.l;
    m = rhs.m;
    return *this;
}

inline DNALength MatchPos::GetLength() const {
    return l;
}

inline DNALength MatchPos::GetX() const {
    return q;
}

inline DNALength MatchPos::GetY() const {
    return t;
}


class ChainedMatchPos : public MatchPos {
private:
    int score;
    ChainedMatchPos *chainPrev;

public:
    ChainedMatchPos(DNALength pt, DNALength pq, DNALength pl, int pm);

    ChainedMatchPos();

    //inline ChainedMatchPos(const ChainedMatchPos &rhs); 

    inline int GetScore();

    int SetScore(int _score); 

    ChainedMatchPos* SetChainPrev(ChainedMatchPos *_chainPrev); 

    inline ChainedMatchPos* GetChainPrev();

    //inline ChainedMatchPos &operator=(const ChainedMatchPos &rhs); 

    friend std::ostream& operator<<(std::ostream & out, ChainedMatchPos &p); 

};

/*
inline ChainedMatchPos::ChainedMatchPos(const ChainedMatchPos &rhs) {
    (*this) = rhs;
}*/

inline ChainedMatchPos* ChainedMatchPos::GetChainPrev() {
    return chainPrev;
}

inline int ChainedMatchPos::GetScore() {
    return score;
}

/*
inline ChainedMatchPos &ChainedMatchPos::operator=(const ChainedMatchPos &rhs) {
    score = rhs.score;
    ((MatchPos&)(*this)) = ((MatchPos&)rhs);
    return *this;
}
*/


template<typename T_MatchPos>
class CompareMatchPos {
    public:
        int operator()(const T_MatchPos &lhs, const T_MatchPos &rhs) const {
            if (lhs.t < rhs.t) 
                return 1;
            else if (lhs.t > rhs.t)
                return 0;
            else { 
                return lhs.q < rhs.q;
            }
        }
};

typedef std::vector<MatchPos> MatchPosList;

template<typename T_MatchPos>
class CompareMatchPosByWeight {
    public:
        int operator()(const T_MatchPos &a, const T_MatchPos &b) const {
            return a.l < b.l;
        }
};

template<typename T_MatchPos>
class CompareMatchPosIndexByWeight {
    public:
        std::vector<T_MatchPos> *list;
        int operator()(const int i, const int j) const {
            return ((*list)[i].w > (*list)[j].w);
        }
};


template<typename T_MatchPos>
class CompareMatchPosIndexByTextPos {
    public:
        std::vector<T_MatchPos> *list;
        int operator()(const int i, const int j) const {
            return (*list)[i].t < (*list)[j].t;
        }
};


template<typename T_MatchPos>
void SortMatchPosList(std::vector<T_MatchPos> &mpl) {
    std::sort(mpl.begin(), mpl.end(), CompareMatchPos<T_MatchPos>());
}

template<typename T_MatchPos>
void SortMatchPosListByWeight(std::vector<T_MatchPos> &mpl) {
    std::sort(mpl.begin(), mpl.end(), CompareMatchPosByWeight<T_MatchPos>());
}

template<typename T_MatchPos>
void SortMatchPosIndexListByWeight(std::vector<T_MatchPos> &mpl, std::vector<int> &indices) {
    CompareMatchPosIndexByWeight<T_MatchPos> cmp;
    cmp.list = &mpl;
    std::sort(indices.begin(), indices.end(), cmp);
}

template<typename T_MatchPos>
void SortMatchPosIndexListByTextPos(std::vector<T_MatchPos> &mpl, std::vector<int> &indices) {
    CompareMatchPosIndexByTextPos<T_MatchPos> cmp;
    cmp.list = &mpl;
    std::sort(indices.begin(), indices.end(), cmp);
}

#endif  // _BLASR_MATCH_POS_HPP_

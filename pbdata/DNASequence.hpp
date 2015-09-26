#ifndef  _BLASR_DNA_SEQUENCE_HPP_
#define  _BLASR_DNA_SEQUENCE_HPP_

#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <cassert>
#include "Types.h"
#include "NucConversion.hpp"
#include "utils.hpp"
#include "libconfig.h"

#ifdef USE_PBBAM
#include "pbbam/BamRecord.h"
#endif


class DNASequence {
public:
    DNALength length;
    Nucleotide *seq;
    int bitsPerNuc;
    bool deleteOnExit;

    inline DNASequence();
    inline ~DNASequence();

    //--- functions ---//
    
    DNALength size();

    inline void CheckBeforeCopyOrReference(const DNASequence & rhs, std::string seqType = "DNASequence"); 

    void TakeOwnership(DNASequence &rhs);

    void Append(const DNASequence &rhs, DNALength appendPos=0);

    DNASequence &Copy(const DNASequence &rhs, DNALength rhsPos=0, DNALength rhsLength=0);

    void ShallowCopy(const DNASequence &rhs);

    DNASequence & Copy(const std::string & rhs);

    int GetStorageSize() const;

    DNASequence &operator=(const DNASequence &rhs);

    DNASequence &operator=(const std::string &rhs);

    void Print(std::ostream &out, int lineLength = 50) const;

    void PrintSeq(std::ostream &out, int lineLength = 50) const;

    void Allocate(DNALength plength);

    void ReferenceSubstring(const DNASequence &rhs, DNALength pos=0, DNALength substrLength=0); 

    DNALength MakeRCCoordinate(DNALength forPos );

    inline void CopyAsRC(DNASequence &rc, DNALength pos=0, DNALength rcLength =0);

    void MakeRC(DNASequence &rc, DNALength pos=0, DNALength rcLength=0);

    void ToTwoBit();

    inline void ToThreeBit();

    void ToFourBit();

    void ConvertThreeBitToAscii();

    void ToAscii();

    void Assign(DNASequence &ref, DNALength start=0, DNALength plength=0);

    void ToLower();

    void ToUpper(); 

    void Concatenate(const Nucleotide *moreSeq, DNALength moreSeqLength); 

    std::string GetTitle() const; 

    void Concatenate(const Nucleotide* moreSeq);

    void Concatenate(DNASequence &seq); 

    int Compare(DNALength pos, DNASequence &rhs, DNALength rhsPos, DNALength length); 

    int LessThanEqual(DNALength pos, DNASequence &rhs, DNALength rhsPos, DNALength length); 

    int Equals(DNASequence &rhs, DNALength rhsPos, DNALength length, DNALength pos=0 ); 

    int LessThan(DNALength pos,  DNASequence &rhs, DNALength rhsPos, DNALength length); 

    void CleanupASCII(); 

    Nucleotide operator[](int i) {
        return seq[i];
    }

    Nucleotide GetNuc(DNALength i) const; 

    DNALength GetRepeatContent() const; 

    void CleanupOnFree();

    virtual void Free(); 

    void Resize(DNALength newLength);

    DNALength GetSeqStorage() const;

#ifdef USE_PBBAM
    /// Copies a BamRecord as a DNASequence.
    DNASequence & Copy(const PacBio::BAM::BamRecord & record);
#endif 
};

inline DNASequence::DNASequence() {
    seq = NULL;
    length = 0;
    bitsPerNuc = 8;
    deleteOnExit = false;
}

inline DNASequence::~DNASequence() {
    DNASequence::Free();
}


// Sanity check:
// If this DNASequence and rhs are pointing to the same seq 
// in memory, make sure this DNASequence's deleteOnExit is false.
// (otherwise, seq in memory will be deleted before being copied)
inline void DNASequence::CheckBeforeCopyOrReference(const DNASequence & rhs, std::string seqType) {

    if (seq == rhs.seq and seq != NULL and deleteOnExit) {
        std::cout << "ERROR, trying to copying a " << seqType << " to itself." << std::endl;
        exit(1);
    }
}

inline void DNASequence::ToThreeBit() {
    DNALength i;
    if (bitsPerNuc != 3) 
        for (i = 0; i < length; i++) { seq[i] = ThreeBit[seq[i]]; }
    bitsPerNuc = 3;
}

inline void DNASequence::CopyAsRC(DNASequence &rc, DNALength pos, DNALength rcLength) {
    // Free rc before copying any data to rc.
    ((DNASequence&)rc).Free();
    //
    // Different way of accounting for position. The position is on
    // the rc strand, not the forward strand.
    //
    if (rcLength == 0) {
        rcLength = length - pos;
    }
    DNALength rcStart = length - (pos + rcLength);
    ((DNASequence&)rc).Resize(rcLength);
    DNALength i;
    for (i = 0; i < rcLength; i++) {
        rc.seq[i] = ReverseComplementNuc[seq[rcStart - 1 + (rcLength - i)]];
    }

    // The reverse complement controls its own memory now.
    rc.deleteOnExit = true;
}


template<typename T>
DNALength ResizeSequence(T &dnaseq, DNALength newLength) {
    assert(newLength > 0);
    ((T&)dnaseq).Free();
    dnaseq.seq = ProtectedNew<Nucleotide>(newLength);
    dnaseq.length = newLength;
    dnaseq.deleteOnExit = true;
    return newLength;
}

#endif // _BLASR_DNA_SEQUENCE_HPP_

#include <cassert>
#include <cstring>
#include "Types.h"
#include "DNASequence.hpp"

DNALength DNASequence::size() {
    return length;
}

void DNASequence::TakeOwnership(DNASequence &rhs) {
    CheckBeforeCopyOrReference(rhs);
    // Free this DNASequence before take owner ship from rhs.
    DNASequence::Free();

    seq = rhs.seq;
    length = rhs.length;
    deleteOnExit = rhs.deleteOnExit;

    rhs.deleteOnExit = false;
}

void DNASequence::Append(const DNASequence &rhs, DNALength appendPos) {
    assert(deleteOnExit); // must have control over seq.
    //
    // Simply append rhs to this seuqence, unless appendPos is nonzero
    // in which case rhs is inserted at attendPos, overwriting this
    // sequence from appendPos to the end.
    //
    Nucleotide *newSeq;
    //
    // Handle the two cases (appendPos == 0 and appendPos > 0)
    // separately in order to handle memory deallocation correctly.
    //
    if (appendPos == 0) {
        DNALength  newSeqLength = length + rhs.length;
        newSeq = ProtectedNew<Nucleotide>(newSeqLength);
        memcpy(newSeq, seq, length);
        memcpy(&newSeq[length], rhs.seq, rhs.length);

        if (length != 0) {
            delete[] seq;
        }
        seq = newSeq;
        length = newSeqLength;
    }
    else {
        if (appendPos + rhs.length < length) {
            memcpy(&seq[appendPos], rhs.seq, rhs.length);
            length = appendPos + rhs.length;
        }
        else {
            DNALength lengthCopy = length;
            length = appendPos;
            DNALength newSeqLength;
            newSeqLength = length + rhs.length;
            newSeq = ProtectedNew<Nucleotide>(newSeqLength);
            memcpy(newSeq, seq, length);
            memcpy(&newSeq[length], rhs.seq, rhs.length);
            if (deleteOnExit and lengthCopy != 0) {
                delete[] seq;
            }
            seq = newSeq;
            length = newSeqLength;
        }
    }
    deleteOnExit = true;
}

// Copy FROM rhs to this DNASequence. 
DNASequence& DNASequence::Copy(const DNASequence &rhs, DNALength rhsPos, DNALength rhsLength) {
    CheckBeforeCopyOrReference(rhs);
    // Free this DNASequence before copying from rhs
    DNASequence::Free();

    //
    // When initializing a vector of DNASequence's, the copy
    // constructor will initialze a list and call this
    // function with a zero-length DNASequence as the rhs to
    // initialize every element in the vector   The check
    // below will fail on zero-length sequences, so add a boundary
    // condition check before that to allow the copy-constructor to
    // work.
    //
    if (rhs.length == 0) {
        seq = NULL;
        length = 0;
        deleteOnExit = true;
        return *this;
    }

    //
    // Silently ignoring this case could lead to problems later on,
    // catastrophically assert here if the input is not valid.
    // In case rhsLength + rhsPos > ULONG_MAX (4294967295), check 
    // both rhsLength and rhsPos, fix bug 21794
    //
    if (not (rhsLength <= rhs.length     && 
             rhsPos    <= rhs.length + 1 &&
             rhsLength + rhsPos <= rhs.length + 2 )) {
        std::cout << "ERROR.  The subsequence to copy is out of bounds." 
            << std::endl
            << "        Failed to copy a subsequence starting at " << rhsPos 
            << std::endl
            << "        with length "<< rhsLength 
            << " from a sequence of length " << rhs.length << "." 
            << std::endl;
        exit(1);
    }

    if (rhsLength == 0) {
        rhsLength = rhs.length - rhsPos;
    }
    if (rhsLength == 0) {
        seq = NULL;
    }
    else {
        seq = ProtectedNew<Nucleotide>(rhsLength);
        memcpy(seq, &rhs.seq[rhsPos], rhsLength);
    }
    length = rhsLength;
    deleteOnExit = true;
    return *this;
}

DNASequence& DNASequence::Copy(const std::string &rhsSeq) {
    // Free this DNASequence before copying from rhs
    DNASequence::Allocate(static_cast<DNALength>(rhsSeq.size()));
    memcpy(seq, rhsSeq.c_str(), length);
    return *this;
}

DNASequence & DNASequence::operator=(const std::string & rhsSeq) {
    return DNASequence::Copy(rhsSeq);
}

void DNASequence::ShallowCopy(const DNASequence &rhs) {
    seq = rhs.seq;
    length = rhs.length;
    deleteOnExit = false;
}

int DNASequence::GetStorageSize() const {
    return (length * sizeof(Nucleotide));
}

DNASequence &DNASequence::operator=(const DNASequence &rhs){ 
    DNASequence::Copy(rhs);
    return *this;
}

//
// synonym for printseq
//
void DNASequence::Print(std::ostream &out, int lineLength) const {
    PrintSeq(out, lineLength);
}

void DNASequence::PrintSeq(std::ostream &out, int lineLength) const {
    if (lineLength == 0) {
        std::string line;
        line.assign((char*)seq, length);
        out << line;
    }
    else {
        //
        // Make sure this isn't 
        assert(lineLength > 0);
        DNALength curPos = 0;
        int curLineLength = lineLength;
        while (curPos < length) {
            if (curPos + curLineLength > length) {
                curLineLength = length - curPos;
            }
            std::string line;
            line.assign((char*) &seq[curPos], curLineLength);
            out << line << std::endl;
            curPos += curLineLength;
        }
    }
}

void DNASequence::Allocate(DNALength plength) {
    DNASequence::Free();
    seq = ProtectedNew<Nucleotide> (plength);
    length = plength;
    deleteOnExit = true;
}

void DNASequence::ReferenceSubstring(const DNASequence &rhs, DNALength pos, DNALength substrLength) {
    CheckBeforeCopyOrReference(rhs); 

    // Free this DNASequence before referencing rhs.
    DNASequence::Free();

    //
    // This makes a reference therefore it should not be deleted.
    //
    assert(pos <= rhs.length && substrLength <= rhs.length);
    if (substrLength == 0) {
        substrLength = rhs.length - pos;
    }
    assert(pos + substrLength <= rhs.length);
    seq = &rhs.seq[pos];
    length = substrLength;
    deleteOnExit = false;
}

DNALength DNASequence::MakeRCCoordinate(DNALength forPos ) {
    return length - forPos - 1;
}

// Create a reverse complement DNAsequence of this DNASequence and assign to rc.
void DNASequence::MakeRC(DNASequence &rc, DNALength pos, DNALength rcLength) {
    if (rcLength == 0) {
        rcLength = length - pos;
    }

    ((DNASequence&)rc).Allocate(rcLength);
    DNALength i;
    for (i = 0; i < rcLength; i++) {
        rc.seq[rcLength - i - 1] = ReverseComplementNuc[seq[i+pos]];
    }
    rc.length = rcLength;
    rc.deleteOnExit = true;
}

void DNASequence::ToTwoBit() {
    DNALength i;
    for (i = 0; i < length; i++) {
        seq[i] = TwoBit[seq[i]];
    }
    bitsPerNuc = 2;
}

void DNASequence::ToFourBit() {
    DNALength i;
    if (bitsPerNuc != 4) 
        for (i = 0; i < length; i++) { seq[i] = FourBit[seq[i]]; }
    bitsPerNuc = 4;
}

void DNASequence::ConvertThreeBitToAscii() {
    DNALength i;
    for (i = 0; i < length; i++ ){
        seq[i] = ThreeBitToAscii[seq[i]];
    }
}

void DNASequence::ToAscii() {
    DNALength i;
    if (bitsPerNuc != 8) {
        for (i = 0; i < length; i++ ){ 
            seq[i] = FourBitToAscii[seq[i]];
        }
        bitsPerNuc = 8;
    }
}

// Copy rhs[start:start+plength] to this DNASequence.
void DNASequence::Assign(DNASequence &ref, DNALength start, DNALength plength) {
    CheckBeforeCopyOrReference(ref);
    // Free this DNASequence before assigning anything
    DNASequence::Free();
    
    if (plength) {
        length = plength;
        seq = ProtectedNew<Nucleotide> (length);
        memcpy(seq, &ref.seq[start], length);
    }
    else if (start) {
        length = ref.length - start;
        seq = ProtectedNew<Nucleotide> (length);
        memcpy(seq, &ref.seq[start], length);
    }
    else {
        ((DNASequence*)this)->Copy(ref);
    }
    deleteOnExit = true;
}

void DNASequence::ToLower() {
    DNALength i;
    for (i = 0; i < length; i++) {
        seq[i] = AllToLower[seq[i]];
    }
}

void DNASequence::ToUpper() {
    DNALength i;
    for (i = 0; i < length; i++) {
        seq[i] = AllToUpper[seq[i]];
    }
}

void DNASequence::Concatenate(const Nucleotide *moreSeq, DNALength moreSeqLength) {
    DNALength prevLength = length;
    length += moreSeqLength;
    Nucleotide *prev = seq;
    seq = ProtectedNew<Nucleotide> (length);
    if (prev != NULL) {
        memcpy(seq, prev, prevLength);
        delete[] prev;
    }
    memcpy((Nucleotide*) &seq[prevLength], moreSeq, moreSeqLength);
    deleteOnExit = true;
}

std::string DNASequence::GetTitle() const {
    return std::string("");
}

void DNASequence::Concatenate(const Nucleotide* moreSeq) {
    DNALength moreSeqLength = strlen((char*) moreSeq);
    Concatenate(moreSeq, moreSeqLength);
}

void DNASequence::Concatenate(DNASequence &seq) {
    Concatenate(seq.seq, seq.length);
}

int DNASequence::Compare(DNALength pos, DNASequence &rhs, DNALength rhsPos, DNALength length) {
    return memcmp(&seq[pos], &rhs.seq[rhsPos], length);
}

int DNASequence::LessThanEqual(DNALength pos, DNASequence &rhs, DNALength rhsPos, DNALength length) {
    int res = Compare(pos, rhs, rhsPos, length);
    if (res <= 0) 
        return 1;
    else
        return 0;
}

int DNASequence::Equals(DNASequence &rhs, DNALength rhsPos, DNALength length, DNALength pos) {
    int res = Compare(pos, rhs, rhsPos, length);
    return res == 0;
}

int DNASequence::LessThan(DNALength pos,  DNASequence &rhs, DNALength rhsPos, DNALength length) {
    int res=  Compare(pos, rhs, rhsPos, length);
    return (res < 0);
}

void DNASequence::CleanupASCII() {
    DNALength i;
    for (i = 0; i < length; i++ ){
        if (ThreeBit[seq[i]] == 255) {
            seq[i] = 'N';
        }
    }
}

Nucleotide DNASequence::GetNuc(DNALength i) const {
    return seq[i];
}

DNALength DNASequence::GetRepeatContent() const {
    DNALength i;
    DNALength nRepeat = 0;
    for (i =0 ; i < length;i++) {
        if (tolower(seq[i]) == seq[i]) { nRepeat++;}
    }
    return nRepeat;
}

void DNASequence::CleanupOnFree() {
    deleteOnExit = true;
}

void DNASequence::Free() {
    if (deleteOnExit == true) {
        // if has control, delete seq 
        // Otherwise, seq in memory is controlled by another object, 
        // and will be deleted later.
        if (seq != NULL) {
            delete[] seq;
        }
    } 
    // Reset seq, length and deleteOnExit
    seq = NULL;
    length = 0;
    deleteOnExit = false;
}

void DNASequence::Resize(DNALength newLength) {
    DNASequence::Free();
    seq = ProtectedNew<Nucleotide>(newLength);
    length = newLength;
    deleteOnExit = true;
}

DNALength DNASequence::GetSeqStorage() const{
    return length;
}

#ifdef USE_PBBAM
DNASequence & DNASequence::Copy(const PacBio::BAM::BamRecord & record) {
    return DNASequence::Copy(record.Sequence());
}
#endif



#include "utils.hpp"
#include "PackedDNASequence.hpp"

const PackedDNAWord PackedDNASequence::NucPosMask[] = {7, 56, 448, 
    3584, 28672, 229376, 
    1835008, 14680064, 117440512, 
    939524096};

const PackedDNAWord PackedDNASequence::NegMask[] = { 4294967288,// 11111111 11111111 11111111 11111000
    4294967239,// 11111111 11111111 11111111 11000111
    4294966847, /// and so on...
    4294963711, 4294938623, 4294737919, 
    4293132287, 4280287231, 4177526783, 
    3355443199};

// Masks from the rightmost pos 'R' all to the end of the left side of
// the word.
const PackedDNAWord PackedDNASequence::MaskRL[] = {1073741823, 1073741816, 1073741760, 
    1073741312, 1073737728, 1073709056,
    1073479680, 1071644672, 1056964608, 
    939524096};

// Masks from starting position 'L' to the rightmost 'R'
const PackedDNAWord PackedDNASequence::MaskLR[] = {7, 63, 511, 4095, 32767, 262143, 2097151, 16777215, 134217727, 1073741823};

/*
 * Count nucleotides by Xor'ing with the complement of the bit pattern. 
 */
const PackedDNAWord PackedDNASequence::xorMask[] = {
    1073741823, // mask A=000 by 111111111111111111111111111111
    920350134,  // mask C=001 by 110110110110110110110110110110
    766958445,  // mask G=010 by 101101101101101101101101101101 
    613566756,  // mask T=011 by 100100100100100100100100100100
    460175067,  // mask N=100 by 011011011011011011011011011011
    153391689,  // mask 001001001001001001001001001001 *unused*, can mask X=110
    306783378,  // mask 010010010010010010010010010010 *unused*, can mask X'=101
    0};         // mask 000000000000000000000000000000 *unused*, mask X''=111



Nucleotide PackedDNASequence::Get(DNALength pos) {
    PackedDNAWord offset = pos% NucsPerWord;
    return (seq[pos/NucsPerWord] >> (3*offset)) & NucMask ;
}

Nucleotide PackedDNASequence::operator[](DNALength pos){
    return Get(pos);
}

PackedDNASequence::PackedDNASequence() {
    nCountInWord = 0;
    nCountNuc = 0;
    length = arrayLength = 0;
    seq = NULL;
}

PackedDNASequence::~PackedDNASequence() {
    nCountInWord = 0;
    nCountNuc = 0;
    length = arrayLength = 0;
    if (seq) {
        delete [] seq; 
        seq = NULL;
    }
}

void PackedDNASequence::Allocate(DNALength numberOfNucleotides) {
    arrayLength = CeilOfFraction(numberOfNucleotides, NucsPerWord);
    length = numberOfNucleotides;
    if (seq) {delete [] seq; seq = NULL;}
    if (arrayLength > 0) {
        seq = ProtectedNew<PackedDNAWord>(arrayLength);
        std::fill(seq, seq + arrayLength, 0);
    }
}

void PackedDNASequence::CreateFromDNASequence(DNASequence &dnaSeq) {
    arrayLength = CeilOfFraction(dnaSeq.length, NucsPerWord);
    length = dnaSeq.length;
    if (seq) {delete [] seq; seq = NULL;}
    if (arrayLength > 0) {
        seq = ProtectedNew<PackedDNAWord>(arrayLength);
        DNALength pos;
        for (pos = 0; pos < dnaSeq.length; pos++) {
            Set(pos, ThreeBit[dnaSeq[pos]]);
        }
    }
}

void PackedDNASequence::Set(DNALength pos, Nucleotide threeBitValue) {
    DNALength wordPos = pos/NucsPerWord;
    DNALength wordOffset = pos%NucsPerWord;
    //
    // Pull the value to update out of memory.
    //
    PackedDNAWord word     = seq[wordPos];

    //
    // Expand the 3 bit value into a whole word.
    //
    PackedDNAWord nuc      = threeBitValue;

    word = word & NegMask[wordOffset];
    nuc  = nuc << (3*(pos%NucsPerWord));
    word = word + nuc;
    //
    // Write back the whole word.
    //
    seq[wordPos] = word;
}


DNALength PackedDNASequence::CountInWord(PackedDNAWord word, PackedDNAWord wordMask, Nucleotide nuc) {

    /*
     * Count the number of times a nucleotide (3-mer) appears in a word.
     * This is done by a series of masking and ands.
     * The packed format is a triplet of bits (b2,b1,b0), where N is a
     * mask bit, and n1n0 specifies the nucleotide.  
     *
     * The sequence CCGN appears as
     * 001001010100
     * 
     * To count C's: perform xor mask with a complement mask that will
     * produce a triplet of 111 on every position where there is a C,
     * followed by an and with a stride-isolation mask.
     *
     *     b0              b1            b2            
     *     001001010100	   001001010100  001001010100  
     * xor 110110110110	   110110110110  110110110110
     *     ------------	   ------------  ------------  
     *     111111100100	   111111100100  111111100100

     * and 001001001001	   010010010010  100100100100  
     *     ------------	   ------------  ------------  
     *     001001000000	   010010000010  000000000100  
     *
     * Shift b1 by 1 to get parity with b0.
     * 010010000010 -> 001001000001
     *                          001001000001
     * And with b0 pattern      001001000000
     * to count masked nucs.    ------------
     *               nuc_count= 001001000000 
     * Shift the b2 by 2 to set up the and-gate.
     * 000000000100 -> 000000000001
     * and ~MASK with nuc_count to get rid of masked columns if the
     * whole word is not being calculated.
     * 111111111110 & 001001000000 = 001001000000 
     * Finally, use 64-bit multiplication bit-hack to count the number
     * of set-bits:
     * CountInWord(001001000000) = 2
     *
     */
    PackedDNAWord w0,w1,w2, w;
    Nucleotide tbn = ThreeBit[nuc];
    PackedDNAWord xorMaskNuc = xorMask[tbn];
    w0 = w1 = w2 = (word ^ xorMaskNuc);
    w0 = (w0)  & Mask0All;
    w1 = ((w1) & Mask1All) >> 1;
    w2 = ((w2) & Mask2All) >> 2;
    // Ideally the architecture will parallelize all these
    // Do a cascaded.
    //		w01 = (w0 & w1);
    //		w12 = (w1 & w2);
    //		w = w01 & w12;
    w = ((w0 & w1) & w2) & wordMask;
    //		++nCountInWord;
    return CountBits(w);
}

DNALength PackedDNASequence::CountNuc(DNALength start, DNALength end, Nucleotide nuc) {
    DNALength startWordIndex, endWordIndex, wordIndex;
    DNALength startInWord, endInWord;
    endInWord = NucsPerWord;
    startInWord = start % NucsPerWord;
    startWordIndex = start / NucsPerWord; 
    endWordIndex   = end   / NucsPerWord;

    // 
    // Process all whole words.
    //
    DNALength nNuc = 0;
    for (wordIndex = startWordIndex; wordIndex < endWordIndex; wordIndex++) {
        endInWord = NucsPerWord;
        nNuc += CountInWord(seq[wordIndex] & MaskRL[startInWord], MaskRL[startInWord], nuc);
        startInWord = 0;
    }
    /* 
     * Look to see if there is extra sequence to process when the seq
     * does not end on a word boundary.
     */
    if (end % NucsPerWord != 0) {
        endInWord = end % NucsPerWord;
        nNuc += CountInWord(seq[wordIndex] & MaskRL[startInWord] & MaskLR[endInWord-1], 
                MaskRL[startInWord] & MaskLR[endInWord-1], 
                nuc);
    }
    //		++nCountNuc;
    return nNuc;
}

void PackedDNASequence::Write(std::ostream &out) {
    out.write((char*)&arrayLength, sizeof(arrayLength));
    out.write((char*)&length, sizeof(length));
    if (arrayLength > 0) {
        out.write((char*)seq, sizeof(PackedDNAWord)*arrayLength);
    }
}

void PackedDNASequence::Read(std::istream &in) {
    in.read((char*)&arrayLength, sizeof(arrayLength));
    in.read((char*)&length, sizeof(length));
    if (seq) {delete [] seq; seq = NULL;}
    if (arrayLength > 0) {
        seq = ProtectedNew<PackedDNAWord>(arrayLength);
        in.read((char*)seq, sizeof(PackedDNAWord)*arrayLength);
    }
}

void PackedDNASequence::PrintUnpacked(std::ostream &out, int lineLength) {
    DNALength p;
    for (p = 0; p < length; p++) {
        out << (char) ThreeBitToAscii[Get(p)];
        if (static_cast<int>(p % lineLength) == lineLength-1) {
            out << std::endl;
        }
    }
    if (p % lineLength != 0) {
        out << std::endl;
    }
}

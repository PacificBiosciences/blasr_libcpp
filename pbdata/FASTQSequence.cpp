#include <cassert>
#include <vector>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdint.h>
#include "Types.h"
#include "NucConversion.hpp"
#include "DNASequence.hpp"
#include "Enumerations.h"
#include "FASTQSequence.hpp"

using namespace std;

//
// Initialize a read with quality probabilities from one with quality values.
//
int FASTQSequence::charToQuality = FASTQ_CHAR_TO_QUALITY;

QVScale FASTQSequence::GetQVScale() const {
    return qvScale;
}

void FASTQSequence::SetQVScale(QVScale qvScaleP) {
    qvScale                   = qvScaleP;
    qual.qvScale              = qvScale;
    deletionQV.qvScale        = qvScale;
    preBaseDeletionQV.qvScale = qvScale;
    insertionQV.qvScale       = qvScale;
    substitutionQV.qvScale    = qvScale;
    mergeQV.qvScale           = qvScale;
}

QualityValueVector<QualityValue>* 
FASTQSequence::GetQVPointerByIndex(int index) {
    if (index == 0) { return &qual; }
    if (index == 1) { return &insertionQV; }
    if (index == 2) { return &deletionQV; }
    if (index == 3) { return &substitutionQV; }
    if (index == 4) { return &mergeQV; }
    return NULL;
}

int FASTQSequence::GetStorageSize() const {
    int total = 0;
    int nQV = 0;
    int nTag =0;
    if (!qual.Empty()) {
        nQV++;
    }
    if (!deletionQV.Empty()) { 
        nQV++;
    }
    if (!preBaseDeletionQV.Empty()) {
        nQV+=4;
    }
    if (!insertionQV.Empty()) {
        nQV++;
    }
    if (!substitutionQV.Empty()) {
        nQV++;
    }
    if (!mergeQV.Empty()) {
        nQV++;
    }
    if (deletionTag != NULL) {
        nTag++;
    }
    if (substitutionTag !=NULL) {
        nTag++;
    }
    total = nQV*sizeof(QualityValue)*length + nTag*sizeof(Nucleotide)*length;
    return total + FASTASequence::GetStorageSize();
}

FASTQSequence::FASTQSequence() : FASTASequence() {
    deletionTag       = NULL;
    substitutionTag   = NULL;

    //
    // For now assume a prior distribution to be the variation of the human genome.
    // FIXME: these were set to 0.001, which ends up being 0 because these priors are integer types
    //        I'm setting these explicitly to 0 to silence the warning and maintain behavior,
    //        but mkinsella recommends revisiting these for potential removal
    //
    deletionQVPrior = 0;
    insertionQVPrior = 0;
    substitutionQVPrior = 0;
    preBaseDeletionQVPrior = 0;
    qvScale = PHRED;
}

QualityValue FASTQSequence::GetDeletionQV(DNALength pos) const {
    assert(pos < ((unsigned int)-1));
    assert(pos < length);
    if (deletionQV.Empty()) {
        return deletionQVPrior;
    }
    else {
        return deletionQV[pos];
    }
}

QualityValue FASTQSequence::GetMergeQV(DNALength pos) const {
    assert(pos < ((unsigned int)-1));
    assert(pos < length);
    if (mergeQV.Empty()) {
        return 0;
    }
    else {
        return mergeQV[pos];
    }
}

Nucleotide FASTQSequence::GetSubstitutionTag(DNALength pos) const {
    if (substitutionTag == NULL) {
        return 'N';
    }
    assert(pos < ((unsigned int)-1));
    assert(pos < length);
    return substitutionTag[pos];
}

Nucleotide FASTQSequence::GetDeletionTag(DNALength pos) const {
    if (deletionTag == NULL) {
        return 'N';
    }
    assert(pos < ((unsigned int)-1));
    assert(pos < length);
    return deletionTag[pos];
}

QualityValue FASTQSequence::GetInsertionQV(DNALength pos) const {
    if (insertionQV.Empty()) {
        return insertionQVPrior;
    }
    assert(pos < ((unsigned int)-1));
    assert(pos < length);
    return insertionQV[pos];
}

QualityValue FASTQSequence::GetSubstitutionQV(DNALength pos) const {
    if (substitutionQV.Empty()) {
        return substitutionQVPrior;
    }
    assert(pos < ((unsigned int)-1));
    assert(pos < length);
    return substitutionQV[pos];
}

QualityValue FASTQSequence::GetPreBaseDeletionQV(DNALength pos, Nucleotide nuc) const {
    if (preBaseDeletionQV.Empty()) {
        return preBaseDeletionQVPrior;
    }
    assert(pos < ((unsigned int)-1));
    assert(pos < length);
    return preBaseDeletionQV[pos*4 + TwoBit[nuc]];
}

void FASTQSequence::ShallowCopy(const FASTQSequence &rhs) {
    CheckBeforeCopyOrReference(rhs, "FASTQSequence");
    FASTQSequence::Free();

    qual.ShallowCopy(rhs.qual, 0, length);
    FASTASequence::ShallowCopy(rhs);
}

void FASTQSequence::ReferenceSubstring(const FASTQSequence &rhs) {
    FASTQSequence::ReferenceSubstring(rhs, 0, rhs.length);
}

void FASTQSequence::ReferenceSubstring(const FASTQSequence &rhs, DNALength pos) {
    FASTQSequence::ReferenceSubstring(rhs, pos, rhs.length - pos);
}

void FASTQSequence::ReferenceSubstring(const FASTQSequence &rhs, DNALength pos, DNALength substrLength) {
    // Sanity check.
    CheckBeforeCopyOrReference(rhs, "FASTQSequence");

    // Free this FASTQSequence before referencing rhs.
    FASTQSequence::Free();

    SetQVScale(rhs.qvScale);
    if (substrLength == 0) {
        substrLength = rhs.length - pos;
    }
    FASTASequence::ReferenceSubstring(rhs,pos,substrLength);
    if (rhs.qual.Empty() == false) {
        qual.ShallowCopy(rhs.qual, pos, substrLength);
    }
    if (rhs.deletionQV.Empty() == false) {
        deletionQV.ShallowCopy(rhs.deletionQV, pos, substrLength);
    }
    if (rhs.mergeQV.Empty() == false) {
        mergeQV.ShallowCopy(rhs.mergeQV, pos, substrLength);
    }
    if (rhs.insertionQV.Empty() == false) {
        insertionQV.ShallowCopy(rhs.insertionQV, pos, substrLength);
    }
    if (rhs.preBaseDeletionQV.Empty() == false ){
        preBaseDeletionQV.ShallowCopy(rhs.preBaseDeletionQV, pos, substrLength);
    }
    if (rhs.deletionTag != NULL) {
        deletionTag = &rhs.deletionTag[pos];
    }
    if (rhs.substitutionTag != NULL) {
        substitutionTag = &rhs.substitutionTag[pos];
    }
    if (rhs.substitutionQV.Empty() == false) {
        substitutionQV.ShallowCopy(rhs.substitutionQV, pos, substrLength);
    }
    deletionQVPrior = rhs.deletionQVPrior;
    insertionQVPrior = rhs.insertionQVPrior;
    substitutionQVPrior = rhs.substitutionQVPrior;
    preBaseDeletionQVPrior = rhs.preBaseDeletionQVPrior;
}

void FASTQSequence::ClearAndNull(QualityValue *value) {
    if (value != NULL) {
        delete[] value;
    }
    value = NULL;
}

void FASTQSequence::CopyQualityValues(const FASTQSequence &rhs) {
    // Make sure QVs and seq are all under control, if seq is referenced
    // while QVs are copied, memory leak can happen. 
    assert(deleteOnExit);

    SetQVScale(rhs.qvScale);
    qual.Copy(rhs.qual, rhs.length);
    deletionQV.Copy(rhs.deletionQV, rhs.length);
    insertionQV.Copy(rhs.insertionQV, rhs.length);
    substitutionQV.Copy(rhs.substitutionQV, rhs.length);
    mergeQV.Copy(rhs.mergeQV, rhs.length);
    //
    // Handle the tags separtely (and verbosely)
    //
    if (rhs.deletionTag) {
        AllocateDeletionTagSpace(rhs.length);
        memcpy(deletionTag, rhs.deletionTag, sizeof(Nucleotide)*rhs.length);
    }
    else { 
        ClearAndNull(deletionTag);
    }

    if (rhs.substitutionTag) {
        AllocateSubstitutionTagSpace(rhs.length);
        memcpy(substitutionTag, rhs.substitutionTag, sizeof(Nucleotide)*rhs.length);
    }
    else {
        ClearAndNull(substitutionTag);
    }
}

void FASTQSequence::AllocateQualitySpace(DNALength qualLength) {
    qual.Allocate(qualLength);
}

void FASTQSequence::AllocateDeletionQVSpace(DNALength qualLength) {
    deletionQV.Allocate(qualLength);
}

void FASTQSequence::AllocateMergeQVSpace(DNALength len) {
    mergeQV.Allocate(len);
}

void FASTQSequence::AllocateDeletionTagSpace(DNALength qualLength) {
    if (deletionTag != NULL) delete[] deletionTag;
    deletionTag = ProtectedNew<Nucleotide>(qualLength);
}

void FASTQSequence::AllocatePreBaseDeletionQVSpace(DNALength qualLength) {
    preBaseDeletionQV.Allocate(qualLength);
}

void FASTQSequence::AllocateInsertionQVSpace(DNALength qualLength) {
    insertionQV.Allocate(qualLength);
}

void FASTQSequence::AllocateSubstitutionQVSpace(DNALength qualLength ){ 
    substitutionQV.Allocate(qualLength);
}

void FASTQSequence::AllocateSubstitutionTagSpace(DNALength qualLength ){ 
    if (substitutionTag != NULL) delete[] substitutionTag;
    substitutionTag = ProtectedNew<Nucleotide>(qualLength);
}

void FASTQSequence::AllocateRichQualityValues(DNALength qualLength) {
    AllocateDeletionQVSpace(qualLength);
    AllocateDeletionTagSpace(qualLength);
    AllocatePreBaseDeletionQVSpace(qualLength);
    AllocateInsertionQVSpace(qualLength);
    AllocateSubstitutionQVSpace(qualLength);
    AllocateSubstitutionTagSpace(qualLength);
    AllocateMergeQVSpace(qualLength);
}

void FASTQSequence::Copy(const FASTQSequence &rhs) {
    CheckBeforeCopyOrReference(rhs, "FASTQSequence");

    // Free *this before copying anything.
    FASTQSequence::Free();

    // Copy FASTASequence from rhs, including seq and title
    FASTASequence::Copy(rhs); 

    assert(deleteOnExit);

    // Copy Quality values from rhs.
    FASTQSequence::CopyQualityValues(rhs);
}

FASTQSequence& FASTQSequence::operator=(const FASTQSequence &rhs) {
    ((FASTQSequence*)this)->Copy(rhs);
    return *this;
}

FASTQSequence::FASTQSequence(const FASTQSequence &rhs) {
    ((FASTQSequence*)this)->Copy(rhs);
}

// Copy rhs to this, including seq, title and QVs.
void FASTQSequence::Assign(FASTQSequence &rhs) {
    CheckBeforeCopyOrReference(rhs);
    FASTQSequence::Free();

    // copy the nucleotide part
    FASTASequence::Assign(rhs);
    // copy the qual part, qual scal is set in CopyQualityValues
    FASTQSequence::CopyQualityValues(rhs);
}

void FASTQSequence::PrintFastq(ostream &out, int lineLength) const {
    PrintSeq(out, lineLength, '@');
    if (lineLength == 0) { 
        out << endl;
    }
    PrintFastqQuality(out, lineLength);
    if (lineLength == 0) {
        out << endl;
    }
}

void FASTQSequence::PrintFastqQuality(ostream &out, int lineLength) const {
    out << "+" << endl;
    PrintAsciiQual(out, lineLength);
}

bool FASTQSequence::GetQVs(const QVIndex & qvIndex, std::vector<uint8_t> & qvs, bool reverse) const {
    qvs.clear();
    uint8_t *  qualPtr;
    int charOffset = charToQuality;
    if (qvIndex == I_QualityValue) {
        qualPtr = qual.data;
    } else if (qvIndex == I_InsertionQV) {
        qualPtr = insertionQV.data;
    } else if (qvIndex == I_DeletionQV) {
        qualPtr = deletionQV.data;
    } else if (qvIndex == I_SubstitutionQV) {
        qualPtr = substitutionQV.data;
    } else if (qvIndex == I_MergeQV) {
        qualPtr = mergeQV.data;
    } else if (qvIndex == I_SubstitutionTag) {
        qualPtr = (uint8_t*)(substitutionTag);
        charOffset = 0;
    } else if (qvIndex == I_DeletionTag) {
        qualPtr = (uint8_t*)(deletionTag);
        charOffset = 0;
    }
    if (qualPtr == NULL) {
        return false;
    }

    qvs.resize(length);
    for (DNALength i = 0; i < length; i++) {
        if (not reverse) { // The same orientation
            qvs[i] = static_cast<uint8_t>(qualPtr[i] + charOffset);
        } else if (qvIndex != I_SubstitutionTag and qvIndex != I_DeletionTag) {
            // Reverse orientation, reverse QVs, except SubstitutionTag and DeletionTag
            qvs[i] = static_cast<uint8_t>(qualPtr[length - i - 1] + charOffset);
        } else { // Reverse and complement SubstitutionTag and DeletionTag
            qvs[i] = static_cast<uint8_t>(ReverseComplementNuc[qualPtr[length - i - 1] + charOffset]);
        }
        //assert(qvs[i] > 32 and qvs[i] < 127);
    }
    return true;
}

QVIndex FASTQSequence::GetQVIndex(const std::string & qvName) const {
    if (qvName == "QualityValue") {
        return I_QualityValue;
    } else if (qvName == "InsertionQV") {
        return I_InsertionQV;
    } else if (qvName == "DeletionQV") {
        return I_DeletionQV;
    } else if (qvName == "SubstitutionQV") {
        return I_SubstitutionQV;
    } else if (qvName == "MergeQV") {
        return I_MergeQV;
    } else if (qvName == "SubstitutionTag") {
        return  I_SubstitutionTag;
    } else if (qvName == "DeletionTag"){
        return I_DeletionTag;
    } else {
        std::cout << "ERROR: unknown Quality Value " << qvName << std::endl;
        assert(false);
    }
}

bool FASTQSequence::GetQVs(const std::string & qvName, std::vector<uint8_t> & qvs, bool reverse) const {
    return GetQVs(GetQVIndex(qvName), qvs, reverse);
}

bool FASTQSequence::GetQVs(const std::string & qvName, std::string & qvsStr, bool reverse) const {
    std::vector<uint8_t> qvs;
    bool OK = GetQVs(qvName, qvs, reverse);
    qvsStr = string(qvs.begin(), qvs.end());
    return OK;
}

void FASTQSequence::PrintAsciiRichQuality(ostream &out, 
        int whichQuality, int lineLength) const {
    vector<uint8_t> qvs;
    bool OK = GetQVs(static_cast<QVIndex>(whichQuality), qvs);
    
    DNALength i;
    if (lineLength == 0) {
        for (i = 0; i < length; i++) {
            if (OK) {
                out << static_cast<char>(qvs[i]);
            }
            else {
                // Fake bad quality
                out << "5";
            }
        }
    }
    else {
        for (i = 0; i < length; i++) {
            if (OK) {
                out << static_cast<char>(qvs[i]);
            }
            else {
                // Fake pretty bad quality.
                out << "5";
            }
            assert(lineLength != 0);
            if (i > 0 and (i+1) % lineLength==0) {
                out << endl;
            }
        }
        if (i == 0 or i % lineLength != 0) {
            out << endl;
        }
    }
}

void FASTQSequence::PrintAsciiQual(ostream &out, int lineLength) const {
    PrintAsciiRichQuality(out, 0, lineLength);
}

void FASTQSequence::PrintQual(ostream &out, int lineLength) const {
    out << ">" << this->title << endl;
    DNALength i;
    for (i = 0; i < length; i++ ){
        out << (int) qual[i];
        if (i > 0 and (i+1) % lineLength == 0)
            out << endl;
        else 
            out << " ";
    }
    if (i == 0 or i % lineLength != 0) {
        out << endl;
    }
}

void FASTQSequence::PrintQualSeq(ostream &out, int lineLength) const {
    FASTASequence::PrintSeq(out, lineLength);
    lineLength /= 4;
    PrintQual(out, lineLength);
}

// Create a reverse complement FASTQSequence of *this and assign to rhs.
void FASTQSequence::MakeRC(FASTQSequence &rc) {
    rc.Free();
    FASTASequence::MakeRC(rc);
    rc.SetQVScale(qvScale);

    if (not qual.Empty()) {
        // QVs are independent of one another. A FASTQSequence can have
        // insertionQV without having QualityValue.
        (static_cast<FASTQSequence*>(&rc))->AllocateQualitySpace(length);
        for (DNALength pos = 0; pos < length; pos++ ){
            rc.qual.data[length - pos - 1] = qual[pos];
        }
    }

    //
    // The read contains rich quality values. Reverse them here.
    //
    if (deletionQV.Empty() == false) {
        (static_cast<FASTQSequence*>(&rc))->AllocateDeletionQVSpace(length);
        for (DNALength pos = 0; pos < length; pos++) {
            rc.deletionQV[length - pos - 1] = deletionQV[pos];
        }
    }

    if (insertionQV.Empty() == false) {
        (static_cast<FASTQSequence*>(&rc))->AllocateInsertionQVSpace(length);
        for (DNALength pos = 0; pos < length; pos++) {
            rc.insertionQV[length - pos - 1] = insertionQV[pos];
        }
    }

    if (substitutionQV.Empty() == false) {
        (static_cast<FASTQSequence*>(&rc))->AllocateSubstitutionQVSpace(length);
        for (DNALength pos = 0; pos < length; pos++) {
            rc.substitutionQV[length - pos - 1] = substitutionQV[pos];
        }
    }

    if (mergeQV.Empty() == false) {
        (static_cast<FASTQSequence*>(&rc))->AllocateMergeQVSpace(length);
        for (DNALength pos = 0; pos < length; pos++) {
            rc.mergeQV[length - pos - 1] = mergeQV[pos];
        }
    }

    if (substitutionTag != NULL) {
        (static_cast<FASTQSequence*>(&rc))->AllocateSubstitutionTagSpace(length);
        for (DNALength pos = 0; pos < length; pos++) {
            rc.substitutionTag[length - pos - 1] = ReverseComplementNuc[substitutionTag[pos]];
        }
    }

    if (deletionTag != NULL) {
        (static_cast<FASTQSequence*>(&rc))->AllocateDeletionTagSpace(length);
        for (DNALength pos = 0; pos < length; pos++) {
            rc.deletionTag[length - pos - 1] = ReverseComplementNuc[deletionTag[pos]];
        }
    }

    if (preBaseDeletionQV.Empty() == false) {
        (static_cast<FASTQSequence*>(&rc))->AllocatePreBaseDeletionQVSpace(length);
        for (DNALength pos = 0; pos < length; pos++) {
            rc.preBaseDeletionQV[length - pos - 1] = preBaseDeletionQV[pos];
        }
    }

    deletionQVPrior = rc.deletionQVPrior;
    insertionQVPrior = rc.insertionQVPrior;
    substitutionQVPrior = rc.substitutionQVPrior;
    preBaseDeletionQVPrior = rc.preBaseDeletionQVPrior;
}

void FASTQSequence::Free() {
    if (deleteOnExit == true) { // Free Quality Values if under control
        qual.Free();
        deletionQV.Free();
        preBaseDeletionQV.Free();
        insertionQV.Free();
        substitutionQV.Free();
        mergeQV.Free();
        if (deletionTag != NULL) {
            delete[] deletionTag;
        }
        if (substitutionTag != NULL) {
            delete[] substitutionTag;
        }
    }
    //Reset deletionTag and substitionTag anyway
    deletionTag = NULL;
    substitutionTag = NULL;

    // Free seq and title, reset deleteOnExit. 
    // Don't call FASTASequence::Free before freeing QVs.
    FASTASequence::Free();
}

void FASTQSequence::LowerCaseMask(int qThreshold) {
    int i;
    if (qual.Empty() == true) return;

   for (i = 0; i < length; i++ ){
        if (qual[i] < qThreshold) {
            seq[i] = tolower(seq[i]);
        }
    }
}

float FASTQSequence::GetAverageQuality() const {
    DNALength p;
    float totalQ;
    if (qual.Empty() == true) { return 0.0; }
    assert(qual.Empty() == false);
    assert(length > 0);
    for (p = 0, totalQ = 0.0; p < length; p++) {
        totalQ += qual[p];
    }
    return totalQ / length;
}

#ifdef USE_PBBAM
void FASTQSequence::Copy(const PacBio::BAM::BamRecord & record) {
    FASTQSequence::Free();

    // Copy title and sequence.
    static_cast<FASTASequence*>(this)->Copy(record);

    // Copy QVs.
    qual.Copy(record.Qualities().Fastq());

    // iq
    if (record.HasInsertionQV()) {
        insertionQV.Copy(record.InsertionQV().Fastq());
    }
    // dq
    if (record.HasDeletionQV()) {
        deletionQV.Copy(record.DeletionQV().Fastq());
    }
    // sq
    if (record.HasSubstitutionQV()) {
        substitutionQV.Copy(record.SubstitutionQV().Fastq());
    }
    // mq
    if (record.HasMergeQV()) {
        mergeQV.Copy(record.MergeQV().Fastq());
    }
    // st
    if (record.HasSubstitutionTag()) {
        std::string qvs = record.SubstitutionTag();
        AllocateSubstitutionTagSpace(static_cast<DNALength>(qvs.size()));
        std::memcpy(substitutionTag, qvs.c_str(), qvs.size() * sizeof(char));
    }
    // dt
    if (record.HasDeletionTag()) {
        std::string qvs = record.DeletionTag();
        AllocateDeletionTagSpace(static_cast<DNALength>(qvs.size()));
        std::memcpy(deletionTag, qvs.c_str(), qvs.size() * sizeof(char));
    }
    // preBaseQVs are not included in BamRecord, and will not be copied.
    if (record.Type() != PacBio::BAM::RecordType::CCS) {
        subreadStart = static_cast<int>(record.QueryStart());
        subreadEnd = static_cast<int>(record.QueryEnd());
    } else {
        subreadStart = 0;
        subreadEnd =  static_cast<int>(record.Sequence().length());;
    }
}
#endif

#ifndef _BLASR_FASTQ_SEQUENCE_HPP_
#define _BLASR_FASTQ_SEQUENCE_HPP_

#include "DNASequence.hpp"
#include "FASTASequence.hpp"
#include "qvs/QualityValue.hpp"
#include "qvs/QualityValueVector.hpp"
#include "matrix/Matrix.hpp"
#include "reads/ZMWGroupEntry.hpp"

enum QVList {InsertionQV=0x1, DeletionQV=0x2, SubstitutionQV=0x4, MergeQV=0x8, SubstitutionTag=0x10, DeletionTag=0x20};
enum QVIndex {I_QualityValue=0, I_InsertionQV=1,I_DeletionQV=2,I_SubstitutionQV=3,I_MergeQV=4, I_SubstitutionTag=5,I_DeletionTag=6};

class FASTQSequence : public FASTASequence {
public:
    static int charToQuality;
    QualityValueVector<QualityValue> qual;
    QualityValueVector<QualityValue> deletionQV;
    QualityValueVector<QualityValue> preBaseDeletionQV;
    QualityValueVector<QualityValue> insertionQV;
    QualityValueVector<QualityValue> substitutionQV;
    QualityValueVector<QualityValue> mergeQV;
    Nucleotide *deletionTag;
    Nucleotide *substitutionTag;
    int subreadStart, subreadEnd;
    QualityValue deletionQVPrior, insertionQVPrior, substitutionQVPrior, preBaseDeletionQVPrior;

    QVScale qvScale;

    QVScale GetQVScale() const; 

    void SetQVScale(QVScale qvScaleP); 

    QualityValueVector<QualityValue>* GetQVPointerByIndex(int index);

    int GetStorageSize() const; 

    FASTQSequence();
    inline ~FASTQSequence();

    QualityValue GetDeletionQV(DNALength pos) const; 

    QualityValue GetMergeQV(DNALength pos) const; 

    Nucleotide GetSubstitutionTag(DNALength pos) const; 

    Nucleotide GetDeletionTag(DNALength pos) const; 

    QualityValue GetInsertionQV(DNALength pos) const; 

    QualityValue GetSubstitutionQV(DNALength pos) const; 

    QualityValue GetPreBaseDeletionQV(DNALength pos, Nucleotide nuc) const; 

    void ShallowCopy(const FASTQSequence &rhs); 

    void ReferenceSubstring(const FASTQSequence &rhs); 

    void ReferenceSubstring(const FASTQSequence &rhs, DNALength pos);

    void ReferenceSubstring(const FASTQSequence &rhs, DNALength pos, DNALength substrLength); 

    void ClearAndNull(QualityValue *value); 

    void CopyQualityValues(const FASTQSequence &rhs);

    void AllocateQualitySpace(DNALength qualLength); 

    void AllocateDeletionQVSpace(DNALength qualLength); 

    void AllocateMergeQVSpace(DNALength len);

    void AllocateDeletionTagSpace(DNALength qualLength); 

    void AllocatePreBaseDeletionQVSpace(DNALength qualLength); 

    void AllocateInsertionQVSpace(DNALength qualLength); 

    void AllocateSubstitutionQVSpace(DNALength qualLength );

    void AllocateSubstitutionTagSpace(DNALength qualLength );

    void AllocateRichQualityValues(DNALength qualLength); 

    void Copy(const FASTQSequence &rhs); 

    FASTQSequence& operator=(const FASTQSequence &rhs); 

    FASTQSequence(const FASTQSequence &rhs); 

    void Assign(FASTQSequence &rhs); 

    void PrintFastq(std::ostream &out, int lineLength=50) const; 

    void PrintFastqQuality(std::ostream &out, int lineLength=50) const; 

    QVIndex GetQVIndex(const std::string & qvName) const;

    /// Get QVs in vector<uint8_t> associated with the given QVIndex.
    /// \returns true if qvs are available, false otherwise
    /// \param [in] qvIndex - enum QVIndex
    /// \param [out] qvs - obtained QVs.
    /// \param [in] reverse - reverse orders of QVs or not
    bool GetQVs(const QVIndex & qvIndex, std::vector<uint8_t> & qvs, bool reverse=false) const;

    /// Get QVs in vector<uint8_t>, given with QV Name.
    /// \returns true if qvs are available, false, otherwise
    /// \param [in] qvName - InsertionQV, DeletionQV, SubstitionQV, MergeQV, SubstitutionTag, DeletionTag
    /// \param [out] qvs - obtians QVs.
    /// \param [in] reverse - reverse orders of QVs or not.
    bool GetQVs(const std::string & qvName, std::vector<uint8_t> & qvs, bool reverse=false) const;

    /// Get QVs in string, given with QV Name.
    /// \returns true if qvs are available, false, otherwise
    /// \param [in] qvName - InsertionQV, DeletionQV, SubstitionQV, MergeQV, SubstitutionTag, DeletionTag
    /// \param [out] qvs - obtians QVs.
    /// \param [in] reverse - reverse order of QVs or not
    bool GetQVs(const std::string & qvName, std::string & qvs, bool reverse=false) const;

    void PrintAsciiRichQuality(std::ostream &out, int whichQuality, int lineLength=50) const;

    void PrintAsciiQual(std::ostream &out, int lineLength=50) const;

    void PrintQual(std::ostream &out, int lineLength = 50) const;

    void PrintQualSeq(std::ostream &out, int lineLength = 50) const;

    void MakeRC(FASTQSequence &rc); 

    void Free(); 

    void LowerCaseMask(int qThreshold); 

    float GetAverageQuality() const; 

#ifdef USE_PBBAM
    /// Copy name, sequence, and QVs from BamRecord.
    void Copy(const PacBio::BAM::BamRecord & record);
#endif 
};

inline FASTQSequence::~FASTQSequence() {
    FASTQSequence::Free();
}

#endif // _BLASR_FASTQ_SEQUENCE_HPP_

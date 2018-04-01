// Author: Mark Chaisson

#ifndef _BLASR_SMRT_SEQUENCE_HPP_
#define _BLASR_SMRT_SEQUENCE_HPP_

#include <cassert>
#include <iostream>
#include <sstream>

#include "Enumerations.h"
#include "FASTQSequence.hpp"
#include "NucConversion.hpp"
#include "Types.h"
#include "reads/RegionTable.hpp"
#include "reads/ZMWGroupEntry.hpp"

class SMRTSequence : public FASTQSequence
{
    friend class HDFZMWReader;
    friend class HDFZMWWriter;
    friend class HDFZMWMetricsWriter;

private:
    enum SnrIndex4Base
    {
        A = 0,
        C = 1,
        G = 2,
        T = 3
    };
    float hqRegionSnr_[4];  // Always saved as 'ACGT'

    DNALength subreadStart_;
    DNALength subreadEnd_;

    // read group id associated with each SMRTSequence
    std::string readGroupId_;

public:
    ZMWGroupEntry zmwData;

    DNALength lowQualityPrefix, lowQualitySuffix;
    int highQualityRegionScore;  // High quality region score in region table.
    float readScore;

    // Whether or not this is originally copied from a BamRecord.
    bool copiedFromBam;

    HalfWord *preBaseFrames;
    HalfWord *widthInFrames;
    //
    // The following are fields that are read in from the pulse file.
    // Because they are not standard in bas.h5 files, these fields
    // should not be preallocated when resizing a SMRTSequence, and
    // memory should be managed separately.  For now, these fields all
    // have the same length as the number of bases, but this could
    // change so that all pulse values are stored in a SMRTSequence.
    //
    HalfWord *meanSignal, *maxSignal, *midSignal;
    float *classifierQV;
    unsigned int *startFrame;
    int *pulseIndex;

public:
    SMRTSequence();

    inline ~SMRTSequence();

    /// \name Sets and gets attributes.
    /// \{
    /// Set HoleNumber.
    /// \returns this SMRTSequence
    SMRTSequence &HoleNumber(UInt holeNumber);

    /// \reutrns HoleNumber
    UInt HoleNumber(void) const;

    /// Set HoleXY
    SMRTSequence &HoleXY(const int x, const int y);

    /// \returns HoleX
    UInt HoleX(void) const;

    /// \returns HoleY
    UInt HoleY(void) const;

    /// Set HoleStatus
    SMRTSequence &HoleStatus(const unsigned char);

    /// \returns HoleStatus
    unsigned char HoleStatus(void) const;

    /// \returns movie name parsed from sequence title
    std::string MovieName(void) const;

    /// \returns start pos of this sequence in coordinate of zmw polymerase sequence
    DNALength SubreadStart(void) const;

    /// Sets subreadStart.
    SMRTSequence &SubreadStart(const DNALength start);

    /// \returns subread end pos of this sequence in coordinate of zmw polymerase sequence
    DNALength SubreadEnd(void) const;

    /// Set subread end pos in coordinate of polymerase sequence.
    SMRTSequence &SubreadEnd(const DNALength end);

    /// A SMRTSequence's this->seq may point to sequence of a whole
    /// polymerase read, but only represents a subread [subreadStart_,  subreadEnd_).
    /// \returns subread length (SubreadEnd() - SubreadStart())
    DNALength SubreadLength(void) const;

    /// \returns read group id for this sequence.
    std::string ReadGroupId(void) const;

    /// Set readGroup Id for this sequence.
    SMRTSequence &ReadGroupId(const std::string &rid);

    /// Access to HQRegion SNRs must be done via public API.
    float HQRegionSnr(const char base) const;

    /// Set HQRegion SNR of base as v.
    SMRTSequence &HQRegionSnr(const char base, float v);

    /// \}

public:
    /// \name Clip subread
    /// \{
    SMRTSequence &Clip(const DNALength subreadStart, const DNALength subreadEnd);
    /// \}

    /// \name Allocate
    /// \{
    /// Allocate space for all possible QVs.
    void Allocate(DNALength length);

    /// Compact allocate space for QVs needed in order to compute alignment score.
    /// \param[in] hasInsertionDeletionQVTag: insertionQV, deletionQV and deletionTag exist
    /// \param[in] hasSubstitutionQVTag: substitutionQV and substitutionTag exist
    void CompactAllocate(const DNALength length, const bool hasInsertionDeletionQVTag,
                         const bool hasSubstitutionQVTag);
    /// \}

    /// Reconstruct this as a polymerase read of a zmw given subreads of the zmw.
    /// Copy only QVs necessary for computing alignment score, ignore irrelavent QVs.
    /// Length of this polymerase read is max {subread.SubreadEnd()}.
    /// Pad sequences and tags of scraps regions with 'N' and QVs of scraps regions with 0s.
    /// \param[in] subreads: subreads of a zmw
    void MadeFromSubreadsAsPolymerase(const std::vector<SMRTSequence> &subreads);

    void SetSubreadTitle(SMRTSequence &subread, DNALength subreadStart, DNALength subreadEnd);

    void SetSubreadBoundaries(SMRTSequence &subread, DNALength subreadStart, DNALength subreadEnd);

    void MakeSubreadAsMasked(SMRTSequence &subread, DNALength subreadStart = 0,
                             int subreadEnd = -1);

    void MakeSubreadAsReference(SMRTSequence &subread, DNALength subreadStart = 0,
                                int subreadEnd = -1);
#ifdef USE_PBBAM
    // Copy native orientated view of record to this->bamRecord
    void MakeNativeOrientedBamRecord(const PacBio::BAM::BamRecord &record);
#endif
    void Copy(const SMRTSequence &rhs);

    void Copy(const SMRTSequence &rhs, DNALength rhsPos, DNALength rhsLength);

    void Print(std::ostream &out) const;

    SMRTSequence &operator=(const SMRTSequence &rhs);

    void Free();

#ifdef USE_PBBAM
public:
    /// \returns if record is a valid bam record.
    /// FIXME: remove this function when pbbam.BamRecord.IsValid is available.
    static bool IsValid(const PacBio::BAM::BamRecord &record);

    // Copy read sequence, title, holeNumber, readGroupId, and QVs
    // (iq, dq, sq, mq, st, dt) from BamRecord to this SMRTSequence.
    // If copyAllQVs is false, also copy all QVs.
    void Copy(const PacBio::BAM::BamRecord &record, bool copyAllQVs = false);

    // Keep track of BamRecord from which this SMRTSequence is
    // originally copied. However, one should NOT assume
    // that this SMRTSequence has the same sequence, title, QVs as
    // the BamRecord, because this SMRTSequence may be created by
    // MakeSubreadAsMasked(...) or MakeRC(...).
    PacBio::BAM::BamRecord bamRecord;
#endif
};

inline SMRTSequence::~SMRTSequence() { SMRTSequence::Free(); }

#endif  // _BLASR_SMRT_SEQUENCE_HPP_

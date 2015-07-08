#ifndef _BLASR_BAM_PRINTER_IMPL_HPP_
#define _BLASR_BAM_PRINTER_IMPL_HPP_

#ifdef USE_PBBAM

#include <algorithm>
#include "utils/SMRTTitle.hpp"
using namespace BAMOutput;
using namespace std;
#include "pbbam/BamRecord.h"
#include "pbbam/BamFile.h"

template<typename T_Sequence>
void AlignmentToBamRecord(T_AlignmentCandidate & alignment, 
        T_Sequence & read, PacBio::BAM::BamRecord & bamRecord,
        AlignmentContext & context, SupplementalQVList & qvList,
        Clipping clipping, bool cigarUseSeqMatch) {
    assert(clipping == SAMOutput::soft or clipping == SAMOutput::subread);

    // Build from scratch if input reads are not from pbbam files.
    // Otherwise, use API provided by pbbam to make alignments from bamRecords.
 
    bool buildFromScratch = false;
    if (dynamic_cast<SMRTSequence*>(&read) == NULL) {
        //not SMRTSequence 
        buildFromScratch = true;
    } else { //is SMRTSequence, but not copied from bam.
        if (not read.copiedFromBam) {
            buildFromScratch = true;
        }
    }

    // build cigar string.
    string cigarString;
    T_Sequence alignedSequence;
    DNALength prefixSoftClip = 0, suffixSoftClip = 0;
    DNALength prefixHardClip = 0, suffixHardClip = 0;
    CreateCIGARString(alignment, read, cigarString, clipping,
                      prefixSoftClip, suffixSoftClip, 
                      prefixHardClip, suffixHardClip,
                      cigarUseSeqMatch);
    SetAlignedSequence(alignment, read, alignedSequence, clipping);
    PacBio::BAM::Cigar cigar = PacBio::BAM::Cigar::FromStdString(cigarString);
 
    // build flag
    uint16_t flag;
    BuildFlag(alignment, context, flag);

    // Get sequence string.
    string seqString;
    seqString.assign((char*)alignedSequence.seq, alignedSequence.length);

    // Get alignment starting position on reference sequence forward strand.
    PacBio::BAM::Position pos = 0;
    PacBio::BAM::Strand strand;
    if (alignment.tStrand == 0) {
        pos = static_cast<PacBio::BAM::Position>(alignment.TAlignStart());
        strand = PacBio::BAM::Strand::FORWARD;
    } else {
        pos = static_cast<PacBio::BAM::Position>(alignment.tLength - (alignment.TAlignStart() + alignment.TEnd()));
        strand = PacBio::BAM::Strand::REVERSE;
    }

    if (buildFromScratch) {
        SMRTTitle smrtTitle(alignment.qName);
        if (smrtTitle.isSMRTTitle) {
            bamRecord.Impl().Name(smrtTitle.ToString());
        } else {
            cout << "ERROR, can not convert non-pacbio reads to pbbam record." << endl;
            exit(-1);
        }
        bamRecord.Impl().SetSequenceAndQualities(seqString, alignedSequence.qual.ToString());
        bamRecord.Impl().CigarData(cigar);
        bamRecord.Impl().Bin(0);
        bamRecord.Impl().InsertSize(0);
        bamRecord.Impl().MapQuality(static_cast<uint8_t>(alignment.mapQV));
        bamRecord.Impl().MatePosition(static_cast<PacBio::BAM::Position>(-1));
        bamRecord.Impl().MateReferenceId(static_cast<int32_t>(-1));
        bamRecord.Impl().Position(pos);
        bamRecord.Impl().ReferenceId(static_cast<int32_t>(alignment.tIndex));

        if (strand == PacBio::BAM::Strand::REVERSE) {
            bamRecord.Impl().SetReverseStrand(true);
        }

        // Add tags required for bax->bam.
        PacBio::BAM::TagCollection tags;
        tags["RG"] = context.readGroupId;
        if (dynamic_cast<CCSSequence*>(&read) == NULL) { // subread
            tags["qs"] = read.subreadStart;
            tags["qe"] = read.subreadEnd;
            /// Number of passes for a subread should always be 1.
            tags["np"] = 1;
        } else { // ccs read
            /// Number of passes for ccs reads.
            tags["np"] = (static_cast<CCSSequence*>(&read))->numPasses;
        }
        tags["zm"] = read.zmwData.holeNumber;

        // Build QV tags.
        // Skip tags not define in BAM specification 3.0. 
        // including XL, XT, XQ, XS, XE, YS, YE.
        //
        // Write out optional quality values.  If qvlist does not
        // have any qv's signaled to print, this is a no-op.
        //
        // First transform characters that are too large to printable ones.
        qvList.FormatQVOptionalFields(alignedSequence);

        // Add QVs to BamRecordImpl.
        string insertionQVs, deletionQVs, substitutionQVs, mergeQVs, substitutionTags, deletionTags;
        bool alnReverse = (alignment.tStrand == 1); // reverse-complement alignment
        // If this is a reverse-complement alignment, bases and QVs of 
        // alignedSequence are reverse(-complement) of the sequence read
        // from bax.h5 file. With PB BAM specification 3.0, QVs will be stored 
        // in BAM Record, therefore, reverse orders of QVs before add QV tags. 
        if (alignedSequence.GetQVs("InsertionQV", insertionQVs, alnReverse)) { 
            tags["iq"] = insertionQVs;
        }
        if (alignedSequence.GetQVs("DeletionQV", deletionQVs, alnReverse)) {
            tags["dq"] = deletionQVs;
        } 
        if (alignedSequence.GetQVs("SubstitutionQV", substitutionQVs, alnReverse)) {
            tags["sq"] = substitutionQVs;
        } 
        if (alignedSequence.GetQVs("MergeQV", mergeQVs, alnReverse)) {
            tags["mq"] = mergeQVs;
        }
        // substitutionTag is not included by default
        if (alignedSequence.GetQVs("DeletionTag", deletionTags, alnReverse)) {
            tags["dt"] = deletionTags;
        }
        bamRecord.Impl().Tags(tags);
    } else {
        // The following code can be used to hard-clip reads, if needed. 
        // PacBio::BAM::Position clipStart = read.bamRecord.QueryStart() + alignment.QAlignStart();
        // PacBio::BAM::Position clipEnd = read.bamRecord.QueryStart() + alignment.QAlignEnd();
        // bamRecord = PacBio::BAM::BamRecord::Clipped(read.bamRecord, 
        //                PacBio::BAM::ClipType::CLIP_TO_QUERY,
        //                clipStart, clipEnd).
        bamRecord = PacBio::BAM::BamRecord::Mapped(read.bamRecord,
                        static_cast<int32_t>(alignment.tIndex),
                        static_cast<PacBio::BAM::Position>(pos),
                        strand, cigar, 
                        static_cast<uint8_t>(alignment.mapQV));
    }

    // Add tags common for bax->bam and bam->bam.
    bamRecord.Impl().AddTag("AS", alignment.score);
    bamRecord.Impl().AddTag("NM", context.editDist);

    // Set Flag 
    bamRecord.Impl().Flag(static_cast<uint32_t>(flag));

}

template<typename T_Sequence>
void BAMOutput::PrintAlignment(T_AlignmentCandidate &alignment, T_Sequence &read,
        PacBio::BAM::BamWriter &bamWriter, AlignmentContext &context, 
        SupplementalQVList & qvList, Clipping clipping, bool cigarUseSeqMatch) {

    PacBio::BAM::BamRecord bamRecord;
    AlignmentToBamRecord(alignment, read, bamRecord, context, qvList, clipping, cigarUseSeqMatch);
    bamWriter.Write(bamRecord);
}
#endif

#endif

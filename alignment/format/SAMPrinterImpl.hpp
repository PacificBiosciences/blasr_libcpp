#include <algorithm> //max
#include <utility> //swap
#include <assert.h> //assert

using namespace SAMOutput; 


template<typename T_Sequence>
void SAMOutput::SetAlignedSequence(T_AlignmentCandidate &alignment, T_Sequence &read,
        T_Sequence &alignedSeq,
        Clipping clipping) {
    //
    // In both no, and hard clipping, the dna sequence that is output
    // solely corresponds to the aligned sequence.
    //
    DNALength clippedReadLength = 0;
    DNALength clippedStartPos   = 0;

    if (clipping == none or clipping == hard) {
        DNALength qStart = alignment.QAlignStart();
        DNALength qEnd   = alignment.QAlignEnd();
        clippedReadLength = qEnd - qStart;
        clippedStartPos   = qStart;
    }
    else if (clipping == soft) {
        clippedReadLength = read.length - read.lowQualityPrefix - read.lowQualitySuffix;
        clippedStartPos = read.lowQualityPrefix;
    }
    else if (clipping == subread) {
        clippedReadLength = read.SubreadLength();
        clippedStartPos = read.SubreadStart();
    }
    else {
        std::cout <<" ERROR! The clipping must be none, hard, subread, or soft when setting the aligned sequence." << std::endl;
        assert(0);
    }

    //
    // Set the aligned sequence according to the clipping boundaries.
    //
    if (alignment.tStrand == 0) {
        alignedSeq.ReferenceSubstring(read, clippedStartPos, clippedReadLength);
    }
    else {
        T_Sequence subSeq;
        subSeq.ReferenceSubstring(read, clippedStartPos, clippedReadLength);
        subSeq.MakeRC(alignedSeq);
        assert(alignedSeq.deleteOnExit);
    }
}

template<typename T_Sequence>
void SAMOutput::SetSoftClip(T_AlignmentCandidate &alignment,
        T_Sequence &read,
        DNALength hardClipPrefix,
        DNALength hardClipSuffix,
        DNALength &softClipPrefix, 
        DNALength &softClipSuffix) {
    DNALength qStart, qEnd;
    qStart = alignment.QAlignStart();
    qEnd   = alignment.QAlignEnd();

    assert(qStart >= hardClipPrefix);
    softClipPrefix = qStart - hardClipPrefix;
    assert(read.length >= hardClipSuffix and read.length >= qEnd and
           qEnd + hardClipSuffix <= read.length);
    softClipSuffix = read.length - hardClipSuffix - qEnd;
}

template<typename T_Sequence>
void SAMOutput::SetHardClip(T_AlignmentCandidate &alignment, 
        T_Sequence &read,
        DNALength &prefixClip,
        DNALength &suffixClip) {
    //
    // Set the hard clipping assuming the read is in the forward
    // direction.
    //
    prefixClip = alignment.QAlignStart();
    suffixClip = read.length - alignment.QAlignEnd();
    if (alignment.tStrand == 1) {
        //
        // If the read is instead reverse, swap the clipping lengths.
        //
        std::swap(prefixClip, suffixClip);
    }
}


//
// Straight forward: create the cigar string allowing some clipping
// The read is provided to give length and hq information.
//
template<typename T_Sequence>
void SAMOutput::CreateCIGARString(T_AlignmentCandidate &alignment,
        T_Sequence &read,
        std::string &cigarString,
        Clipping clipping,
        DNALength & prefixSoftClip, DNALength & suffixSoftClip, 
        DNALength & prefixHardClip, DNALength & suffixHardClip,
        bool cigarUseSeqMatch, const bool allowAdjacentIndels) {

    cigarString = "";
    // All cigarString use the no clipping core
    std::vector<int> opSize;
    std::vector<char> opChar;
    CreateNoClippingCigarOps(alignment, opSize, opChar, cigarUseSeqMatch, allowAdjacentIndels);

    // Clipping needs to be added

    if (clipping == hard) {
      SetHardClip(alignment, read, prefixHardClip, suffixHardClip);
      if (prefixHardClip > 0) {
        opSize.insert(opSize.begin(), prefixHardClip);
        opChar.insert(opChar.begin(), 'H');
      }
      if (suffixHardClip > 0) {
        opSize.push_back(suffixHardClip);
        opChar.push_back('H');
      }
      prefixSoftClip = 0;
      suffixSoftClip = 0;
    }
    if (clipping == soft or clipping == subread) {
      //
      // Even if clipping is soft, the hard clipping removes the 
      // low quality regions
      //
      if (clipping == soft) {
          prefixHardClip = read.lowQualityPrefix;
          suffixHardClip = read.lowQualitySuffix;
      }
      else if (clipping == subread) {
          prefixHardClip = std::max((DNALength) read.SubreadStart(), read.lowQualityPrefix);
          suffixHardClip = std::max((DNALength)(read.length - read.SubreadEnd()), read.lowQualitySuffix);
      }

      SetSoftClip(alignment, read, prefixHardClip, suffixHardClip, prefixSoftClip, suffixSoftClip);

      if (alignment.tStrand == 1) {
        std::swap(prefixHardClip, suffixHardClip);
        std::swap(prefixSoftClip, suffixSoftClip);
      }

      //
      // Insert the hard and soft clipping so that they are in the
      // order H then S if both exist.
      //
      if (prefixSoftClip > 0) {
        opSize.insert(opSize.begin(), prefixSoftClip);
        opChar.insert(opChar.begin(), 'S');
      }
      if (prefixHardClip > 0) {
        opSize.insert(opSize.begin(), prefixHardClip);
        opChar.insert(opChar.begin(), 'H');
      }
      
      //
      // Append the hard and soft clipping so they are in the order S
      // then H. 
      //
      if (suffixSoftClip > 0) {
        opSize.push_back(suffixSoftClip);
        opChar.push_back('S');
      }
      if (suffixHardClip > 0) {
        opSize.push_back(suffixHardClip);
        opChar.push_back('H');
      }
    }

    CigarOpsToString(opSize, opChar, cigarString);
}

template<typename T_Sequence>
void SAMOutput::PrintAlignment(T_AlignmentCandidate &alignment,
        T_Sequence &read,
        std::ostream &samFile,
        AlignmentContext &context,
        SupplementalQVList & qvList,
        Clipping clipping,
        bool cigarUseSeqMatch,
        const bool allowAdjacentIndels) {

    std::string cigarString;
    uint16_t flag;
    T_Sequence alignedSequence;
    DNALength prefixSoftClip = 0, suffixSoftClip = 0;
    DNALength prefixHardClip = 0, suffixHardClip = 0;

    CreateCIGARString(alignment, read, cigarString, clipping, prefixSoftClip, suffixSoftClip, prefixHardClip, suffixHardClip, cigarUseSeqMatch, allowAdjacentIndels);
    SetAlignedSequence(alignment, read, alignedSequence, clipping);
    BuildFlag(alignment, context, flag);
    samFile << alignment.qName << "\t" 
            << flag << "\t" 
            << alignment.tName << "\t";   // RNAME
    if (alignment.tStrand == 0) {
      samFile << alignment.TAlignStart() + 1 << "\t"; 
      // POS, add 1 to get 1 based coordinate system
    }
    else {
      samFile << alignment.tLength - (alignment.TAlignStart() + alignment.TEnd()) + 1 << "\t"; // includes - 1 for rev-comp,  +1 for one-based
    }
    samFile << (int) alignment.mapQV << "\t"// MAPQ
            << cigarString << "\t"; // CIGAR
      
      //
      // Determine RNEXT

    std::string rNext;
    rNext = "*";
    /*
    if (context.hasNextSubreadPos == false) {
      rNext = "*";
    }
    else {
      if (context.rNext == alignment.tName) {
        rNext = "=";
      }
      else {
        rNext = context.rNext;
      }
    }
    */
    samFile << rNext << "\t"; // RNEXT
    
    DNALength nextSubreadPos = 0;
    /*
    if (context.hasNextSubreadPos) {
      nextSubreadPos = context.nextSubreadPos + 1;
      }*/
    samFile << nextSubreadPos << "\t"; // RNEXT, add 1 for 1 based
                                           // indexing

    //DNALength tLen = alignment.GenomicTEnd() - alignment.GenomicTBegin();
    //SAM v1.5, tLen is set as 0 for single-segment template
    samFile << 0 << "\t"; // TLEN
    // Print the sequence on one line, and suppress printing the
    // newline (by setting the line length to alignedSequence.length
    (static_cast<DNASequence*>(&alignedSequence))->PrintSeq(samFile, 0);  // SEQ
    samFile << "\t";
    if (alignedSequence.qual.data != NULL && qvList.useqv == 0) {
        alignedSequence.PrintAsciiQual(samFile, 0);  // QUAL
    }
    else {
      samFile <<"*";
    }
    samFile << "\t";
    //
    // Add optional fields
    //
    samFile << "RG:Z:" << context.readGroupId << "\t";
    samFile << "AS:i:" << alignment.score << "\t";

    //
    // "RG" read group Id
    // "AS" alignment score
    // "XS" read alignment start position without counting previous soft clips (1 based) 
    // "XE" read alignment end position without counting previous soft clips (1 based) 
    // "XL" aligned read length 
    // "XQ" query sequence length
    // "XT" # of continues reads, always 1 for blasr 
    // "NM" edit distance 
    // "FI" read alignment start position (1 based) 
    //
    
    DNALength qAlignStart = alignment.QAlignStart();
    DNALength qAlignEnd = alignment.QAlignEnd();

    if (clipping == none) {
      samFile << "XS:i:" << qAlignStart + 1 << "\t";
      samFile << "XE:i:" << qAlignEnd + 1 << "\t";
    }
    else if (clipping == hard or clipping == soft or clipping == subread) {
        DNALength xs = prefixHardClip;
        DNALength xe = read.length - suffixHardClip;
        if (alignment.tStrand == 1) {
            xs = suffixHardClip;
            xe = read.length - prefixHardClip;
        }
        samFile << "XS:i:" << xs + 1 << "\t"; // add 1 for 1-based indexing in sam
        assert(read.length - suffixHardClip == prefixHardClip + alignedSequence.length);
        samFile << "XE:i:" << xe + 1 << "\t";
    }
    samFile << "YS:i:" << read.SubreadStart() << "\t";
    samFile << "YE:i:" << read.SubreadEnd() << "\t";
    samFile << "ZM:i:" << read.HoleNumber() << "\t";
    samFile << "XL:i:" << alignment.qAlignedSeq.length << "\t";
    samFile << "XT:i:1\t"; // reads are allways continuous reads, not
                        // referenced based circular consensus when
                        // output by blasr.
    samFile << "NM:i:" << context.editDist << "\t";
    samFile << "FI:i:" << alignment.qAlignedSeqPos + 1;
    // Add query sequence length
    samFile << "\t" << "XQ:i:" << alignment.qLength;

    //
	// Write out optional quality values.  If qvlist does not 
	// have any qv's signaled to print, this is a no-op.
	//
	// First transform characters that are too large to printable ones.
	qvList.FormatQVOptionalFields(alignedSequence);
	qvList.PrintQVOptionalFields(alignedSequence, samFile);

    samFile << std::endl;
}

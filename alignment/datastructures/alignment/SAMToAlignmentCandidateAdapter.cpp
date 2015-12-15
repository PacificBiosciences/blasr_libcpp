#include <algorithm>
#include "utils/SMRTTitle.hpp"
#include "SAMToAlignmentCandidateAdapter.hpp"

void InitializeCandidateFromSAM(SAMAlignment &sam,
                                AlignmentCandidate<> &candidate) {
  candidate.qName = sam.qName;
  candidate.tName = sam.rName;
  candidate.tPos  = sam.pos;
  candidate.mapQV = sam.mapQV;
}

size_t ProcessGap(std::vector<int> &lengths, std::vector<char> &ops, 
               size_t &opIndex, size_t opEnd, blasr::GapList &gaps, 
               int &qAdvance, int &tAdvance) {
  //
  // Default is no gap.
  //
  qAdvance = tAdvance = 0;
  gaps.clear();
  if (opIndex >= lengths.size()) {
    return opIndex;
  }

  //
  // Look to see if the suffix of operations stars with a gap.
  //
  if (ops[opIndex] != 'I' and
      ops[opIndex] != 'D') {
    return opIndex;
  }

  while (opIndex < opEnd and (ops[opIndex] == 'I' or ops[opIndex] == 'D')) {
    if (ops[opIndex] == 'I') {
      gaps.push_back(blasr::Gap(blasr::Gap::Target, lengths[opIndex]));
      qAdvance += lengths[opIndex];
    }
    else {
      assert(ops[opIndex] == 'D');
      gaps.push_back(blasr::Gap(blasr::Gap::Query, lengths[opIndex]));
      tAdvance += lengths[opIndex];
    }
    opIndex++;
  }
  return opIndex;
}

int AdvancePastClipping(std::vector<int> &lengths,
                        std::vector<char> &ops,
                        size_t &opIndex,
                        int &numSoftClipped) {
  int numClipped = 0;
  numSoftClipped = 0;
  while (opIndex < lengths.size() and IsClipping(ops[opIndex])) {
    numClipped += lengths[opIndex];
    if (ops[opIndex] == 'S') {
      numSoftClipped+= lengths[opIndex];
    }

    ++opIndex;
  }
  return numClipped;
}

int AdvancePastSkipped(std::vector<int> &lengths,
                       std::vector<char> &ops,
                       size_t &opIndex) {
  int numSkipped = 0;
  while (opIndex < lengths.size() and IsSkipped(ops[opIndex])) {
    numSkipped += lengths[opIndex];
    opIndex++;
  }
  return numSkipped;
}

int ProcessMatch(std::vector<int> &lengths, std::vector<char> &ops,
                 size_t &opIndex, size_t opEnd) {

  //
  // Make sure this starts on some sort of blockn
  //

  int blockLength = 0;
  while (opIndex < opEnd and (IsAssignChar(ops[opIndex]))) {
    blockLength += lengths[opIndex];
    opIndex++;
  }
  return blockLength;
}

void CIGAROpsToBlocks(std::vector<int> &lengths, std::vector<char> &ops,
                      size_t cigarStart,
                      size_t cigarEnd,
                      AlignmentCandidate<> &aln) {
  cigarStart = 0;
  cigarEnd   = lengths.size();
  CIGAROpsToBlocks(lengths, ops, cigarStart, cigarEnd, aln);
}

int AdvancePosToAlignmentEnd(std::vector<char> &ops, size_t &pos) {
  size_t start = pos;
  while (pos < ops.size() and ops[pos] != 'N' and !IsClipping(ops[pos])) {
    pos++;
  }
  return pos - start;
}

int GetAlignedQueryLengthByCIGARSum(std::vector<char> &ops, std::vector<int> &lengths) {
  size_t i;
  for (i = 0; i < ops.size(); i++) {
    if (ops[i] != 'S' && ops[i] != 'H') {
      break;
    }
  }
  int queryLength = 0;
  for (; i < ops.size() && ops[i] != 'S' && ops[i] != 'H'; i++) {
    if (IsAssignChar(ops[i]) or ops[i] == 'I' or ops[i] == 'N') {
      queryLength += lengths[i];
    }
  }
  return queryLength;
}

int GetAlignedReferenceLengthByCIGARSum(std::vector<char> &ops, std::vector<int> &lengths) {
  size_t i;
  for (i = 0; i < ops.size(); i++) {
    if (ops[i] != 'S' && ops[i] != 'H') {
      break;
    }
  }
  int refLength = 0;
  for (; i < ops.size() && ops[i] != 'S' && ops[i] != 'H'; i++) {
    if (IsAssignChar(ops[i]) or ops[i] == 'D' or ops[i] == 'N') {
      refLength += lengths[i];
    }
  }
  return refLength;
}

void CIGAROpsToBlocks(std::vector<int> &lengths, std::vector<char> &ops,
                      size_t &cigarPos,
                      size_t &cigarEnd,
                      int &qPos, int &tPos,
                      AlignmentCandidate<> &aln) {

  DNALength qStart = qPos, tStart = tPos;
  assert(cigarPos >= cigarEnd or !IsClipping(ops[cigarPos]));

  //
  // Advance past any skipped portion.
  //
  int numSkipped = AdvancePastSkipped(lengths, ops, cigarPos);
  tPos += numSkipped;

  //
  // Process the gaps before the first match.
  //

  //
  // If there is nothing left, just bail.
  blasr::GapList gap;

  cigarEnd = cigarPos;
  AdvancePosToAlignmentEnd(ops, cigarEnd);
  if (cigarPos >= cigarEnd) {
    return;
  }
  
  // 
  // Process any gap that the aligner produces before the first match
  // begins.
  //
  int qAdvance, tAdvance;
  ProcessGap(lengths, ops, cigarPos, cigarEnd, gap, qAdvance, tAdvance);
  aln.gaps.push_back(gap);
  qPos += qAdvance;
  tPos += tAdvance;
  //
  // Now add gaps.
  //
  while(cigarPos < cigarEnd) {
    //
    // The next operation must be a match.
    //
    int matchLength = ProcessMatch(lengths, ops, cigarPos, cigarEnd);
    blasr::Block b;
    b.qPos = qPos - qStart;
    b.tPos = tPos - tStart;
    b.length = matchLength;
    aln.blocks.push_back(b);
    qPos += b.length;
    tPos += b.length;

    ProcessGap(lengths, ops, cigarPos, cigarEnd, gap, qAdvance, tAdvance);
    aln.gaps.push_back(gap);
    tPos += tAdvance;
    qPos += qAdvance;
  }
  
}

void SAMAlignmentsToCandidates(SAMAlignment &sam,
                               std::vector<FASTASequence> &referenceSequences,
                               std::map<std::string,int> & refNameToRefListIndex,
                               std::vector<AlignmentCandidate<> > &candidates, 
                               bool parseSmrtTitle,
                               bool keepRefAsForward,
                               bool copyQVs) {
  //
  // First determine how many alignments there are from CIGAR string.
  //
  std::vector<int> lengths;
  std::vector<char> ops;
  sam.cigar.Vectorize(lengths, ops);

  DNASequence querySeq;
  // For now just reference the query sequence.
  querySeq.deleteOnExit = false;
  querySeq.seq = (Nucleotide*) sam.seq.c_str();
  querySeq.length = sam.seq.size();

  DNALength samTEnd = 0;
  DNALength samTStart = sam.pos - 1;

  std::vector<std::string> optionalQVs;
  if (copyQVs) {
      sam.CopyQVs(&optionalQVs);    
  }
  if (keepRefAsForward == false and IsReverseComplement(sam.flag)) {
    ReverseAlignmentOperations(lengths, ops);
    DNASequence rcQuerySeq;
    querySeq.CopyAsRC(rcQuerySeq);
    //
    // Zero out the query seq so that the string memory is not
    // deleted.
    //
    querySeq.seq = NULL;
    querySeq.length = 0;
    querySeq = rcQuerySeq;
    rcQuerySeq.Free();
    samTEnd = GetAlignedReferenceLengthByCIGARSum(ops, lengths);
    
    // We also need to reverse any optional QVs
    if (copyQVs) {
      for(size_t i=0; i<optionalQVs.size(); i++) {
        std::reverse(optionalQVs[i].begin(), optionalQVs[i].end());
      }
    }
  }


  if (ops.size() == 0) {
    return;
  }
  
  //
  // Advance past any clipping.  This advances in both query and
  // reference position.
  //
  size_t cigarPos = 0;
  int qPos = 0; 
  int tPos = 0;

  DNALength queryPosOffset = 0;
  if (parseSmrtTitle) {
    //
    // The aligned sequence is really a subread of a full
    // sequence. The position of the aligments start at 0, the
    // beginning of the query sequence, but in the sam file, they
    // may appear as subreads, and are offset from the start of the
    // subread.  By convention, the subread coordinates are embedded
    // in the title of the query, if it is a smrtTitle. 
    // Two types of smrtTitle are supported:
    // movie/zmw/start_end
    // movie/zmw/start_end/start2_end2
    SMRTTitle stitle = SMRTTitle(sam.qName);

    if (not stitle.isSMRTTitle) {
      std::cout << "ERROR. Could not parse title " << sam.qName << std::endl;
      exit(1);
    }
    queryPosOffset = stitle.start;
  }
  else if (sam.xs) {
    queryPosOffset += sam.xs - 1;
  }


  while (cigarPos < lengths.size()) {
    //
    // Sequence clipping becomes offsets into the q/t alignedSeqPos
    //


    int numSoftClipped;
    AdvancePastClipping(lengths, ops, cigarPos, numSoftClipped);

    //
    // End loop now.
    //
    if (cigarPos >= lengths.size()) {
      break;
    }
    qPos += numSoftClipped;

    //
    // Skipped sequences are just advances in the tPos.
    //
    int numSkipped = AdvancePastSkipped(lengths, ops, cigarPos);
    tPos += numSkipped;

    if (cigarPos >= lengths.size()) {
      break;
    }


    AlignmentCandidate<> alignment;
    //
    // The aligned sequence must start at a match therefore the tpos
    // and qpos are 0.
    //
    alignment.qPos = 0;
    alignment.tPos = 0;

    // qAlignStart is the start of the alignment relative to the sequence in the SAM file.
    DNALength qAlignStart = qPos;
    // tAlignStart is the start of the alignment in the genome.
    DNALength tAlignStart = tPos;
    
    size_t cigarEnd = cigarPos;
    AdvancePosToAlignmentEnd(ops, cigarEnd);

    CIGAROpsToBlocks(lengths, ops,          
                     cigarPos, cigarEnd,
                     qPos, tPos,
                     alignment);

    GetAlignedQueryLengthByCIGARSum(ops, lengths);
    GetAlignedReferenceLengthByCIGARSum(ops, lengths);

    alignment.qAlignedSeqLength = qPos - qAlignStart;
    alignment.tAlignedSeqLength = tPos - tAlignStart;

    //
    // Assign candidate sequences.
    //
    // First, the query sequence is straight from the SAM line.
    ((DNASequence*)&alignment.qAlignedSeq)->Copy(querySeq, qAlignStart, alignment.qAlignedSeqLength);
    if (copyQVs) {
      alignment.ReadOptionalQVs(optionalQVs, qAlignStart, alignment.qAlignedSeqLength);
    }
    
    // The SAM Alignments a
    alignment.qStrand = IsReverseComplement(sam.flag);
    alignment.tStrand = 0;
    alignment.mapQV   = sam.mapQV;

    //
    // Assign the offsets into the original sequence where the
    // subsequence starts.
    //

    alignment.qAlignedSeqPos = queryPosOffset + qAlignStart;    
    alignment.tAlignedSeqPos = samTStart + tAlignStart;
    
    if (sam.rName == "*") {
      //
      // No reference, do not add the alignment to the list of
      // candidates.
      //
      continue;
    }
    else {
      int refIndex;
      if (refNameToRefListIndex.find(sam.rName) == refNameToRefListIndex.end()) {
        std::cout <<" ERROR.  SAM Reference " << sam.rName << " is not found in the list of reference contigs." << std::endl;
        exit(1);
      }
      
      refIndex = refNameToRefListIndex[sam.rName];
     
      alignment.tLength = referenceSequences[refIndex].length;
      alignment.qLength = sam.seq.size(); 
      alignment.qName = sam.qName;
      alignment.tName = sam.rName;


      if (keepRefAsForward == false and alignment.qStrand == 1) {

        //
        // Now that the reference sequence has been copied, if it is
        // on the reverse strand, make the reverse complement for
        // proper printing.
        //
        alignment.tAlignedSeqPos = samTStart + (samTEnd - tAlignStart - alignment.tAlignedSeqLength);
		if (alignment.tAlignedSeqLength > referenceSequences[refIndex].length ||
			alignment.tAlignedSeqPos    > referenceSequences[refIndex].length ||
			alignment.tAlignedSeqLength + alignment.tAlignedSeqPos > referenceSequences[refIndex].length + 2) {
            //alignment.tAlignedSeqPos is 1 based and unsigned.
			std::cout << "WARNING. The mapping of read " << alignment.qName  
				 << " to reference "      << alignment.tName 
                 << " is out of bounds."  << std::endl
                 << "         StartPos (" << alignment.tAlignedSeqPos  
                 << ") + AlnLength (" << alignment.tAlignedSeqLength 
                 << ") > RefLength (" << referenceSequences[refIndex].length
                 << ") + 2 "          << std::endl;
            continue;
		}
        ((DNASequence*)&alignment.tAlignedSeq)->Copy(referenceSequences[refIndex], alignment.tAlignedSeqPos, alignment.tAlignedSeqLength);             
        alignment.tAlignedSeq.ReverseComplementSelf();
        // either ref or read is defined as being in the forward
        // orientation.  Here, since refAsForward is false, the read
        // is forward.  Since the read is forward, the aligned
        // sequences are stored as the reverse complement of the read
        // and the references.
        //
        alignment.tStrand = 1;
        alignment.qStrand = 0;
      }
      else {
        if (alignment.tAlignedSeqLength > referenceSequences[refIndex].length ||
			alignment.tAlignedSeqPos    > referenceSequences[refIndex].length ||
			alignment.tAlignedSeqLength + alignment.tAlignedSeqPos > referenceSequences[refIndex].length + 2) {
            //alignment.tAlignedSeqPos is 1 based and unsigned. 
			std::cout << "WARNING. The mapping of read " << alignment.qName  
				 << " to reference "      << alignment.tName 
                 << " is out of bounds."  << std::endl
                 << "         StartPos (" << alignment.tAlignedSeqPos  
                 << ") + AlnLength (" << alignment.tAlignedSeqLength 
                 << ") > RefLength (" << referenceSequences[refIndex].length
                 << ") + 2 "          << std::endl;
            continue;
		}
        ((DNASequence*)&alignment.tAlignedSeq)->Copy(referenceSequences[refIndex], 
                                                     alignment.tAlignedSeqPos, 
                                                     alignment.tAlignedSeqLength);
      }
    }

    if (alignment.blocks.size() > 0) {
      candidates.push_back(alignment);
    }
  }
  if (candidates.size() > 0 and keepRefAsForward == false and candidates[0].tStrand == 1) {
    std::reverse(candidates.begin(), candidates.end());
  }
  querySeq.Free();
}

void ReverseAlignmentOperations(std::vector<int> &lengths, std::vector<char> &ops) {
  reverse(lengths.begin(), lengths.end());
  reverse(ops.begin(), ops.end());
}

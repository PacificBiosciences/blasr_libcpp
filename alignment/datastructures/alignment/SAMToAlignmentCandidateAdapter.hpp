#ifndef _BLASR_SAM_TO_ALIGNMENT_CANDIDATE_ADAPTER_HPP_
#define _BLASR_SAM_TO_ALIGNMENT_CANDIDATE_ADAPTER_HPP_

#include <map>

#include "sam/SAMAlignment.hpp"
#include "sam/SAMFlag.h"
#include "datastructures/alignment/AlignmentCandidate.hpp"

void InitializeCandidateFromSAM(SAMAlignment &sam,
                                AlignmentCandidate<> &candidate);

//
// Change a list of gap operations into a gap structure between
// blocks. 
//
size_t ProcessGap(std::vector<int> &lengths, std::vector<char> &ops, 
               size_t &opIndex, size_t opEnd, blasr::GapList &gaps, 
               int &qAdvance, int &tAdvance);

inline bool IsClipping(char c) {
  return (c == 'S' or c == 'H');
}

inline int IsSkipped(char c) {
  return c == 'N';
}
//
//
// Returns true if the corresponding characters are assigned to each
// other in the alignment.  This is either a match or a mismatch.
//

inline bool IsAssignChar(char c) {
  return (c == 'M' or c == 'X' or c == '=');
}

int AdvancePastClipping(std::vector<int> &lengths,
                        std::vector<char> &ops,
                        size_t &opIndex,
                        int &numSoftClipped);

int AdvancePastSkipped(std::vector<int> &lengths,
                       std::vector<char> &ops,
                       size_t &opIndex);



int ProcessMatch(std::vector<int> &lengths, std::vector<char> &ops,
                 int &opIndex, int opEnd);
                 

void CIGAROpsToBlocks(std::vector<int> &lengths, std::vector<char> &ops,
                      size_t cigarStart,
                      size_t cigarEnd,
                      AlignmentCandidate<> &aln);
  
int AdvancePosToAlignmentEnd(std::vector<char> &ops, int &pos);

int GetAlignedQueryLengthByCIGARSum(std::vector<char> &ops, std::vector<int> &lengths);

int GetAlignedReferenceLengthByCIGARSum(std::vector<char> &ops, std::vector<int> &lengths);

void CIGAROpsToBlocks(std::vector<int> &lengths, std::vector<char> &ops,
                      size_t &cigarPos,
                      size_t &cigarEnd,
                      int &qPos, int &tPos,
                      AlignmentCandidate<> &aln);

void ReverseAlignmentOperations(std::vector<int> &lengths, std::vector<char> &ops);

void SAMAlignmentsToCandidates(SAMAlignment &sam,
                               std::vector<FASTASequence> &referenceSequences,
                               std::map<std::string,int> & refNameToRefListIndex,
                               std::vector<AlignmentCandidate<> > &candidates, 
                               bool parseSmrtTitle = false,
                               bool keepRefAsForward = true,
                               bool copyQVs = false);
#endif

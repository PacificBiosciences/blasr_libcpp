#ifndef _BLASR_FIND_MAX_INTERVAL_HPP_
#define _BLASR_FIND_MAX_INTERVAL_HPP_

#include <semaphore.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include "LongestIncreasingSubsequence.hpp"
#include "GlobalChain.hpp"
#include "BasicEndpoint.hpp"
#include "../../datastructures/anchoring/WeightedInterval.hpp"
#include "../../datastructures/anchoring/MatchPos.hpp"
#include "../../datastructures/anchoring/ClusterList.hpp"
#include "../../statistics/VarianceAccumulator.hpp"

unsigned int NumRemainingBases(DNALength curPos, DNALength intervalLength);

class IntervalSearchParameters {
public:
    bool  advanceHalf;
    int   globalChainType;
    float maxPValue;
    float aboveCategoryPValue;
    bool  warp;
    bool  fastMaxInterval; 
    bool  aggressiveIntervalCut;
    int   verbosity;
    float ddPValueThreshold;
    IntervalSearchParameters();
};


template<typename T_Sequence, typename T_AnchorList>
class DefaultWeightFunction {
public:
	float operator()(T_Sequence &text, T_Sequence &read, 
        T_AnchorList matchPosList);
};


template<typename T_Pos>
class MatchPosQueryOrderFunctor {
public:
	int operator()(T_Pos &pos);
};


template<typename T_MatchList>
void PrintLIS(T_MatchList &matchList, DNALength curPos,  
    DNALength curGenomePos, DNALength nextGenomePos, DNALength clp, 
    DNALength cle);

template<typename T_MatchList, typename T_SequenceDB>
void FilterMatchesAsLIMSTemplateSquare(T_MatchList &matches, 
    DNALength queryLength, DNALength limsTemplateLength, 
    T_SequenceDB &seqDB);
        
template<typename T_MatchList, typename T_SequenceBoundaryDB>
void AdvanceIndexToPastInterval(T_MatchList &pos, DNALength nPos,
    DNALength intervalLength, DNALength contigLength,
    T_SequenceBoundaryDB &SeqBoundary,
    DNALength startIndex, DNALength startIntervalBoundary,
    DNALength &index, DNALength &indexIntervalBoundary);

template<typename T_MatchList>  
int RemoveZeroLengthAnchors(T_MatchList &matchList);

template<typename T_MatchList>
int RemoveOverlappingAnchors(T_MatchList &matchList);

template<typename T_MatchList>
int SumAnchors(T_MatchList &pos, int start, int end);

template<typename T_MatchList,
         typename T_SequenceBoundaryDB>
void StoreLargestIntervals(
    T_MatchList &pos, 
    // End search for intervals at boundary positions
    // stored in seqBoundaries
	T_SequenceBoundaryDB & ContigStartPos,
	// parameters
	// How many values to search through for a max set.
	DNALength intervalLength,  
	// How many sets to keep track of
	int minSize,
	vector<DNALength> &start,
	vector<DNALength> &end);

template<typename T_MatchList,
         typename T_PValueFunction, 
         typename T_WeightFunction,
         typename T_SequenceBoundaryDB,
         typename T_ReferenceSequence,
         typename T_Sequence>
int FindMaxIncreasingInterval(
        // Input
        // readDir is used to indicate if the interval that is being stored is 
        // in the forward or reverse strand.  This is important later when 
        // refining alignments so that the correct sequence is aligned back 
        // to the reference.
        int readDir, 
        T_MatchList &pos, 
        // How many values to search through for a max set.
        DNALength intervalLength,  
        // How many sets to keep track of
        VectorIndex nBest, 
        // End search for intervals at boundary positions
        // stored in seqBoundaries
        T_SequenceBoundaryDB & ContigStartPos,
        // First rand intervals by their p-value
        T_PValueFunction &MatchPValueFunction,  
        // When ranking intervals, sum over weights determined by MatchWeightFunction
        T_WeightFunction &MatchWeightFunction,  
        // Output.
        // The increasing interval coordinates, 
        // in order by queue weight.
        WeightedIntervalSet &intervalQueue, 
        T_ReferenceSequence &reference, 
        T_Sequence &query,
        IntervalSearchParameters &params,
        std::vector<BasicEndpoint<ChainedMatchPos> > *chainEndpointBuffer,
        ClusterList &clusterList,
        VarianceAccumulator<float> &accumPValue, 
        VarianceAccumulator<float> &accumWeight,
        VarianceAccumulator<float> &accumNumAnchorBases);

template<typename T_MatchList,
         typename T_PValueFunction, 
         typename T_WeightFunction,
         typename T_SequenceBoundaryDB,
         typename T_ReferenceSequence,
         typename T_Sequence>
int FastFindMaxIncreasingInterval(
        int readDir, 
        T_MatchList &pos, 
        DNALength intervalLength,  
        VectorIndex nBest, 
        T_SequenceBoundaryDB & ContigStartPos,
        T_PValueFunction &MatchPValueFunction,  
        T_WeightFunction &MatchWeightFunction,  
        WeightedIntervalSet &intervalQueue, 
        T_ReferenceSequence &reference, 
        T_Sequence &query,
        IntervalSearchParameters &params,
        std::vector<BasicEndpoint<ChainedMatchPos> > *chainEndpointBuffer,
        ClusterList &clusterList,
        VarianceAccumulator<float> &accumPValue, 
        VarianceAccumulator<float> &accumWeight);

template<typename T_MatchList,
         typename T_PValueFunction, 
         typename T_WeightFunction,
         typename T_SequenceBoundaryDB,
         typename T_ReferenceSequence,
         typename T_Sequence>
int ExhaustiveFindMaxIncreasingInterval(
        int readDir, 
        T_MatchList &pos, 
        DNALength intervalLength,  
        VectorIndex nBest, 
        T_SequenceBoundaryDB & ContigStartPos,
        T_PValueFunction &MatchPValueFunction,  
        T_WeightFunction &MatchWeightFunction,  
        WeightedIntervalSet &intervalQueue, 
        T_ReferenceSequence &reference, 
        T_Sequence &query,
        IntervalSearchParameters &params,
        std::vector<BasicEndpoint<ChainedMatchPos> > *chainEndpointBuffer,
        ClusterList &clusterList,
        VarianceAccumulator<float> &accumPValue, 
        VarianceAccumulator<float> &accumWeight);

#include "FindMaxIntervalImpl.hpp"
#endif

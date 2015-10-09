#ifndef _BLASR_REGION_UTILS_HPP_
#define _BLASR_REGION_UTILS_HPP_

#include <algorithm>
#include <cmath>
#include "SMRTSequence.hpp"
#include "statistics/StatUtils.hpp"
#include "reads/ReadInterval.hpp"
#include "reads/RegionTable.hpp"

bool LookupHQRegion(int holeNumber, RegionTable &regionTable, 
    int &start, int &end, int &score);

template<typename T_Sequence>
bool MaskRead(T_Sequence &fastaRead, ZMWGroupEntry &zmwData,
    RegionTable &regionTable);


template<typename T_Sequence>
bool GetReadTrimCoordinates(T_Sequence &fastaRead,
	ZMWGroupEntry &zmwData,	RegionTable &regionTable,
	DNALength &readStart ,DNALength &readEnd, int &score);


// Given a vecotr of ReadInterval objects and their corresponding 
// directions, intersect each object with an interval 
// [hqStart, hqEnd), if there is no intersection or the intersected
// interval is less than minIntervalLength, remove this object and
// their corresponding directions; otherwise, replace this object 
// with the intersected interval and keep their directions. 
// Return index of the (left-most) longest subread interval in the
// updated vector.
int GetHighQualitySubreadsIntervals(
    std::vector<ReadInterval> & subreadIntervals, 
    std::vector<int> & subreadDirections, 
    int hqStart, int hqEnd, int minIntervalLength = 0); 

// Given a vector of subreads and a vector of adapters, return
// indices of all full-pass subreads.
std::vector<int> GetFullPassSubreadIndices(
    std::vector<ReadInterval> & subreadIntervals,
    std::vector<ReadInterval> & adapterIntervals);

// Given a vector of subreads and a vector of adapters, return
// index of the (left-most) longest subread which has both
// adapters before & after itself. If no full-pass subreads are
// available, return -1.
int GetLongestFullSubreadIndex(
    std::vector<ReadInterval> & subreadIntervals,
    std::vector<ReadInterval> & adapterIntervals);

// Given a vector of subreads and a vector of adapters, return
// index of the median length subread which has both
// adapters before & after itself. If no full-pass subreads are
// available, return -1.
int GetMedianLengthFullSubreadIndex(
    std::vector<ReadInterval> & subreadIntervals,
    std::vector<ReadInterval> & adapterIntervals);


// Given a vector of subreads and a vector of adapters, return
// index of the typical fullpass subread which can represent subreads
// of this zmw.
// * if there is no fullpass subread, return -1;
// * if number of fullpass subreads is less than 4, return index of the 
//   left-most longest subread
// * if number of fullpass subreads is greater than or equal 4, 
//   * if length of the longest read does not exceed 
//      meanLength + 1.96 * deviationLength
//     then, return index of the longest left-most subread 
//   * otherwise, return index of the second longest left-most subread
int GetTypicalFullSubreadIndex(
    std::vector<ReadInterval> & subreadIntervals,
    std::vector<ReadInterval> & adapterIntervals);

// Create a vector of n directions consisting of interleaved 0 and 1s.
void CreateDirections(std::vector<int> & directions, const int & n);

// Flop all directions in the given vector, if flop is true.
void UpdateDirections(std::vector<int> & directions, bool flop = false);

#include "utils/RegionUtilsImpl.hpp"

#endif

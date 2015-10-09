/*
 * =====================================================================================
 *
 *       Filename:  RegionUtils_gtest.cpp
 *
 *    Description:  Test alignment/utils/RegionUtils.hpp
 *
 *        Version:  1.0
 *        Created:  11/29/2012 05:11:56 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */


#include "gtest/gtest.h"
#include "utils/RegionUtils.hpp"
#include "reads/RegionTable.hpp"
#include "HDFRegionTableReader.hpp"
#include "pbdata/testdata.h"

class RegionUtilTestFixture: public testing::Test {
public:
    RegionUtilTestFixture() {
        fileName = baxFile1; 
        reader = new HDFRegionTableReader();
        int rev = reader->Initialize(fileName);
        EXPECT_TRUE(rev);
        reader->ReadTable(regionTable);
        reader->Close();
    }

    void SetUp() {
        /* The region table for hole number 14798 is as follows.
         * type start end  score
         *   1   0    712   -1 
         *   1   760  2040  -1 
         *   1   2098 3452  -1
         *   0   712  760   937
         *   0   2040 2098  741
         *   2   0    3424  819
         *   where type 1 = insertion, 0 = adapter, 2 = HQRegion
         * */
        hqStart = 0;
        hqEnd   = 3424;
        hqScore = 819;
        holeNumber = 14798;

        int srs[29][2] = {{95, 353},
              {404, 955},
              {998, 1243},
              {1289, 1529},
              {1576, 1824},
              {1869, 2132},
              {2178, 2444},
              {2486, 2727},
              {2767, 3004},
              {3047, 3286},
              {3329, 3537},
              {3583, 3818},
              {3862, 4096},
              {4139, 4404},
              {4447, 4673},
              {4717, 4956},
              {5004, 5228},
              {5275, 5508},
              {5552, 5790},
              {5835, 6072},
              {6113, 6338},
              {6389, 6604},
              {6647, 7145},
              {7189, 7421},
              {7468, 7693},
              {7733, 7951},
              {7989, 8197},
              {8243, 8482},
              {8522, 9244}};

        int rgn[59][5] = {{9, 1, 95, 353, -1},
            {9, 1, 404, 955, -1},
            {9, 1, 998, 1243, -1},
            {9, 1, 1289, 1529, -1},
            {9, 1, 1576, 1824, -1},
            {9, 1, 1869, 2132, -1},
            {9, 1, 2178, 2444, -1},
            {9, 1, 2486, 2727, -1},
            {9, 1, 2767, 3004, -1},
            {9, 1, 3047, 3286, -1},
            {9, 1, 3329, 3537, -1},
            {9, 1, 3583, 3818, -1},
            {9, 1, 3862, 4096, -1},
            {9, 1, 4139, 4404, -1},
            {9, 1, 4447, 4673, -1},
            {9, 1, 4717, 4956, -1},
            {9, 1, 5004, 5228, -1},
            {9, 1, 5275, 5508, -1},
            {9, 1, 5552, 5790, -1},
            {9, 1, 5835, 6072, -1},
            {9, 1, 6113, 6338, -1},
            {9, 1, 6389, 6604, -1},
            {9, 1, 6647, 7145, -1},
            {9, 1, 7189, 7421, -1},
            {9, 1, 7468, 7693, -1},
            {9, 1, 7733, 7951, -1},
            {9, 1, 7989, 8197, -1},
            {9, 1, 8243, 8482, -1},
            {9, 1, 8522, 9244, -1},
            {9, 0, 41, 95, 722},
            {9, 0, 353, 404, 784},
            {9, 0, 955, 998, 697},
            {9, 0, 1243, 1289, 804},
            {9, 0, 1529, 1576, 744},
            {9, 0, 1824, 1869, 844},
            {9, 0, 2132, 2178, 847},
            {9, 0, 2444, 2486, 666},
            {9, 0, 2727, 2767, 824},
            {9, 0, 3004, 3047, 883},
            {9, 0, 3286, 3329, 860},
            {9, 0, 3537, 3583, 891},
            {9, 0, 3818, 3862, 795},
            {9, 0, 4096, 4139, 906},
            {9, 0, 4404, 4447, 906},
            {9, 0, 4673, 4717, 886},
            {9, 0, 4956, 5004, 708},
            {9, 0, 5228, 5275, 957},
            {9, 0, 5508, 5552, 886},
            {9, 0, 5790, 5835, 733},
            {9, 0, 6072, 6113, 804},
            {9, 0, 6338, 6389, 666},
            {9, 0, 6604, 6647, 813},
            {9, 0, 7145, 7189, 795},
            {9, 0, 7421, 7468, 829},
            {9, 0, 7693, 7733, 774},
            {9, 0, 7951, 7989, 657},
            {9, 0, 8197, 8243, 760},
            {9, 0, 8482, 8522, 675},
            {9, 2, 0, 9244, 834}};

        for (int i = 0; i < 29; i++) {
            subreadIntervals.push_back(ReadInterval(srs[i][0], srs[i][1]));
        }

        for (int i = 0; i < 59; i++) {
            adapterIntervals.push_back(ReadInterval(rgn[i][2], rgn[i][3], rgn[i][4]));
        }
    }

    void TearDown() {
    }

    ~RegionUtilTestFixture() {
        delete reader;
    }

    HDFRegionTableReader * reader;
    RegionTable regionTable;
    string fileName;
    int hqStart, hqEnd, hqScore, holeNumber;

    // Data for testing Get*FullSubreadIndex.
    vector<ReadInterval> subreadIntervals;
    vector<ReadInterval> adapterIntervals;
};

TEST_F(RegionUtilTestFixture, LookupHQRegion) {
    int start, end, score;
    bool rev = LookupHQRegion(holeNumber, regionTable, start, end, score);
    EXPECT_EQ(rev, true);
    EXPECT_EQ(start, hqStart);
    EXPECT_EQ(end, hqEnd);
    EXPECT_EQ(score, hqScore);
}

TEST_F(RegionUtilTestFixture, GetHighQulitySubreadsIntervals) {
    vector<ReadInterval> intervals;
    intervals.push_back(ReadInterval(0, 712));
    intervals.push_back(ReadInterval(760, 2040));
    intervals.push_back(ReadInterval(2098, 3452));

    vector<int> directions;
    directions.push_back(0);
    directions.push_back(1);
    directions.push_back(0);

    int indx = GetHighQualitySubreadsIntervals(intervals, directions, hqStart, hqEnd);
    EXPECT_EQ(intervals.size(), 3);
    EXPECT_EQ(indx, 2);
    int starts [3] = {0, 760, 2098};
    int ends   [3] = {712, 2040, 3424};
    int ds     [3] = {0, 1, 0};
    for(int i=0; i < 3; i++) {
        EXPECT_EQ(intervals[i].start, starts[i]);
        EXPECT_EQ(intervals[i].end  , ends[i]  );
        EXPECT_EQ(directions[i]     , ds[i]    );
    }

    indx = GetHighQualitySubreadsIntervals(intervals, directions, hqStart, hqEnd, 800);
    EXPECT_EQ(intervals.size(), 2);
    // The first interval and its direction has been removed as the length is less
    // than 800.
    for(int i=0; i < 2; i++) {
        EXPECT_EQ(intervals[i].start, starts[i+1]);
        EXPECT_EQ(intervals[i].end  , ends[i+1]  );
    }
}

TEST_F(RegionUtilTestFixture, GetFullPassSubreadIndices ) {
    vector<int> vi = GetFullPassSubreadIndices(subreadIntervals, adapterIntervals);
    // vi = 0, ..., 27
    EXPECT_EQ(vi.size(), 28);
    for(int i=0; i < 28; i++) {
        EXPECT_EQ(vi[i], i);
    }
}

TEST_F(RegionUtilTestFixture, GetMedianLengthFullSubreadIndex) {
    int idx = GetMedianLengthFullSubreadIndex(subreadIntervals, adapterIntervals);
    EXPECT_EQ(idx, 18); // Median length subread (5552, 5790) 
}

TEST_F(RegionUtilTestFixture, GetLongestFullSubreadIndex) {
    int idx = GetLongestFullSubreadIndex(subreadIntervals, adapterIntervals);
    EXPECT_EQ(idx, 1); 
    // The longest full pass subread (404, 955) 
}

TEST_F(RegionUtilTestFixture, GetTypicalFullSubreadIndex) {
    int idx = GetTypicalFullSubreadIndex(subreadIntervals, adapterIntervals);
    EXPECT_EQ(idx, 22); 
    // Typical = the second longest full pass subread (6647, 7145)
}

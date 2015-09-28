/*
 * =====================================================================================
 *
 *       Filename:  FragmentCCSIterator_gtest.cpp
 *
 *    Description:  Test alignment/files/FragmentCCSIterator.hpp
 *
 *        Version:  1.0
 *        Created:  11/29/2012 04:51:29 PM
 *       Revision:  08/20/2014
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#include <algorithm>
#include <stdio.h>
#include <string.h>
#include "files/FragmentCCSIterator.hpp"
#include "reads/RegionTable.hpp"
#include "HDFRegionTableReader.hpp"
#include "pbdata/testdata.h"
#include <gtest/gtest.h>
using namespace std;

class FragmentCCSIteratorTestFixture: public testing::Test {
public:
    FragmentCCSIteratorTestFixture() {
    }

    void SetUp() {
        fileName = baxFile1;
        reader = new HDFRegionTableReader();
        ccs = new CCSSequence();
        rgn = new RegionTable();

        int rev = reader->Initialize(fileName);
        EXPECT_TRUE(rev);
        reader->ReadTable(*rgn);
        reader->Close();

        rgn->SortTableByHoleNumber();
    }

    void TearDown() {
        if (reader) delete reader;
        if (ccs) delete ccs;
        if (rgn) delete rgn;
    }

    ~FragmentCCSIteratorTestFixture() {
    }

    string fileName;
    HDFRegionTableReader * reader;
    CCSSequence * ccs;
    RegionTable * rgn;
    FragmentCCSIterator it;
};


TEST_F(FragmentCCSIteratorTestFixture, Initialize) {
    // void Initialize(CCSSequence *_seqPtr, RegionTable *_regionTablePtr) {
    ccs->HoleNumber(10);
    ccs->unrolledRead.Allocate(7000);
    it.Initialize(ccs, rgn);

    int numPasses = it.GetNumPasses();
    EXPECT_EQ(numPasses, 7);
/*
 * The region table of zmw 10 is:
 *
    (52,0): 10, 1, 0, 443, -1,
    (53,0): 10, 1, 487, 1168, -1,
    (54,0): 10, 1, 1213, 1907, -1,
    (55,0): 10, 1, 1956, 2619, -1,
    (56,0): 10, 1, 2668, 3423, -1,
    (57,0): 10, 1, 3474, 4205, -1,
    (58,0): 10, 1, 4256, 5949, -1,
    (59,0): 10, 1, 5997, 6161, -1,
    (60,0): 10, 0, 443, 487, 863,
    (61,0): 10, 0, 1168, 1213, 822,
    (62,0): 10, 0, 1907, 1956, 836,
    (63,0): 10, 0, 2619, 2668, 693,
    (64,0): 10, 0, 3423, 3474, 764,
    (65,0): 10, 0, 4205, 4256, 862,
    (66,0): 10, 0, 5949, 5997, 812,
    (67,0): 10, 2, 0, 4920, 788,
 *
 */
    int exp_directions[7] = {0, 1, 0, 1, 0, 1, 0};
    int exp_starts[7] = {0, 487, 1213, 1956, 2668, 3474, 4256};
    int exp_ends[7]  = {443, 1168, 1907, 2619, 3423, 4205, 4920};

    int passDirection, start, numBases;
    for(int i=0; i < numPasses; i++) {
        it.GetNext(passDirection, start, numBases);
        EXPECT_EQ(passDirection, exp_directions[i]);
        EXPECT_EQ(start, exp_starts[i]);
        EXPECT_EQ(start + numBases, exp_ends[i]);
    }
}

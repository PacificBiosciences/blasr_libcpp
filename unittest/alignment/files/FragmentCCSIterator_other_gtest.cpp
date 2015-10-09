/*
 * =====================================================================================
 *
 *       Filename:  CCSIterator_gtest.cpp
 *
 *    Description:  Test alignment/files/CCSIterator.hpp
 *
 *        Version:  1.0
 *        Created:  11/29/2012 04:51:02 PM
 *       Revision:  08/20/2014 
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#include "gtest/gtest.h"
#define private public
#include "files/CCSIterator.hpp"
#include "files/FragmentCCSIterator.hpp"
#include "reads/RegionTable.hpp"
#include "HDFRegionTableReader.hpp"
#include "pbdata/testdata.h"

using namespace std;

static const UInt HOLENUMBER = 76772;

// Adater - 0, Insert - 1, HQRegion - 2 
static const std::vector<RegionType> TYPES = {Adapter, Insert, HQRegion};
static const std::vector<std::string> TYPESTRS = {"Adapter", "Insert", "HQRegion"};

static const std::vector<RegionAnnotation> INSERTS = {
    RegionAnnotation(HOLENUMBER, 1, 0,   253,  -1),
    RegionAnnotation(HOLENUMBER, 1, 301, 678,  -1),
    RegionAnnotation(HOLENUMBER, 1, 724, 1101, -1),
    RegionAnnotation(HOLENUMBER, 1, 1150, 1534, -1),
    RegionAnnotation(HOLENUMBER, 1, 1575, 1956, -1),
    RegionAnnotation(HOLENUMBER, 1, 1999, 2379, -1),
    RegionAnnotation(HOLENUMBER, 1, 2417, 2803, -1),
    RegionAnnotation(HOLENUMBER, 1, 2852, 3245, -1),
    RegionAnnotation(HOLENUMBER, 1, 3287, 3727, -1),
    RegionAnnotation(HOLENUMBER, 1, 3778, 4176, -1),
    RegionAnnotation(HOLENUMBER, 1, 4221, 4618, -1),
    RegionAnnotation(HOLENUMBER, 1, 4661, 4862, -1)
};

static const std::vector<RegionAnnotation> ADAPTERS = {
    RegionAnnotation(HOLENUMBER, 0, 253, 301, 854),
    RegionAnnotation(HOLENUMBER, 0, 678, 724, 978),
    RegionAnnotation(HOLENUMBER, 0, 1101, 1150, 897),
    RegionAnnotation(HOLENUMBER, 0, 1534, 1575, 804),
    RegionAnnotation(HOLENUMBER, 0, 1956, 1999, 930),
    RegionAnnotation(HOLENUMBER, 0, 2379, 2417, 736),
    RegionAnnotation(HOLENUMBER, 0, 2803, 2852, 918),
    RegionAnnotation(HOLENUMBER, 0, 3245, 3287, 928),
    RegionAnnotation(HOLENUMBER, 0, 3727, 3778, 784),
    RegionAnnotation(HOLENUMBER, 0, 4176, 4221, 911),
    RegionAnnotation(HOLENUMBER, 0, 4618, 4661, 767)
};

static const std::vector<RegionAnnotation> HQREGION = {
    RegionAnnotation(HOLENUMBER, 2, 0, 4861, 865)
};

static const DNALength EXPECTED_HQSTART = 0;

static const DNALength EXPECTED_HQEND   = 4861;

static const DNALength EXPECTED_SCORE   = 865;

static const DNALength WHOLE_LENGTH = 5000;

static const int EXPECTED_NUM_SUBREADS = 12;


TEST(CCSFragmentIterator, Constructor) {
    std::vector<RegionAnnotation> regions = INSERTS;
    regions.insert(regions.end(), HQREGION.begin(), HQREGION.end());
    regions.insert(regions.end(), ADAPTERS.begin(), ADAPTERS.end());

    CCSSequence ccs;
    ccs.HoleNumber(HOLENUMBER);
    ccs.Allocate(WHOLE_LENGTH);
    ccs.unrolledRead.Allocate(WHOLE_LENGTH);

    RegionTable table;
    table.ConstructTable(regions, TYPESTRS);

    FragmentCCSIterator it;
    it.Initialize(&ccs, &table);

    EXPECT_EQ(it.subreadIntervals.size(), EXPECTED_NUM_SUBREADS);

    EXPECT_EQ(it.subreadIntervals[0], ReadInterval(0, 253, 865));

    EXPECT_EQ(it.subreadIntervals[EXPECTED_NUM_SUBREADS-1], ReadInterval(4661, 4861, 865));
}


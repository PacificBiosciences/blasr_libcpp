/*
 * ==================================================================
 *
 *       Filename:  RegionAnnotations_gtest.cpp
 *
 *    Description:  Test pbdata/reads/RegionAnnotations.hpp
 *
 *        Version:  1.0
 *        Created:  09/27/2015 03:54:55 PM
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * ==================================================================
 */
#include "gtest/gtest.h"
#define private public
#include "reads/ReadInterval.hpp"
#include "reads/RegionAnnotations.hpp"
#include <algorithm>
#include <iostream>

using namespace std;

static const UInt HOLENUMBER = 1720;

// Adater - 0, Insert - 1, HQRegion - 2 
static const std::vector<RegionType> TYPES = {Adapter, Insert, HQRegion};
static const std::vector<std::string> TYPESTRS = {"Adapter", "Insert", "HQRegion"};

static const std::vector<RegionAnnotation> REGIONS = {
    RegionAnnotation(HOLENUMBER, 2, 50,  900, 900),// hqregion
    RegionAnnotation(HOLENUMBER, 1, 700, 999, -1), // insert 
    RegionAnnotation(HOLENUMBER, 0, 649, 700, 700),// adapter
    RegionAnnotation(HOLENUMBER, 1, 300, 650, -1), // insert
    RegionAnnotation(HOLENUMBER, 0, 249, 329, 800),// adapter
    RegionAnnotation(HOLENUMBER, 1, 0,   250, -1)  // insert
};

static const std::vector<RegionAnnotation> EXPECTED_HQREGIONS = {
    RegionAnnotation(HOLENUMBER, 2, 50,  900, 900) // hqregion
};

static const DNALength EXPECTED_HQSTART = 50;

static const DNALength EXPECTED_HQEND   = 900;

static const DNALength WHOLE_LENGTH = 1200;

static const std::vector<RegionAnnotation> EXPECTED_ADAPTERS = {
    RegionAnnotation(HOLENUMBER, 0, 249, 329, 800),// adapter
    RegionAnnotation(HOLENUMBER, 0, 649, 700, 700) // adapter
};

static const std::vector<RegionAnnotation> EXPECTED_INSERTS = {
    RegionAnnotation(HOLENUMBER, 1, 0,   250, -1),// insert
    RegionAnnotation(HOLENUMBER, 1, 300, 650, -1),// insert
    RegionAnnotation(HOLENUMBER, 1, 700, 999, -1) // insert 
};

static const std::vector<ReadInterval> EXPECTED_SUBREAD_INTERVALS_BYADAPTER_NOHQ = {
    ReadInterval(0,   249, 0),   // by Adapter, subread score unknown.
    ReadInterval(329, 649, 0), 
    ReadInterval(700, 1200, 0)
};

static const std::vector<ReadInterval> EXPECTED_SUBREAD_INTERVALS_BYADAPTER_HQ = {
    ReadInterval(50,  249, 900), // by HQ, subread score = HQRegion score
    ReadInterval(329, 649, 900), 
    ReadInterval(700, 900, 900)
};

static const std::vector<ReadInterval> EXPECTED_SUBREAD_INTERVALS_NOHQ = {
    ReadInterval(0,   250, -1),  // not by adapter, not by HQ, use the original score.
    ReadInterval(300, 650, -1), 
    ReadInterval(700, 999, -1)
};

static const std::vector<ReadInterval> EXPECTED_SUBREAD_INTERVALS_HQ = {
    ReadInterval(50,  250, 900),  // by HQ, subread score = HQRegion score
    ReadInterval(300, 650, 900), 
    ReadInterval(700, 900, 900)
};

static const std::vector<ReadInterval> EXPECTED_ADAPTER_INTERVALS = {
    ReadInterval(249,  329, 800), 
    ReadInterval(649, 700, 700)
};

static const std::vector<RegionAnnotation> REGIONS_SORTED_BY_POS = {
    RegionAnnotation(HOLENUMBER, 1, 0,   250, -1), // insert
    RegionAnnotation(HOLENUMBER, 2, 50,  900, 900),// hqregion
    RegionAnnotation(HOLENUMBER, 0, 249, 329, 800),// adapter
    RegionAnnotation(HOLENUMBER, 1, 300, 650, -1), // insert
    RegionAnnotation(HOLENUMBER, 0, 649, 700, 700),// adapter
    RegionAnnotation(HOLENUMBER, 1, 700, 999, -1)  // insert 
};

static const std::vector<RegionAnnotation> REGIONS_SORTED_BY_TYPE = {
    RegionAnnotation(HOLENUMBER, 0, 249, 329, 800),// adapter
    RegionAnnotation(HOLENUMBER, 0, 649, 700, 700),// adapter
    RegionAnnotation(HOLENUMBER, 1, 0,   250, -1), // insert
    RegionAnnotation(HOLENUMBER, 1, 300, 650, -1), // insert
    RegionAnnotation(HOLENUMBER, 1, 700, 999, -1), // insert 
    RegionAnnotation(HOLENUMBER, 2, 50,  900, 900) // hqregion
};

/*
TEST(RegionAnnotationTest, Sort_By_Pos) {
    std::vector<RegionAnnotation> ras = REGIONS;
    std::sort(ras.begin(), ras.end());
    EXPECT_EQ(ras, REGIONS_SORTED_BY_POS);
}

TEST(RegionAnnotationTest, Sort_By_Type) {
    std::vector<RegionAnnotation> ras = REGIONS;
    std::sort(ras.begin(), ras.end(), compare_region_annotation_by_type);
    EXPECT_EQ(ras, REGIONS_SORTED_BY_TYPE);
}

TEST(RegionAnnotationsTest, Constructor) {
    RegionAnnotations ras(HOLENUMBER, REGIONS, TYPES);
    EXPECT_EQ(ras.table_, REGIONS_SORTED_BY_TYPE);
    EXPECT_EQ(ras.HoleNumber(), HOLENUMBER);
}

TEST(RegionAnnotationsTest, RegionAnnotationsOfType) {
    RegionAnnotations ras(HOLENUMBER, REGIONS, TYPES);
    EXPECT_EQ(ras.Adapters(), EXPECTED_ADAPTERS); 
    EXPECT_EQ(ras.HQRegions(), EXPECTED_HQREGIONS); 
    EXPECT_EQ(ras.Inserts(), EXPECTED_INSERTS); 
    EXPECT_EQ(ras.HQStart(), EXPECTED_HQSTART);
    EXPECT_EQ(ras.HQEnd(), EXPECTED_HQEND);
}

TEST(RegionAnnotationsTest, SubreadIntervals) {
    RegionAnnotations ras(HOLENUMBER, REGIONS, TYPES);
    vector<ReadInterval> ris = ras.SubreadIntervals(WHOLE_LENGTH, true, false);
    EXPECT_EQ(ris, EXPECTED_SUBREAD_INTERVALS_BYADAPTER_NOHQ);

    ris = ras.SubreadIntervals(WHOLE_LENGTH, true, true);
    EXPECT_EQ(ris, EXPECTED_SUBREAD_INTERVALS_BYADAPTER_HQ);

    ris = ras.SubreadIntervals(WHOLE_LENGTH, false, false);
    EXPECT_EQ(ris, EXPECTED_SUBREAD_INTERVALS_NOHQ);

    ris = ras.SubreadIntervals(WHOLE_LENGTH, false, true);
    EXPECT_EQ(ris, EXPECTED_SUBREAD_INTERVALS_HQ);
}

TEST(RegionAnnotationsTest, AdapterIntervals) {
    RegionAnnotations ras(HOLENUMBER, REGIONS, TYPES);
    EXPECT_EQ(ras.AdapterIntervals(), EXPECTED_ADAPTER_INTERVALS);
}
*/

TEST(RegionAnnotationsTest, SubreadIntervals_2) {
    std::vector<RegionAnnotation> regions({
            RegionAnnotation(HOLENUMBER, 0, 0, 112, -1)// adapter, no insert, no hq
            });
    RegionAnnotations ras(HOLENUMBER, regions, TYPES);

    vector<ReadInterval> ris = ras.SubreadIntervals(WHOLE_LENGTH, true, false);
    EXPECT_EQ(ris.size(), 1); // (112, WHOLE_LENGTH, -1)
    EXPECT_EQ(ris[0].start, 112);
    EXPECT_EQ(ris[0].end, WHOLE_LENGTH);

    // no insert, no hq && require adapter, require hq
    ris = ras.SubreadIntervals(WHOLE_LENGTH, true, true); 
    EXPECT_EQ(ris.size(), 0);

    // no require adapter, no require hq
    ris = ras.SubreadIntervals(WHOLE_LENGTH, false, false);// no insert
    EXPECT_EQ(ris.size(), 0); 

    // no require adapter, require hq
    ris = ras.SubreadIntervals(WHOLE_LENGTH, false, true); // no hq
    EXPECT_EQ(ris.size(), 0);
}

TEST(RegionAnnotationsTest, SubreadIntervals_3) {
    std::vector<RegionAnnotation> regions({
            RegionAnnotation(HOLENUMBER, 1, 0, 170, -1),// insert
            RegionAnnotation(HOLENUMBER, 2, 0, 0, 0)   //  hq length = 0
            });
    RegionAnnotations ras(HOLENUMBER, regions, TYPES);

    // require adapter, no require hq
    vector<ReadInterval> ris = ras.SubreadIntervals(WHOLE_LENGTH, true, false);
    EXPECT_EQ(ris.size(), 0);

    // require adapter, require hq
    ris = ras.SubreadIntervals(WHOLE_LENGTH, true, true);
    EXPECT_EQ(ris.size(), 0);

    // no require adapter, no require hq
    ris = ras.SubreadIntervals(WHOLE_LENGTH, false, false);
    EXPECT_EQ(ris.size(), 1); 
    EXPECT_EQ(ris[0], ReadInterval(0, 170, -1)); 

    // no require adapter, require hq
    ris = ras.SubreadIntervals(WHOLE_LENGTH, false, true);
    EXPECT_EQ(ris.size(), 0);
}

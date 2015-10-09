/*
 * ==================================================================
 *
 *       Filename:  RegionTypeMap_gtest.cpp
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
#include "reads/RegionTypeMap.hpp"

using namespace std;


// Adater - 0, Insert - 1, HQRegion - 2 
const vector<RegionType> TYPES = {Adapter, Insert, HQRegion};

TEST(RegionTypeMapTest, ToString) {
    EXPECT_EQ(RegionTypeMap::ToString(Adapter),  "Adapter");
    EXPECT_EQ(RegionTypeMap::ToString(HQRegion), "HQRegion");
    EXPECT_EQ(RegionTypeMap::ToString(Insert),   "Insert");
}

TEST(RegionTypeMapTest, ToRegionType) {
    EXPECT_EQ(RegionTypeMap::ToRegionType("Adapter"),  Adapter);
    EXPECT_EQ(RegionTypeMap::ToRegionType("HQRegion"), HQRegion);
    EXPECT_EQ(RegionTypeMap::ToRegionType("Insert"),   Insert);
}

TEST(RegionTypeMapTest, ToIndex) {
    // In most bas.h5 files, order of region types:
    std::vector<std::string> typeStrs = {"Insert", "Adapter", "HQRegion"};

    EXPECT_EQ(RegionTypeMap::ToIndex(Insert,   typeStrs), 0);
    EXPECT_EQ(RegionTypeMap::ToIndex(Adapter,  typeStrs), 1);
    EXPECT_EQ(RegionTypeMap::ToIndex(HQRegion, typeStrs), 2);

    EXPECT_EQ(RegionTypeMap::ToIndex("Insert",   typeStrs), 0);
    EXPECT_EQ(RegionTypeMap::ToIndex("Adapter",  typeStrs), 1);
    EXPECT_EQ(RegionTypeMap::ToIndex("HQRegion", typeStrs), 2);

    // Test given a different region type order.
    typeStrs = {"Insert", "HQRegion", "Adapter", "BarCode"};

    EXPECT_EQ(RegionTypeMap::ToIndex(Insert,   typeStrs), 0);
    EXPECT_EQ(RegionTypeMap::ToIndex(HQRegion, typeStrs), 1);
    EXPECT_EQ(RegionTypeMap::ToIndex(Adapter,  typeStrs), 2);

    EXPECT_EQ(RegionTypeMap::ToIndex("Insert",   typeStrs), 0);
    EXPECT_EQ(RegionTypeMap::ToIndex("HQRegion", typeStrs), 1);
    EXPECT_EQ(RegionTypeMap::ToIndex("Adapter",  typeStrs), 2);
}

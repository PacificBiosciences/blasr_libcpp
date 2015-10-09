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
#include "files/CCSIterator.hpp"
#include "reads/RegionTable.hpp"
#include "HDFRegionTableReader.hpp"
#include "pbdata/testdata.h"

using namespace std;

class CCSIteratorTestFixture: public testing::Test {
public:
    CCSIteratorTestFixture() {
    }

    void SetUp() {
        fileName = baxFile2; 
        reader = new HDFRegionTableReader();
        ccs = new CCSSequence();
        rgn = new RegionTable();

        int rev = reader->Initialize(fileName);
        EXPECT_TRUE(rev);
        reader->ReadTable(*rgn);
        reader->Close();
    }

    void TearDown() {
        if (reader) delete reader;
        if (ccs) delete ccs;
        if (rgn) delete rgn;
    }

    ~CCSIteratorTestFixture() {
    }

    string fileName;
    HDFRegionTableReader * reader;
    CCSSequence * ccs;
    RegionTable * rgn;
    CCSIterator it;
};


TEST_F(CCSIteratorTestFixture, Initialize) {
 
}

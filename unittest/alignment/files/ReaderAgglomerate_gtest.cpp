/*
 * =====================================================================================
 *
 *       Filename:  ReaderAgglomerate_gtest.cpp
 *
 *    Description:  Test alignment/files/ReaderAgglomerate.hpp
 *
 *        Version:  1.0
 *        Created:  02/25/2015 04:51:29 PM
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
#include "files/ReaderAgglomerate.hpp"
#include "pbdata/testdata.h"
#include <gtest/gtest.h>
using namespace std;

class ReaderAgglomerateTest: public testing::Test {
public:
    void SetUp() {
        reader = new ReaderAgglomerate();
    }

    void TearDown() {
        if (reader) {
            delete reader;
            reader = NULL;
        }
    }
    ReaderAgglomerate * reader;
};

#define INIT_READER(INFILE) \
    fn = INFILE; \
    reader->SetReadFileName(fn); \
    EXPECT_EQ(reader->Initialize(), 1);

TEST_F(ReaderAgglomerateTest, Initialize) {
    string fn;
    INIT_READER(fastaFile1)
    reader->Close();

    INIT_READER(baxFile1)
    reader->Close();

    INIT_READER(baxFile3)
    reader->Close();
}


#define SIMPLE(INFILE, SEQ) \
    INIT_READER(INFILE) \
    for(int i = 0; i < 10; i++) {\
        reader->GetNext(SEQ);\
    }\
    reader->Close();

TEST_F(ReaderAgglomerateTest, Simple) {
    string fn;

    // Fasta
    FASTASequence fastaSeq;
    SIMPLE(fastaFile1, fastaSeq)
    
    // Bax
    SMRTSequence smrtSeq;
    SIMPLE(baxFile1, smrtSeq)
    SIMPLE(baxFile3, smrtSeq)
}

#define GET_MOVIE_NAME(INFILE) \
    INIT_READER(INFILE) \
    reader->GetMovieName(movieName); \
    reader->Close();

TEST_F(ReaderAgglomerateTest, GetMovieName) {
    string fn;
    string movieName = "";
    // Fasta
    GET_MOVIE_NAME(fastaFile1)
    EXPECT_EQ(movieName, "/mnt/secondary-siv/testdata/BlasrTestData/utest/data/read.fasta");

    // Bax
    GET_MOVIE_NAME(baxFile1)
    EXPECT_EQ(movieName, "m130220_114643_42129_c100471902550000001823071906131347_s1_p0");

    GET_MOVIE_NAME(baxFile3)
    EXPECT_EQ(movieName, "m150223_190837_42175_c100735112550000001823160806051530_s1_p0");
}

#define GET_CHEMISTRY_TRIPLE(INFILE, A, B, C) \
    INIT_READER(INFILE) \
    bindingKit = sequencingKit = version; \
    reader->GetChemistryTriple(bindingKit, sequencingKit, version); \
    reader->Close(); \
    EXPECT_EQ(bindingKit, A); \
    EXPECT_EQ(sequencingKit, B); \
    EXPECT_EQ(version, C);

TEST_F(ReaderAgglomerateTest, GetChemistryTriple) {
    string fn;
    string bindingKit, sequencingKit, version;

    GET_CHEMISTRY_TRIPLE(baxFile3, "100356300", "100356200", "2.3")
} 

TEST_F(ReaderAgglomerateTest, ReadFromBam) {
    string fn (bamFile1);
    reader->SetReadFileName(fn);
    EXPECT_EQ(reader->Initialize(), 1);

    SMRTSequence seq;
    int ret, count=0;
    while (ret = reader->GetNext(seq) and ret != 0) { 
        count++;
    }

    EXPECT_EQ(count, 116);

    reader->Close();
}

TEST_F(ReaderAgglomerateTest, ReadsFromBam) {
    string fn (bamFile1);
    reader->SetReadFileName(fn);
    EXPECT_EQ(reader->Initialize(), 1);

    vector<SMRTSequence> seqs;
    int ret, count=0;
    while (ret = reader->GetNext(seqs) and ret != 0) { 
        count+ = seqs.size();
    }

    EXPECT_EQ(count, 116);

    reader->Close();
}

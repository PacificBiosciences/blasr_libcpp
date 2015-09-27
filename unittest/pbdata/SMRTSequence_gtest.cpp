/*
 * =====================================================================================
 *
 *       Filename:  SMRTSequenceSequence_gtest.cpp
 *
 *    Description:  Test pbdata/SMRTSequenceSequence.h
 *
 *        Version:  1.0
 *        Created:  10/29/2012 05:17:04 PM
 *       Revision:  08/19/2014 
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#include "gtest/gtest.h"
#include "SMRTSequence.hpp"

using namespace std;

Nucleotide seqnt[] = "ATATGGGGATTAGGGGATA"; 
const string seqst("ATATGGGGATTAGGGGATA"); 

class SMRTSequenceTest: public:: testing:: Test{
public:
    void SetUp() {
        smrt.seq = seqnt;
        int len = sizeof(seqnt) / sizeof(Nucleotide) - 1;
        smrt.length = len; 
        smrt.HoleNumber(1);
        smrt.SubreadStart(0);
        smrt.SubreadEnd  (19);
        smrt.AllocateDeletionQVSpace(len);
        for(int i=0; i < 19; i ++) {
            smrt.deletionQV[i] = i;
        }
    }
    void TearDown() {
        smrt.Free();
    }
    SMRTSequence smrt;
};


TEST_F(SMRTSequenceTest, Print) {
   stringstream ss;
   smrt.Print(ss);
   ASSERT_EQ(ss.str(), (string("SMRTSequence for zmw 1, [0, 19)\n")
                        + seqst + "\n"));
}


TEST_F(SMRTSequenceTest, GetDeletionQV) {
    for(int i = 0; i < smrt.length; i ++){
        ASSERT_EQ(smrt.GetDeletionQV(i), i);
    }
}



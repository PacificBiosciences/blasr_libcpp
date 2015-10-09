/*
 * =====================================================================================
 *
 *       Filename:  CCSSequence_gtest.cpp
 *
 *    Description:  Test pbdata/CCSSequence.hpp
 *
 *        Version:  1.0
 *        Created:  10/29/2012 05:17:04 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#include "gtest/gtest.h"
#include "CCSSequence.hpp"

using namespace std;

Nucleotide sr0[] = "ATATGGGGATTAGGGGATA"; 
Nucleotide sr1[] = "TAATCCCGTAATCCCGGTAT"; //rc = TAATCCCGTAATCCCGGTAT
Nucleotide sr2[] = "ATAGGGGGATTAGGGGATTCA"; 
Nucleotide adapter[] = "CCC";
Nucleotide unrolled_seq[] = "ATATGGGGATTAGGGGATACCCTAATCCCGTAATCCCGGTATCCCATAGGGGGATTAGGGGATTCA"; 
int sz0 = 19;
int sz1 = 20;
int sz2 = 21;
int adaptersz = 3;
int unrolledsz = sz0 + sz1 + sz2 + adaptersz * 2;
const int numSubreads = 3;


class CCSSequenceTest: public:: testing:: Test{
public:
    void CreateSMRTSequence(SMRTSequence & smrt, Nucleotide* seq,
            int holeNumber, int start, int end) {
        int size = end - start;
        smrt.seq = new Nucleotide[size];
        memcpy(smrt.seq, seq, size*sizeof(Nucleotide));
        smrt.length = size;
        smrt.deleteOnExit = false;

        smrt.HoleNumber  (holeNumber);
        smrt.SubreadStart(start);
        smrt.SubreadEnd  (end);

        stringstream ss;
    }

    void SetUp() {
        //Create a vector of subreads.
        subreads.resize(numSubreads);
        int s = 0; 
        CreateSMRTSequence(subreads[0], &sr0[0], 1, s, s + sz0); 
        s += sz0 + adaptersz; 
        CreateSMRTSequence(subreads[1], &sr1[0], 1, s, s + sz1); 
        s += sz1 + adaptersz; 
        CreateSMRTSequence(subreads[2], &sr2[0], 1, s, s + sz2); 

        //Create ccs
        stringstream ss;
        subreads[0].Print(ss);

        CreateSMRTSequence(ccs.unrolledRead, &unrolled_seq[0], 1, 0, unrolledsz);
        ccs.numPasses = numSubreads;
        CreateSMRTSequence((SMRTSequence&)ccs, &sr0[0], 1, 0, sz0); 
        ccs.numConsensusBases = sz0;
        ccs.passStartBase.resize(numSubreads);
        ccs.passNumBases.resize(numSubreads);
        ccs.passDirection.resize(numSubreads);
        s = 0;
        for(int i=0; i < ccs.numPasses; i++) {
            ccs.passStartBase[i] = subreads[i].SubreadStart();
            ccs.passDirection[i] = (i%2==0)?(0):(1);
            ccs.passNumBases[i] = subreads[i].length;
        }
    }

    void TearDown() {
        ccs.Free();
        for(int i=0; i<numSubreads; i++) {
            subreads[i].Free();
        }
    }

    CCSSequence ccs;
    vector<SMRTSequence> subreads;
};


TEST_F(CCSSequenceTest, Print) {
    stringstream ss, ss1;
    ccs.Print(ss);
    ccs.unrolledRead.Print(ss1);
    ASSERT_EQ(ss.str(), (string(
        "SMRTSequence for zmw 1, [0, 19)\nATATGGGGATTAGGGGATA\n")));
    ASSERT_EQ(ss1.str(), (string(
        "SMRTSequence for zmw 1, [0, 66)\nATATGGGGATTAGGGGATACCCTAATCCCGTAATCCCGGTATCCCATAGG\nGGGATTAGGGGATTCA\n")));
}

TEST_F(CCSSequenceTest, Explode) {
    vector<SMRTSequence> exploded_subreads;

    ccs.Explode(exploded_subreads);
    for(int i = 0; i < numSubreads; i++){
        stringstream ss, ss1;
        subreads[i].PrintSeq(ss);
        exploded_subreads[i].PrintSeq(ss1);
        ASSERT_EQ(ss.str(), ss1.str());
    }
}



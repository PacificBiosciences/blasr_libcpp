/*
 * =====================================================================================
 *
 *       Filename:  SAMHeaderPrinter_gtest.cpp
 *
 *    Description:  Test alignment/format/SAMHeaderPrinter.hpp
 *
 *        Version:  1.0
 *        Created:  03/24/2015 04:51:29 PM
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#define private public
#define protected public

#include <algorithm>
#include <stdio.h>
#include <string.h>
#include "format/SAMPrinter.hpp"
#include <gtest/gtest.h>
using namespace std;

TEST(SAMPrinterTest, AddMatchBlockCigarOps) {
    DNASequence qSeq; qSeq.Copy("XXXXXAAAAAGGGGGCCCCC");
    DNASequence tSeq; tSeq.Copy(     "AAAAANGGGGCCCCC");
    blasr::Block b; 
    b.qPos = 5;
    b.tPos = 0;
    b.length=15;
    vector<int> opSize;
    vector<char> opChar;

    DNALength qSeqPos = 0;
    DNALength tSeqPos = 0;

    const vector<int>  expOpSize = {5  , 1  , 9  };
    const vector<char> expOpChar = {'=', 'X', '='};

    AddMatchBlockCigarOps(qSeq, tSeq, b, qSeqPos, tSeqPos, opSize, opChar);

    EXPECT_EQ(opSize, expOpSize);
    EXPECT_EQ(opChar, expOpChar);
}

std::string merge_indels(const std::vector<int>& opSize, 
                         const std::vector<char> & opChar) {
    std::vector<int> opSize_ = opSize;
    std::vector<char> opChar_ = opChar;
    SAMOutput::MergeAdjacentIndels(opSize_, opChar_, 'X');
    std::string ret = "";
    SAMOutput::CigarOpsToString(opSize_, opChar_, ret);
    return ret;
}

TEST(SAMPrinterTest, MergeAdjacentIndels) {
    std::vector<int>  opSize({10,  1,   5,   7,  6});
    std::vector<char> opChar({'=', '=', '=', 'X', 'X'});

    EXPECT_EQ(merge_indels(opSize, opChar), "16=13X");

    opChar = std::vector<char>({'I', 'D', 'D', '=', 'I'});
    EXPECT_EQ(merge_indels(opSize, opChar), "6X4I7=6I");

    opChar = std::vector<char>({'I', 'D', 'D', 'I', 'I'});
    EXPECT_EQ(merge_indels(opSize, opChar), "6X17I");

    opChar = std::vector<char>({'=', 'D', 'D', 'I', 'I'});
    EXPECT_EQ(merge_indels(opSize, opChar), "10=6X7I");

    opChar = std::vector<char>({'I', 'D', '=', 'I', 'D'});
    EXPECT_EQ(merge_indels(opSize, opChar), "1X9I5=6X1I");

    opSize = std::vector<int >({1,  1 });
    opChar = std::vector<char>({'I','D'});
    EXPECT_EQ(merge_indels(opSize, opChar), "1X");

    opSize = std::vector<int >({1,  10 });
    opChar = std::vector<char>({'I','D'});
    EXPECT_EQ(merge_indels(opSize, opChar), "1X9D");
}

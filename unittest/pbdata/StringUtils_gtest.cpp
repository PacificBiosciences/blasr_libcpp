/*
 * ==================================================================
 *
 *       Filename:  StringUtils_gtest.cpp
 *
 *    Description:  Test pbdata/StringUtils.hpp
 *
 *        Version:  1.0
 *        Created:  03/28/2015 03:54:55 PM
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * ==================================================================
 */
#include "gtest/gtest.h"
#include "reads/ReadType.hpp"
#include "StringUtils.hpp"

using namespace std;


TEST(StringUtilTest, MakeReadGroupId) {

    string movieName = "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0";
    ReadType::ReadTypeEnum readType = ReadType::SUBREAD;
    string expectedReadGroupId = "b89a4406";
    EXPECT_EQ(MakeReadGroupId(movieName, readType), expectedReadGroupId);

    movieName = "movie32";
    readType = ReadType::CCS;
    expectedReadGroupId = "f5b4ffb6";
    EXPECT_EQ(MakeReadGroupId(movieName, readType), expectedReadGroupId);
}

TEST(StringUtilTest, Splice) {
    vector<string> tokens;

    Splice("movie/zmw/0_1", "/", tokens);
    vector<string> exp = {"movie", "zmw", "0_1"};
    EXPECT_EQ(tokens, exp); 

    string test = "abc,ef,12,4";
    Splice(test, ",", tokens);
    exp = vector<string>{"abc", "ef", "12", "4"};
    EXPECT_EQ(tokens, exp); 
    
    Splice(test, "ef,", tokens);
    exp = vector<string>{"abc,", "12,4"};
    EXPECT_EQ(tokens, exp);

    Splice("", ",", tokens);
    exp = vector<string>{""};
    EXPECT_EQ(tokens, exp); 

    Splice(",", ",", tokens);
    exp = vector<string>{"", ""};
    EXPECT_EQ(tokens, exp); 
    
    Splice(",abc,", ",", tokens);
    exp = vector<string>{"", "abc", ""};
    EXPECT_EQ(tokens, exp); 

    Splice("abc,", ",", tokens);
    exp = vector<string>{"abc", ""};
    EXPECT_EQ(tokens, exp); 

    Splice(",abc", ",",  tokens);
    exp = vector<string>{"", "abc"};
    EXPECT_EQ(tokens, exp); 

    Splice("abc", "abc",  tokens);
    exp = vector<string>{"", ""};
    EXPECT_EQ(tokens, exp); 

    Splice("a\tb\tc", "\t",  tokens);
    exp = vector<string>{"a", "b", "c"};
    EXPECT_EQ(tokens, exp); 
}

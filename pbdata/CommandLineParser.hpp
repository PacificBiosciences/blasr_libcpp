#ifndef _BLASR_COMMAND_LINE_PARSER_HPP_
#define _BLASR_COMMAND_LINE_PARSER_HPP_

#include <vector>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <assert.h>
#include <sstream>
#include "StringUtils.hpp"

class CommandLineParser {
public:
    enum ErrorValue {
        CLGood,
        CLBadOption,
        CLMissingOption,
        CLMissingValue,
        CLInvalidInteger,
        CLInvalidPositiveInteger,
        CLInvalidNonNegativeInteger,
        CLInvalidFloat,
        CLInvalidPositiveFloat,
        CLInvalidNonNegativeFloat };


    enum OptionType {
        Flag,
        Integer,
        PositiveInteger, // > 0
        NonNegativeInteger, // >= 0
        IntegerList,
        Float,
        PositiveFloat,  // > 0
        NonNegativeFloat, // >= 0
        String,
        StringList };

    std::vector<bool*> boolValues;
    std::vector<int*> intValues;
    std::vector<float*> floatValues;
    std::vector<std::string*> stringValues;
    std::vector<std::vector<std::string> *>stringListValues;
    std::vector<std::vector<int> *> intListValues;
    std::vector<int*> flagList;
    std::vector<std::string> optionList;
    std::vector<OptionType>    optionTypeList;
    std::vector<int>    optionValueIndexList;
    std::vector<std::string> descriptions;
    std::vector<char>   optionRequired;
    std::vector<char>   optionUsed;
    std::vector<char>   named;

    std::string programName;
    std::string programSummary;
    std::string conciseHelp;
    std::string verboseHelp;
    std::string helpString;
    std::string examples;
    std::string version;

    int lineLength;
    int numUnnamedOptions;
    bool specialVersionFlag;

    CommandLineParser();

    void SetProgramSummary(std::string summaryp);

    void SetHelp(std::string _help);

    void SetConciseHelp(std::string _conciseHelp) ;

    void SetProgramName(std::string namep);

    void SetVersion(std::string versionp);

    void SetVerboseHelp(std::string helpp);

    void SetExamples(std::string examplesp);

    void RegisterPreviousFlagsAsHidden();

    void RegisterVersionFlag(bool *value);

    void RegisterFlagOption(std::string option, bool *value, 
        std::string description, bool required=false);

    void RegisterIntOption(std::string option, int *value, 
        std::string description, OptionType type, 
        bool required=false);

    void RegisterFloatOption(std::string option, float *value, 
        std::string description, OptionType type, bool required=false);

    void RegisterStringOption(std::string option, std::string *value, 
        std::string description, bool required=false);

    void RegisterStringListOption(std::string option, 
        std::vector<std::string> *value, std::string description, 
        bool required=false);

    void RegisterIntListOption(std::string option, std::vector<int> *value, 
        std::string description, bool required=false);

    int IsOption(char *str);

    int IsInteger(char *str);

    int IsFloat(char *str);

    int FindOption(char *option);

    static void CommandLineToString(int argc, char* argv[], 
        std::string& commandLine);

    int ParseCommandLine(int argc, char* argv[], bool isProgramNameOnlyAllowed=false);

    int ParseCommandLine(int argc, char* argv[], 
        std::vector<std::string> &unflaggedValues, bool isProgramNameOnlyAllowed=false);

    ErrorValue ParseOption(int optionIndex,
        int &argi, int argc, char *argv[]);

    void PrintErrorMessage(ErrorValue ev, char *option);

    ErrorValue ParseFlag(int optionValueIndex);

    ErrorValue ParseInteger(int optionValueIndex, 
        int &argi, int argc, char *argv[]);

    ErrorValue ParsePositiveInteger(int optionValueIndex, 
        int &argi, int argc, char *argv[]);

    ErrorValue ParseNonNegativeInteger(int optionValueIndex, 
        int &argi, int argc, char *argv[]);

    ErrorValue ParseFloat(int optionValueIndex, 
        int &argi, int argc, char *argv[]);

    ErrorValue ParsePositiveFloat(int optionValueIndex, 
        int &argi, int argc, char *argv[]);

    ErrorValue ParseNonNegativeFloat(int optionValueIndex, 
        int &argi, int argc, char *argv[]);

    ErrorValue ParseString(int optionValueIndex, int &argi, 
        int argc, char *argv[]);

    bool IsValuedOption(OptionType optType);

    ErrorValue ParseIntList(int optionValueIndex, int &argi, 
        int argc, char *argv[]);

    ErrorValue ParseStringList(int optionValueIndex, int &argi, 
        int argc, char *argv[]);

    void PrintVersion();

    void PrintUsage();

    int GetNextWordLength(std::string &text, int pos);

    void PrintIndentedText(std::ostream &out, std::string &text, 
        int allLineIndent, int lineLength = 80, int firstLineIndent=0);

    unsigned int GetMaxOptionLength();

    ErrorValue PrintErrorOnMissingOptions();
};

#endif

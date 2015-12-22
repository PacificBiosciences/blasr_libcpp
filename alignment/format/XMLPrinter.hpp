#ifndef _BLASR_XML_PRINTER_HPP_
#define _BLASR_XML_PRINTER_HPP_

#include "utils/SimpleXMLUtils.hpp"

namespace XMLOutput{

template<typename T_Alignment, typename T_Sequence> 
void Print(T_Alignment &alignment, 
    T_Sequence &query, T_Sequence &text, std::ostream &out,
    int qPrintStart = 0,
    int tPrintStart = 0,
    int maxPrintLength = 50) {
    (void)(qPrintStart); (void)(tPrintStart); (void)(maxPrintLength);
    /*
     * Sample alignment:
     *
     <hit name="x15_y33_1220208-0008_m081205_152444_Uni_p2_b15" unalignedLength="1301" start="1051" end="1016" strand="-" targetStart="1" targetEnd="44" targetStrand="+">
     <zScore value="-6.091"/>
     <nInsert value="1" percent="2.86" />
     <nDelete value="9" percent="25.71" />
     <nMismatch value="1" percent="2.86" />
     <nCorrect value="24" percent="68.57" />
     <alignment><query>
     AG--CGTTCC-TATGG-TG-GGGTCGTTA-ACT---GTCGCCAG
     </query><target>
     AGCCCG-TCCTTATGGTTGAGGGTTGTTACACTTCGGTCGCCAG
     </target></alignment>
     </hit>
     */
    char strand[2] = {'+', '-'};
    std::string tAlignStr, alignStr, qAlignStr;
    CreateAlignmentStrings(alignment, query.seq, text.seq, tAlignStr, alignStr, qAlignStr);
    int alignLength = tAlignStr.size();
    if (alignLength == 0) {
        alignLength = 1; // Make sure there are no divide by zero.
        alignment.nIns = 0;
        alignment.nDel = 0;
        alignment.nMismatch = 0;
        alignment.nMatch = 0;
    }
    out << BeginDataEntry(std::string("hit"),
            CreateKeywordValuePair(std::string("name"), alignment.qName) +
            CreateKeywordValuePair(std::string("unalignedLength"), alignment.qLength) +
            CreateKeywordValuePair(std::string("start"), alignment.qPos) + 
            CreateKeywordValuePair(std::string("end"), alignment.qPos + alignment.qAlignLength) +
            CreateKeywordValuePair(std::string("strand"), strand[alignment.qStrand]) + 
            CreateKeywordValuePair(std::string("targetStart"), alignment.tPos) +													 
            CreateKeywordValuePair(std::string("targetEnd"), alignment.tPos + alignment.tAlignLength) +
            CreateKeywordValuePair(std::string("targetStrand"), strand[alignment.tStrand])) 
        << std::endl;
    out << CreateDataEntry(std::string("zScore"),
            CreateKeywordValuePair(std::string("value"), alignment.zScore)) 
        << std::endl;
    out << CreateDataEntry(std::string("nInsert"),
            CreateKeywordValuePair(std::string("value"), alignment.nIns) + " " +
            CreateKeywordValuePair(std::string("percent"), alignment.nIns*0.5/alignLength)) 
        << std::endl;
    out << CreateDataEntry(std::string("nDelete"),
            CreateKeywordValuePair(std::string("value"), alignment.nDel) + " " + 
            CreateKeywordValuePair(std::string("percent"), alignment.nDel*0.5/alignLength)) 
        << std::endl;
    out << CreateDataEntry(std::string("nMismatch"),
            CreateKeywordValuePair(std::string("value"), alignment.nMismatch) + " " + 
            CreateKeywordValuePair(std::string("percent"), alignment.nMismatch*0.5/alignLength)) 
        << std::endl;
    out << CreateDataEntry(std::string("nCorrect"),
            CreateKeywordValuePair(std::string("value"), alignment.nMatch) + " " +
            CreateKeywordValuePair(std::string("percent"), alignment.nMatch*0.5/alignLength)) 
        << std::endl;

    out << CreateStartEntry(std::string("alignment"), std::string("")) 
        << CreateStartEntry(std::string("query"), std::string("")) 
        << std::endl;
    out << qAlignStr << std::endl;
    out << CreateEndEntry(std::string("query")) << std::endl;
    out << CreateStartEntry(std::string("target"), std::string("")) 
        << std::endl;
    out << tAlignStr << std::endl;
    out << CreateEndEntry(std::string("target")) << std::endl;
    out << CreateEndEntry(std::string("alignment")) << std::endl;
    out << CreateEndEntry(std::string("hit")) << std::endl;
}
}
#endif

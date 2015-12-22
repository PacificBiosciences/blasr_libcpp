#ifndef _BLASR_SAM_READER_HPP_
#define _BLASR_SAM_READER_HPP_

#include "SAMKeywordValuePair.hpp"
#include "ReadGroup.hpp"
#include "ReferenceSequence.hpp"

#include "AlignmentSet.hpp"
#include "../StringUtils.hpp"
#include "../utils.hpp"

template<typename T_ReferenceSequence, typename T_ReadGroup, typename T_SAMAlignment>
class SAMReader {
  public:
  int lineNumber;
  std::ifstream samFile;
  std::istream *samFilePtr;
  bool Initialize(std::string samFileName);

  void Close();

  enum LineType {Blank, HSHeader, HSSequence, HSReadGroup, HSProgram, HSComment, Alignment, Error};

  void GetLine(std::istream &in, std::string &line);

  bool LineTypeIsHeader(LineType lineType);

  bool PeekLineIsHeader(std::istream &in);

  LineType GetLineType(std::string &line);
 
  void StoreKVPairs(std::string line, std::vector<SAMKeywordValuePair> &kvPairs);
 
  void StoreHeader(std::vector<SAMKeywordValuePair> &kvPairs,
                  AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments);

  void StoreReferenceSequence(std::vector<SAMKeywordValuePair> &kvPairs,
                              AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments);

  void StoreReadGroup(std::vector<SAMKeywordValuePair> &kvPairs,
                      AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments);

  void StoreAlignment(std::string & line,
                      AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments);
    
  // Not implemented
  void StoreProgram(std::vector<SAMKeywordValuePair> &kvPairs,
                    AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments ) {(void)(kvPairs); (void)(alignments);}

  void Read(std::string samFileName, AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments);

  std::vector<std::string> ReadHeader(AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments);
 
  void Read(AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments);

  bool GetNextAlignment(SAMAlignment &alignment);
};

#include "SAMReaderImpl.hpp"

#endif

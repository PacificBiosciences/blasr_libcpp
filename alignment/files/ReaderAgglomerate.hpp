#ifndef _BLASR_READER_AGGLOMERATE_HPP_
#define _BLASR_READER_AGGLOMERATE_HPP_

#include <cstdlib>

#include "Enumerations.h"
#include "reads/ReadType.hpp"
#include "files/BaseSequenceIO.hpp"

#include "FASTAReader.hpp"
#include "FASTQReader.hpp"
#include "CCSSequence.hpp"
#include "SMRTSequence.hpp"
#include "StringUtils.hpp"

#include "HDFBasReader.hpp"
#include "HDFCCSReader.hpp"

#ifdef USE_PBBAM
#include "pbbam/BamFile.h"
#include "pbbam/EntireFileQuery.h"
#include "pbbam/GroupQuery.h"
#include "pbbam/BamRecord.h"
#endif

class ReaderAgglomerate : public BaseSequenceIO {
  FASTAReader fastaReader;
  FASTQReader fastqReader;
  int readQuality;
  int stride;
  int start;
  float subsample;
  bool useRegionTable;
  bool ignoreCCS;
  ReadType::ReadTypeEnum readType;

public:
  //
  // Create interfaces for reading hdf 
  //
  T_HDFBasReader<SMRTSequence>  hdfBasReader;
  HDFCCSReader<CCSSequence>     hdfCcsReader;
  vector<SMRTSequence>          readBuffer;
  vector<CCSSequence>           ccsBuffer;
  string readGroupId;

public:
  void SetToUpper();

  void InitializeParameters();
  ReaderAgglomerate();

  ReaderAgglomerate(float _subsample);

  ReaderAgglomerate(int _stride);

  ReaderAgglomerate(int _start, int _stride);

  void GetMovieName(string &movieName);

  /// Get BindingKit, SequencingKit and Base Caller Version from h5.
  ///
  /// /param [out] sequencingKit - sequencingKit from 
  /// /ScanData/RunInfo/SequencingKit.
  ///
  /// /param [out] bindingKit - BindingKit from 
  /// /ScanData/RunInfo/BindingKit.
  ///
  /// /param [out] baseCallerVersion - Base Caller Version
  /// from /PulseData/BaseCalls/ChangeListID.
  ///
  void GetChemistryTriple(string & bindingKit, 
                          string & sequencingKit, 
                          string & baseCallerVersion);
 
  bool FileHasZMWInformation();

  void SkipReadQuality();

  void IgnoreCCS();

  void UseCCS();

  int Initialize(string &pFileName);

  bool SetReadFileName(string &pFileName);

  int Initialize(FileType &pFileType, string &pFileName);

  bool HasRegionTable();

  int Initialize();

  ReaderAgglomerate &operator=(ReaderAgglomerate &rhs);

  bool Subsample(float rate);

  // Set read type to SUBREAD, CCS, or UNKNOWN.
  void SetReadType(const ReadType::ReadTypeEnum & readType_);
  
  // returns read type, SUBREAD, CCS, or UNKNOWN
  ReadType::ReadTypeEnum GetReadType();

public:
  int GetNext(FASTASequence &seq);
  int GetNext(FASTQSequence &seq);
  int GetNext(SMRTSequence &seq);
  int GetNext(CCSSequence &seq);
  int GetNext(vector<SMRTSequence> & reads);

  template<typename T_Sequence>
      int GetNext(T_Sequence & seq, int & randNum);

  int GetNextBases(SMRTSequence & seq, bool readQVs);

  int Advance(int nSteps);

  void Close();

#ifdef USE_PBBAM
public:
  // Define reader to fetch sequences from bam.
  PacBio::BAM::BamFile * bamFilePtr;
  PacBio::BAM::EntireFileQuery * entireFileQueryPtr;
  PacBio::BAM::EntireFileQuery::iterator bamIterator;
  PacBio::BAM::QNameQuery * zmwGroupQueryPtr;
  PacBio::BAM::QNameQuery::iterator zmwGroupIterator;
#endif
};


template<typename T_Sequence>
int ReadChunkByNReads(ReaderAgglomerate &reader, vector<T_Sequence> &reads, int maxNReads);

template<typename T_Sequence>
int ReadChunkBySize (ReaderAgglomerate &reader, vector<T_Sequence> &reads, int maxMemorySize);

#include "files/ReaderAgglomerateImpl.hpp"

#endif

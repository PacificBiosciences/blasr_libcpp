#ifndef _BLASR_AMOS_AFG_BAS_WRITER_HPP_
#define _BLASR_AMOS_AFG_BAS_WRITER_HPP_

#include <cmath>
#include <fstream>
#include <sstream>
#include "SMRTSequence.hpp"

class AfgBasWriter {
    // TODO correct file type? outFile
    // TODO change SMRTSequence to FASTQ sequence? then maybe fasta2afg would be 
    // easier to write
    std::string afgFileName;
    std::ofstream afgOut;
    bool firstRecord;
    int recordCount;
    int defaultQuality;
    // Use 122 = '{' - 1 because '{' and '}' are the message delimiters in AMOS.
    static const unsigned char maxAfgQuality = 122; 
    static const unsigned char charToQuality = 48;
    static const unsigned char minAfgQuality = charToQuality + 1; 
    static const int lineLength = 80;

  private:

    unsigned char pacbioQVtoPhredQV(unsigned char pacbio);
    
  public:
    // Assume that closing the afg file must be done
    // manually and not in a destructor.
    ~AfgBasWriter();

    AfgBasWriter();

    void Initialize(std::string _afgFileName);

    void Close();

    int Write(SMRTSequence &seq);

    void SetDefaultQuality(int _defaultQuality);

  private:

    void WriteHeader();

    void WriteOpen(void);

    void WriteIdentifier(SMRTSequence &seq);

    void WriteClose();

    void WriteBases(SMRTSequence &seq);

    void WriteQualities(SMRTSequence &seq);
};

#endif

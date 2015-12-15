#ifndef _BLASR_FASTA_READER_HPP_
#define _BLASR_FASTA_READER_HPP_
#include <stdint.h>
#include <string>
#include "FASTASequence.hpp"
#include "metagenome/SequenceIndexDatabase.hpp"

class FASTAReader {
protected:
    GenomeLength fileSize;
    int fileDes;
    char* filePtr;
    GenomeLength curPos;
    int padding;
    char endOfReadDelim;
    char readStartDelim;
    bool doToUpper;
    unsigned char *convMat;
    //
    // Quick check to see how much to read.
    //
    void SetFileSize(); 

    void ReadTitle(GenomeLength &p, char *&title, int &titleLength); 

public:
    bool computeMD5;
    std::string curReadMD5;

    void Init(); 

    FASTAReader(); 

    FASTAReader(std::string &fileName); 

    void SetSpacePadding(int _padding); 

    void SetToUpper(); 
    
    //
    // Synonym for Init() for consistency.
    //
    int Initialize(std::string &seqInName); 

    int Init(std::string &seqInName, int passive=0); 

    void AdvanceToTitleStart(GenomeLength &p, char delim='>'); 

    void CheckValidTitleStart(GenomeLength &p, char delim='>'); 

    GenomeLength ReadAllSequencesIntoOne(FASTASequence &seq, SequenceIndexDatabase<FASTASequence> *seqDBPtr=NULL); 

    void ReadTitle(GenomeLength &p, FASTASequence & seq);

    int GetNext(FASTASequence &seq); 
    /*
       Advance to the read nSeq forward.

input: nSeq, the number of sequences to skip.
output: 
returns 1 if after advancing nSeq sequences, the file pointer is pointing to 
a new sequence.
0 otherwise. 
A return value of 0 will signal that the file is done being processed if it is 
iterting over reads.
*/
    int Advance(int nSeq); 

    int CriticalGetNext(FASTASequence &seq); 
    
    int ConcatenateNext(FASTASequence &cur); 

    void Close(); 

    void ReadAllSequences(vector<FASTASequence> &sequences); 

}; 


#endif // _BLASR_FASTA_READER_HPP_

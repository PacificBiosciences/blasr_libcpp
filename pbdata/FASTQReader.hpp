#ifndef _BLASR_FASTQ_READER_HPP_
#define _BLASR_FASTQ_READER_HPP_

#include "FASTASequence.hpp"
#include "FASTAReader.hpp"
#include "FASTQSequence.hpp"
#include "qvs/QualityValue.hpp"

class FASTQReader : public FASTAReader {
public:

    FASTQReader();

    int GetNext(FASTASequence &seq); // return 1 if suceed.

    unsigned char phredQVtoPacbioQV(unsigned char phredQV);

    int GetNext(FASTQSequence &seq); 

    int Advance(int nSteps); 
};


#endif // _BLASR_FASTQ_READER_HPP_

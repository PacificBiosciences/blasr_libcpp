#ifndef _BLASR_CCS_SEQUENCE_HPP_
#define _BLASR_CCS_SEQUENCE_HPP_

#include "Enumerations.h"
#include "SMRTSequence.hpp"
#include "VectorUtils.hpp"


//
// A CCS Sequence is both a SMRTSequence itself, and contains a list of SMRTSequences.
//
class CCSSequence : public SMRTSequence {
private:
    void UNDEFINED(const std::string & msg) const {
        std::cout << msg  << " is not defined in CCSSequence." << std::endl; exit(1);
    }

public:
	UInt numPasses;
	UInt numConsensusBases;
	std::vector<DNALength> passStartPulse, passNumPulses, passStartBase, passNumBases;
	std::vector<Byte> passDirection;
	std::vector<Byte>      adapterHitBefore, adapterHitAfter, adapterHitConfidence;
	//
	// The CCS Sequence originates from a full length read.  That read
	// is stored here for reference later on.  The ccs read is stored in
	// the inherited fields from SMRTSequence so that it may be worked
	// with as if it were a normal non-ccs sequence.
	//
	SMRTSequence      unrolledRead;

public:
    inline ~CCSSequence();
	void Free(); 

    UInt HoleNumber(void) const;

    CCSSequence & HoleNumber(const UInt holeNumber);

    DNALength SubreadStart(void) const {UNDEFINED("SubreadStart");}

    DNALength SubreadEnd(void) const {UNDEFINED("SubreadEnd");}

    DNALength SubreadLength(void) const {UNDEFINED("SubreadLength");}

	int GetStorageSize(); 

    /// \name
    /// \{
        /// In the first iteration, Explode simply pulls the subreads out
        /// that are used in the ccs.   Eventually, it will pull out all
        /// high-quality subreads.
        ///
	void Explode(std::vector<SMRTSequence> &subreads); 
    /// \}
};

inline CCSSequence::~CCSSequence() {
    CCSSequence::Free();
}

#endif // _BLASR_CCS_SEQUENCE_HPP_

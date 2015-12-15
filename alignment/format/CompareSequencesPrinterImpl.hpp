#ifndef _BLASR_COMPARE_SEQUENCES_PRINTER_IMPL_HPP_
#define _BLASR_COMPARE_SEQUENCES_PRINTER_IMPL_HPP_

template<typename T_Alignment, typename T_QuerySequence, typename T_TargetSequence>
void CompareSequencesOutput::Print(T_Alignment &alignment, 
    T_QuerySequence &qseq, T_TargetSequence &tseq, 
    std::ostream &out, bool refForward) {

	std::string queryStr, alignStr, textStr;
	CreateAlignmentStrings(alignment, qseq, tseq, textStr, alignStr, queryStr);

	if (refForward == false) {
		if (alignment.qStrand == 1 and alignment.tStrand == 0) {
			DNALength alignedSeqToEnd = 0;
			//			DNALength alignedTSeqToEnd = 0;
			if (alignment.blocks.size() > 0) {
				// First compute the offset of the reverse of the substring that was aligned.
				
				alignedSeqToEnd = alignment.qLength - (alignment.qAlignedSeqPos + alignment.qAlignedSeq.length);
				DNALength alignEndToSubstrEnd = alignment.qAlignedSeq.length - (alignment.qPos + alignment.blocks[alignment.blocks.size()-1].qPos + alignment.blocks[alignment.blocks.size()-1].length);
				alignment.qPos = alignEndToSubstrEnd;
			}
			alignment.qAlignedSeqPos = alignedSeqToEnd;
			alignment.qStrand = 0;
			alignment.tStrand = 1;
							
		}	
	}
	
	PrintCompareSequencesAlignmentStats(alignment, out);
	// change the spaces in the align string to *s for easy parsing of alignment
	VectorIndex i;
	for (i = 0; i < alignStr.size(); i++ ) { 
		if (alignStr[i] == ' ') alignStr[i] = '*';
	}

	if (refForward == false and alignment.tStrand == 1) {
		//
		// Build reverse complement strings.
		//
        std::string queryStrRC, alignStrRC, textStrRC;
		queryStrRC.resize(queryStr.size());
		alignStrRC.resize(alignStr.size());
		textStrRC.resize(alignStr.size());
		
		DNALength pos;
		DNALength alignStringLength = alignStr.size();
		for (pos = 0; pos < alignStringLength; pos++ ) {
			if (queryStr[pos] != '-') {
				queryStrRC[alignStringLength-pos-1] = ReverseComplementNuc[static_cast<int>(queryStr[pos])];
			}
			else {
				queryStrRC[alignStringLength-pos-1] = '-';
			}
			alignStrRC[alignStringLength-pos-1] = alignStr[pos];
			
			if (textStr[pos] != '-') {
				textStrRC[alignStringLength-pos-1] = ReverseComplementNuc[static_cast<int>(textStr[pos])];
			}
			else {
				textStrRC[alignStringLength-pos-1] = '-';
			}
		}
		queryStr = queryStrRC;
		alignStr = alignStrRC;
		textStr  = textStrRC;
	}
					
    // Headers of m5 format are: 
    // qName qSeqLength qStart qEnd qStrand 
    // tName tSeqLength tStart tEnd tStrand
    // score numMatch numMismatch numIns numDel
    // mapQV qAlignedSeq matchPattern tAlignedSeq
	out << queryStr << " " << alignStr << " " << textStr << std::endl;
}

#endif

/*
 * =====================================================================================
 *
 *       Filename:  FilterCriteria.h
 *
 *    Description:  Criteria for filtering alignment hits. 
 *
 *        Version:  1.0
 *        Created:  04/14/2015 11:33:00 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#include "FilterCriteria.hpp"

constexpr float Score::errorunit;

Score::Score(const float & value, const ScoreSign & sign)
    : _value(value)
    , _sign(sign)
{ }

Score::Score(const Score & another) 
    : _value(another._value)
    , _sign(another._sign)
{ }

Score::~Score() { }

float Score::Value() const { return _value; }

ScoreSign Score::Sign() const { return _sign; }

bool Score::operator == (const Score & another) const {
    return (another.Value() - this->_value < errorunit and 
            another.Value() - this->_value > -errorunit);
}

bool Score::BetterThan(const Score & another) const {
    if (this->_value == another.Value()) return false;
    if (this->_sign == ScoreSign::POSITIVE) return (this->_value > another.Value());
    else return (this->_value < another.Value());
}

bool Score::BetterThanOrEqual(const Score & another) const {
    return BetterThan(another) or ((*this) == another);
}

bool Score::WorseThan(const Score & another) const {
    return (not BetterThanOrEqual(another));
}

/// -----------------------------------------
/// HitPolicy 
/// -----------------------------------------
HitPolicy::HitPolicy(const std::string & hitPolicyStr, const ScoreSign & sign) {
    // Create a multiple hit policy.
    std::string str(hitPolicyStr);
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);

    if (str == "RANDOM")    { _hit = HitPolicyEnum::RANDOM; }
    else if (str == "ALL")  { _hit = HitPolicyEnum::ALL; }
    else if (str == "ALLBEST") { _hit = HitPolicyEnum::ALLBEST; }
    else if (str == "RANDOMBEST") { _hit = HitPolicyEnum::RANDOMBEST; }
    else if (str == "LEFTMOST") { _hit = HitPolicyEnum::LEFTMOST; }
    else {
        std::cout <<"ERROR, the specified multiple hit policy "
                  << hitPolicyStr <<" is not supported." << std::endl;
        exit(1);
    }

    _sign = sign;
}

HitPolicy::~HitPolicy(void) { }

ScoreSign HitPolicy::Sign() const {
    return _sign;
}

bool HitPolicy::IsRandom() const { return _hit == HitPolicyEnum::RANDOM; }

bool HitPolicy::IsRandombest() const { return _hit == HitPolicyEnum::RANDOMBEST; }

bool HitPolicy::IsAll() const {return _hit == HitPolicyEnum::ALL; }

bool HitPolicy::IsAllbest() const {return _hit == HitPolicyEnum::ALLBEST; }

bool HitPolicy::IsLeftmost() const {return _hit == HitPolicyEnum::LEFTMOST; }

const std::string HitPolicy::ToString() const {
    switch (_hit) {
        case (HitPolicyEnum::RANDOM)    : return "random";
        case (HitPolicyEnum::RANDOMBEST): return "randombest";
        case (HitPolicyEnum::ALL)       : return "all";
        case (HitPolicyEnum::ALLBEST)   : return "allbest";
        case (HitPolicyEnum::LEFTMOST)  : return "leftmost";
        default: {
            assert(false);
        }
    }
}

std::ostream & HitPolicy::operator << (std::ostream & os) {
    return os << this->ToString();
}

const std::string HitPolicy::Help(const std::string & pad) const{
    std::stringstream ss;
    ss << "(" << this->ToString() 
       << ") Specify a policy to treat multiple hits from "
       << "[all, allbest, random, randombest, leftmost]\n" << pad
       << "  all       : report all alignments.\n" << pad
       << "  allbest   : report all equally top scoring alignments.\n"  << pad
       << "  random    : report a random alignment.\n" << pad
       << "  randombest: report a random alignment from multiple "
       << "equally top scoring alignments.\n" << pad
       << "  leftmost  : report an alignment which has the best alignment" 
       << "score and has the smallest mapping coordinate in any reference.";
    return ss.str();
}

std::vector<T_AlignmentCandidate*> 
HitPolicy::Apply(const std::vector<T_AlignmentCandidate*> alnPtrs, 
                 const bool & createRand,
                 const int  & passedRand) const {
    // Shallow copy pointers.
    if (alnPtrs.empty() or IsAll()) return alnPtrs;

    int rint = createRand?rand():passedRand;
    if (IsRandom()) {
        return std::vector<T_AlignmentCandidate*>({alnPtrs[rint%alnPtrs.size()]});
    }

    std::vector<T_AlignmentCandidate*> ret = alnPtrs;
    // Sort AlignmentCandidate pointers according to score and target start position.
    ScoreSign sign = this->Sign();
    if (sign == ScoreSign::POSITIVE) {
        std::cout << "ERROR: ScoreSign POSITIVE not supported yet." << std::endl;
        exit(1);
    }
    std::sort(ret.begin(), ret.end(), SortAlignmentPointersByScore()); 

    // Get the best alignments whose alignment scores are the best.
    // Assume that input alignments share the same query name and
    // are sorted by score and tPos asscendingly: worst, ...., best
    int bestScore = ret[0]->score;
    ret.erase(std::remove_if(ret.begin(), ret.end(),
              [&bestScore](const T_AlignmentCandidate* x)->bool{return x->score != bestScore;}),
              ret.end());

    if (IsAllbest()) {
        return ret;
    } else if (IsRandombest()) {
        return std::vector<T_AlignmentCandidate*>({ret[rint%ret.size()]});
    } else if (IsLeftmost()) {
        return std::vector<T_AlignmentCandidate*>({ret[0]});
    } else {
        assert("Unsupported hit policy" == 0);
    }
}

/*
std::vector<AlignmentCandidate<>> HitPolicy::Apply(std::vector<AlignmentCandidate<>> & records) const {
    if (records.empty() or IsALL()) return;
    sort(records.begin(), records.end(), compareByQNameScoreTStart);
    assert(records.size() > 0);

    if (IsRandom) {
        //records = vector<AlignmentCandidate>({records[
    }
}
*/

#ifdef USE_PBBAM
bool HitPolicy::compareByQNameScoreTStart(const PacBio::BAM::BamRecord & a, const PacBio::BAM::BamRecord & b) const {
    assert(a.Impl().HasTag(AS) and b.Impl().HasTag(AS));
    const int aScore = a.Impl().TagValue(AS).ToInt32();
    const int bScore = b.Impl().TagValue(AS).ToInt32();
    if (a.FullName() == b.FullName()) {
        if (aScore == bScore)
            return a.ReferenceStart() < b.ReferenceEnd();
        return Score(aScore, _sign).WorseThan(Score(bScore, _sign));
    }
    return (a.FullName() < b.FullName());
}

bool HitPolicy::compareByScoreTStart(const PacBio::BAM::BamRecord & a, const PacBio::BAM::BamRecord & b) const {
    assert(a.Impl().HasTag(AS) and b.Impl().HasTag(AS));
    assert(a.FullName() == b.FullName());
    const int aScore = a.Impl().TagValue(AS).ToInt32();
    const int bScore = b.Impl().TagValue(AS).ToInt32();
    if (aScore == bScore)
        return a.ReferenceStart() < b.ReferenceEnd();
    return Score(aScore, _sign).WorseThan(Score(bScore, _sign));
}

std::vector<PacBio::BAM::BamRecord> 
HitPolicy::Apply(const std::vector<PacBio::BAM::BamRecord> & records,
                 const bool & createRand,
                 const int  & passedRand) const {
    if (records.empty() or IsAll()) return records;

    int rint = createRand?rand():passedRand;
    //cout << "FilterCriteria " << ", " << records[0].FullName() << ", " << rint << endl;
    if (IsRandom()) {
        return std::vector<PacBio::BAM::BamRecord>({records[rint%records.size()]});
    }

    std::vector<PacBio::BAM::BamRecord> ret = records;
    // Sort bam records according to score and target start position.
    ScoreSign sign = this->Sign();
    std::sort(ret.begin(), ret.end(), 
              [&sign](const PacBio::BAM::BamRecord & a, 
                      const PacBio::BAM::BamRecord & b)->bool {
                      assert(a.Impl().HasTag(AS) and b.Impl().HasTag(AS));
                      assert(a.FullName() == b.FullName());
                      const int aScore = a.Impl().TagValue(AS).ToInt32();
                      const int bScore = b.Impl().TagValue(AS).ToInt32();
                      if (aScore == bScore) return a.ReferenceStart() < b.ReferenceEnd();
                      else return Score(aScore, sign).WorseThan(Score(bScore, sign));
              });

    // Get the best alignments whose alignment scores are the best.
    // Assume that input alignments share the same query name and
    // are sorted by score and tPos asscendingly: worst, ...., best
    int bestScore = ret[0].Impl().TagValue(AS).ToInt32();
    ret.erase(std::remove_if(ret.begin(), ret.end(),
              [&bestScore](const PacBio::BAM::BamRecord & x)->bool{return x.Impl().TagValue(AS).ToInt32() != bestScore;}),
              ret.end());

    if (IsAllbest()) {
        return ret;
    } else if (IsRandombest()) {
        return std::vector<PacBio::BAM::BamRecord>({ret[rint%ret.size()]});
    } else if (IsLeftmost()) {
        return std::vector<PacBio::BAM::BamRecord>({ret[0]});
    } else {
        assert("Unsupported hit policy" == 0);
    }
}
#endif

/// -----------------------------------------
/// FilterCriteria
/// -----------------------------------------
FilterCriteria::FilterCriteria(const DNALength & minAlnLength,
                               const float & minPctSimilarity,
                               const float & minPctAccuracy,
                               const bool & useScore,
                               const Score & score)
    : _minAlnLength(minAlnLength)
    , _minPctSimilarity(minPctSimilarity)
    , _minPctAccuracy(minPctAccuracy)
    , _useScore(useScore)
    , _scoreCutoff(score)
    , _verbose(false)
{ }

FilterCriteria::~FilterCriteria(void) { }

/*
const std::string FilterCriteria::Help(const std::string & pad) {
    return std::string (pad + minAlnLengthHelp() +
                        pad + minPctSimilarityHelp() +
                        pad + minPctAccuracyHelp() +
                        pad + scoreSignHelp() + 
                        pad + scoreCutoffHelp());
}*/

const std::string FilterCriteria::MinAlnLengthHelp() {
    return std::string("(") + std::to_string(_minAlnLength) + ") " +
           "Report alignments only if their lengths are greater than minAlnLength.";
}

const std::string FilterCriteria::MinPctSimilarityHelp() {
    return std::string("(") + std::to_string(static_cast<int>(_minPctSimilarity)) + ") " +
           "Report alignments only if their percentage similairty is greater than minPctSimilarity.";
}

const std::string FilterCriteria::MinPctAccuracyHelp() {
    return std::string("(") + std::to_string(static_cast<int>(_minPctAccuracy)) + ") " + 
           "Report alignments only if their percentage accuray is greater than minAccuracy.";
}

const std::string FilterCriteria::ScoreSignHelp() {
    return "(-1) Whether higher or lower scores are better. -1: lower is better; 1: higher is better.";
}

const std::string FilterCriteria::ScoreCutoffHelp() {
    return "(INF) Report alignments only if their scores are no worse than score cut off.";
}

bool FilterCriteria::MakeSane(std::string & errMsg) const{
    if (_minPctSimilarity > 100 or _minPctSimilarity < 0) {
        errMsg = "ERROR, minimum similarity not in [0, 100].";
        return false;
    }
    if (_minPctAccuracy > 100 or _minPctAccuracy < 0) {
        errMsg = "ERROR, minimum accuracy not in [0, 100].";
        return false;
    }
    return true;
}

void FilterCriteria::Verbose(bool verbose) { _verbose = verbose; }

/*
template<typename T_TSequence, typename T_QSequence>
bool FilterCriteria::Satisfy(AlignmentCandidate<T_TSequence, T_QSequence> & a) const;
//bool FilterCriteria::Satisfy(AlignmentCandidate<> & a) const {
    float pctAccuracy = 100 * a.nMatch / static_cast<float>(a.nMismatch
                        + a.nMatch + a.nIns + a.nDel);
    Score s(a.score, _scoreCutoff.Sign());
    return Satisfy(a.qAlignedSeqLength, a.pctSimilarity, pctAccuracy, s);
}
*/

bool FilterCriteria::Satisfy(const DNALength & alnLength, 
                             const float & pctSimilarity, 
                             const float & pctAccuracy,
                             const Score & score) const {
    if (alnLength < _minAlnLength) {
        if (_verbose > 0) 
            std::cout << "Alignment length " << alnLength << " is too short." << std::endl;
        return false;
    } 
    if (pctSimilarity < _minPctSimilarity) {
        if (_verbose > 0) 
            std::cout << "Percentage similarity " << pctSimilarity << " is too low." << std::endl;
        return false;
    }
    if (pctAccuracy < _minPctAccuracy) {
        if (_verbose > 0)
            std::cout << "Percentage accuracy " << pctAccuracy << " is too low." << std::endl;
        return false;
    }
    if (_useScore and not score.BetterThanOrEqual(_scoreCutoff)) {
        if (_verbose) 
            std::cout << "Alignment score " << score.Value() << " worse than cut off." << std::endl;
        return false;
    }

    return true;
}

#ifdef USE_PBBAM
bool FilterCriteria::Satisfy(const PacBio::BAM::BamRecord & record) const {
    assert(record.IsMapped() and record.Impl().HasTag(AS));
    DNALength alnLength = static_cast<DNALength>(record.Sequence(PacBio::BAM::Orientation::NATIVE, true, true).size());
    PacBio::BAM::Cigar cigar = record.CigarData();

    uint32_t nMatch = 0, nMismatch = 0, nIns = 0, nDel = 0;
    for(auto op: cigar) {
        uint32_t n = op.Length();
        switch (op.Type()) {
            case PacBio::BAM::CigarOperationType::SEQUENCE_MATCH: nMatch += n; break;
            case PacBio::BAM::CigarOperationType::SEQUENCE_MISMATCH: nMismatch += n; break;
            case PacBio::BAM::CigarOperationType::ALIGNMENT_MATCH: nMismatch += n; break;
            case PacBio::BAM::CigarOperationType::INSERTION: nIns += n; break;
            case PacBio::BAM::CigarOperationType::DELETION: nDel += n; break;
            default: break;
        }
    }
    //TODO: Use = X instead of M in blasr.

    float pctSimilarity = 100 * nMatch / static_cast<float>(alnLength);
    float pctAccuracy = 100 * nMatch / static_cast<float>(nMismatch + nMatch + nIns + nDel);
    Score s(static_cast<float>(record.Impl().TagValue("AS").ToInt32()), ScoreSign::NEGATIVE);
    return Satisfy(alnLength, pctSimilarity, pctAccuracy, s);
}
#endif

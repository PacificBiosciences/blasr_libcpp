#include <assert.h>
#include "SAMHeaderPrinter.hpp"

const std::string SAMVERSION("1.5");
const std::string PBBAMVERSION("3.0.1");
const std::string PACBIOPL("PACBIO");

std::vector<SAMHeaderItem> MakeSAMHeaderItems(const std::string & fromString){
    std::vector<SAMHeaderItem> items;

    std::vector<std::string> vs;
    Splice(fromString, ";", vs);
    std::vector<std::string>::iterator it;
    for (it = vs.begin(); it != vs.end(); it++) {
        items.push_back(SAMHeaderItem(*it));
    }
    return items;
}

// SAMHeaderItem
SAMHeaderItem::SAMHeaderItem(const std::string & fromString) {
    std::size_t pos = fromString.find("=");
    if (pos != std::string::npos) {
        _key = fromString.substr(0, pos);
        _val = fromString.substr(pos + 1);
    }
}

std::string SAMHeaderItem::ToString() {
    std::stringstream ss;
    if (_key != "") ss << _key << "=" << _val;
    return ss.str();
}

// SAMHeaderTag
SAMHeaderTag::SAMHeaderTag(const std::string & fromString) {
    size_t pos = fromString.find(":");
    if (pos != string::npos) {
        _tagName = fromString.substr(0, pos);
        string tagValue = fromString.substr(pos + 1);
        if (tagValue.find("=") != std::string::npos) {
            AddItems(tagValue);
        } else {
            _tagValue = tagValue;
        }
    } else {
        cout << "Unable to parse SAM/BAM header" << fromString << endl;
        exit(1);
    }
}

std::string SAMHeaderTag::ToString() {
    std::stringstream ss;
    if (_tagName != "") {
        ss << _tagName << ":";
        if (_tagValue != "") {
            ss << _tagValue;
        } else {
            std::vector<SAMHeaderItem>::iterator it;
            for(it = _tagItems.begin(); it != _tagItems.end(); it++) {
                if (it != _tagItems.begin() and (*it).ToString() != "") {
                    ss << ";";
                }
                ss << (*it).ToString();
            }
        }
    }
    return ss.str();
}

std::string SAMHeaderTag::TagName() const {
    return _tagName;
}

void SAMHeaderTag::AddItem(SAMHeaderItem & item) {
    _tagItems.push_back(item);
}

void SAMHeaderTag::AddItem(const std::string & fromString) {
    _tagItems.push_back(SAMHeaderItem(fromString));
}

void SAMHeaderTag::AddItems(const std::string & fromString) {
    // SAM Header Items: key=value;key=value
    std::vector<SAMHeaderItem> items = MakeSAMHeaderItems(fromString);
    _tagItems.insert(_tagItems.begin(), items.begin(), items.end());
}

// SAMHeaderGroup
SAMHeaderGroup::SAMHeaderGroup(const std::string & fromString) {
    if (fromString == "" || fromString[0] != '@') return;
    std::vector<std::string> vs;
    Splice(fromString.substr(1), "\t", vs);
    if (vs.size() >= 1) {
        std::vector<std::string>::iterator it = vs.begin();
        _groupName = (*it);
        it++;
        for(; it != vs.end(); it++) {
            _tags.push_back(SAMHeaderTag(*it));
        }
    }
}

std::string SAMHeaderGroup::ToString() {
    std::stringstream ss;
    ss << "@" << _groupName;
    SAMHeaderTags::iterator it;
    for (it = _tags.begin(); it != _tags.end(); it++) {
        ss << "\t" << (*it).ToString();
    }
    return ss.str();
}

std::ostream & operator << (std::ostream& os, SAMHeaderGroup & g) {
    return os << g.ToString();
}

/// returns true if this SAMHeaderGroup contains a tag tagName.
bool SAMHeaderGroup::HasTag(std::string tagName) {
    for(SAMHeaderTags::iterator it = _tags.begin(); it != _tags.end(); it++) {
        if ((*it).TagName() == tagName) {
            return true;
        } 
    }
    return false;
}

/// returns value of a tag if it exists, otherwise, return an empty string. 
std::string SAMHeaderGroup::Tag(const std::string & tagName) {
    for(SAMHeaderTags::iterator it = _tags.begin(); it != _tags.end(); it++) {
        if ((*it).TagName() == tagName) {
            return (*it).ToString();
        }
    }
    return "";
}

// SAMHeaderGroupWithID
std::string SAMHeaderGroupWithID::ID() const {return _id;}

SAMHeaderGroupWithID::SAMHeaderGroupWithID(const std::string & fromString)
: SAMHeaderGroup(fromString) {
    assert (HasTag("ID"));
    if (not HasTag("ID")) {
        assert("ERROR! SAM Header read/program group must has ID tag." == 0);
    }
    _id = Tag("ID");
}

// SAMHeaderRG
SAMHeaderRG::SAMHeaderRG(const std::string & id, const std::string & pl, 
            const std::string & pu, const std::vector<SAMHeaderItem> & dsItems) {
    _groupName = "RG";
    _id = id;
    _tags.push_back(SAMHeaderTag("ID", id));
    _tags.push_back(SAMHeaderTag("PU", pu));
    _tags.push_back(SAMHeaderTag("PL", pl));
    _tags.push_back(SAMHeaderTag("DS", dsItems));
}

SAMHeaderRG::SAMHeaderRG(const std::string & fromString)
:SAMHeaderGroupWithID(fromString) {
    if (_groupName != "RG") {
        assert("ERROR! SAM Header read group should start with @RG" == 0);
    }
}

// SAMHeaderPG
SAMHeaderPG::SAMHeaderPG(const std::string & id, const std::string & progName, 
        const std::string & progVersion, const std::string & commandLine) {
    _groupName = "PG";
    _id = id;
    _tags.push_back(SAMHeaderTag("ID", id));
    _tags.push_back(SAMHeaderTag("PN", progName));
    _tags.push_back(SAMHeaderTag("VN", progVersion));
    _tags.push_back(SAMHeaderTag("CL", commandLine));
}

SAMHeaderPG::SAMHeaderPG(const std::string & fromString)
: SAMHeaderGroupWithID(fromString) {
    if (_groupName != "PG") {
        assert("ERROR! SAM Header program group must start with @PG" == 0);
    }
}

// SAMHeaderSQ
SAMHeaderSQ::SAMHeaderSQ(const std::string & sn, const DNALength & ln, const std::string & md5)
: SAMHeaderSQ(sn, std::to_string(ln), md5) {}

SAMHeaderSQ::SAMHeaderSQ(const std::string & sn, const string & ln, const std::string & md5)
: SAMHeaderGroup() {
    _groupName = "SQ";
    _tags.push_back(SAMHeaderTag("SN", sn));
    _tags.push_back(SAMHeaderTag("LN", ln));
    _tags.push_back(SAMHeaderTag("M5", md5));
}

SAMHeaderSQ::SAMHeaderSQ(const std::string & fromString)
: SAMHeaderGroup(fromString) {
    if (_groupName != "SQ") {
        assert("ERROR! SAM Header soring order must start with @SO" == 0);
    }
}

// SAMHeaderPrinter
std::string SAMHeaderPrinter::ToString() {
    std::stringstream ss;
    ss << _hd.ToString()
        << _sqs.ToString()
        << _rgs.ToString() 
        << _pgs.ToString()
        << _cos.ToString();
        return ss.str();
}
/*
/// Add a SAM Header @RG entry, including the following tags:
/// ID (identifier), PL (platform), PU (platform unit), 
/// DS (description, which may have many items)
SAMHeaderPrinter & SAMHeaderPrinter::AddRG(std::string & id, std::string & pl, std::string & pu, std::string & ds) {
    _rgs.Add(SAMHeaderRG(id, pl, pu, ds));
    return *this;
}*/

/// Add a SAM Header @RG entry from string, which may contain arbitary tags.
SAMHeaderPrinter & SAMHeaderPrinter::AddRG(const std::string & fromString) {
    _rgs.Add(SAMHeaderRG(fromString)); 
    return *this;
}

/// Add a SAM Header @PG entry, including the following tags:
/// ID, progName, progVersion, commandLine. @PG must have unique ID.
SAMHeaderPrinter & SAMHeaderPrinter::AddPG(std::string & id, std::string & progName, std::string & progVersion, std::string & commandLine) {
    _pgs.Add(SAMHeaderPG(id, progName, progVersion, commandLine));
    return *this;
}

/// Add a SAM Header @PG entry from string, which may contain arbitary tags.
SAMHeaderPrinter & SAMHeaderPrinter::AddPG(const std::string & fromString) {
    _pgs.Add(SAMHeaderPG(fromString));
    return *this;
}

/// Add a SAM Header @CO entry from string.
SAMHeaderPrinter & SAMHeaderPrinter::AddCO(const std::string & fromString) {
    _cos.Add(fromString);
    return *this;
}

/// \returns SAMHeaderGroup @HD
SAMHeaderGroup SAMHeaderPrinter::MakeHD(const std::string & sortingOrder) {
    std::stringstream ss;
    ss << "@HD" << "\t" 
        << "VN:" << SAMVERSION << "\t" 
        << "SO:" << sortingOrder << "\t"
        << "pb:" << PBBAMVERSION << std::endl;
    return SAMHeaderGroup(ss.str());
}

/// \returns SAMHeaderSQs, SAM Header @SQ groups
SAMHeaderSQs SAMHeaderPrinter::MakeSQs(SequenceIndexDatabase<FASTASequence> & seqdb) {
    SAMHeaderSQs sqs;
    for (int i = 0; i < seqdb.nSeqPos-1; i++) {
        std::string md5 = "";
        if (static_cast<int>(seqdb.md5.size()) == seqdb.nSeqPos - 1) {
            md5 = seqdb.md5[i];
        }
        string seqName;
        seqdb.GetName(i, seqName);
        sqs.Add(SAMHeaderSQ(seqName, seqdb.GetLengthOfSeq(i), md5));
    }
    return sqs;
}

SAMHeaderPrinter::SAMHeaderPrinter(const std::string & sortingOrder,
        SequenceIndexDatabase<FASTASequence> & seqdb,
        const std::vector<std::string> & readsFiles,
        const ReadType::ReadTypeEnum & readType,
        const SupplementalQVList & samQVs,
        const std::string & progName, 
        const std::string & progVersion,
        const std::string & commandLine)  
    :_sortingOrder(sortingOrder), _readsFiles(readsFiles), _seqdb(seqdb) {
    if (readsFiles.size() == 0) {
        assert("Must specify input reads files" == 0);
    }
    // all read types supported, no check needed.

    // Determine fileType based on extension of the first read file.
    std::string rf = readsFiles[0];
    BaseSequenceIO::DetermineFileTypeByExtension(rf, fileType);

    // Make @HD, @SQ, @RG 
    _hd  = MakeHD(sortingOrder);
    _sqs = MakeSQs(seqdb);
    _rgs = MakeRGs(readsFiles, readType, samQVs);

    // Make PGs and COs 
    _pgs = MakePGs(readsFiles, progName, progVersion, commandLine);
    _cos = MakeCOs(readsFiles);
}

/// \param[in] readsFiles - incoming reads files in BAM or other formats
/// \param[in] readType - read type, must be either SUBREAD or CCS or UNKNOWN
/// \param[in] samQVs - SupplementalQVList, an object that handles which 
///                      QVs to print in SAM/BAM file.
SAMHeaderRGs SAMHeaderPrinter::MakeRGs(const std::vector<std::string> & readsFiles,
        const ReadType::ReadTypeEnum & readType,
        const SupplementalQVList & samQVs) {
    SAMHeaderRGs rgs;

    if (fileType != FileType::PBBAM and fileType != FileType::PBDATASET) {
        ReaderAgglomerate * reader = new ReaderAgglomerate();
        assert(reader != nullptr);
        std::vector<std::string>::const_iterator rfit;
        for(rfit = readsFiles.begin(); rfit != readsFiles.end(); rfit++) {
            std::string rf(*rfit);
            reader->SetReadFileName(rf);
            reader->SetReadType(readType);
            reader->Initialize();

            // Get movie name from ReaderAgglomerate
            std::string movieName;
            reader->GetMovieName(movieName);

            string bindingKit, sequencingKit, baseCallerVersion;
            reader->GetChemistryTriple(bindingKit, sequencingKit, baseCallerVersion);
            reader->Close();

            std::vector<SAMHeaderItem> dsItems;
            dsItems.push_back(SAMHeaderItem("READTYPE", ReadType::ToString(readType)));
            dsItems.push_back(SAMHeaderItem("BINDINGKIT", bindingKit));
            dsItems.push_back(SAMHeaderItem("SEQUENCINGKIT", sequencingKit));
            dsItems.push_back(SAMHeaderItem("BASECALLERVERSION", baseCallerVersion));

            // Add QVs, e.g., InsertionQV=iq;DeletionQV=dq;...
            if (samQVs.useqv) {
                for (int i = 0; i < samQVs.nTags; i++) {
                    if (samQVs.useqv & (1 << i)) {
                        dsItems.push_back(SAMHeaderItem(samQVs.qvNames[i], samQVs.qvTags[i]));
                    }
                }
            }
            rgs.Add(SAMHeaderRG(reader->readGroupId, PACBIOPL, movieName, dsItems));
        }
        delete reader;
    } else {
#ifdef USE_PBBAM
        if (fileType == FileType::PBBAM) {
        // TODO: use Derek's API to merge bamHeaders from different files when 
        // it is in place. Use the following code for now. 
        std::vector<std::string>::const_iterator rfit;
        for(rfit = readsFiles.begin(); rfit != readsFiles.end(); rfit++) {
            try {
                PacBio::BAM::BamFile bamFile(*rfit);
                PacBio::BAM::BamHeader header = bamFile.Header();
                // Get read groups from bam header.
                std::vector<PacBio::BAM::ReadGroupInfo> vrgs = header.ReadGroups();
                std::vector<PacBio::BAM::ReadGroupInfo>::iterator rgit;
                for (rgit = vrgs.begin(); rgit != vrgs.end(); rgit++) {
                    rgs.Add(SAMHeaderRG((*rgit).ToSam()));
                }
            } catch (std::exception e) {
                cout << "ERROR, unable to open bam file " << (*rfit) << endl;
                exit(1);
            }
        }
        } else if (fileType == FileType::PBDATASET) {
            for (auto xmlfn: readsFiles) {
                for (auto bamFile: PacBio::BAM::DataSet(xmlfn).BamFiles()) {
                    for (auto rg: bamFile.Header().ReadGroups())
                        rgs.Add(SAMHeaderRG(rg.ToSam()));
                }
            }
        }
#else
        REQUEST_PBBAM_ERROR();
#endif
    }
    return rgs;
}


SAMHeaderPGs SAMHeaderPrinter::MakePGs(const std::vector<std::string> & readsFiles, 
        const std::string & progName, const std::string & progVersion, 
        const std::string & commandLine) {
    SAMHeaderPGs pgs;

    // program id, unique identifier for @PG lines;
    int prog_id = 1;

    if (fileType != FileType::PBBAM) {
        // Reads files are not in BAM format, no other @PG lines.
    } else {
#ifdef USE_PBBAM
        // Reads files are in BAM format, keep all @PG lines from Bam files.
        // TODO: use Derek's API to merge bamHeaders from different files when 
        // it is in place. Use the following code for now. 
        std::vector<std::string>::const_iterator rfit;
        for(rfit = readsFiles.begin(); rfit != readsFiles.end(); rfit++) {
            PacBio::BAM::BamFile bamFile(*rfit);
            PacBio::BAM::BamHeader bamHeader = bamFile.Header();
            std::vector<PacBio::BAM::ProgramInfo> progs = bamHeader.Programs();
            std::vector<PacBio::BAM::ProgramInfo>::iterator it;
            for (it = progs.begin(); it != progs.end(); it++) {
                pgs.Add(SAMHeaderPG((*it).ToSam()));
                prog_id++; // Increase prog_id;
            }
        }
#else
        REQUEST_PBBAM_ERROR();
#endif
    }

    pgs.Add(SAMHeaderPG(std::to_string(prog_id), progName, progVersion, commandLine));
    return pgs;
}

SAMHeaderCOs SAMHeaderPrinter::MakeCOs(const std::vector<std::string> & readsFiles) {
    SAMHeaderCOs cos;

    if (fileType == FileType::PBBAM) {
#ifdef USE_PBBAM
        // Reads files are in BAM format, keep all @CO lines from Bam files.
        std::vector<std::string>::const_iterator rfit;
        for(rfit = readsFiles.begin(); rfit != readsFiles.end(); rfit++) {
            PacBio::BAM::BamFile bamFile(*rfit);
            PacBio::BAM::BamHeader bamHeader = bamFile.Header();
            cos.Append(bamHeader.Comments());
        }
#else
        REQUEST_PBBAM_ERROR();
#endif
    } // reads files are not in BAM format, no comments.
    return cos;
}

#ifdef USE_PBBAM
PacBio::BAM::BamHeader SAMHeaderPrinter::ToBamHeader() {
    std::string headerStr = ToString();
    return PacBio::BAM::BamHeader(headerStr);
}
#endif


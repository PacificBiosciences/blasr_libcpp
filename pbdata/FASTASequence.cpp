#include <stdlib.h>
#include "FASTASequence.hpp"

using namespace std;

FASTASequence::FASTASequence() : DNASequence() {
    title = NULL;
    titleLength = 0;
    deleteTitleOnExit = false; 
    // If deleteTitleOnExist is false, whether to delete title
    // or not should depend on deleteOnExit; otherwise, delete title
    // regardless of deleteOnExit.
}

void FASTASequence::PrintSeq(ostream &out, int lineLength, char delim) const {
    out << delim;
    if (title) out << title;
    out << endl;
    static_cast<const DNASequence*>(this)->PrintSeq(out, lineLength); 
}

int FASTASequence::GetStorageSize() const {
    if (!title) 
        return DNASequence::GetStorageSize();
    return strlen(title) + DNASequence::GetStorageSize();
}

string FASTASequence::GetName() const {
    string name;
    int i;
    for (i = 0; i < titleLength; i++) {
        if (title[i] != ' ' and
                title[i] != '\t' and
                title[i] != '\n' and
                title[i] != '\r') {
            name.push_back(title[i]);
        }
        else {
            break;
        }
    }
    return name;
}

void FASTASequence::ShallowCopy(const FASTASequence &rhs) {
    CheckBeforeCopyOrReference(rhs, "FASTASequence");
    FASTASequence::Free();

    static_cast<DNASequence*>(this)->ShallowCopy(rhs);

    title = rhs.title;
    titleLength = rhs.titleLength;
    deleteTitleOnExit = false;
}

string FASTASequence::GetTitle() const {
    return string(title);
}

// Delete title if this FASTASequence is under control or 
// only title is under control.
void FASTASequence::DeleteTitle() {
    if (deleteOnExit or deleteTitleOnExit) {
        if (title != NULL) {
            delete[] title;
        }
    } // otherwise, title is controlled by another obj
    title = NULL;
    titleLength = 0;
    deleteTitleOnExit = false;
}

void FASTASequence::CopyTitle(const char* str, int strlen) {
    FASTASequence::DeleteTitle();

    // No segfault when str is NULL;
    if (str == NULL) {
        title = NULL;
        titleLength = 0;
    } else {
        title = new char[strlen+1];
        memcpy(title, str, strlen);
        titleLength = strlen;
        title[titleLength] = '\0';
    }

    // In some cases, (e.g., when ReferenceSubstring and CopyTitle
    // are called together), this Sequence may only have control over
    // title but not seq.
    deleteTitleOnExit = true;
}

void FASTASequence::CopyTitle(string str) {
    FASTASequence::CopyTitle(str.c_str(), str.size());
}

void FASTASequence::GetFASTATitle(string& fastaTitle) const {
    // look for the first space, and return the string until there.
    int i;
    for (i = 0; i < titleLength; i++ ){
        if (title[i] == ' ' or
                title[i] == '\t') {
            break;
        }
    }
    fastaTitle.assign(title, i);
}

// Copy rhs.seq[readStart:readEnd] to this Sequence.
void FASTASequence::CopySubsequence(FASTASequence &rhs, int readStart, int readEnd) {
    CheckBeforeCopyOrReference(rhs, "FASTASequence");
    
    // Free before copying anything 
    FASTASequence::Free();

    if (readEnd == -1) {
        readEnd = rhs.length;
    }

    if (readEnd > readStart) {
        length = readEnd - readStart;
        DNASequence::Copy(rhs, readStart, length);
    }
    else {
        seq = NULL;
        length = 0;
        deleteOnExit = true;
    }
    FASTASequence::CopyTitle(rhs.title);
}

void FASTASequence::AppendToTitle(string str) {
    int newLength = titleLength + str.size() + 1;
    if (newLength == 0) {
        DeleteTitle();
        return;
    }

    char *tmpTitle = new char[newLength];
    memcpy(tmpTitle, title, titleLength);
    memcpy(&tmpTitle[titleLength], str.c_str(), str.size());
    tmpTitle[newLength-1] = '\0';
    delete[] title;
    title = tmpTitle;
    titleLength = newLength;
    deleteTitleOnExit = true;
}

void FASTASequence::Assign(FASTASequence &rhs) {
    *this = (FASTASequence&)rhs;
}

// Create a reverse complement FASTASequence of *this and assign to rhs.
void FASTASequence::MakeRC(FASTASequence &rhs, DNALength rhsPos, DNALength rhsLength) {
    rhs.Free();
    DNASequence::MakeRC((DNASequence&) rhs, rhsPos, rhsLength);
    if (title != NULL) {
        (static_cast<FASTASequence*>(&rhs))->CopyTitle(title);
    }
}

void FASTASequence::ReverseComplementSelf() {
    DNALength i;
    for (i = 0; i < length/2 + length % 2; i++) {
        char c = seq[i];
        seq[i] = ReverseComplementNuc[seq[length - i - 1]];
        seq[length - i - 1] = ReverseComplementNuc[static_cast<int>(c)];
    }
}

void FASTASequence::operator=(const FASTASequence &rhs) {
    CheckBeforeCopyOrReference(rhs, "FASTASequence");

    // Free before copying anything 
    FASTASequence::Free();

    // Copy seq from rhs
    ((DNASequence*)this)->Copy((DNASequence&)rhs);

    assert(deleteOnExit);

    // Copy title from rhs
    FASTASequence::CopyTitle(rhs.title, rhs.titleLength);

    assert(deleteOnExit);
}

void FASTASequence::Copy(const std::string & rhsTitle, const std::string & rhsSeq) {
    this->Copy(rhsSeq);
    this->CopyTitle(rhsTitle);
}

void FASTASequence::Copy(const std::string & rhsSeq) {
    (static_cast<DNASequence*>(this))->Copy(rhsSeq);
}

void FASTASequence::Copy(const FASTASequence &rhs) {
    *this = (FASTASequence&)rhs;
}

#ifdef USE_PBBAM
void FASTASequence::Copy(const PacBio::BAM::BamRecord & record) {
    FASTASequence::Copy(record.Impl().Name(), record.Sequence());
}
#endif

void FASTASequence::Free() {
    // Delete title if title is under control, reset deleteTitleOnExit.
    FASTASequence::DeleteTitle();

    // Delete seq if under control, reset deleteOnExit.
    // Don't call Free() before calling DeleteTitle().
    DNASequence::Free();
}

template DNALength ResizeSequence<FASTASequence>(FASTASequence &, DNALength);

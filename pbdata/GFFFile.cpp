#include "GFFFile.hpp"

GFFEntry::GFFEntry(std::string & _name, std::string & _source,
                  std::string & _type, UInt & _start, UInt & _end,
                  float & _score, char & _strand,
                  std::string & _frame, std::string _attributes) {
    name = _name; 
    source = _source;
    type = _type;
    start = _start;
    end = _end;
    score = _score;
    strand = _strand;
    frame = _frame;
    attributes = _attributes;
}

void GFFFile::ReadAll(std::string & gffFileName) {
  std::fstream gffIn;
  CrucialOpen(gffFileName, gffIn, std::ios::in);
  while(gffIn) {
    std::string line;
    getline(gffIn, line);
    std::stringstream linestrm(line);
    std::string name, source, type;
    UInt start, end;
    char strand;
    float score;
    std::string frame, attributes;
    // A sample record in adapterGffFile:
    // ref000001   .   adapter 10955   10999   0.00    +   .   xxxx
    linestrm >> name >> source >> type 
             >> start >> end >> score 
             >> strand >> frame >> attributes;
    entries.push_back(GFFEntry(
      name, source, type, start, end, 
      score, strand, frame, attributes));
    }
    gffIn.close();
}

#ifndef _SIMULATOR_CONTEXT_SET_HPP_
#define _SIMULATOR_CONTEXT_SET_HPP_

#include <iostream>
#include <map>

#include <alignment/simulator/ContextSample.hpp>

typedef std::map<std::string, ContextSample *> T_ContextSampleMap;

class ContextSampleMap : public T_ContextSampleMap
{
public:
    int contextLength;

    ContextSampleMap();
    void Write(std::ofstream &out);
    void Read(std::ifstream &in);
};

#endif

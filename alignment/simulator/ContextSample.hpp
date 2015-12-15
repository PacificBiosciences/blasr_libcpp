#ifndef _SIMULATOR_CONTEXT_SAMPLE_HPP_
#define _SIMULATOR_CONTEXT_SAMPLE_HPP_

#include<string>
#include<vector>
#include<iostream>
#include "simulator/QualitySample.hpp"
#include "statistics/StatUtils.hpp"

class ContextSample {
public:
    std::vector<QualitySample> samples;
    size_t minSamples;
    size_t maxSamples;
    int reachedMinSamples;

    ContextSample();

    size_t GetNSamples();

    int AppendSample(SMRTSequence &seq, DNALength pos);

    QualitySample* GetRandomQualitySample();

    void Write(std::ofstream &out);

    void Read(std::ifstream &in);
};

#endif

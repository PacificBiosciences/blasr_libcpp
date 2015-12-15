#include "ContextSample.hpp"

ContextSample::ContextSample() {
    minSamples = 0;
    maxSamples = 0;
    reachedMinSamples = 0;
}

size_t ContextSample::GetNSamples() {
    return samples.size();
}

int ContextSample::AppendSample(SMRTSequence &seq, DNALength pos) {
    // 
    // Check to see if the space has been sampled enough.
    //
    if (maxSamples != 0 and samples.size() >= maxSamples)
        return 0;

    samples.resize(samples.size()+1);
    samples[samples.size()-1].CreateFromRead(seq,pos);
    if (samples.size() > minSamples and reachedMinSamples == 0) {
        reachedMinSamples = 1;
        return 1;
    }
    else {
        return 0;
    }
}

QualitySample* ContextSample::GetRandomQualitySample() {
    int sampleIndex = RandomInt(samples.size());
    return &samples[sampleIndex];
}

void ContextSample::Write(std::ofstream &out) {
    int sampleSize = samples.size();
    out.write((char*) &sampleSize, sizeof(sampleSize));
    for (size_t s = 0; s < samples.size(); s++ ){
        samples[s].Write(out);
    }
}

void ContextSample::Read(std::ifstream &in) {
    int s;
    int sampleSize;
    in.read((char*)&sampleSize, sizeof(sampleSize));
    samples.resize(sampleSize);
    for (s = 0; s < sampleSize; s++ ){
        samples[s].Read(in);
    }
}

#include "simulator/LengthHistogram.hpp"

int LengthHistogram::Read(std::string &inName) {
    std::ifstream in;
    CrucialOpen(inName, in, std::ios::in);
    return Read(in);
}

int LengthHistogram::Read(std::ifstream &in) {
    while(in) {
        int length, count;
        in >> length;
        in >> count;
        lengthHistogram.data.push_back(length);
        if (lengthHistogram.cdf.size() == 0) {
            lengthHistogram.cdf.push_back(count);
        }
        else {
            lengthHistogram.cdf.push_back(lengthHistogram.cdf[lengthHistogram.cdf.size()-1] + count);
        }
    }
    return 1;
}

void LengthHistogram::GetRandomLength(int &length) {
    lengthHistogram.SelectRandomValue(length);
}

void LengthHistogram::BuildFromAlignmentLengths(std::vector<int> &lengths) {
    sort(lengths.begin(), lengths.end());
    for (size_t f = 0, i = 1; i < lengths.size(); i++) {
        if (lengths[i] != lengths[f]) {
            lengthHistogram.data.push_back(lengths[f]);
            lengthHistogram.cdf.push_back(i);
            f = i;
        }
    }
    if (lengths.size() != 0) {
        lengthHistogram.data.push_back(lengths[lengths.size()-1]);
        lengthHistogram.cdf.push_back(lengths.size());
    }
    /* Tests:
     * indices                0 1 2 3 4  5  6
     * lengths:               1 3 5 9 10 10 11
     * lengthHistogram.data:  1 3 5 9 10 11
     * lengthHistogram.cdf :  1 2 3 4 6  7
     *
     * lengths:               1 3 5 9 10 11
     * lengthHistogram.data:  1 3 5 9 10 11
     * lengthHistogram.cdf :  1 2 3 4 5  6
     *
     * lengths:               10
     * lengthHistogram.data:  10
     * lengthHistogram.cdf :  1 
     */ 
}

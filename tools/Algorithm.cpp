#include "Algorithm.h"

void Algorithm::Permutation(int max, int numSamples, int numTimes, std::vector<std::vector<int>>& samples)
{
    auto rng = std::default_random_engine {};
    std::vector<int> indices(max,0);
    std::iota(indices.begin(), indices.end(),0);

    samples.clear();
    samples.resize(numTimes);

    for (int i=0;i<numTimes;i++)
    {
        std::shuffle(std::begin(indices), std::end(indices), rng);
        std::vector<int> newInd(numSamples);
        for (int j=0;j<numSamples;j++)
        {
            newInd[j] = indices[j];
        }
        samples[i] = newInd;
    }
}

void Algorithm::Permutation(int max, int numSamples, std::vector<int>& samples)
{
    std::default_random_engine rng(time(0));
    std::vector<int> indices(max,0);
    std::iota(indices.begin(), indices.end(),0);

    samples.clear();
    samples.resize(numSamples);

    std::shuffle(std::begin(indices), std::end(indices), rng);
    for (int j=0;j<numSamples;j++)
    {
        samples[j] = indices[j];
    }
}
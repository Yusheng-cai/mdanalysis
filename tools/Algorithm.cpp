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

void Algorithm::find_all_cycles(const std::vector<std::vector<int>>& Neighbor_list, 
int cycle_length, std::vector<std::vector<int>>& all_cycles)
{
    #pragma omp parallel
    {
        std::vector<std::vector<int>> c;
        #pragma omp for 
        for (int i=0;i<Neighbor_list.size();i++){
            std::vector<int> temp_path;
            find_cycle(i, temp_path, Neighbor_list, cycle_length, c);
        }

        #pragma omp critical
        {
            all_cycles.insert(all_cycles.end(), c.begin(), c.end());
        }
    }
}

void Algorithm::find_cycle(int start, std::vector<int>& path, 
                    const std::vector<std::vector<int>>& neighbor_list, 
                    int cycle_length, std::vector<std::vector<int>>& all_cylcles)
{
    path.push_back(start);

    if (path.size() == cycle_length){
        if (Algorithm::contain(neighbor_list[start], path[0])){
            all_cylcles.push_back(path);
        }
        path.pop_back();
        return;
    }

    for (int neighbor : neighbor_list[start]){
        if (! Algorithm::contain(path, neighbor)){
            find_cycle(neighbor, path, neighbor_list, cycle_length, all_cylcles);
        }
    }

    path.pop_back();
}
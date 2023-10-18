#include "IdentifyIceNearSubstrate.hpp"
#include "SimulationState.h"

namespace CalculationRegistry
{
    registry_<IdentifyIceNearSubstrate> registerIINS("IdentifyIceNearSubstrate");
}

IdentifyIceNearSubstrate::IdentifyIceNearSubstrate(const CalculationInput& input)
: Calculation(input)
{
    pack_.ReadString("Substrate", ParameterPack::KeyType::Required, substrate_name_);
    pack_.ReadString("index_file", ParameterPack::KeyType::Required, index_file_);
    pack_.ReadNumber("radius", ParameterPack::KeyType::Optional, radius_);
    pack_.ReadNumber("threshold", ParameterPack::KeyType::Optional, threshold_);
    radius_sq_ = radius_ * radius_;

    // add the atom groups
    addAtomgroup(substrate_name_);

    readFile(index_file_);

    // initialize cell grid
    cell_ = cellptr(new CellGrid(simstate_,radius_));

    // resize neighbor waters
    const auto& substrate_pos = getAtomGroup(substrate_name_);
    neighbor_water_.resize(substrate_pos.getNumAtoms(), 0);
    neighbor_water_iter_.resize(substrate_pos.getNumAtoms(),0);


    // register outputs
    registerOutputFunction("NumberOfNeighborWater", [this](std::string fname) -> void {this -> printNumNeighborWater(fname);});
    registerPerIterOutputFunction("NumberOfNeighborWater", [this](std::ofstream& ofs) -> void {this -> printNumNeighborWaterPerIter(ofs);});
    registerPerIterOutputFunction("CrossThresholdIndex", [this](std::ofstream& ofs) -> void {this -> printAtomsCrossThreshold(ofs);});
}

void IdentifyIceNearSubstrate::calculate(){
    std::vector<Real3> total_pos = simstate_.getTotalAtomPos();
    const auto& substrate = getAtomGroup(substrate_name_);
    const auto& substrate_pos = substrate.getAtomPositions();
    std::vector<std::vector<int>> substrate_cell_indices = cell_->calculateIndices(substrate_pos);
    neighbor_water_iter_.clear();
    neighbor_water_iter_.resize(substrate_pos.size(),0.0);
    int step = simstate_.getFrameNumber();

    // #pragma omp parallel
    // {
    //     std::vector<Real> local_vec(substrate_pos.size(),0);

    //     #pragma omp for
    //     for (int i=0;i<indices_[step].size();i++){
    //         int water_ind = indices_[step][i];

    //         for (int j =0; j<substrate_pos.size();j++){
    //             Real3 dist;
    //             Real distsq;
    //             simstate_.getSimulationBox().calculateDistance(water_pos[water_ind], substrate_pos[j], dist, distsq);
    //             if (distsq <= radius_sq_){
    //                 local_vec[j] += 1;
    //             }
    //         }
    //     }

    //     #pragma omp critical
    //     {
    //         neighbor_water_ = neighbor_water_ + local_vec;
    //     }
    // }

    #pragma omp parallel
    {
        std::vector<Real> local_vec(substrate_pos.size(), 0);

        #pragma omp for
        for (int i=0;i<indices_[step].size();i++){
            int water_ind = indices_[step][i];
            int water_intind = cell_->getCellGridIntIndex(total_pos[water_ind]);
            std::vector<int> neighbor_index = cell_->getNeighborIndex(total_pos[water_ind]);

            for (auto& vec_neigh : substrate_cell_indices){
                for (int ind : vec_neigh){
                    Real3 dist;
                    Real dist_sq;
                    simstate_.getSimulationBox().calculateDistance(total_pos[water_ind], substrate_pos[ind], dist, dist_sq);

                    if (dist_sq <= radius_sq_){
                        local_vec[ind] += 1;
                    }
                }
            }
        }

        #pragma omp critical
        {
            neighbor_water_ = neighbor_water_ + local_vec;
            neighbor_water_iter_ = neighbor_water_iter_ + local_vec;
        }
    }

    cross_threshold_index_.clear();
    const auto& substrate_atoms = substrate.getAtoms();
    for (int i=0;i<substrate_atoms.size();i++){
        if (neighbor_water_iter_[i] >= threshold_){
            cross_threshold_index_.push_back(substrate_atoms[i].index + 1);
        }
    }
}

void IdentifyIceNearSubstrate::printAtomsCrossThreshold(std::ofstream& ofs){
    int time = simstate_.getTime();

    ofs << time << " ";
    for (int i=0;i<cross_threshold_index_.size();i++){
        ofs << cross_threshold_index_[i] << " ";
    }
    ofs << "\n";
}

void IdentifyIceNearSubstrate::update(){
    cell_->update();
}

void IdentifyIceNearSubstrate::finishCalculate(){
    int numFrames = simstate_.getTotalFrames();
    neighbor_water_ = neighbor_water_ /  numFrames;
}

void IdentifyIceNearSubstrate::printNumNeighborWater(std::string fname){
    std::ofstream ofs;
    ofs.open(fname);

    for (int i=0;i<neighbor_water_.size();i++){
        ofs << neighbor_water_[i] << "\n";
    }

    ofs.close();
}

void IdentifyIceNearSubstrate::printNumNeighborWaterPerIter(std::ofstream& ofs){
    for (int i=0;i<neighbor_water_iter_.size();i++){
        ofs << neighbor_water_iter_[i] << " ";
    }
    ofs << "\n";
}

void IdentifyIceNearSubstrate::readFile(const std::string& filename){
    std::ifstream ifs;
    ifs.open(filename);
    ASSERT((ifs.is_open()), "File " << filename << " cannot be opened.");

    std::stringstream ss;
    std::string sentence;
    while(std::getline(ifs, sentence)){
        ss.str(sentence);
        std::vector<int> local_index;
        int index;
        // get rid of the first one
        ss >> index;
        while (ss >> index){
            local_index.push_back(index-1);
        }
        indices_.push_back(local_index);
        ss.clear();
    }
    ifs.close(); 

    ASSERT((simstate_.getTotalFrames() <= indices_.size()), "The index file has less frames than the simulation.");
}

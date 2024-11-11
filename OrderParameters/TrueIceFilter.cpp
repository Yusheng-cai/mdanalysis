#include "TrueIceFilter.hpp"
#include "SimulationState.h"

namespace CalculationRegistry
{
    registry_<TrueIceFilter> registerTrueIceFilter("TrueIceFilter");
}

TrueIceFilter::TrueIceFilter(const CalculationInput& input)
: Calculation(input)
{
    // read in the solid like atoms index
    pack_.ReadString("solid_like_atoms", ParameterPack::KeyType::Optional, filename_);
    StringTools::ParseIndexFile(filename_, ice_like_indices_total_);

    // read atom group
    pack_.ReadString("atomgroup", ParameterPack::KeyType::Required, atomgroup_name_);
    addAtomgroup(atomgroup_name_);

    // register per iter output functions
    registerPerIterOutputFunction("ice_like_atoms", [this](std::ofstream& ofs) -> void {this -> printIceLikeAtoms(ofs);});
    registerPerIterOutputFunction("added_ice", [this](std::ofstream& ofs) -> void {this -> printAddedIce(ofs);});

    // read surface atom group
    pack_.ReadVectorString("surfaceAtomGroups", ParameterPack::KeyType::Optional, surface_atomgroups_);
    for (auto s : surface_atomgroups_){
        addAtomgroup(s);
    }

    // whether we are performing surface correction for number of ice like atoms 
    pack_.ReadNumber("surface_cutoff", ParameterPack::KeyType::Optional, surface_cutoff_);
    pack_.ReadNumber("ice_cutoff", ParameterPack::KeyType::Optional, ice_cutoff_);
    pack_.ReadNumber("surface_threshold", ParameterPack::KeyType::Optional, surface_threshold_);
    pack_.ReadNumber("ice_threshold", ParameterPack::KeyType::Optional, ice_threshold_);
    ice_cutoff_sq_ = ice_cutoff_ * ice_cutoff_;
    surface_cutoff_sq_ = surface_cutoff_ * surface_cutoff_;

    cell_ice_corr_ = cellptr(new CellGrid(simstate_, ice_cutoff_));
    cell_surface_corr_ = cellptr(new CellGrid(simstate_, surface_cutoff_));
}

void TrueIceFilter::calculate(){
    // obtain the atom groups
    const auto& pos = getAtomGroup(atomgroup_name_).getAtomPositions();
    const auto& ag  = getAtomGroup(atomgroup_name_);
    int num_atoms   = ag.getNumAtoms();

    // get the frame number
    int framenum = simstate_.getFrameNumber();

    // here it is in global index
    ice_like_indices_ = ice_like_indices_total_[framenum];
    is_ice_like_.clear();
    is_ice_like_.resize(num_atoms, 0);
    num_ice_like_atoms_before_correction_ = ice_like_indices_.size();

    // convert vector to a map --> water index is in atom group index
    auto map_num = Algorithm::vectorToUnorderedMap(ice_like_indices_total_[framenum]);
    for (int i=0;i<num_atoms;i++){
        int global_index = ag.AtomGroupIndices2GlobalIndices(i) + 1;
        if (! Algorithm::IsInUnorderedMap(map_num, global_index)){
            water_indices_.push_back(i);
        }
    }

    // fill is ice like
    for (int i=0;i<ice_like_indices_.size();i++){
        int index              = ice_like_indices_[i];
        int ag_index           = ag.GlobalIndices2AtomGroupIndices(index-1); 
        is_ice_like_[ag_index] = 1;
    }

    // correct ice like if necessary
    CorrectIceLikeAtomsBasedOnSurface();
}

void TrueIceFilter::CorrectIceLikeAtomsBasedOnSurface(){
    // get all atom positions
    const auto& ag = getAtomGroup(atomgroup_name_);
    const auto& pos= ag.getAtomPositions();
    std::cout << "before correction ice indices = " << ice_like_indices_.size() << "\n";

    int total_atoms = ice_like_indices_.size() + water_indices_.size();
    std::cout << "total_atoms = " << total_atoms << std::endl;
    addedIce_.clear();

    // make cell list for all the water atoms 
    auto ice_cell_corr_list = cell_ice_corr_->calculateIndices(pos);
    std::vector<std::vector<std::vector<int>>> surface_cell_corr_list;
    for (auto s : surface_atomgroups_){
        const auto& posA = getAtomGroup(s).getAtomPositions();
        surface_cell_corr_list.push_back(cell_surface_corr_->calculateIndices(posA));
    }

    int iter=1;
    while (true){ 
        // obtain water indices  --> find how many ice like indices are around it 
        std::vector<int> ice_neighbors(water_indices_.size(),0);
        std::vector<int> surface_neighbors(water_indices_.size(),0);

        #pragma omp parallel for
        for (int i=0;i<water_indices_.size();i++){
            int index = water_indices_[i];

            std::vector<int> ice_neighbor_indices = cell_ice_corr_->getNeighborIndex(pos[index]);
            std::vector<int> surface_neighbor_indices = cell_surface_corr_->getNeighborIndex(pos[index]);

            // check all the neighbor cells 
            int ice_neighbor=0;
            int surface_neighbor=0;

            // no CELL grid --> for debug purposes
            // check the ice 
            // for (int j=0;j<pos.size();j++){
            //     if ((is_ice_like_[j]) && (j != index)){
            //         // check the distance 
            //         Real3 distance;
            //         Real distsq;
            //         simstate_.getSimulationBox().calculateDistance(pos[index], pos[j], distance, distsq);
            //         if (distsq <= ice_cutoff_sq_){
            //             ice_neighbor += 1;
            //         }
            //     }
            // }

            for (int neighbor_cell_ind : ice_neighbor_indices){
                for (int neighbor_ind : ice_cell_corr_list[neighbor_cell_ind]){
                    // only perform calculation if it's ice like
                    if (is_ice_like_[neighbor_ind]){
                        // check the distance 
                        Real3 distance;
                        Real distsq;
                        simstate_.getSimulationBox().calculateDistance(pos[index], pos[neighbor_ind], distance, distsq);
                        if (distsq <= ice_cutoff_sq_){
                            ice_neighbor += 1;
                        }
                    }
                }
            }

            // no CellGrid --> for debug purposes
            // for (int j=0;j<surface_atomgroups_.size();j++){
            //     const auto& surface_pos = getAtomGroup(surface_atomgroups_[j]).getAtomPositions();
            //     for (int k=0;k<surface_pos.size();k++){
            //         // check the distance
            //         Real3 distance;
            //         Real distsq;
            //         simstate_.getSimulationBox().calculateDistance(pos[index], surface_pos[k], distance, distsq);
            //         if (distsq <= surface_cutoff_sq_){
            //             surface_neighbor += 1;
            //         }
            //     }
            // }


            for (int j=0;j<surface_atomgroups_.size();j++){
                const auto& surface_pos = getAtomGroup(surface_atomgroups_[j]).getAtomPositions();
                for (int neighbor_cell_ind : surface_neighbor_indices){
                    for (int neighbor_ind : surface_cell_corr_list[j][neighbor_cell_ind]){
                        // check the distance
                        Real3 distance;
                        Real distsq;
                        simstate_.getSimulationBox().calculateDistance(pos[index], surface_pos[neighbor_ind], distance, distsq);
                        if (distsq <= surface_cutoff_sq_){
                            surface_neighbor += 1;
                        }
                    }
                }
            }

            surface_neighbors[i] = surface_neighbor;
            ice_neighbors[i] = ice_neighbor;
        }

        // set up the new water indices --> iterate over existing water indices as no ice should turn into water 
        std::vector<int> new_water_indices;
        for (int i=0;i<water_indices_.size();i++){
            int index = water_indices_[i];
            // if the water atom meets the criteria
            if ((surface_neighbors[i] >= surface_threshold_) && (ice_neighbors[i] >= ice_threshold_)){
                ASSERT((is_ice_like_[index] == 0), "There must be something wrong.");
                // change this atom to be ice like
                is_ice_like_[index] = 1;
                ice_like_indices_.push_back(ag.AtomGroupIndices2GlobalIndices(index));
                addedIce_.push_back(ag.AtomGroupIndices2GlobalIndices(index));
            }
            else{
                new_water_indices.push_back(index);
            }
        }

        iter++;

        // break if no new water indices has been added
        if (water_indices_.size() == new_water_indices.size()){
            break;
        }

        water_indices_.clear(); 
        water_indices_.insert(water_indices_.end(), new_water_indices.begin(), new_water_indices.end());
    }

    ASSERT((water_indices_.size() + ice_like_indices_.size() == total_atoms), "After correction, the total number of atoms changed.");
    std::cout << "After correction = " << ice_like_indices_.size() << "\n";
}

void TrueIceFilter::printIceLikeAtoms(std::ofstream& ofs){
    int time = simstate_.getTime();
    ofs << time << " ";

    const auto& ag = getAtomGroup(atomgroup_name_);
    for (int i=0;i<ice_like_indices_.size();i++){
        ofs << ice_like_indices_[i] + 1 << " ";
    }
    ofs << "\n";
}

void TrueIceFilter::printAddedIce(std::ofstream& ofs){
    int time = simstate_.getTime();
    ofs << time << " ";

    const auto& ag = getAtomGroup(atomgroup_name_);
    for (int i=0;i<addedIce_.size();i++){
        ofs << addedIce_[i] + 1 << " ";
    }
    ofs << "\n";
}

void TrueIceFilter::finishCalculate()
{

}

void TrueIceFilter::update()
{
    cell_ice_corr_->update();
    cell_surface_corr_->update();

    // clear water and ice_like_indices
    water_indices_.clear();
    ice_like_indices_.clear();
}
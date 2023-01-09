#include "IdentifySolventNearSubstrate.hpp"

namespace CalculationRegistry
{
    registry_<IdentifySolventNearSubstrate> registerISNS("IdentifySolventNearSubstrate");
}

IdentifySolventNearSubstrate::IdentifySolventNearSubstrate(const CalculationInput& input)
: Calculation(input)
{
    pack_.ReadString("solvent_name", ParameterPack::KeyType::Required, solvent_name_);
    pack_.ReadString("substrate_name", ParameterPack::KeyType::Required, substrate_name_);
    pack_.ReadNumber("radius", ParameterPack::KeyType::Required, radius_);
    radius_sq_ = radius_ * radius_;

    // add the atom groups
    addAtomgroup(solvent_name_);
    addAtomgroup(substrate_name_);
    num_substrate_atoms_ = getAtomGroup(substrate_name_).getNumAtoms();
    num_atoms_in_sphere_.resize(num_substrate_atoms_,0);
    PerIter_num_atoms_in_sphere_.resize(num_substrate_atoms_,0);

    // initialize the cell grid
    cell_ = cellptr(new CellGrid(simstate_, radius_));

    // register output functions
    registerPerIterOutputFunction("AtomsInSphere", [this](std::ofstream& ofs) -> void {this -> printPerIterAtomsInSphere(ofs);});
}


void IdentifySolventNearSubstrate::calculate()
{
    const auto& solvent_pos = getAtomGroup(solvent_name_).getAtomPositions();
    const auto& substrate_pos = getAtomGroup(substrate_name_).getAtomPositions();

    std::vector<std::vector<int>> solvent_indices = cell_->calculateIndices(solvent_pos);
    std::fill(PerIter_num_atoms_in_sphere_.begin(), PerIter_num_atoms_in_sphere_.end(),0.0);


    for (int i=0;i<substrate_pos.size();i++){
        std::vector<int> neighbor_index = cell_->getNeighborIndex(substrate_pos[i]);
        for (auto neigh : neighbor_index){
            for (auto ind : solvent_indices[neigh]){
                Real distsq;
                Real3 dist;
                simstate_.getSimulationBox().calculateDistance(substrate_pos[i], solvent_pos[ind], dist, distsq);

                if (distsq <= radius_sq_){
                    PerIter_num_atoms_in_sphere_[i] += 1;
                    num_atoms_in_sphere_[i] += 1;
                }
            }
        }
    }
}

void IdentifySolventNearSubstrate::printPerIterAtomsInSphere(std::ofstream& ofs){
    int frameNumber = simstate_.getFrameNumber();

    for (int i=0;i<PerIter_num_atoms_in_sphere_.size();i++){
        ofs << PerIter_num_atoms_in_sphere_[i] << " ";
    }
    ofs << "\n";
}

void IdentifySolventNearSubstrate::update(){

}

void IdentifySolventNearSubstrate::finishCalculate(){
    int numTotalFrames = simstate_.getTotalFrames();

    num_atoms_in_sphere_ = num_atoms_in_sphere_ / numTotalFrames;
}
#include "Dimer.h"
#include "SimulationState.h"

namespace CalculationRegistry
{
    registry_<Dimer> registerDimer("Dimer");
}

Dimer::Dimer(const CalculationInput& input) : Calculation(input)
{
    pack_.ReadString("residue", ParameterPack::KeyType::Required, resname_);
    initializeResidueGroup(resname_);

    pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional, headindex_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, tailindex_);
    headindex_--;
    tailindex_--;

    pack_.ReadNumber("alignment_cutoff", ParameterPack::KeyType::Optional, alignment_cutoff_);
    pack_.ReadNumber("distance_cutoff", ParameterPack::KeyType::Optional, distance_cutoff_);
    pack_.ReadNumber("distance_cutoff_B1B2", ParameterPack::KeyType::Optional, distance_cutoff_B1B2_);
    pack_.ReadNumber("numtbins", ParameterPack::KeyType::Optional, numtbins_);

    orientation_dimer_.resize(numtbins_,0);
    orientation_monomer_.resize(numtbins_,0);

    // initialize bin
    bin_ = Binptr(new Bin({{-1,1}}, numtbins_));

    initializeProbeVolumes();
    initializeNotInProbeVolumes();

    registerPerIterOutputFunction("DimerIndices", [this](std::ofstream& ofs)->void{this->printDimerIndices(ofs);});
    registerPerIterOutputFunction("nonDimerIndices", [this](std::ofstream& ofs)->void{this->printNonDimerIndices(ofs);});
    registerOutputFunction("OrientationDistribution", [this](std::string name)->void{this->printOrientation(name);});
    registerOutputFileOutputs("numDimers", [this]()->Real{return this->getNumDimer();});

    COMIndicesB1_ = COMIndicesB1_-1;
    COMIndicesB2_ = COMIndicesB2_-1;
}

void Dimer::update(){
    dimer_vec_.clear();
    non_dimer_vec_.clear();
}

void Dimer::printDimerIndices(std::ofstream& ofs){
    const auto& resgroup = getResidueGroup(resname_).getResidues();
    int numStep = simstate_.getStep();
    ofs << numStep << " ";

    for (int i=0;i<dimer_vec_.size();i++){
        int ind = dimer_vec_[i];
        for (int j=0;j<resgroup[ind].atoms_.size();j++){
            ofs << resgroup[ind].atoms_[j].atomNumber_ << " ";
        }
    }
    ofs << "\n";
}

void Dimer::printNonDimerIndices(std::ofstream& ofs){
    const auto& resgroup = getResidueGroup(resname_).getResidues();
    int numStep = simstate_.getStep();
    ofs << numStep << " ";

    for (int i=0;i<non_dimer_vec_.size();i++){
        int ind = non_dimer_vec_[i];
        for (int j=0;j<resgroup[ind].atoms_.size();j++){
            ofs << resgroup[ind].atoms_[j].atomNumber_  << " ";
        }
    }
    ofs << "\n";
}

void Dimer::calculate(){
    const auto& resgroup = getResidueGroup(resname_).getResidues();
    COM_.clear();
    COM_.resize(resgroup.size());
    COMB1_.clear();
    COMB1_.resize(resgroup.size());
    COMB2_.clear();
    COMB2_.resize(resgroup.size());
    uij_.clear();
    uij_.resize(resgroup.size());
    num_dimers_ = 0;
    std::vector<int> is_in_PV;

    // calculate which COM are inside the probe volume
    for (int i=0;i<resgroup.size();i++){
        COM_[i] = calcCOM(resgroup[i]);
        COMB1_[i] = calcCOM(resgroup[i], COMIndicesB1_);
        COMB2_[i] = calcCOM(resgroup[i], COMIndicesB2_);

        if (isInPV(COM_[i])){
            is_in_PV.push_back(i);
        }

        // calculate the uij 
        Real3 headatom = resgroup[i].atoms_[headindex_].positions_;
        Real3 tailatom = resgroup[i].atoms_[tailindex_].positions_;
        Real3 vec_dist;
        Real dist_sq;

        simstate_.getSimulationBox().calculateDistance(headatom, tailatom, vec_dist, dist_sq);

        // normalize the vector
        LinAlg3x3::normalize(vec_dist);

        // store the vector 
        uij_[i] = vec_dist;
    }

    std::map<int, bool> mapIndicesToDone;

    // calculate the pair of residues inside the probe volume 
    #pragma omp parallel
    {
        std::vector<int> local_dimer_vec;

        #pragma omp for
        for (int i=0;i<is_in_PV.size();i++){
            // get the com
            Real3 com = COM_[is_in_PV[i]];
            Real3 comb1 = COMB1_[is_in_PV[i]];
            Real3 comb2 = COMB2_[is_in_PV[i]];
            bool not_dimer=true;

            // iterate over the rest of the COM
            for (int j=0;j<is_in_PV.size();j++){
                if (j != i){
                    Real3 com2 = COM_[is_in_PV[j]];
                    Real3 comb12 = COMB1_[is_in_PV[j]];
                    Real3 comb22 = COMB2_[is_in_PV[j]];

                    Real distsq, dist, distsqb1b2, distb1b2, distsqb2b1, distb2b1; 
                    Real3 vec_dist;

                    // calculate the distance between CN's 
                    simstate_.getSimulationBox().calculateDistance(com , com2, vec_dist, distsq);
                    simstate_.getSimulationBox().calculateDistance(comb1, comb22, vec_dist, distsqb1b2);
                    simstate_.getSimulationBox().calculateDistance(comb2, comb12, vec_dist, distsqb2b1);
                    dist = std::sqrt(distsq);
                    distb2b1 = std::sqrt(distsqb2b1);
                    distb1b2 = std::sqrt(distsqb1b2);

                    // dot product 
                    Real dot = LinAlg3x3::DotProduct(uij_[is_in_PV[i]], uij_[is_in_PV[j]]);

                    if ((dot <= alignment_cutoff_ && dist <= distance_cutoff_ ) \
                     || (dot<= alignment_cutoff_ && distb1b2 <= distance_cutoff_B1B2_) \
                     || (dot <= alignment_cutoff_ && distb2b1 <= distance_cutoff_B1B2_)){
                        local_dimer_vec.push_back(is_in_PV[i]);
                        local_dimer_vec.push_back(is_in_PV[j]);
                        bool a = true;
                        num_dimers_ += 1;
                        not_dimer = false;
                    }
                }
            }
        }

        #pragma omp critical
        {
            dimer_vec_.insert(dimer_vec_.end(), local_dimer_vec.begin(), local_dimer_vec.end());
        }
    }

    Algorithm::unique(dimer_vec_);
    for (int i : dimer_vec_){
        mapIndicesToDone.insert(std::make_pair(i,true));
    }
    

    for (int i=0;i<is_in_PV.size();i++){
        bool a;
        if (! Algorithm::IsInMap(mapIndicesToDone, is_in_PV[i],a)){
            non_dimer_vec_.push_back(is_in_PV[i]);
        }
    }

    Algorithm::unique(dimer_vec_);
    Algorithm::unique(non_dimer_vec_);

    for (int i=0;i<dimer_vec_.size();i++){
        int index = dimer_vec_[i];
        Real dot  = uij_[index][2];
        if (bin_->isInRange(dot)){
            int ind = bin_->findBin(dot);
            orientation_dimer_[ind] += 1;
        }
    }

    for (int i=0;i<non_dimer_vec_.size();i++){
        int index = non_dimer_vec_[i];
        Real dot  = uij_[index][2];
        if (bin_->isInRange(dot)){
            int ind = bin_->findBin(dot);
            orientation_monomer_[ind] += 1;
        }
    }

    std::cout << "Number of dimers = " << dimer_vec_.size() << " Total atoms = " << is_in_PV.size() << "\n";
    std::cout << "Non Dimer = " << non_dimer_vec_.size() << "\n";
}

void Dimer::printOrientation(std::string name){
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<numtbins_;i++){
        ofs << bin_->getCenterLocationOfBin(i) << " " << orientation_dimer_[i] << " " << orientation_monomer_[i] << "\n";
    }

    ofs.close();
}

void Dimer::finishCalculate(){
    int numFrames = simstate_.getTotalFrames();
    orientation_dimer_ = orientation_dimer_ / numFrames;
    orientation_monomer_ = orientation_monomer_ / numFrames;
}
#include "SlabDimer.hpp"
#include "SimulationState.h"

namespace CalculationRegistry
{
    registry_<SlabDimer> registerSlabDimer("SlabDimer");
}

SlabDimer::SlabDimer(const CalculationInput& input):
Calculation(input)
{
    // initialize the bins first --> zbin and thetabin
    auto zbinPack = input.pack_.findParamPack("zbin", ParameterPack::KeyType::Optional);
    pack_.ReadString("direction", ParameterPack::KeyType::Optional, direction_);

    if(zbinPack != nullptr){
        zBin_    = Binptr(new Bin(*zbinPack));
        numzbins_= zBin_->getNumbins();
    }
    else
    {
        pack_.ReadNumber("numzbins", ParameterPack::KeyType::Required, numzbins_);
        pack_.ReadNumber("above", ParameterPack::KeyType::Optional, above_);
        pack_.ReadNumber("below", ParameterPack::KeyType::Optional, below_);
        zBin_    = Binptr(new Bin());
        usingMinMax_ = true;
    }
    BinLocation_.resize(numzbins_,0);

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

    tBin_ = Binptr(new Bin({{-1,1}}, numtbins_));

    // initialize direction
    auto it = MapNameToDirection.find(direction_);
    ASSERT((it != MapNameToDirection.end()), "The direction " << direction_ << " is not recognized.");
    index_ = it -> second;

    // get the COM indices 
    COMIndicesB1_ = COMIndicesB1_-1;
    COMIndicesB2_ = COMIndicesB2_-1;
    num_dimers_.clear();
    num_dimers_.resize(numzbins_,0);
    ratio_dimer_.clear();
    ratio_dimer_.resize(numzbins_,0);
    orientation_dimer_.resize(numzbins_, std::vector<Real>(numtbins_,0.0));
    orientation_monomer_.resize(numzbins_, std::vector<Real>(numtbins_,0.0));

    registerOutputFunction("ratio_dimer", [this](std::string name)->void{this -> printDimerRatio(name);});
    registerOutputFunction("dimer_orientation", [this](std::string name)->void{this->printDimerOrientation(name);});
    registerOutputFunction("monomer_orientation", [this](std::string name)->void{this->printMonomerOrientation(name);});
}


void SlabDimer::binUsingMinMax()
{
    Real slight_shift=1e-3;

    std::vector<Real> zdir;
    for (int i=0;i<COM_.size();i++){
        if ((COM_[i][index_] > above_) && (COM_[i][index_] < below_)){
            zdir.push_back(COM_[i][index_]);
        }
    }

    Real min   = Algorithm::min(zdir) - slight_shift;
    Real max   = Algorithm::max(zdir) + slight_shift;

    Range range = {{min, max}};

    zBin_->update(range, numzbins_);

    // sum over the bin location 
    for (int i=0;i<numzbins_;i++){
        BinLocation_[i] += zBin_ -> getCenterLocationOfBin(i);
    }
    std::cout << "Min = " << min << ", Max = " << max << "\n";
}



void SlabDimer::calculate(){
    const auto& res = getResidueGroup(resname_).getResidues();
    COM_.clear();
    COM_.resize(res.size());
    COMB1_.clear();
    COMB1_.resize(res.size());
    COMB2_.clear();
    COMB2_.resize(res.size());
    uij_.clear();
    uij_.resize(res.size());
    BinIndices_.clear();
    BinIndices_.resize(numzbins_);
    DimerIndices_.clear();
    DimerIndices_.resize(numzbins_, {});
    MonomerIndices_.clear();
    MonomerIndices_.resize(numzbins_,{});
    numCOM_.clear();
    numCOM_.resize(numzbins_,0);

    // obtain the center of mass of each of the residues
    #pragma omp parallel for 
    for (int i=0;i<res.size();i++){
        COM_[i] = calcCOM(res[i]);
        COMB1_[i] = calcCOM(res[i], COMIndicesB1_);
        COMB2_[i] = calcCOM(res[i], COMIndicesB2_);

        // calculate the uij 
        Real3 headatom = res[i].atoms_[headindex_].positions_;
        Real3 tailatom = res[i].atoms_[tailindex_].positions_;
        Real3 vec_dist;
        Real dist_sq;

        simstate_.getSimulationBox().calculateDistance(headatom, tailatom, vec_dist, dist_sq);

        // normalize the vector
        LinAlg3x3::normalize(vec_dist);

        // store the vector 
        uij_[i] = vec_dist;
    }

    if (usingMinMax_){
        binUsingMinMax();
    }

    // let's then check which bin each COM falls into 
    for (int i=0;i<COM_.size();i++){
        if (zBin_->isInRange(COM_[i][index_])){
            int binnum = zBin_->findBin(COM_[i][index_]);
            BinIndices_[binnum].push_back(i);
            numCOM_[binnum] += 1;
        }
    }

    // for each bin, perform the dimer calculation
    // calculate the pair of residues inside the probe volume 
    for (int i=0;i<numzbins_;i++){
        #pragma omp parallel
        {
            std::vector<int> local_dimer_vec;
            std::vector<int> local_monomer_vec;

            #pragma omp for
            for (int j=0;j<BinIndices_[i].size();j++){
                // get the com
                int ind = BinIndices_[i][j];
                Real3 com = COM_[ind];
                Real3 comb1 = COMB1_[ind];
                Real3 comb2 = COMB2_[ind];
                bool isDimer = false;

                // iterate over the rest of the COM
                for (int k=0;k<res.size();k++){
                    int ind2 = k;
                    if (ind2 != ind){
                        Real3 com2 = COM_[ind2];
                        Real3 comb12 = COMB1_[ind2];
                        Real3 comb22 = COMB2_[ind2];

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
                        Real dot = LinAlg3x3::DotProduct(uij_[ind], \
                                                        uij_[ind2]);

                        if ((dot <= alignment_cutoff_ && dist <= distance_cutoff_ ) \
                        || (dot<= alignment_cutoff_ && distb1b2 <= distance_cutoff_B1B2_) \
                        || (dot <= alignment_cutoff_ && distb2b1 <= distance_cutoff_B1B2_)){
                            local_dimer_vec.push_back(ind);
                            isDimer = true;
                            break;
                        }
                    }
                }
                if (! isDimer){
                    local_monomer_vec.push_back(ind);
                }
            }

            #pragma omp critical
            {
                DimerIndices_[i].insert(DimerIndices_[i].end(),\
                                        local_dimer_vec.begin(), \
                                        local_dimer_vec.end());
                MonomerIndices_[i].insert(MonomerIndices_[i].end(),\
                                        local_monomer_vec.begin(), \
                                        local_monomer_vec.end());
            }
        }
        Algorithm::unique(DimerIndices_[i]);
        Algorithm::unique(MonomerIndices_[i]);
        num_dimers_[i] = DimerIndices_[i].size();
    }

    for (int i=0;i<numzbins_;i++){
        ratio_dimer_[i] += (Real)num_dimers_[i] /  (Real)numCOM_[i];
    }

    // get the orientational information for dimers and monomers
    for (int i=0;i<numzbins_;i++){
        for (int j=0;j<DimerIndices_[i].size();j++){
            int ind = DimerIndices_[i][j];
            if (tBin_->isInRange(uij_[ind][2])){
                int binNum = tBin_->findBin(uij_[ind][2]);
                orientation_dimer_[i][binNum] += 1;
            }
        }

        for (int j=0;j<MonomerIndices_[i].size();j++){
            int ind = MonomerIndices_[i][j];
            if (tBin_->isInRange(uij_[ind][2])){
                int binNum = tBin_->findBin(uij_[ind][2]);
                orientation_monomer_[i][binNum] += 1;
            }
        }
    }
}

void SlabDimer::printDimerRatio(std::string name){
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<ratio_dimer_.size();i++){
        ofs << BinLocation_[i] << " " << ratio_dimer_[i] << "\n";
    }

    ofs.close();
}

void SlabDimer::printDimerOrientation(std::string name){
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<numzbins_;i++){
        for (int j=0;j<numtbins_;j++){
            ofs << i << " " << j << " " << orientation_dimer_[i][j] << "\n";
        }
    }

    ofs.close();
}

void SlabDimer::printMonomerOrientation(std::string name){
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<numzbins_;i++){
        for (int j=0;j<numtbins_;j++){
            ofs << i << " " << j << " " << orientation_monomer_[i][j] << "\n";
        }
    }


    ofs.close();
}

void SlabDimer::finishCalculate(){
    int numFrames = simstate_.getTotalFrames();

    BinLocation_ = BinLocation_ / numFrames;
    ratio_dimer_ = ratio_dimer_ / numFrames;

    orientation_dimer_ = orientation_dimer_ / numFrames;
    orientation_monomer_ = orientation_monomer_ / numFrames;
}

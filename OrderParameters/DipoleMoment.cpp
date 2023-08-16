#include "DipoleMoment.hpp"

namespace CalculationRegistry
{
    registry_<DipoleMoment> registerDipoleMoment("DipoleMoment");
}

DipoleMoment::DipoleMoment(const CalculationInput& input)
: Calculation(input)
{
    // read in the residue
    pack_.ReadString("residue", ParameterPack::KeyType::Required, resname_);
    pack_.ReadNumber("numbins", ParameterPack::KeyType::Optional, numbins_);
    bool read_dipole_index =  \
            pack_.ReadVectorNumber("dipole_index", ParameterPack::KeyType::Optional, dipole_index_);

    // initialize bin
    bin_ = Binptr(new Bin(range_, numbins_));
    azi_bin_ = Binptr(new Bin({{-Constants::PI, Constants::PI}}, numbins_));
    azimuthal_angle_dist_.resize(numbins_,0.0);
    zenithal_angle_dist_.resize(numbins_,0.0);

    initializeResidueGroup(resname_);

    if (! read_dipole_index){
        auto res = getResidueGroup(resname_).getResidues()[0];
        int numatoms = res.atoms_.size();
        dipole_index_.clear();
        dipole_index_.resize(numatoms);
        std::iota(dipole_index_.begin(), dipole_index_.end(), 0);
    }

    registerOutputFunction("distribution", \
    [this](std::string name) -> void {this -> printAngleDistribution(name);});
};

void DipoleMoment::calculate(){
    const auto& res = getResidueGroup(resname_).getResidues();

    mu_.clear();
    mu_.resize(res.size());

    azimuthal_angle_.clear();
    azimuthal_angle_.resize(res.size());

    zenithal_angle_.clear();
    zenithal_angle_.resize(res.size());

    for (int i=0;i<res.size();i++){
        auto& r = res[i];

        Real3 mu_i = {{0,0,0}};
        mu_[i] =  CalculationTools::CalculateResidueDipole(r, simstate_,dipole_index_);

        // Real3 dist;
        // Real distsq;
        // simstate_.getSimulationBox().calculateDistance(r.atoms_[0].positions_, r.atoms_[1].positions_,\
        //                                                         dist, distsq);
        // LinAlg3x3::normalize(dist);
        // std::cout << "CN vec = " << dist << "\n";
        

        // calculate the azimuthal angle
        azimuthal_angle_[i] = std::atan2(mu_[i][1], mu_[i][0]);

        if (azi_bin_->isInRange(azimuthal_angle_[i])){
            int numBin = azi_bin_->findBin(azimuthal_angle_[i]);
            azimuthal_angle_dist_[numBin] += 1;
        }

        // calculate the zenithal angle 
        zenithal_angle_[i] = mu_[i][2];

        if (bin_->isInRange(zenithal_angle_[i])){
            int numBin  = bin_->findBin(zenithal_angle_[i]);
            zenithal_angle_dist_[numBin] += 1;
        }
    }
}

void DipoleMoment::finishCalculate(){
    int num_totalframes = simstate_.getTotalFrames();

    for (int i=0;i<numbins_;i++){
        azimuthal_angle_dist_[i] = azimuthal_angle_dist_[i] / num_totalframes;
        zenithal_angle_dist_[i]  = zenithal_angle_dist_[i] / num_totalframes;
    }
}

void DipoleMoment::printAngleDistribution(std::string name){
    std::ofstream ofs;
    ofs.open(name);


    ofs << "# azi_angle(radian) azimuthal_distribution  zen_angle(cos) zenithal_distribution\n";

    for (int i=0;i<numbins_;i++){
        ofs << azi_bin_->getCenterLocationOfBin(i) <<  " " << azimuthal_angle_dist_[i] << " " \
        << bin_->getCenterLocationOfBin(i) << " " << zenithal_angle_dist_[i] << "\n";
    }

    ofs.close();

}

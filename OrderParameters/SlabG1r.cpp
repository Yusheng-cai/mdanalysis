#include "SlabG1r.h"

namespace CalculationRegistry
{
    registry_<SlabG1r> registerSlabG1r("SlabG1r");
}

SlabG1r::SlabG1r(const CalculationInput& input)
: Calculation(input)
{
    // read in information about the residue
    pack_.ReadString("residue", ParameterPack::KeyType::Required, residue_);
    initializeResidueGroup(residue_);

    // read in the distance COM 
    pack_.ReadVectorNumber("distanceCOM", ParameterPack::KeyType::Required, distanceCOM_);

    // initialize z bins
    pack_.ReadNumber("numzbins", ParameterPack::KeyType::Required, numzbins_);
    pack_.ReadNumber("above", ParameterPack::KeyType::Optional, above_);
    zbin_ = binptr(new Bin());
    IndicesZbin_.resize(numzbins_);
    zBinLocation_.resize(numzbins_);

    // initialize t bins
    pack_.ReadNumber("numtbins", ParameterPack::KeyType::Optional, numtbins_);
    tbin_ = binptr(new Bin({{-1,1}}, numtbins_));

    // initialize r bin
    auto rBinPack = pack_.findParamPack("rbin", ParameterPack::KeyType::Required);
    rbin_ = binptr(new Bin(*rBinPack));
    numrbins_ = rbin_->getNumbins();

    // resize g1r
    G1r_.resize(numzbins_, std::vector<Real>(numrbins_, 0.0));
    NumberResidueG1r_.resize(numzbins_, std::vector<Real>(numrbins_,0.0));
    G1r_3d_.resize(numzbins_, std::vector<std::vector<Real>>(numrbins_, std::vector<Real>(numtbins_, 0.0)));

    // head tail business
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional, headindex_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, tailindex_);
    headindex_--;
    tailindex_--;

    // perform g1(r) calculation within each z bin 
    pack_.Readbool("within_zbin", ParameterPack::KeyType::Optional, g1r_within_zbin_);

    // calculate the volume
    volume_.resize(rbin_->getNumbins(),0.0);
    Real dr = rbin_->getStep();
    Real dz = zbin_->getStep();
    if (! g1r_within_zbin_){
        volume_[0] = 4.0/3.0*Constants::PI*std::pow(dr,3.0);
        for (int i=1;i<numrbins_;i++){
            volume_[i] = 4 * Constants::PI * std::pow(rbin_ ->getLeftLocationOfBin(i),2) * dr;
        }
    }
    else{
        Real dz = zbin_->getStep();
        // pi * r^2 * L
        volume_[0] = Constants::PI * std::pow(dr, 2.0) * dz;

        for (int i=0;i<numrbins_;i++){
            volume_[i] = 2 * Constants::PI * rbin_->getLeftLocationOfBin(i) * dz * dr;
        }
    }

    // register outputs
    registerOutputFunction("SlabG1r", [this](std::string name) -> void { this -> printSlabG1r(name);});
    registerOutputFunction("SlabG1r_3d", [this](std::string name) -> void {this -> printSlabG1r3d(name);});
}

void SlabG1r::binUsingMinMax()
{
    Real slight_shift=1e-3;

    std::vector<Real> zdir;
    for (int i=0;i<COM_.size();i++){
        if (COM_[i][index_] > above_){
            zdir.push_back(COM_[i][index_]);
        }
    }

    auto minit = std::min_element(zdir.begin(), zdir.end());
    auto maxit = std::max_element(zdir.begin(), zdir.end());

    Real min   = *minit - slight_shift;
    Real max   = *maxit + slight_shift;

    Range range = {{min, max}};

    zbin_->update(range, numzbins_);

    // sum over the bin location 
    for (int i=0;i<numzbins_;i++)
    {
        zBinLocation_[i] += zbin_ -> getCenterLocationOfBin(i);
    }
    std::cout << "Min = " << min << ", Max = " << max << "\n";

    if (g1r_within_zbin_){
        Real dz = zbin_->getStep();
        Real dr = rbin_->getStep();
        // pi * r^2 * L
        volume_[0] = Constants::PI * std::pow(dr, 2.0) * dz;

        for (int i=1;i<numrbins_;i++){
            volume_[i] = 2 * Constants::PI * rbin_->getLeftLocationOfBin(i) * dz * dr;
        }
    }
}


void SlabG1r::calculate()
{
    auto& res = getResidueGroup(residue_).getResidues();
    uij_.clear();
    uij_.resize(res.size(),{});
    std::vector<std::vector<int>> IndicesPerZbin(numzbins_);

    std::vector<std::vector<std::vector<Real>>> G1r_3d_perIter(numzbins_, std::vector<std::vector<Real>>(numrbins_, std::vector<Real>(numtbins_,0.0)));
    std::vector<Real> InsideNum(numzbins_);

    // calculate COM
    #pragma omp parallel for 
    for (int i=0;i<res.size();i++)
    {
        COM_[i] = calcCOM(res[i]);

        Real3 director;
        Real sq_dist;
        auto headpos = res[i].atoms_[headindex_].positions_;
        auto tailpos = res[i].atoms_[tailindex_].positions_;

        simstate_.getSimulationBox().calculateDistance(headpos, tailpos, director, sq_dist);

        LinAlg3x3::normalize(director);

        uij_[i] = director;
    }

    // bin using the min max of COM 
    binUsingMinMax();

    // find which residues fall into each of the z bins 
    for (int i=0;i<res.size();i++){
        if (zbin_->isInRange(COM_[i][index_])){
            int index = zbin_->findBin(COM_[i][index_]);
            IndicesPerZbin[index].push_back(i);
        }
    }

    for (int i=0;i<numzbins_;i++){
        InsideNum[i] = IndicesPerZbin[i].size();
    }

    // split calculation depending on if we are calculating g1(r) within each z bin or not 
    if (! g1r_within_zbin_){
        #pragma omp parallel
        {
            std::vector<std::vector<Real>> G1r_local_(numzbins_, std::vector<Real>(numrbins_,0.0));
            std::vector<std::vector<Real>> Numres_local_(numzbins_, std::vector<Real>(numrbins_,0.0));
            std::vector<std::vector<std::vector<Real>>> G1r_3d_local(numzbins_, std::vector<std::vector<Real>>(numrbins_, std::vector<Real>(numtbins_,0.0)));

            #pragma omp for
            for (int i=0;i<res.size();i++){
                for (int j=0;j<numzbins_;j++){
                    auto& indices = IndicesPerZbin[j];
                    for (int ind : indices){
                        if (ind != i){
                            // calculate distance squared
                            Real dist_sq;
                            Real3 dist;
                            simstate_.getSimulationBox().calculateDistance(COM_[ind], COM_[i], dist, dist_sq);
                            Real r = std::sqrt(dist_sq);

                            if (rbin_->isInRange(r)){
                                int index = rbin_->findBin(r);
                                Real dotproduct = LinAlg3x3::DotProduct(uij_[ind], uij_[i]);
                                G1r_local_[j][index] += dotproduct;
                                Numres_local_[j][index] += 1;

                                if (tbin_->isInRange(dotproduct)){
                                    int tind = tbin_->findBin(dotproduct);
                                    G1r_3d_local[j][index][tind] += 1;
                                }
                            }
                        }
                    }
                }
            }

            #pragma omp critical
            for (int i=0;i<numzbins_;i++){
                for (int j=0;j<numrbins_;j++){
                    NumberResidueG1r_[i][j] += Numres_local_[i][j];
                    G1r_[i][j] += G1r_local_[i][j];
                    for (int k=0;k<numtbins_;k++){
                        G1r_3d_perIter[i][j][k] += G1r_3d_local[i][j][k];
                    }
                }
            }
        }
    }
    else{
        #pragma omp parallel
        {
            std::vector<std::vector<Real>> G1r_local_(numzbins_, std::vector<Real>(numrbins_,0.0));
            std::vector<std::vector<Real>> Numres_local_(numzbins_, std::vector<Real>(numrbins_,0.0));
            std::vector<std::vector<std::vector<Real>>> G1r_3d_local(numzbins_, std::vector<std::vector<Real>>(numrbins_, std::vector<Real>(numtbins_,0.0)));

            #pragma omp for
            for (int i=0;i<numzbins_;i++){
                auto indices = IndicesPerZbin[i];
                int numRes = indices.size();
                for (int j=0;j<numRes;j++){
                    for (int k=j+1;k<numRes;k++){
                        int ind1 = indices[j];
                        int ind2 = indices[k];

                        // calculate distance squared
                        Real dist_sq;
                        Real3 dist;
                        simstate_.getSimulationBox().calculateDistance(COM_[ind1], COM_[ind2], dist, dist_sq);
                        Real r = std::sqrt(dist_sq);

                        if (rbin_->isInRange(r)){
                            int index = rbin_->findBin(r);
                            Real dotproduct = LinAlg3x3::DotProduct(uij_[ind1], uij_[ind2]);
                            G1r_local_[i][index] += dotproduct;
                            Numres_local_[i][index] += 1;

                            if (tbin_->isInRange(dotproduct)){
                                int tind = tbin_->findBin(dotproduct);
                                G1r_3d_local[i][index][tind] += 1;
                            }
                        }
                    }
                }
            }

            #pragma omp critical
            for (int i=0;i<numzbins_;i++){
                for (int j=0;j<numrbins_;j++){
                    NumberResidueG1r_[i][j] += Numres_local_[i][j];
                    G1r_[i][j] += G1r_local_[i][j];
                    for (int k=0;k<numtbins_;k++){
                        G1r_3d_perIter[i][j][k] += G1r_3d_local[i][j][k];
                    }
                }
            }
        }
    }

    for (int i=0;i<numzbins_;i++){
        for (int j=0;j<numrbins_;j++){
            for (int k=0;k<numtbins_;k++){
                if (InsideNum[i]!=0){
                    G1r_3d_[i][j][k] += (G1r_3d_perIter[i][j][k] / (InsideNum[i] * volume_[j]));
                }
            }
        }
    }
}

void SlabG1r::finishCalculate()
{
    int numFrames = simstate_.getTotalFrames();
    // first process the g1r 3d 
    for (int i=0;i<numzbins_;i++){
        for (int j=0;j<numrbins_;j++){
            G1r_3d_[i][j] = G1r_3d_[i][j] / numFrames;
        }
    }

    // then process the number of residues
    for (int i=0;i<numzbins_;i++){
        for (int j=0;j<numrbins_;j++){
            if (NumberResidueG1r_[i][j] != 0) G1r_[i][j] = G1r_[i][j] / NumberResidueG1r_[i][j];
            else G1r_[i][j] = 0.0;
        }
    }
}

void SlabG1r::printSlabG1r(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<numzbins_;i++){
        for (int j=0;j<numrbins_;j++){
            ofs << G1r_[i][j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}

void SlabG1r::printSlabG1r3d(std::string name){
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<numzbins_;i++){
        for (int j=0;j<numrbins_;j++){
            for (int k=0;k<numtbins_;k++){
                ofs << i << " " << j << " " << k << " " << G1r_3d_[i][j][k] <<"\n";
            }
        }
    }

    ofs.close();
}
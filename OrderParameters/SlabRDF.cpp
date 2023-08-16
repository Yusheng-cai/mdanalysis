#include "SlabRDF.h"

namespace CalculationRegistry
{
    registry_<SlabRDF> registerSlabRDF("SlabRDF");
}

SlabRDF::SlabRDF(const CalculationInput& input) 
: Calculation(input){
    pack_.ReadString("residue", ParameterPack::KeyType::Required, resname_);
    initializeResidueGroup(resname_);
    auto zbin = pack_.findParamPack("zbin", ParameterPack::KeyType::Optional);

    if (zbin != nullptr){
        zbin_ = binptr(new Bin(*zbin));
        numzbins_ = zbin_->getNumbins();
        dz_ = zbin_->getStep();
    }
    else{
        pack_.ReadNumber("numzbins", ParameterPack::KeyType::Required, numzbins_);
        pack_.ReadNumber("above", ParameterPack::KeyType::Optional, above_);
        binUsingMinMax_=true;
        zbin_ = binptr(new Bin());
    }
    zBinLocation_.resize(numzbins_,0.0);

    auto rbin = pack_.findParamPack("rbin", ParameterPack::KeyType::Required);
    rbin_ = binptr(new Bin(*rbin));
    numrbins_ = rbin_->getNumbins();

    slabRDF_.resize(numzbins_, std::vector<Real>(numrbins_, 0.0));
    slabN_.resize(numzbins_,0.0);

    pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional, headindex_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, tailindex_);
    headindex_--;
    tailindex_--;

    pack_.ReadNumber("direction", ParameterPack::KeyType::Optional, direction_);

    pack_.Readbool("within_z", ParameterPack::KeyType::Optional, within_z_);

    // calculate the volume
    volume_.resize(numrbins_,0.0);
    Real dr = rbin_->getStep();
    if (! within_z_){
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

    // register
    registerOutputFunction("RDF", [this](std::string name) -> void {this -> printSlabRDF(name);});
}

void SlabRDF::binUsingMinMax()
{
    Real slight_shift=1e-3;

    std::vector<Real> zdir;
    for (int i=0;i<COM_.size();i++){
        if (COM_[i][direction_] > above_){
            zdir.push_back(COM_[i][direction_]);
        }
    }

    auto minit = std::min_element(zdir.begin(), zdir.end());
    auto maxit = std::max_element(zdir.begin(), zdir.end());

    Real min   = *minit - slight_shift;
    Real max   = *maxit + slight_shift;

    Range range = {{min, max}};

    zbin_->update(range, numzbins_);

    // sum over the bin location 
    for (int i=0;i<numzbins_;i++){
        zBinLocation_[i] += zbin_ -> getCenterLocationOfBin(i);
    }
    std::cout << "Min = " << min << ", Max = " << max << "\n";

    if (within_z_){
        Real dz = zbin_->getStep();
        Real dr = rbin_->getStep();
        // pi * r^2 * L
        volume_[0] = Constants::PI * std::pow(dr, 2.0) * dz;

        for (int i=1;i<numrbins_;i++){
            volume_[i] = 2 * Constants::PI * rbin_->getLeftLocationOfBin(i) * dz * dr;
        }
    }
}

void SlabRDF::calculate(){
    auto& res = getResidueGroup(resname_).getResidues();
    std::vector<std::vector<int>> IndicesPerZbin(numzbins_);
    std::vector<Real> slabN_iter_(numzbins_,0.0);
    std::vector<Real> density(numzbins_,0.0);

    // calculate COM
    for (int i=0;i<res.size();i++){
        COM_[i] = calcCOM(res[i]);
    }

    if (binUsingMinMax_){
        binUsingMinMax();
    }

    for (int i=0;i<COM_.size();i++){
        if (zbin_->isInRange(COM_[i][direction_])){
            int index = zbin_->findBin(COM_[i][direction_]);
            IndicesPerZbin[index].push_back(i);
        }
    }

    for (int i=0;i<numzbins_;i++){
        slabN_[i] += IndicesPerZbin[i].size();
        slabN_iter_[i] = IndicesPerZbin[i].size();

        // Real xy_area = 1;
        // Real3 sim_box = simstate_.getSimulationBox().getSides();
        // for (int j=0;j<3;j++){
        //     if (j != direction_){
        //         xy_area *= sim_box[j];
        //     }
        // }
        // Real v = xy_area * zbin_->getStep();

        // density[i] = slabN_iter_[i] / v;
    }

    std::vector<std::vector<Real>> slabRDF_iter(numzbins_, std::vector<Real>(numrbins_,0.0));

    if (! within_z_){
        #pragma omp parallel
        {
            std::vector<std::vector<Real>> slabRDF_local(numzbins_, std::vector<Real>(numrbins_,0.0));
            #pragma omp for
            for (int i=0;i<COM_.size();i++){
                for (int j=0;j<numzbins_;j++){
                    for (int k=0;k<IndicesPerZbin[j].size();k++){
                        int index = IndicesPerZbin[j][k];
                        if (i != index){
                            Real3 distance;
                            Real dist_sq;
                            simstate_.getSimulationBox().calculateDistance(COM_[i], COM_[index],distance, dist_sq);
                            Real dist = std::sqrt(dist_sq);

                            if (rbin_->isInRange(dist)){
                                int binNum = rbin_->findBin(dist);
                                slabRDF_local[j][binNum] += 1;
                            }
                        }
                    }
                }
            }

            #pragma omp critical
            {
                for (int i=0;i<numzbins_;i++){
                    for (int j=0;j<numrbins_;j++){
                        slabRDF_iter[i][j] += slabRDF_local[i][j];
                    }
                }
            }
        }
    }
    else{
        #pragma omp parallel
        {
            std::vector<std::vector<Real>> slabRDF_local(numzbins_, std::vector<Real>(numrbins_,0.0));
            #pragma omp for
            for (int i=0;i<numzbins_;i++){
                int numRes = IndicesPerZbin[i].size();
                for (int j=0;j<numRes;j++){
                    for (int k=j+1;k<numRes;k++){
                        int ind1 = IndicesPerZbin[i][j];
                        int ind2 = IndicesPerZbin[i][k];

                        Real3 distance;
                        Real dist_sq;
                        simstate_.getSimulationBox().calculateDistance(COM_[ind1], COM_[ind2], distance, dist_sq);
                        Real dist = std::sqrt(dist_sq);

                        if (rbin_->isInRange(dist)){
                            int binNum = rbin_->findBin(dist);
                            slabRDF_local[i][binNum] += 1;
                        }
                    }
                }
            }

            #pragma omp critical
            {
                for (int i=0;i<numzbins_;i++){
                    for (int j=0;j<numrbins_;j++){
                        slabRDF_iter[i][j] += slabRDF_local[i][j];
                    }
                }
            }

        }
    }

    for (int i=0;i<numzbins_;i++){
        for (int j=0;j<numrbins_;j++){
            slabRDF_[i][j] += (slabRDF_iter[i][j] / (slabN_iter_[i] * volume_[j])); 
        }
    }
}

void SlabRDF::printSlabRDF(std::string name){
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<numzbins_;i++){
        for (int j=0;j<numrbins_;j++){
            ofs << i << " " << j << " " << zBinLocation_[i] << " " << rbin_->getCenterLocationOfBin(j) << " " << slabRDF_[i][j] << "\n";
        }
    }

    ofs.close();
}

void SlabRDF::finishCalculate(){
    int numframes = simstate_.getTotalFrames();

    for (int i=0;i<numzbins_;i++){
        for (int j=0;j<numrbins_;j++){
            slabRDF_[i][j] /= numframes;
        }
    }

    if (binUsingMinMax_){
        zBinLocation_ = zBinLocation_ / numframes;
    }
    else{
        for (int i=0;i<numzbins_;i++){
            zBinLocation_[i] = zbin_->getCenterLocationOfBin(i);
        }
    }
}
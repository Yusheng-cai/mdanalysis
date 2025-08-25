#include "DipoleDistributionRDF.hpp"

namespace CalculationRegistry
{
    registry_<DipoleDistribution> registerDipoleDistribution("DipoleDistribution");
}

DipoleDistribution::DipoleDistribution(const CalculationInput& input)
: Calculation(input)
{
    // initialize the probe volumes
    initializeProbeVolumes();
    initializeNotInProbeVolumes();

    // check for probe volume initialization
    ASSERT((probevolumes_.size() > 0 || NotInprobevolumes_.size() > 0), "No probe volume passed in");

    // intialize the bins 
    auto binPack = input.pack_.findParamPack("rbin", ParameterPack::KeyType::Required);
    bin_ = binptr(new Bin(*binPack));
    numbins_ = bin_ ->getNumbins();

    // read inputs
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional,headindex_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, tailindex_);
    headindex_--;
    tailindex_--;

    // read 1 residues 
    input.pack_.ReadString("residue", ParameterPack::KeyType::Required, residueName_);
    initializeResidueGroup(residueName_, "COMIndices", COMIndices_, COM_);
    auto& res = getResidueGroup(residueName_);
    numresidues_ = res.getResidues().size();
    numatoms_ = res[0].atoms_.size();

    // costheta bins is always assumed to go from 0 to 1
    input.pack_.ReadNumber("numtbins", ParameterPack::KeyType::Optional, numtbins_);
    tbin_ = binptr(new Bin(trange_, numtbins_));

    // initialize histogram 2d which contains (R, t)
    histogram2d_.resize(numbins_, std::vector<Real>(numtbins_,0.0));

    registerOutputFunction("histogram2d", [this](std::string name) -> void {this ->printHistogram2d(name);});
}

void DipoleDistribution::calculate(){
    const auto& res = getResidueGroup(residueName_).getResidues(); 
    uij_.clear();
    uij_.resize(res.size());

    // find the COM of 1 & 2
    #pragma omp parallel
    {
        #pragma omp for
        for (int i=0;i<numresidues_;i++)
        {
            COM_[i] = CalculationTools::getCOM(res[i], simstate_, COMIndices_); 

            Real3 headpos = res[i].atoms_[headindex_].positions_;
            Real3 tailpos = res[i].atoms_[tailindex_].positions_;
            Real3 distance;
            Real distance_sq;
            simstate_.getSimulationBox().calculateDistance(headpos, tailpos, distance, distance_sq);

            LinAlg3x3::normalize(distance);

            uij_[i] = distance;
        }
    }

    // find the distance between pair of COM 
    // (numresidues, numbins)
    std::vector<std::vector<Real3>> global_director(numresidues_, std::vector<Real3>(numbins_));
    std::vector<std::vector<Real>> num_molecules(numresidues_, std::vector<Real>(numbins_));
    #pragma omp parallel
    {
        // local histogram for theta bins
        std::vector<std::vector<Real3>> local_director(numresidues_, std::vector<Real3>(numbins_));
        std::vector<std::vector<Real>> local_molecules(numresidues_, std::vector<Real>(numbins_));

        // find distances inside PV
        #pragma omp for
        for (int i=0;i<numresidues_;i++){
            for (int j=0;j<numresidues_;j++){
                if (i != j){
                    // find the distance between ith and jth COM 
                    Real3 distance;
                    Real val;

                    simstate_.getSimulationBox().calculateDistance(COM_[i], COM_[j], distance, val);
                    val = std::sqrt(val);

                    if (bin_ -> isInRange(val)){
                        int rbinnum = bin_ -> findBin(val);
                        for (int k=rbinnum; k < bin_->getNumbins(); k++){
                            local_director[i][k] = local_director[i][k] + uij_[j];
                            local_molecules[i][k] += 1;
                        }
                    }
                }
            }
        }


        #pragma omp critical
        for (int i=0;i<numresidues_;i++){
            for (int j=0;j<numbins_;j++){
                global_director[i][j] = global_director[i][j] + local_director[i][j];
                num_molecules[i][j]   = num_molecules[i][j] + local_molecules[i][j];
            }
        }
    }
    std::cout << "Done 2d." << std::endl;

    for (int i=0;i<global_director.size();i++){
        for (int j=0;j<global_director[i].size();j++){
            global_director[i][j] = global_director[i][j] / num_molecules[i][j];
        }
    }

    // now we bin it 
    // first iterate over the R bins
    for (int j=0;j<numbins_;j++){
        for (int i=0;i < numresidues_;i++){
            // ith molecule, jth R 
            Real magnitude = LinAlg3x3::norm(global_director[i][j]);
            if (tbin_->isInRange(magnitude)){
                int tbin_num = tbin_->findBin(magnitude);
                histogram2d_[j][tbin_num] += 1;
            }
        }
    }
}

void DipoleDistribution::finishCalculate(){
    int numframes = simstate_.getTotalFrames();
    std::cout << "num frames = " << numframes << std::endl;

    // Total cosine theta divided by total histogram
    for (int i=0;i<histogram2d_.size();i++){
        histogram2d_[i] = histogram2d_[i] / numframes;
    }
}

void DipoleDistribution::printHistogram2d(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<numbins_;i++){
        for (int j=0;j<numtbins_;j++){
            ofs << i << "\t" << j << "\t" << bin_ -> getCenterLocationOfBin(i) << "\t" << tbin_ -> getCenterLocationOfBin(j) \
            << "\t" << histogram2d_[i][j] << "\n";
        }
    }

    ofs.close();
}
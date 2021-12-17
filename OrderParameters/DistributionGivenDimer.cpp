#include "DistributionGivenDimer.h"


namespace CalculationRegistry
{
    registry_<DistributionGivenDimer> registerDistributionGivenDimer("DistGivenDimer");
}

DistributionGivenDimer::DistributionGivenDimer(const CalculationInput& input)
:Calculation(input)
{
    pack_.ReadString("residue", ParameterPack::KeyType::Required, resName_);
    initializeResidueGroup(resName_);
    auto& res = getResidueGroup(resName_).getResidues();
    numres_ = res.size();
    AngleWithSurface_.resize(numres_,0.0);
    COM_.resize(numres_);
    DimerPerResidue_.resize(numres_);
    uij_.resize(numres_);

    // read head index and tail index
    pack_.ReadNumber("heaindex", ParameterPack::KeyType::Optional, headindex_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, tailindex_);
    headindex_--;
    tailindex_--;

    // read in the binning information
    auto binPack = pack_.findParamPack("bin", ParameterPack::KeyType::Required);
    bin_ = binptr(new Bin(*binPack));
    numbins_ = bin_ -> getNumbins();
    histogram_.resize(numbins_,0.0);

    // read in the surface normal
    pack_.ReadArrayNumber("surfacenormal", ParameterPack::KeyType::Optional, surfaceNormal_);

    // read in the cosine theta max and r max 
    pack_.ReadNumber("cosmax", ParameterPack::KeyType::Required, cosmax_);
    pack_.ReadNumber("rmax", ParameterPack::KeyType::Required, rmax_);

    initializeProbeVolumes();
    initializeNotInProbeVolumes();

    registerOutputFunction("histogram", [this](std::string name) -> void {this->printHistogram(name);});
}

void DistributionGivenDimer::calculate()
{
    auto& res = getResidueGroup(resName_).getResidues();
    std::vector<int> InsideIndices;
    std::cout << "Got to before first loop" << std::endl;
    for (int i=0;i<res.size();i++)
    {
        COM_[i] = CalculationTools::getCOM(res[i], simstate_, COMIndices_);

        Real3 headpos = res[i].atoms_[headindex_].positions_;
        Real3 tailpos = res[i].atoms_[tailindex_].positions_;

        Real3 distance;
        Real distsq;
        simstate_.getSimulationBox().calculateDistance(headpos, tailpos, distance, distsq);
        
        LinAlg3x3::normalize(distance);

        AngleWithSurface_[i] = LinAlg3x3::DotProduct(distance, surfaceNormal_);
        uij_[i] = distance;
    }
    std::cout << "Got here." << std::endl;

    // find the inside indices 
    InsideIndices = InsidePVIndices(COM_);
    int size = InsideIndices.size();
    std::vector<std::vector<Real>> pairDistances(size, std::vector<Real>(size,0.0));
    std::vector<std::vector<Real>> costhetaPair(size, std::vector<Real>(size,0.0));

    #pragma omp parallel for
    for (int i=0;i<size;i++)
    {
        for (int j=i+1;j<size;j++)
        {
            int index1 = InsideIndices[i];
            int index2 = InsideIndices[j];

            Real3 distance;
            Real distsq;

            simstate_.getSimulationBox().calculateDistance(COM_[index1], COM_[index2], distance, distsq);

            pairDistances[i][j] = std::sqrt(distsq);

            Real costheta = LinAlg3x3::DotProduct(uij_[index1], uij_[index2]);
            costhetaPair[i][j] = costheta;
        }
    }

    for (int i=0;i<size;i++)
    {
        int index1 = InsideIndices[i];
        int count = 0;
        for (int j=i+1;j<size;j++)
        {
            int index2 = InsideIndices[j];
            if (pairDistances[i][j] <= rmax_ && costhetaPair[i][j] <= cosmax_)
            {
                Real angle1 = AngleWithSurface_[index1];
                Real angle2 = AngleWithSurface_[index2];

                if(count ==0)
                {
                    if (bin_ ->isInRange(angle1))
                    {
                        int binnum1 = bin_ -> findBin(angle1);
                        histogram_[binnum1] += 1;
                    }
                }

                if (bin_ -> isInRange(angle2))
                {
                    int binnum2 = bin_ -> findBin(angle2);
                    histogram_[binnum2] += 1;
                }

                count += 1;
            }
        }
    }
}

void DistributionGivenDimer::finishCalculate()
{
    int numFrames = simstate_.getTotalFrames();

    for (int i=0;i<histogram_.size();i++)
    {
        histogram_[i] /= numFrames;
    }
}

void DistributionGivenDimer::printHistogram(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<histogram_.size();i++)
    {
        ofs << bin_ -> getCenterLocationOfBin(i) << "\t" << histogram_[i] << "\n";
    }

    ofs.close();
}

#include "grcostorientation.h"

namespace CalculationRegistry
{
    registry_<grcostorientation> registergrcostorientation("grcostorientation");
}

grcostorientation::grcostorientation(const CalculationInput& input)
:Calculation(input)
{
    // read head and tail index 
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional, headindex_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional, tailindex_);
    headindex_--;
    tailindex_--;

    // initialize the r bin 
    auto rbinPack = pack_.findParamPack("rbin", ParameterPack::KeyType::Required);
    rbin_ = binptr(new Bin(*rbinPack));
    numrbins_ = rbin_ -> getNumbins();

    // initialize the costheta bin
    pack_.ReadNumber("numtbin", ParameterPack::KeyType::Optional, numtbins_);
    tbin_ = binptr(new Bin(trange_, numtbins_));

    // initialize the 2d histogram
    histogram2d_.resize(numrbins_, std::vector<Real>(numtbins_,0.0));
    histogram2dPerIter_.resize(numrbins_, std::vector<Real>(numtbins_,0.0));

    // initialize the residues
    pack_.ReadString("residue", ParameterPack::KeyType::Required, resName_);
    initializeResidueGroup(resName_);
    auto& res = getResidueGroup(resName_);
    COM_.resize(res.getResidues().size());
    AngleWithSurface_.resize(res.getResidues().size(),0.0);
    uij_.resize(res.getResidues().size());

    // initialize probe volumes
    initializeProbeVolumes();
    initializeNotInProbeVolumes();

    // register outputs
    registerOutputFunction("histogram2d", [this](std::string name) -> void {this -> printHistogram2d(name);});
}

void grcostorientation::calculate()
{
    auto& res = getResidueGroup(resName_).getResidues();
    histogram2dPerIter_.clear();
    histogram2dPerIter_.resize(numrbins_, std::vector<Real>(numtbins_,0.0));

    std::vector<std::vector<Real>> numhistogram(numrbins_, std::vector<Real>(numtbins_,0.0));

    std::vector<int> InsideIndices;

    for (int i=0;i<res.size();i++)
    {
        COM_[i] = CalculationTools::getCOM(res[i], simstate_, COMIndices_);

        Real3 headpos = res[i].atoms_[headindex_].positions_;
        Real3 tailpos = res[i].atoms_[tailindex_].positions_;

        Real3 distance;
        Real distsq;
        simstate_.getSimulationBox().calculateDistance(headpos, tailpos, distance, distsq);
        
        LinAlg3x3::normalize(distance);
        uij_[i] = distance;

        AngleWithSurface_[i] = LinAlg3x3::DotProduct(distance, surfaceNormal_);
    }

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

    // now we start binning and write to data 
    for (int i=0;i<size;i++)
    {
        // find the angle that the molecule forms with the surface 
        int index  = InsideIndices[i]; 
        Real angle = AngleWithSurface_[index];

        int binangle = tbin_ -> findBin(angle); 

        for (int j=i+1;j<size;j++)
        {
            Real pD = pairDistances[i][j];

            if (rbin_->isInRange(pD))
            {
                int binR = rbin_ -> findBin(pD);

                histogram2dPerIter_[binR][binangle] += costhetaPair[i][j];
                numhistogram[binR][binangle] += 1;
            }
        }
    }

    for (int i=0;i<numrbins_;i++)
    {
        for (int j=0;j<numtbins_;j++)
        {
            if (numhistogram[i][j] != 0)
            {
                histogram2dPerIter_[i][j] /= numhistogram[i][j];
                histogram2d_[i][j] += histogram2dPerIter_[i][j];
            }
        }
    }
}

void grcostorientation::printHistogram2d(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<numrbins_;i++)
    {
        for (int j=0;j<numtbins_;j++)
        {
            ofs << i << "\t" << j << "\t" << rbin_->getCenterLocationOfBin(i) << "\t" << tbin_->getCenterLocationOfBin(j) << \
            "\t" << histogram2d_[i][j] << "\n";
        }
    }

    ofs.close();
}

void grcostorientation::finishCalculate()
{
    int numframes = simstate_.getTotalFrames();

    for (int i=0;i<numrbins_;i++)
    {
        for (int j=0;j<numtbins_;j++)
        {
            histogram2d_[i][j] /= numframes;
        }
    }
}
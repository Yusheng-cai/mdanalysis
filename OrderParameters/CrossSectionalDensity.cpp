#include "CrossSectionalDensity.h"

CrossSectionalDensity::CrossSectionalDensity(const CalculationInput& input)
: Calculation(input)
{
    auto xbinInput = pack_.findParamPack("xbin", ParameterPack::KeyType::Required);
    xbin_          = binptr(new Bin(*xbinInput));
    xbinsize_      = xbin_->getNumbins();

    auto ybinInput = pack_.findParamPack("ybin", ParameterPack::KeyType::Required);
    ybin_          = binptr(new Bin(*ybinInput));
    ybinsize_      = ybin_->getNumbins();

    histogram_.resize(xbinsize_, std::vector<Real>(ybinsize_,0.0));

    // read and initialize the residue group
    pack_.ReadString("residue", ParameterPack::KeyType::Required, resname_);
    initializeResidueGroup(resname_);

    // register the printing functions
    registerOutputFunction("histogram", [this](std::string name) -> void { this -> printHistogram(name);});
}


void CrossSectionalDensity::calculate()
{
    // calculate the COM
    auto res = getResidueGroup(resname_).getResidues();

    #pragma omp parallel for 
    for (int i=0;i<res.size();i++)
    {
        COM_[i] = calcCOM(res[i]);
    }

    // bin the COM along x and y and store it in histogram 
    for (int i=0;i<COM_.size();i++)
    {
        if (xbin_->isInRange(COM_[i][0]) && ybin_ -> isInRange(COM_[i][1]))
        {
            int xindex = xbin_->findBin(COM_[i][0]);
            int yindex = ybin_->findBin(COM_[i][1]);
            histogram_[xindex][yindex] += 1;
        }
    }
}

void CrossSectionalDensity::finishCalculate()
{
    int numframes = simstate_.getTotalFrames();

    for (int i=0;i<xbinsize_;i++)
    {
        for (int j=0;j<ybinsize_;j++)
        {
            histogram_[i][j]  = histogram_[i][j] / numframes;
        }
    }
}

void CrossSectionalDensity::printHistogram(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<xbinsize_;i++)
    {
        for (int j=0;j<ybinsize_;j++)
        {
            ofs << histogram_[i][j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}
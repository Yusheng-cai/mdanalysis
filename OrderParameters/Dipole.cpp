#include "Dipole.h"

Dipole::Dipole(const CalculationInput& input)
:Calculation(input)
{
    pack_.ReadString("residue", ParameterPack::KeyType::Required, residueName_);
    pack_.ReadArrayNumber("direction", ParameterPack::KeyType::Optional, direction_);

    // add the corresponding residuegroup
    addResidueGroup(residueName_);

    // resize the dipole direction vector
    auto& res = getResidueGroup(residueName_).getResidues();
    dipoledirection_.resize(res.size());

    initializeBins();
}

void Dipole::initializeBins()
{
    auto binPack = pack_.findParamPack("bin", ParameterPack::KeyType::Required);

    bin_ = Binptr(new Bin(*binPack));

    histogram_.resize(bin_->getNumbins(), 0.0);
}

void Dipole::calculate()
{
    // get the residue group
    auto& res = getResidueGroup(residueName_).getResidues();
    int residueSize = res.size();

    std::vector<Real> costheta(residueSize,0.0);
    
    // calculate the dipoles
    for (int i=0;i<residueSize;i++)
    {
        int numAtoms = res[i].atoms_.size();
        auto& atoms  = res[i].atoms_;
        Real3 dipole = {};
        for (int j=0;j<numAtoms;j++)
        {
            for(int k=0;k<3;k++)
            {
                dipole[k] += atoms[j].charge_* atoms[j].positions_[k];
            }
        }

        Real dot_product = LinAlg3x3::DotProduct(dipole, direction_);
        costheta[i] += dot_product;
    }

    // bin the data 
    for (int i=0;i<costheta.size();i++)
    {
        if (bin_-> isInRange(costheta[i]))
        {
            int binNum = bin_ -> findBin(costheta[i]);
            histogram_[binNum] += 1;
        }
    }
}

void Dipole::finishCalculate()
{
    int numSteps = simstate_.getTotalFrames();

    #ifdef MY_DEBUG
    std::cout << "Number of steps = " << numSteps << std::endl;
    #endif  

    for (int i=0;i<histogram_.size();i++)
    {
        histogram_[i] /= numSteps;
    }
}

void Dipole::printHistogram(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    ofs_ << "# positions  Number\n";
    for (int i=0;i<histogram_.size();i++)
    {
        ofs_ << bin_->getCenterLocationOfBin(i) << "\t" << histogram_[i] << "\n"; 
    }

    ofs_.close();
}
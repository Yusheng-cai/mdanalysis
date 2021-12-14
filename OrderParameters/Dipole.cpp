#include "Dipole.h"

namespace CalculationRegistry
{
    registry_<Dipole> registerDipole("dipole");
}

Dipole::Dipole(const CalculationInput& input)
:Calculation(input)
{
    pack_.ReadString("residue", ParameterPack::KeyType::Required, residueName_);
    pack_.ReadArrayNumber("direction", ParameterPack::KeyType::Optional, direction_);

    registerOutputFunction("costheta", [this](std::string name) -> void {this -> printHistogram(name);});
    registerOutputFunction("costhetasquared", [this](std::string name) -> void {this -> printHistogramsquared(name);});

    // add the corresponding residuegroup
    addResidueGroup(residueName_);

    // resize the dipole direction vector
    auto& res = getResidueGroup(residueName_).getResidues();
    dipoledirection_.resize(res.size());

    initializeBins();
    initializeAtomIndices();
}

void Dipole::initializeAtomIndices()
{
    auto& res = getResidueGroup(residueName_).getResidues();
    int atomNum = res[0].atoms_.size();
    Atomindices_.resize(atomNum);

    std::iota(Atomindices_.begin(), Atomindices_.end(), 1);

    pack_.ReadVectorNumber("atomIndices", ParameterPack::KeyType::Optional, Atomindices_);

    for (int i=0;i<Atomindices_.size();i++)
    {
        Atomindices_[i] -= 1;
    }
}

void Dipole::initializeBins()
{
    // initialize histogram and bin for cosine theta 
    auto binPack = pack_.findParamPack("bin", ParameterPack::KeyType::Required);
    bin_ = Binptr(new Bin(*binPack));
    histogram_.resize(bin_->getNumbins(), 0.0);

    // initialize histogram and bin for cosine squared theta
    pack_.ReadNumber("numbinssquared", ParameterPack::KeyType::Optional, numsquaredbin_);
    Range range = {{0,1}}; 
    binsquared_ = Binptr(new Bin(range, numsquaredbin_));
    histogramsquared_.resize(numsquaredbin_, 0.0);
}

void Dipole::calculate()
{
    // get the residue group
    auto& res = getResidueGroup(residueName_).getResidues();
    int residueSize = res.size();

    std::vector<Real> costheta(residueSize,0.0);
    std::vector<Real> costhetasquared(residueSize,0.0);
    
    // calculate the dipoles
    for (int i=0;i<residueSize;i++)
    {
        auto& atoms  = res[i].atoms_;
        Real3 dipole = {};

        Real3 COM = CalculationTools::getCOM(res[i], simstate_, Atomindices_);

        #ifdef MY_DEBUG
        std::cout << "COM for the " << i  << "th molecule is " << COM[0] << " " << COM[1] << " " << COM[2] << std::endl;
        #endif 

        for (int j=0;j<Atomindices_.size();j++)
        {
            int id_ = Atomindices_[j];
            Real3 dist;
            Real distsq;
            simstate_.getSimulationBox().calculateDistance(atoms[id_].positions_, COM, dist, distsq);

            #ifdef MY_DEBUG
            std::cout << "distance of " << id_ << " atom with respect to COM is " << dist[0] << " " << dist[1] << " " << dist[2] << std::endl;
            std::cout << "Atom charge for " << id_ << " atom is " << atoms[id_].charge_ << std::endl;
            #endif 

            // LinAlg3x3::normalize(dist); 

            for(int k=0;k<3;k++)
            {
                dipole[k] += atoms[id_].charge_* dist[k];
            }
        }
        // LinAlg3x3::normalize(dipole);
        Real dot_product = LinAlg3x3::DotProduct(dipole, direction_);
        costheta[i] = dot_product;
        costhetasquared[i] = dot_product * dot_product;

        #ifdef MY_DEBUG
        std::cout << "dipole for the " << i << "th molecule is " << dipole[0] << " " << dipole[1] << " " << dipole[2] << std::endl;
        std::cout << "costheta for the " << i << "th molecule is " << dot_product << std::endl;
        #endif
    }

    // bin the data 
    for (int i=0;i<costheta.size();i++)
    {
        if (bin_-> isInRange(costheta[i]))
        {
            int binNum = bin_ -> findBin(costheta[i]);
            histogram_[binNum] += 1;
        }

        if (binsquared_ -> isInRange(costhetasquared[i]))
        {
            int binNum = binsquared_ -> findBin(costhetasquared[i]);
            histogramsquared_[binNum] += 1;
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

    for (int i=0;i<histogramsquared_.size();i++)
    {
        histogramsquared_[i] /= numSteps;
    }
}

void Dipole::printHistogram(std::string name)
{
    std::ofstream ofs_;
    ofs_ << std::fixed << std::setprecision(precision_);
    ofs_.open(name);

    ofs_ << "# positions  Number\n";
    for (int i=0;i<histogram_.size();i++)
    {
        ofs_ << bin_->getCenterLocationOfBin(i) << "\t" << histogram_[i] << "\n"; 
    }

    ofs_.close();
}

void Dipole::printHistogramsquared(std::string name)
{
    std::ofstream ofs_;
    ofs_ << std::fixed << std::setprecision(precision_);
    ofs_.open(name);

    ofs_ << "# positions  Number\n";
    for (int i=0;i<histogramsquared_.size();i++)
    {
        ofs_ << binsquared_->getCenterLocationOfBin(i) << "\t" << histogramsquared_[i] << "\n";
    }

    ofs_.close();
}
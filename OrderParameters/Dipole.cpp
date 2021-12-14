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
    histogram_.resize(bin_->getNumbins(), {});

    // initialize histogram and bin for cosine squared theta
    pack_.ReadNumber("numbinssquared", ParameterPack::KeyType::Optional, numsquaredbin_);
    Range range = {{0,1}}; 
    binsquared_ = Binptr(new Bin(range, numsquaredbin_));
    histogramsquared_.resize(numsquaredbin_, {});
}

void Dipole::calculate()
{
    // get the residue group
    auto& res = getResidueGroup(residueName_).getResidues();
    int residueSize = res.size();

    std::vector<Real3> dipole_total(residueSize,{0,0,0});
    std::vector<Real3> dipolesquared(residueSize,{0,0,0});
    
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

            for(int k=0;k<3;k++)
            {
                dipole[k] += atoms[id_].charge_* dist[k];
            }
        }
        LinAlg3x3::normalize(dipole);
        dipole_total[i] = dipole;
        Real3 dipoles;
        for (int j=0;j<3;j++)
        {
            dipoles[j] = dipole[j] * dipole[j];
        }
        dipolesquared[i] = dipoles;

        #ifdef MY_DEBUG
        std::cout << "dipole for the " << i << "th molecule is " << dipole[0] << " " << dipole[1] << " " << dipole[2] << std::endl;
        #endif
    }

    // bin the data 
    for (int i=0;i<dipole_total.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            if (bin_-> isInRange(dipole_total[i][j]))
            {
                int binNum = bin_ -> findBin(dipole_total[i][j]);
                histogram_[binNum][j] += 1;
            }

            if (binsquared_ -> isInRange(dipolesquared[i][j]))
            {
                int binNum = binsquared_ -> findBin(dipolesquared[i][j]);
                histogramsquared_[binNum][j] += 1;
            }
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
        for (int j=0;j<3;j++)
        {
            histogram_[i][j] /= numSteps;
        }
    }

    for (int i=0;i<histogramsquared_.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            histogramsquared_[i][j] /= numSteps;
        }
    }
}

void Dipole::printHistogram(std::string name)
{
    std::ofstream ofs_;
    ofs_ << std::fixed << std::setprecision(precision_);
    ofs_.open(name);

    ofs_ << "# positions\tNumberx\tNumbery\tNumberz\n";
    for (int i=0;i<histogram_.size();i++)
    {
        ofs_ << bin_->getCenterLocationOfBin(i) << "\t";
        for (int j=0;j<3;j++)
        {
            ofs_ << histogram_[i][j] << "\t";
        }
        ofs_ << "\n";
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
        ofs_ << binsquared_->getCenterLocationOfBin(i) << "\t"; 
        for (int j=0;j<3;j++)
        {
            ofs_ << histogramsquared_[i][j] << "\t";
        }
        ofs_ << "\n";
    }

    ofs_.close();
}
#include "SlabAtomicProperty.h"

namespace CalculationRegistry
{
    registry_<SlabAtomicProperty> registerSlabAtomicProperty("SlabAtomicProperty");
}

SlabAtomicProperty::SlabAtomicProperty(const CalculationInput& input)
: Calculation(input)
{
    pack_.ReadString("residue", ParameterPack::KeyType::Required, residueName_);
    pack_.ReadString("property", ParameterPack::KeyType::Required, property_);
    pack_.ReadNumber("direction", ParameterPack::KeyType::Optional, direction_);

    auto binPack = pack_.findParamPack("bin", ParameterPack::KeyType::Required);
    bin_ = Binptr(new Bin(*binPack));
    numbins_ = bin_->getNumbins();
    dz_ = bin_->getStep();
    AtomicProperty_.resize(numbins_, 0.0);

    // add the residue group
    initializeResidueGroup(residueName_);
    
    // register the calculation functions 
    registerCalcFunc("density", [this](std::vector<Real3>& COM) -> std::vector<Real> {return this -> CalculateDensity(COM);});
    registerOutputFunction("property", [this](std::string name) -> void {this -> printProperty(name);});
}

void SlabAtomicProperty::registerCalcFunc(std::string name, calcfunc func)
{
    auto it  = MapNameToCalculationFunc_.find(name);
    ASSERT((it == MapNameToCalculationFunc_.end()), "The function " << name << " is registered more than once.");
    MapNameToCalculationFunc_.insert(std::make_pair(name, func));
}

void SlabAtomicProperty::calculate()
{
    auto& res = getResidueGroup(residueName_).getResidues();
    int numres= res.size();
    COM_.clear();
    COM_.resize(numres);

    #pragma omp parallel for
    for (int i=0;i<numres;i++)
    {
        COM_[i] = calcCOM(res[i]);
    }

    auto it = MapNameToCalculationFunc_.find(property_);
    ASSERT((it != MapNameToCalculationFunc_.end()), "The property with name " << property_ << " is not found.");
    std::vector<Real> tempProperty = it -> second(COM_);

    for (int i=0;i<numbins_;i++)
    {
        AtomicProperty_[i] = AtomicProperty_[i] + tempProperty[i];
    }
}

void SlabAtomicProperty::printProperty(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << "# Bin  Property\n";
    for (int i=0;i<numbins_;i++)
    {
        ofs << i+1 << " " << AtomicProperty_[i] << "\n";
    }

    ofs.close();
}

void SlabAtomicProperty::finishCalculate()
{
    int numframes = simstate_.getTotalFrames();

    for (int i=0;i<numbins_;i++)
    {
        AtomicProperty_[i] = AtomicProperty_[i] / numframes;
    }
}

std::vector<SlabAtomicProperty::Real> SlabAtomicProperty::CalculateDensity(std::vector<Real3>& COM)
{
    std::vector<Real> Density(numbins_, 0.0);
    auto box = simstate_.getSimulationBox().getBox();
    Real area = 1.0;
    for (int i=0;i<3;i++)
    {
        if (i != direction_)
        {
            area *= box[i][i];
        }
    }

    #pragma omp parallel
    {
        std::vector<Real> localDensity(numbins_, 0.0);

        #pragma omp for
        for (int i=0;i<COM_.size();i++)
        {
            if (bin_->isInRange(COM_[i][direction_]))
            {
                int binIndex = bin_->findBin(COM_[i][direction_]);
                localDensity[binIndex] += 1;
            }
        }

        #pragma omp critical
        {
            for (int i =0;i<numbins_;i++)
            {
                Density[i] += localDensity[i];
            }
        }
    }

    for (int i=0;i<Density.size();i++)
    {
        Density[i] = Density[i] / (area * dz_);
    }

    return Density;
}
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

    auto binPack = pack_.findParamPack("bin", ParameterPack::KeyType::Optional);
    if(binPack != nullptr)
    {
        bin_    = Binptr(new Bin(*binPack));
        numbins_= bin_->getNumbins();
        dz_     = bin_->getStep();
    }
    else
    {
        pack_.ReadNumber("numzbins", ParameterPack::KeyType::Required, numbins_);
        pack_.ReadNumber("above", ParameterPack::KeyType::Required, above_);
        bin_    = Binptr(new Bin());
        usingMinMax_ = true;
        dz_     = 0.0;
    }

    ResidueLocationPerBin_.resize(numbins_,0.0);
    AtomicProperty_.resize(numbins_, 0.0);

    // add the residue group
    initializeResidueGroup(residueName_);
    
    // register the calculation functions 
    registerCalcFunc("density", [this](void) -> std::vector<Real> {return this -> CalculateDensity();});
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

    auto it = MapNameToCalculationFunc_.find(property_);
    ASSERT((it != MapNameToCalculationFunc_.end()), "The property with name " << property_ << " is not found.");
    std::vector<Real> tempProperty = it -> second();

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
        ofs << ResidueLocationPerBin_[i] << " " << AtomicProperty_[i] << "\n";
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

    if (usingMinMax_)
    {
        for (int i=0;i<numbins_;i++)
        {
            ResidueLocationPerBin_[i] = ResidueLocationPerBin_[i] / numframes;
        }
    }
    else
    {
        for (int i=0;i<numbins_;i++)
        {
            ResidueLocationPerBin_[i] = bin_->getCenterLocationOfBin(i);
        }
    }
}

void SlabAtomicProperty::binUsingMinMax(const std::vector<Real3>& positions)
{
    Real slight_shift=1e-3;
    std::vector<Real> ZdirectionNum;

    for (int i=0;i<positions.size();i++){
        Real val = positions[i][direction_];

        if (val > above_){
            ZdirectionNum.push_back(val);
        }
    }

    auto maxit = std::max_element(ZdirectionNum.begin(), ZdirectionNum.end());
    auto minit = std::min_element(ZdirectionNum.begin(), ZdirectionNum.end());

    Real max = *maxit + slight_shift;
    Real min = *minit - slight_shift;
    Range2 range = {{min, max}};

    std::cout << "Max = " << max << " Min = " << min << std::endl;

    bin_ -> update(range, numbins_);
    dz_ = bin_->getStep();

    for (int i=0;i<numbins_;i++){
        ResidueLocationPerBin_[i] += bin_->getCenterLocationOfBin(i);
    }
}

std::vector<SlabAtomicProperty::Real> SlabAtomicProperty::CalculateDensity()
{
    auto& Res = getResidueGroup(residueName_).getResidues();

    std::vector<Real> Density(numbins_, 0.0);
    auto box = simstate_.getSimulationBox().getBox();
    Real area = 1.0;
    for (int i=0;i<3;i++){
        if (i != direction_){
            area *= box[i][i];
        }
    }

    std::vector<Real3> positions;
    for (int i=0;i<Res.size();i++){
        for (int j=0;j<Res[i].atoms_.size();j++){
            positions.push_back(Res[i].atoms_[j].positions_);
        }
    }

    if (usingMinMax_){
        binUsingMinMax(positions);
    }

    #pragma omp parallel
    {
        std::vector<Real> localDensity(numbins_, 0.0);

        #pragma omp for
        for (int i=0;i<positions.size();i++){
            if (bin_->isInRange(positions[i][direction_])){
                int binIndex = bin_->findBin(positions[i][direction_]);
                localDensity[binIndex] += 1;
            }
        }

        #pragma omp critical
        {
            for (int i =0;i<numbins_;i++){
                Density[i] += localDensity[i];
            }
        }
    }

    for (int i=0;i<Density.size();i++){
        Density[i] = Density[i] / (area * dz_);
    }

    return Density;
}
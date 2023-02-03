#include "SlabOrientation.h"

namespace CalculationRegistry
{
    static registry_<SlabOrientation> registerSlabOrientation("SlabOrientation");
}

SlabOrientation::SlabOrientation(const CalculationInput& input)
:Calculation(input)
{
    // initialize the bins first --> zbin and thetabin
    auto zbinPack = input.pack_.findParamPack("zbin", ParameterPack::KeyType::Optional);

    // default to be 30
    pack_.ReadNumber("numtbins", ParameterPack::KeyType::Optional, numtbins_);
    costBin_ = Binptr(new Bin({{-1,1}}, numtbins_));
    costsquaredBin_ = Binptr(new Bin({{0,1}}, numtbins_));

    if(zbinPack != nullptr){
        zBin_    = Binptr(new Bin(*zbinPack));
        numzbins_= zBin_->getNumbins();
    }
    else
    {
        pack_.ReadNumber("numzbins", ParameterPack::KeyType::Required, numzbins_);
        pack_.ReadNumber("above", ParameterPack::KeyType::Required, above_);
        zBin_    = Binptr(new Bin());
        usingMinMax_ = true;
    }

    // register the output functions
    RegisterOutputs();

    // read the inputs 
    ReadInputs();

    // initialize residue group
    initializeResidueGroup(residueGroupName_);
    auto& res = getResidueGroup(residueGroupName_);

    // size histogram 2d      
    histogram2d_.resize(numzbins_, std::vector<Real>(numtbins_,0.0));
    histogram2d_squared_.resize(numzbins_, std::vector<Real>(numtbins_,0.0));

    // resize the molecular director size
    uij_.resize(res.size());

    // resize the number of residues to size of residue
    numResiduePerBin_.resize(numzbins_,0.0);
    ResidueLocationPerBin_.resize(numzbins_,0.0);
}

void SlabOrientation::RegisterOutputs()
{
    registerOutputFunction("PCosthetaZ", [this](std::string name)->void {this -> printHistogram(name);});
    registerOutputFunction("PCosthetaZSquared", [this](std::string name) -> void {this -> printHistogramSquared(name);});
    registerOutputFunction("ResiduePerBin", [this](std::string name) -> void {this -> printNumResidue(name);});
}

void SlabOrientation::ReadInputs()
{
    pack_.ReadString("residue", ParameterPack::KeyType::Required, residueGroupName_);
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Optional,headIndex_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Optional,tailIndex_);
    pack_.ReadString("direction", ParameterPack::KeyType::Optional, direction_);
    pack_.ReadArrayNumber("array", ParameterPack::KeyType::Optional, arr_);
    pack_.ReadNumber("precision", ParameterPack::KeyType::Optional, precision_);

    // correct for the indexing on head and tail index 
    headIndex_--;
    tailIndex_--;

    // find the direction
    directionIndex_ = MapdirectionToIndex_.find(direction_) -> second;
}

void SlabOrientation::calculate()
{
    // obtain the residue group
    auto& res = getResidueGroup(residueGroupName_).getResidues();

    // find all the COM of the residues in the system
    #pragma omp parallel for 
    for (int i=0;i<res.size();i++)
    {
        COM_[i] = CalculationTools::getCOM(res[i], simstate_, COMIndices_);

        Real3 distance;
        Real dist_sq;
        Real3 headPos = res[i].atoms_[headIndex_].positions_;
        Real3 tailPos = res[i].atoms_[tailIndex_].positions_;
        simstate_.getSimulationBox().calculateDistance(headPos, tailPos, distance, dist_sq);
        LinAlg3x3::normalize(distance);
        uij_[i] = distance;
    }

    // bin using min max of the molecules if needed 
    if (usingMinMax_)
    {
        binUsingMinMax();
    }

    // find the distances between all pairs of residues
    for (int i=0;i<res.size();i++)
    {
        Real dist = COM_[i][directionIndex_];
        Real cost = LinAlg3x3::DotProduct(uij_[i], arr_);
        if (zBin_->isInRange(dist) && costBin_->isInRange(cost))
        {  
            int zbinNum = zBin_->findBin(dist);
            int tbinNum = costBin_->findBin(cost);
            int tsquaredbinNum = costsquaredBin_->findBin(cost*cost);

            histogram2d_[zbinNum][tbinNum] += 1;
            numResiduePerBin_[zbinNum] += 1;

            histogram2d_squared_[zbinNum][tsquaredbinNum] += 1;
        }
    }
}

void SlabOrientation::binUsingMinMax()
{
    Real slight_shift=1e-3;
    std::vector<Real> ZdirectionNum;

    for (int i=0;i<COM_.size();i++){
        Real val = COM_[i][directionIndex_];

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

    zBin_ -> update(range, numzbins_);

    for (int i=0;i<numzbins_;i++)
    {
        ResidueLocationPerBin_[i] += zBin_->getCenterLocationOfBin(i);
    }
}

void SlabOrientation::finishCalculate()
{
    int Numframes = simstate_.getTotalFrames();

    // first normalize by the number of frames 
    for (int i=0;i<numResiduePerBin_.size();i++)
    {
        numResiduePerBin_[i] /= Numframes;
        ResidueLocationPerBin_[i] /= Numframes;
    }

    // normalize the entire histogram such as \sum \sum P(z, cos(\beta)) = 1
    // first divide P(z, cos(\beta)) by numframes
    for (int i=0;i<histogram2d_.size();i++)
    {
        for (int j=0;j<histogram2d_[0].size();j++)
        {
            histogram2d_[i][j] /= Numframes;
            histogram2d_squared_[i][j] /= Numframes;
        }
    }


    // calculate anchoring energy 
    Real3 sides = simstate_.getSimulationBox().getSides(); 
    Real AreaXY = 1;
    for (int i=0;i<3;i++)
    {
        if (i != directionIndex_)
        {
            AreaXY *= sides[i];
        }
    }
    
    // Calculate normalized 2d histogram
    std::vector<std::vector<Real>> PcostZ(numzbins_, std::vector<Real>(numtbins_,0.0));
    BWcostZ_.resize(numzbins_, std::vector<Real>(numtbins_,0.0));
    Real sum = 0.0;
    for (int i=0;i<histogram2d_.size();i++){
        for (int j=0;j<histogram2d_[0].size();j++){
            sum += histogram2d_[i][j];
        }
    }

    for (int i=0;i<PcostZ.size();i++){
        for (int j=0;j<PcostZ[i].size();j++){
            PcostZ[i][j] = histogram2d_[i][j] / sum;
            BWcostZ_[i][j] = 0;
            if (PcostZ[i][j] != 0){
                BWcostZ_[i][j] = -std::log(PcostZ[i][j])/AreaXY * numResiduePerBin_[i];
            }
        }
    }

    for (int i=0;i<BWcostZ_.size();i++)
    {
        Real mini= *std::min_element(BWcostZ_[i].begin(), BWcostZ_[i].end());
        for (int j=0;j<BWcostZ_[i].size();j++)
        {
            BWcostZ_[i][j] = BWcostZ_[i][j] - mini;
        }
    }
}

void SlabOrientation::printHistogram(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << std::fixed << std::setprecision(precision_);

    for (int i=0;i<histogram2d_.size();i++)
    {
        for (int j=0;j<histogram2d_[i].size();j++)
        {
            ofs << i << "\t" << j << "\t" << histogram2d_[i][j] << "\t" << BWcostZ_[i][j] << "\n";
        }
    }
    ofs.close();
}

void SlabOrientation::printHistogramSquared(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << std::fixed << std::setprecision(precision_);

    for (int i=0;i<histogram2d_squared_.size();i++)
    {
        for (int j=0;j<histogram2d_squared_[i].size();j++)
        {
            ofs << i << "\t" << j << "\t" << ResidueLocationPerBin_[i] << "\t" << costBin_->getCenterLocationOfBin(j) \
            << "\t" << histogram2d_squared_[i][j] << "\n";
        }
    }
}

void SlabOrientation::printNumResidue(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << std::fixed << std::setprecision(precision_);

    for (int i=0;i<numResiduePerBin_.size();i++)
    {
        ofs << ResidueLocationPerBin_[i] << " " << numResiduePerBin_[i] << "\n";
    }
    ofs.close();
}

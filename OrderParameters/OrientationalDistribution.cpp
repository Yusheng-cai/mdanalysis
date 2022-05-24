#include "OrientationalDistribution.h"

namespace CalculationRegistry
{
    registry_<OrientationalDistribution> registerOrientationalDistribution("OrientationalDistribution");
}

OrientationalDistribution::OrientationalDistribution(const CalculationInput& input)
:Calculation(input)
{
    pack_.ReadString("residue", ParameterPack::KeyType::Required, ResidueGroupName_);
    pack_.ReadNumber("headindex", ParameterPack::KeyType::Required,HeadIndex_);
    pack_.ReadNumber("tailindex", ParameterPack::KeyType::Required,TailIndex_);
    pack_.ReadNumber("numbins", ParameterPack::KeyType::Required, NumBins_);
    pack_.ReadArrayNumber("array", ParameterPack::KeyType::Optional, arr_);
    pack_.Readbool("usedirector", ParameterPack::KeyType::Optional, useDirector_);
    pack_.ReadNumber("xbins", ParameterPack::KeyType::Optional, numxbin_);
    pack_.ReadNumber("ybins", ParameterPack::KeyType::Optional, numybin_);
    HeadIndex_--;
    TailIndex_--;

    // resize the xy histogram
    PcostDistributionXY_.resize(numxbin_, std::vector<Real>(numybin_,0.0));
    HisotgramXY_.resize(numxbin_, std::vector<Real>(numybin_,0.0));
    xbin_ = Binptr(new Bin());
    ybin_ = Binptr(new Bin());

    // initialize Residue
    initializeResidueGroup(ResidueGroupName_);

    // initialize the probe volumes
    initializeProbeVolumes();
    initializeNotInProbeVolumes();

    ASSERT((probevolumes_.size() > 0 || NotInprobevolumes_.size() > 0), "Probe volume list must be provided");

    // initialize the bins 
    CosThetaBin_ = Binptr(new Bin(CosThetaRange_, NumBins_));
    CosThetaSquaredBin_ = Binptr(new Bin(CosThetaSquaredRange_, NumBins_));
    PCosTheta_.resize(NumBins_,0.0);
    PCosThetaSquared_.resize(NumBins_,0.0);

    // register the output functions 
    registerOutputs();
    registerOutputfile();
}

void OrientationalDistribution::registerOutputfile()
{
    registerOutputFileOutputs("costheta", [this](void)-> Real {return this -> getAvgCostheta();});
    registerOutputFileOutputs("costhetasquared", [this](void)-> Real {return this -> getAvgCosthetasquared();});
}


void OrientationalDistribution::registerOutputs()
{
    registerOutputFunction("Distribution", [this](std::string name) -> void {this -> PrintDistribution(name);});
    registerPerIterOutputFunction("costhetasquared_betafactor", [this](std::ofstream& ofs) -> void {this -> PrintCosthetasquared_betafactors(ofs);});
    registerPerIterOutputFunction("ResidueAngles", [this](std::ofstream& ofs)-> void {this -> printResidueAngles(ofs);});
    registerOutputFunction("DistributionXY", [this](std::string name) -> void {this->PrintDistributionXY(name);});
}

void OrientationalDistribution::finishCalculate()
{
    int numFrames = simstate_.getTotalFrames();
    Real CosThetaTot = 0.0;
    Real CosSquaredThetaTot = 0.0;

    for (int i=0;i<NumBins_;i++)
    {
        PCosTheta_[i] = PCosTheta_[i]/numFrames;
        PCosThetaSquared_[i] = PCosThetaSquared_[i] / numFrames;

        CosThetaTot += PCosTheta_[i];
        CosSquaredThetaTot += PCosThetaSquared_[i];
    }

    for (int i=0;i<numxbin_;i++)
    {
        for (int j=0;j<numybin_;j++)
        {
            if (HisotgramXY_[i][j] != 0)
            {
                PcostDistributionXY_[i][j] = PcostDistributionXY_[i][j] / HisotgramXY_[i][j];
            }
        }
    }
}

void OrientationalDistribution::PrintDistribution(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << "# Bin\tcostheta\tcos2theta\tP(costheta)\tP(cos2theta)\n";

    for (int i=0;i<NumBins_;i++)
    {
        ofs <<CosThetaBin_->getCenterLocationOfBin(i) << " " << CosThetaSquaredBin_->getCenterLocationOfBin(i) << " " \
        << PCosTheta_[i] << " " << PCosThetaSquared_[i] << "\n";
    }
    ofs.close();
}

void OrientationalDistribution::update()
{
    auto& res    = getResidueGroup(ResidueGroupName_).getResidues();
    int AtomSize = simstate_.getTotalNumberAtoms();
    COM_.clear();
    COM_.resize(res.size());
    uij_.clear();
    uij_.resize(res.size());
    Real3 BoxSides = simstate_.getSimulationBox().getSides();

    #pragma omp parallel for 
    for (int i=0;i<res.size();i++)
    {
        COM_[i] = calcCOM(res[i]);

        Real3 HeadPos = res[i].atoms_[HeadIndex_].positions_;
        Real3 TailPos = res[i].atoms_[TailIndex_].positions_;
        Real3 distance;
        Real dist_sq;

        simstate_.getSimulationBox().calculateDistance(HeadPos, TailPos,distance, dist_sq);
        LinAlg3x3::normalize(distance);

        uij_[i] = distance;
    }

    // find atoms within probe volume
    AtomIndices_ = InsidePVIndices(COM_);

    if (useDirector_)
    {
        Matrix Qtensor = {};
        #pragma omp parallel
        {
            Matrix Qlocal = {};
            for (int i=0;i<AtomIndices_.size();i++)
            {
                int index = AtomIndices_[i];
                Matrix Q  = LinAlg3x3::LocalQtensor(uij_[index]);
                LinAlg3x3::matrix_accum_inplace(Qlocal, Q);
            }

            #pragma omp critical
            {
                LinAlg3x3::matrix_accum_inplace(Qtensor, Qlocal);
            }
        }
        LinAlg3x3::matrix_mult_inplace(Qtensor, 1.0/(2.0 * (Real)AtomIndices_.size()));
        auto res = LinAlg3x3::OrderEigenSolver(Qtensor);

        for (int i=0;i<3;i++)
        {
            arr_[i] = res.second[i][0];
        }
    }

    // update the x and y bins
    xbin_->update({{0,BoxSides[0]}}, numxbin_);
    ybin_->update({{0,BoxSides[1]}}, numybin_);
}

void OrientationalDistribution::calculate()
{
    auto& res    = getResidueGroup(ResidueGroupName_).getResidues();
    int AtomSize = simstate_.getTotalNumberAtoms();
    costhetasquared_betafactor_.clear();
    costhetasquared_betafactor_.resize(AtomSize,0.0);
    MapIndexToAngle_.clear();
    AvgCostheta_ = 0.0;
    AvgCosthetasquared_= 0.0;

    // iterate over the atom indices
    #pragma omp parallel
    {
        std::vector<Real> localcostheta(NumBins_, 0.0);
        std::vector<Real> localcosthetasquared(NumBins_, 0.0);
        Real localavgcostheta = 0.0;
        Real localavgcosthetasquared = 0.0;

        // local map from residue index to its cos(theta) value with array
        std::map<int, Real> localMap;

        // local xy histogram
        std::vector<std::vector<Real>> localHistogram(numxbin_, std::vector<Real>(numybin_,0.0));
        std::vector<std::vector<Real>> localcosthetaxy(numxbin_, std::vector<Real>(numybin_,0.0));

        #pragma omp for
        for (int i=0;i<AtomIndices_.size();i++)
        {
            int k = AtomIndices_[i];
            Real cost = LinAlg3x3::DotProduct(uij_[k], arr_);
            Real costsq = cost * cost;

            // map index to the cos theta of the angle 
            localMap.insert(std::make_pair(k, cost));

            localavgcostheta += cost;
            localavgcosthetasquared += costsq;

            int CosThetaBinNum = CosThetaBin_->findBin(cost);
            int CosThetaSquaredBinNum = CosThetaSquaredBin_->findBin(costsq);

            localcostheta[CosThetaBinNum] += 1;
            localcosthetasquared[CosThetaSquaredBinNum] += 1;

            for (int j=0;j<res[k].atoms_.size();j++)
            {
                int index = res[k].atoms_[j].atomNumber_-1;
                costhetasquared_betafactor_[index] = costsq;
            }

            // find the x and y index 
            int xindex = xbin_->findBin(COM_[k][0]);
            int yindex = ybin_->findBin(COM_[k][1]);
            localHistogram[xindex][yindex] += 1;
            localcosthetaxy[xindex][yindex]+= cost;
        }

        #pragma omp critical
        {
            for (int i=0;i<NumBins_;i++)
            {
                PCosTheta_[i] = PCosTheta_[i] + localcostheta[i];
                PCosThetaSquared_[i] = PCosThetaSquared_[i] + localcosthetasquared[i];
            }

            AvgCosthetasquared_ = AvgCosthetasquared_ + localavgcosthetasquared;
            AvgCostheta_ = AvgCostheta_ + localavgcostheta;

            MapIndexToAngle_.insert(localMap.begin(), localMap.end());

            for (int i=0;i<numxbin_;i++)
            {
                for (int j=0;j<numybin_;j++)
                {
                    HisotgramXY_[i][j] += localHistogram[i][j];
                    PcostDistributionXY_[i][j] += localcosthetaxy[i][j];
                }
            }
        }
    }

    AvgCosthetasquared_ = AvgCosthetasquared_ * 1.0/AtomIndices_.size();
    AvgCostheta_ = AvgCostheta_ * 1.0/AtomIndices_.size();
}

void OrientationalDistribution::PrintCosthetasquared_betafactors(std::ofstream& ofs)
{
    int totalatoms = simstate_.getTotalNumberAtoms();
    int framenum   = simstate_.getFrameNumber();

    ofs << framenum << " ";
    for (int i=0;i<totalatoms;i++)
    {
        ofs << costhetasquared_betafactor_[i] << " ";
    }
    ofs << "\n";
}

void OrientationalDistribution::printResidueAngles(std::ofstream& ofs)
{
    int framenum = simstate_.getFrameNumber();

    ofs << "# Frame " << framenum << "\n"; 

    for (auto it = MapIndexToAngle_.begin(); it != MapIndexToAngle_.end(); it ++)
    {
        ofs << it -> first << " " << it -> second << "\n";
    }
}

void OrientationalDistribution::PrintDistributionXY(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<numxbin_;i++)
    {
        for (int j=0;j<numybin_;j++)
        {
            ofs << PcostDistributionXY_[i][j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}
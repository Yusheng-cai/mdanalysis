#include "RDFresidue.h"

namespace CalculationRegistry
{
    registry_<RDFresidue> registerRDFresidue("RDFresidue");
}

RDFresidue::RDFresidue(const CalculationInput& input)
:Calculation(input)
{
    input.pack_.ReadString("residue1", ParameterPack::KeyType::Required, resname1_);
    input.pack_.ReadString("residue2", ParameterPack::KeyType::Required, resname2_);

    // initialize bins
    auto binPack = input.pack_.findParamPack("bin", ParameterPack::KeyType::Required);
    bins_ = Binptr(new Bin(*binPack));

    // initialize the volume 
    volume_.resize(bins_->getNumbins());
    Real dr = bins_ -> getStep();
    volume_[0] = 4.0/3.0*Constants::PI*std::pow(dr,3.0);
    for (int i=1;i<bins_->getNumbins();i++)
    {
        Real r = bins_ -> getLeftLocationOfBin(i);
        volume_[i] = 4.0*Constants::PI*std::pow(r,2.0)*dr;
    }

    // add the residue groups
    initializeResidueGroup(resname1_, "COMIndices1", COMIndices1_, COM1_);
    initializeResidueGroup(resname2_, "COMIndices2", COMIndices2_, COM2_);

    // resize rdf vector to be same lengths as number of bins
    rdf_.resize(bins_->getNumbins());

    // register outputs
    registerOutputFunction("unnormalizedRDF", [this](std::string name) -> void {this -> printRDFUnnormalized(name);});
    registerOutputFunction("RDF", [this](std::string name) -> void {this -> printRDF(name);});
    registerOutputFunction("numneighbors", [this](std::string name) -> void {this -> printNumNeighbors(name);});
}

void RDFresidue::calculate()
{
    const auto& res1 = getResidueGroup(resname1_).getResidues();
    const auto& res2 = getResidueGroup(resname2_).getResidues();

    int N = res1.size();

    // obtain the density (num Residues/ nm3)
    Real rho = res1.size()/simstate_.getSimulationBox().getVolume();
    rho_ += rho;

    COM1_.clear();
    COM1_.resize(res1.size());
    COM2_.clear();
    COM2_.resize(res2.size());

    #pragma omp parallel
    {
        #pragma omp for
        for (int i=0;i<COM1_.size();i++)
        {
            Real3 C1 = CalculationTools::getCOM(res1[i], simstate_, COMIndices1_);

            COM1_[i] = C1;
        }

        #pragma omp for
        for (int i=0;i<COM2_.size();i++)
        {
            Real3 C2 = CalculationTools::getCOM(res2[i], simstate_, COMIndices2_);

            COM2_[i] = C2;
        }
    }


    std::vector<int> NumCountsPerBin(bins_->getNumbins(), 0);
    #pragma omp parallel
    {
        std::vector<int> NumCountPerBinLocal(bins_ -> getNumbins(),0);
        #pragma omp for
        for (int i=0;i<COM1_.size();i++)
        {
            for (int j=0;j<COM2_.size();j++)
            {
                Real3 distance;
                Real sq_dist;
                simstate_.getSimulationBox().calculateDistance(COM1_[i], COM2_[j],distance, sq_dist);

                Real dist = std::sqrt(sq_dist);

                if (dist != 0)
                {
                    if (bins_-> isInRange(dist))
                    {
                        int binnum = bins_ -> findBin(dist);
                        NumCountPerBinLocal[binnum] += 1;
                    }
                }
            }
        }

        #pragma omp critical 
        for (int i=0;i<NumCountsPerBin.size();i++)
        {
            NumCountsPerBin[i] += NumCountPerBinLocal[i];
        }
    }

    for (int i=0;i<rdf_.size();i++)
    {
        rdf_[i] = rdf_[i] + NumCountsPerBin[i];
    }

}

void RDFresidue::finishCalculate()
{
    // The differential r is the step in the bins
    Real dr = bins_->getStep();
    int numBins = bins_->getNumbins();

    // normalize the numCountsPerBin By the number of frames first
    int numFrames_ = simstate_.getTotalFrames();

    for (int i=0;i<rdf_.size();i++)
    {
        rdf_[i] = rdf_[i]/numFrames_;
    }

    rho_ /= numFrames_;
}

void RDFresidue::printNumNeighbors(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    std::vector<Real> temp(rdf_.size(),0.0);
    const auto& res1 = getResidueGroup(resname1_).getResidues();
    int N = res1.size();

    for (int i=0;i<rdf_.size();i++)
    {
        temp[i] = rdf_[i]/N; 
    }

    for (int i=0;i<rdf_.size();i++)
    {
        ofs << bins_->getLeftLocationOfBin(i) << "\t" << temp[i] << "\n";
    }

    ofs.close();
}

void RDFresidue::printRDFUnnormalized(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    std::vector<Real> rdftemp(rdf_.size(),0.0);
    const auto& res1 = getResidueGroup(resname1_);
    int N = res1.size();

    for (int i=1;i<rdf_.size();i++)
    {
        rdftemp[i] = rdf_[i] / (volume_[i] * N); 
    }

    for (int i=0;i<rdftemp.size();i++)
    {
        ofs << bins_->getCenterLocationOfBin(i) << "\t" << rdftemp[i] << "\n";
    }

    ofs.close();
}

void RDFresidue::printRDF(std::string name)
{
    std::ofstream ofs; 
    ofs.open(name);

    std::vector<Real> rdftemp(rdf_.size(),0.0);
    const auto& res1 = getResidueGroup(resname1_);
    int N = res1.size();

    for (int i=1;i<rdf_.size();i++)
    {
        rdftemp[i] = rdf_[i] / (volume_[i] * N * rho_);
    }

    for (int i=0;i<rdftemp.size();i++)
    {
        ofs << bins_ -> getCenterLocationOfBin(i) << "\t" << rdftemp[i] << "\n";
    }

    ofs.close();
}
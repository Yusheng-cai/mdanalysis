#include "RDF2d.h"

RDF2d::RDF2d(const CalculationInput& input)
: Calculation(input)
{
    pack_.ReadString("residue1", ParameterPack::KeyType::Required, resname1_);
    pack_.ReadString("residue2", ParameterPack::KeyType::Required, resname2_);
    pack_.ReadArrayNumber("positionIndex", ParameterPack::KeyType::Optional, position_index_);

    // initialize bins
    auto binPack = input.pack_.findParamPack("bin", ParameterPack::KeyType::Required);
    bins_ = Binptr(new Bin(*binPack));

    // initialize the volume 
    Area_.resize(bins_->getNumbins());
    Real dr = bins_ -> getStep();
    Area_[0] = 2.0 * Constants::PI * dr * dr;

    for (int i=1;i<bins_->getNumbins();i++)
    {
        Real r = bins_ -> getLeftLocationOfBin(i);
        Area_[i] = 2.0*Constants::PI*r*dr;
    }

    // add the residue groups
    initializeResidueGroup(resname1_, "COMIndices1", COMIndices1_, COM1_);
    initializeResidueGroup(resname2_, "COMIndices2", COMIndices2_, COM2_);

    // resize rdf vector to be same lengths as number of bins
    rdf_.resize(bins_->getNumbins());
}

void RDF2d::calculate()
{

}
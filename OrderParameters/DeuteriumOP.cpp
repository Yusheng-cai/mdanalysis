#include "DeuteriumOP.h"

namespace CalculationRegistry
{
    registry_<DeuteriumOP> registerDeuteriumOP("DeuteriumOP");
}

DeuteriumOP::DeuteriumOP(const CalculationInput& input)
: Calculation(input)
{
    // read in the carbon indices that user inputs
    pack_.ReadVectorNumber("CarbonIndices", ParameterPack::KeyType::Required, CarbonIndices_);
    numCarbons_ = CarbonIndices_.size()-2; // exclude the terminal and first carbon 
    avgSCD_.resize(numCarbons_,0.0);
    avgSZZ_.resize(numCarbons_,0.0);
    avgSCC_.resize(numCarbons_+1,0.0);

    // read in residue
    pack_.ReadString("residue", ParameterPack::KeyType::Required, residue_);
    initializeResidueGroup(residue_);

    // read in surface normal, optional
    pack_.ReadArrayNumber("SurfaceNormal", ParameterPack::KeyType::Optional, SurfaceNormal_);

    // register output functions
    registerOutputFunction("DeuteriumOP", [this](std::string name)->void {this -> printDeuteriumOP(name);});
}

void DeuteriumOP::calculate()
{
    auto& res = getResidueGroup(residue_).getResidues();

    for (int i=0;i<res.size();i++)
    {
        auto& r = res[i];
        for (int j=1;j<CarbonIndices_.size()-1;j++)
        {
            // the carbon below it 
            Real3 posBelow = r.atoms_[j-1].positions_;
            Real3 posAbove = r.atoms_[j+1].positions_;

            // find out the z axis 
            Real3 zaxis = LinAlg3x3::vec_sub(posAbove, posBelow);
            LinAlg3x3::normalize(zaxis);

            // find the rotation matrix that rotates from zaxis to [0,0,1]
            Matrix InvRotMat = LinAlg3x3::GetRotationMatrix({{0,0,1}}, zaxis);

            // rotate y axis to the frame of reference
            Real3 yaxis = LinAlg3x3::MatrixDotVector(InvRotMat, {{0,1,0}});
            Real3 xaxis = LinAlg3x3::MatrixDotVector(InvRotMat, {{1,0,0}});

            // calculate Scd
            Real cosy = LinAlg3x3::DotProduct(yaxis, SurfaceNormal_);
            Real cosx = LinAlg3x3::DotProduct(xaxis, SurfaceNormal_);
            Real cosz = LinAlg3x3::DotProduct(zaxis, SurfaceNormal_);

            // Scd
            Real Sxx  = 1.5 * cosx * cosx  - 0.5;
            Real Syy  = 1.5 * cosy * cosy  - 0.5;
            Real Szz  = 1.5 * cosz * cosz  - 0.5;
            Real Scd  = 2.0/3.0 * Sxx + 1.0/3.0 * Syy;

            avgSCD_[j-1] += Scd;
            avgSZZ_[j-1] += Szz;
        }

        for (int j=1;j<CarbonIndices_.size();j++)
        {
            // the first carbon 
            Real3 FirstCarbon = r.atoms_[0].positions_;
            Real3 ThisCarbon  = r.atoms_[j].positions_;

            Real3 dir   = LinAlg3x3::vec_sub(ThisCarbon, FirstCarbon);
            LinAlg3x3::normalize(dir);

            Real cosz   = LinAlg3x3::DotProduct(dir, SurfaceNormal_);

            // Scc
            Real Scc    = 1.5 * cosz * cosz - 0.5;

            avgSCC_[j-1] += Scc;
        }
    }
}

void DeuteriumOP::finishCalculate()
{
    int numFrames = simstate_.getTotalFrames();
    auto& res     = getResidueGroup(residue_);
    Real ressize   = (Real)res.size();

    for (int i=0;i<numCarbons_;i++)
    {
        avgSCD_[i] = avgSCD_[i] / (numFrames * ressize);
        avgSZZ_[i] = avgSZZ_[i] / (numFrames * ressize);
    }

    for (int i=0;i<numCarbons_+1;i++)
    {
        avgSCC_[i] = avgSCC_[i] / (numFrames * ressize);
    }
}

void DeuteriumOP::printDeuteriumOP(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << "# NumCarbon  Scd  Scdc  Scc\n";
    for (int i=0;i<numCarbons_;i++)
    {
        ofs << i+1 << " " << avgSCD_[i] << " " << avgSZZ_[i] << " " << avgSCC_[i] << "\n";
    }

    ofs << numCarbons_+1 << "\t" << "\t" << "\t" << avgSCC_[numCarbons_] << "\n";

    ofs.close();
}

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

            // Scd
            Real Sxx  = 1.5 * cosx * cosx  - 0.5;
            Real Syy  = 1.5 * cosy * cosy  - 0.5;
            Real Scd  = 2.0/3.0 * Sxx + 1.0/3.0 * Syy;

            avgSCD_[j-1] += Scd;
        }
    }
}

void DeuteriumOP::finishCalculate()
{
    int numFrames = simstate_.getTotalFrames();

    for (int i=0;i<numCarbons_;i++)
    {
        avgSCD_[i] = avgSCD_[i] / numFrames;
    }
}

void DeuteriumOP::printDeuteriumOP(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << "# NumCarbon  DeuteriumOP\n";
    for (int i=0;i<numCarbons_;i++)
    {
        ofs << i+1 << " " << avgSCD_[i] << "\n";
    }

    ofs.close();
}

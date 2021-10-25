#include "ProbeVolumeTiltedCylinder.h"

namespace ProbeVolumes
{
    static const registry_<ProbeVolumeTiltedCylinder> register_box("tilted_cylinder");
}

ProbeVolumeTiltedCylinder::ProbeVolumeTiltedCylinder(ProbeVolumeInput& input)
: ProbeVolume(input)
{
    input.ParamPack.ReadNumber("radius", ParameterPack::KeyType::Required, radius_);

    if (isDynamic())
    {
        addDynamicAtomGroup(dynamicAtomName_);

        auto& DynamicAG = getDynamicAtomGroup(dynamicAtomName_);

        ASSERT((DynamicAG.getNumAtoms() == 2), "The dynamic atom for dynamic cylinder is specified by 2 atoms.");
    }

    // if not dynamic, then we read in the base & top
    if (! isDynamic())
    {
        input.ParamPack.ReadArrayNumber("base", ParameterPack::KeyType::Required, base_);
        input.ParamPack.ReadArrayNumber("top", ParameterPack::KeyType::Required, top_);

        Real dist_ = 0.0;
        for (int i=0;i<3;i++)
        {
            Zvector_[i] = top_[i] - base_[i];
            dist_ += Zvector_[i] * Zvector_[i];
        }
        dist_ = std::sqrt(dist_);

        rotationMat_ = LinAlg3x3::GetRotationMatrix(Zvector_, refVector);

        cylinder_ -> setGeometry(radius_, dist_, ac_, sigma_);
    }
}

void ProbeVolumeTiltedCylinder::setGeometry()
{
    if (isDynamic())
    {
        auto& DynamicAG = getDynamicAtomGroup(dynamicAtomName_);

        const auto& Atoms = DynamicAG.getAtoms();

        base_ = Atoms[0].position;
        top_  = Atoms[1].position;

        // calculate the pbc corrected distance between top and base 
        getSimulationBox().calculateDistance(base_, top_, Zvector_, distance_);

        LinAlg3x3::normalize(Zvector_);

        // make the rotation matrix 
        rotationMat_ = LinAlg3x3::GetRotationMatrix(Zvector_, refVector);

        cylinder_ -> setGeometry(radius_, distance_, ac_, sigma_);
    }
}

ProbeVolumeOutput ProbeVolumeTiltedCylinder::calculate(const Real3& x) const
{
    Real3 dist_;
    Real d_;
    getSimulationBox().calculateDistance(x, base_, dist_, d_);

    Real3 rotatedVec = LinAlg3x3::MatrixDotVector(rotationMat_, dist_);

    auto output = cylinder_ -> calculate(rotatedVec);

    return output;
}

void ProbeVolumeTiltedCylinder::update()
{
    setGeometry();
}
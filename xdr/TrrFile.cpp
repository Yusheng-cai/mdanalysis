#include "TrrFile.h"

namespace XdrFiles
{
    static const registry_<TrrFile> register_trrfile("trr");
}

TrrFile::TrrFile()
:XdrWrapper()
{};

void TrrFile::readNextFrame()
{
    ASSERT((isOpen()), "The file is not opened.");

    auto& position = frame_.accessPositions();
    ASSERT((position.size() == natoms_), "The position vector does not have the correct number of atoms.");
    auto& velocities = frame_.accessVelocities();
    ASSERT((velocities.size() == natoms_), "The velocity vector does not have the correct number of atoms.");
    auto& forces = frame_.accessForces();
    ASSERT((forces.size() == natoms_), "The force vector does not have the correct number of atoms.");

    rvec* position_ptr = (rvec*)position.data();
    rvec* velocities_ptr = (rvec*)velocities.data(); 
    rvec* forces_ptr = (rvec*)forces.data();

    matrix box;
    int step;
    Frame::Real time, lambda; 
    int has_prop;
    read_trr(file_, natoms_, &step, &time, &lambda, box, position_ptr, velocities_ptr, forces_ptr, &has_prop);

    frame_.setTime(time);
    frame_.setStep(step);
    Frame::Matrix box_;
    for (int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            box_[i][j] = box[i][j];
        }
    }
    frame_.setBox(box_);
}

void TrrFile::readNumAtoms()
{
    ASSERT((isOpen()), "The file is not open.");

    int sucess = read_trr_natoms(const_cast<char*>(name_.c_str()),&natoms_);
    ASSERT((sucess == exdrOK), "Reading natoms failed.");
}
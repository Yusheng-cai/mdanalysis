#include "SimulationState.h"

Calculation::Calculation(const CalculationInput& input)
:simstate_(input.simstate_), pack_(input.pack_)
{
    pack_.ReadVectorString("outputs", ParameterPack::KeyType::Optional, vectorOutputs_);
    pack_.ReadVectorString("outputNames", ParameterPack::KeyType::Optional, vectorOutputNames_);

    ASSERT((vectorOutputs_.size() == vectorOutputNames_.size()), "The size of the vector outputs does not agree with the \
    size of names.");

    pack_.ReadVectorString("perIteroutputs", ParameterPack::KeyType::Optional, perIteroutputs_);
    pack_.ReadVectorString("perIteroutputNames", ParameterPack::KeyType::Optional, perIteroutputNames_);

    ASSERT((perIteroutputs_.size() == perIteroutputNames_.size()), "The number of output files must match the number of outputs.");
    pack_.ReadNumber("precision", ParameterPack::KeyType::Optional, precision_);
    pack_.ReadString("name", ParameterPack::KeyType::Optional, name_);

    pack_.ReadString("COMmode", ParameterPack::KeyType::Optional, COM_mode_);

    for (int i=0;i<perIteroutputs_.size();i++)
    {
        ofsVector_.push_back(ofsptr(new std::ofstream));
    }

    for (int i=0;i<perIteroutputs_.size();i++)
    {
        ofsVector_[i]->open(perIteroutputNames_[i]);
    }
}

Calculation::Real3 Calculation::calcCOM(const Molecule::residue& residues)
{
    if (COM_mode_ == "mass")
    {
        return CalculationTools::getCOM(residues, simstate_, COMIndices_);
    }
    else if (COM_mode_ == "charge")
    {
        return CalculationTools::getCOC(residues, simstate_, COMIndices_);
    }
    else
    {
        ASSERT((true == false), "The mode " << COM_mode_ << " is not yet available.");
    }
}


void Calculation::initializeNotInProbeVolumes()
{
    pack_.ReadVectorString("NotInprobevolumes", ParameterPack::KeyType::Optional, NotInprobevolumeNames_);

    for (int i=0;i<NotInprobevolumeNames_.size();i++)
    {
        auto& pv = simstate_.getProbeVolume(NotInprobevolumeNames_[i]);
        NotInprobevolumes_.push_back(&pv);
    }

}
void Calculation::initializeProbeVolumes()
{
    // read in the names of the probevolumes
    pack_.ReadVectorString("probevolumes", ParameterPack::KeyType::Optional, probevolumeNames_);

    for (int i=0;i<probevolumeNames_.size();i++)
    {
        auto& pv = simstate_.getProbeVolume(probevolumeNames_[i]);
        probevolumes_.push_back(&pv);
    }
}


void Calculation::registerOutputFileOutputs(std::string name, OutputValue::ValueFunction func)
{
    OutputValue val(name, func);

    output_.insert(std::make_pair(name, val));
}

void Calculation::addAtomgroup(std::string name)
{
    // Check if the AtomGroupName is not in the map
    auto it = MapAtomGroupNameToIndex_.find(name);
    ASSERT((it == MapAtomGroupNameToIndex_.end()), "The AtomGroup with name " << name << " is registered twice.");

    // Get the AtomGroup reference from simulation state
    AtomGroup& ag = simstate_.getAtomGroup(name);

    int index = AtomGroups_.size();
    AtomGroups_.push_back(&ag);

    MapAtomGroupNameToIndex_.insert(std::make_pair(name, index));
}

void Calculation::addResidueGroup(std::string name)
{
    // Check if the ResidueGroupName is not in the map
    auto it = MapResidueGroupNameToIndex_.find(name);
    ASSERT((it == MapResidueGroupNameToIndex_.end()), "The AtomGroup with name " << name << " is registered twice.");

    // Get the ResidueGroup reference from simulation state
    ResidueGroup& rg = simstate_.getResidueGroup(name);

    int index = ResidueGroups_.size();
    ResidueGroups_.push_back(&rg);

    MapResidueGroupNameToIndex_.insert(std::make_pair(name, index));
}

const ResidueGroup& Calculation::getResidueGroup(std::string name) const 
{
    auto it  = MapResidueGroupNameToIndex_.find(name);

    ASSERT((it != MapResidueGroupNameToIndex_.end()), "The residue with name " << name << " is not registered.");

    return *ResidueGroups_[it->second]; 
}

void Calculation::initializeResidueGroup(const std::string& residueName)
{
    ASSERT((! residueName.empty()), "The residue name is not provided.");

    // add the residue group
    addResidueGroup(residueName);

    // initialize the COM indice
    auto& res = getResidueGroup(residueName).getResidues();
    COMIndices_.resize(res[0].atoms_.size());
    std::iota(COMIndices_.begin(), COMIndices_.end(), 1);

    // COMIndices are in 1-based counting 
    pack_.ReadVectorNumber("COMIndices", ParameterPack::KeyType::Optional, COMIndices_);
    for (int i=0;i<COMIndices_.size();i++)
    {
        COMIndices_[i] -= 1;
    }

    // resize the COM 
    COM_.resize(res.size());
}

void Calculation::printOutput()
{
    for (int i=0;i<vectorOutputs_.size();i++)
    {
        getOutputByName(vectorOutputs_[i])(vectorOutputNames_[i]);
    }

    closeAllOutputPerIter();
}

void Calculation::registerOutputFunction(std::string name, outputFunc func)
{
    auto it = MapNameToOutputFunction_.find(name);

    ASSERT((it == MapNameToOutputFunction_.end()), "The output with name " << name << " is already registered in calculation.");

    MapNameToOutputFunction_.insert(std::make_pair(name, func));
}

void Calculation::registerPerIterOutputFunction(std::string name, perIteroutputFunc func)
{
    auto it = MapNameToPerIterOutput_.find(name);

    ASSERT((it == MapNameToPerIterOutput_.end()), "The per iter output with name " << name << " is already registered in calculation.");

    MapNameToPerIterOutput_.insert(std::make_pair(name, func));
}

Calculation::perIteroutputFunc& Calculation::getIterOutputByName(std::string name)
{
    auto it = MapNameToPerIterOutput_.find(name);

    ASSERT((it != MapNameToPerIterOutput_.end()), "The per iter output with name " << name << " is not registered.");

    return it -> second;
}


Calculation::outputFunc& Calculation::getOutputByName(std::string name)
{
    auto it = MapNameToOutputFunction_.find(name);

    ASSERT((it != MapNameToOutputFunction_.end()), "The output with name " << name << " is not found within calculation");

    return it -> second;
}

void Calculation::printOutputOnStep()
{
    for (int i=0;i<perIteroutputs_.size();i++)
    {
        std::string name = perIteroutputs_[i];
        getIterOutputByName(name)(*ofsVector_[i]);
    }
}

void Calculation::closeAllOutputPerIter()
{
    for (int i=0;i<ofsVector_.size();i++)
    {
        ofsVector_[i]->close();
    }
}

CalculationTools::Real3 CalculationTools::getCOM(const Molecule::residue& residueGroup, const SimulationState& simstate, \
std::vector<int>& indices_)
{
    auto& simbox = simstate.getSimulationBox();

    // For COM calculation, for each residue, we shift the atoms with respect to the first atom
    // obtain the position of the first atom in the residue group
    int index0 = indices_[0];
    auto& pos1 = residueGroup.atoms_[index0].positions_;

    // Total mass of the atoms of interest
    Real massTot = 0;
    Real3 COM_pos = {{0,0,0}};
    

    // iterate over the indices of interest
    for (int j=0;j<indices_.size();j++)
    {
        Real3 distance;
        Real distsq;

        int index = indices_[j];
        Real3 shiftWRTatom1 = simbox.calculateShift(residueGroup.atoms_[index].positions_, pos1);

        Real3 shiftedPos;

        for (int k=0;k<3;k++)
        {
            shiftedPos[k] = residueGroup.atoms_[index].positions_[k] + shiftWRTatom1[k];
        }

        for (int k=0;k<3;k++)
        {
            COM_pos[k] += shiftedPos[k] * residueGroup.atoms_[index].mass_;
        }

        massTot += residueGroup.atoms_[index].mass_;
    }

    for (int j=0;j<3;j++)
    {
        COM_pos[j] = COM_pos[j]/massTot;
    }

    // check if it is outside of the box, if it is, then move it inside the box
    Real3 shift = simbox.calculateShift(COM_pos, simbox.getCenter());

    for (int j=0;j<3;j++)
    {
        COM_pos[j] += shift[j];
    }

    return COM_pos;
}

CalculationTools::Real3 CalculationTools::getCOC(const Molecule::residue& residueGroup, const SimulationState& simstate, \
std::vector<int>& indices_)
{
    auto& simbox = simstate.getSimulationBox();

    // For COM calculation, for each residue, we shift the atoms with respect to the first atom
    // obtain the position of the first atom in the residue group
    int index0 = indices_[0];
    auto& pos1 = residueGroup.atoms_[index0].positions_;

    // Total mass of the atoms of interest
    Real chargeTot = 0;
    Real3 COM_pos = {{0,0,0}};
    

    // iterate over the indices of interest
    for (int j=0;j<indices_.size();j++)
    {
        Real3 distance;
        Real distsq;

        int index = indices_[j];
        Real3 shiftWRTatom1 = simbox.calculateShift(residueGroup.atoms_[index].positions_, pos1);

        Real3 shiftedPos;

        for (int k=0;k<3;k++)
        {
            shiftedPos[k] = residueGroup.atoms_[index].positions_[k] + shiftWRTatom1[k];
        }

        for (int k=0;k<3;k++)
        {
            COM_pos[k] += shiftedPos[k] * residueGroup.atoms_[index].charge_;
        }
    }

    // check if it is outside of the box, if it is, then move it inside the box
    Real3 shift = simbox.calculateShift(COM_pos, simbox.getCenter());

    for (int j=0;j<3;j++)
    {
        COM_pos[j] += shift[j];
    }

    return COM_pos;
}

std::vector<int> Calculation::InsidePVIndices(std::vector<Real3>& pos)
{
    std::vector<int> insideindices;
    OpenMP::OpenMP_buffer<std::vector<int>> insideindicesbuffer;
    insideindicesbuffer.set_master_object(insideindices);

    // first find which set of COM are inside of PV
    #pragma omp parallel
    {
        auto& buffer_ = insideindicesbuffer.access_buffer_by_id();
        buffer_.clear();

        #pragma omp for
        for (int i=0;i<pos.size();i++)
        {
            bool inPV = false;
            for (auto pv : NotInprobevolumes_)
            {
                auto out = pv -> calculate(pos[i]);

                if (out.hx_ == 1)
                {
                    // break breaks out of the closest enclosing for loop
                    inPV = true;
                    break;
                }
            }

            if (! inPV)
            {
                for (auto pv : probevolumes_)
                {
                    auto out = pv -> calculate(pos[i]);
                    if (out.hx_ == 1)
                    {
                        buffer_.push_back(i);
                        break;
                    }
                }
            }
        }
    }

    int size=insideindices.size();
    for (auto it = insideindicesbuffer.beginworker(); it != insideindicesbuffer.endworker(); it ++)
    {
        size += it -> size();
    }

    insideindices.reserve(size);

    for (auto it = insideindicesbuffer.beginworker(); it != insideindicesbuffer.endworker(); it ++)
    {
        insideindices.insert(insideindices.end(), it -> begin(), it -> end());
    }

    return insideindices;
}

std::vector<int> Calculation::InsidePVIndices(std::vector<Real3>& pos, std::vector<int>& outsideIndices)
{
    std::vector<int> insideindices;
    outsideIndices.clear();

    OpenMP::OpenMP_buffer<std::vector<int>> insideindicesbuffer;
    OpenMP::OpenMP_buffer<std::vector<int>> outsideindicesbuffer;

    insideindicesbuffer.set_master_object(insideindices);
    outsideindicesbuffer.set_master_object(outsideIndices);

    // first find which set of COM are inside of PV
    #pragma omp parallel
    {
        auto& buffer_ = insideindicesbuffer.access_buffer_by_id();
        auto& outsidebuffer = outsideindicesbuffer.access_buffer_by_id();
        buffer_.clear();
        outsidebuffer.clear();

        #pragma omp for
        for (int i=0;i<pos.size();i++)
        {
            bool inPV = false;
            for (auto pv : NotInprobevolumes_)
            {
                auto out = pv -> calculate(pos[i]);

                if (out.hx_ == 1)
                {
                    // break breaks out of the closest enclosing for loop
                    inPV = true;
                    break;
                }
            }

            if (! inPV)
            {
                bool foundinpv=false;
                for (auto pv : probevolumes_)
                {
                    auto out = pv -> calculate(pos[i]);
                    if (out.hx_ == 1)
                    {
                        buffer_.push_back(i);
                        foundinpv = true;
                        break;
                    }
                }

                if (! foundinpv)
                {
                    outsidebuffer.push_back(i);
                }
            }
        }
    }

    int size=insideindices.size();
    int outsidesize = outsideIndices.size();

    for (auto it = insideindicesbuffer.beginworker(); it != insideindicesbuffer.endworker(); it ++)
    {
        size += it -> size();
    }

    for (auto it = outsideindicesbuffer.beginworker(); it != outsideindicesbuffer.endworker(); it ++)
    {
        outsidesize += it -> size();
    }

    insideindices.reserve(size);
    outsideIndices.reserve(outsidesize);

    for (auto it = insideindicesbuffer.beginworker(); it != insideindicesbuffer.endworker(); it ++)
    {
        insideindices.insert(insideindices.end(), it -> begin(), it -> end());
    }

    for (auto it = outsideindicesbuffer.beginworker(); it != outsideindicesbuffer.endworker(); it ++)
    {
        outsideIndices.insert(outsideIndices.end(), it -> begin(), it -> end());
    }

    outsideindicesbuffer.clearMasterObject();

    return insideindices;
}



CalculationTools::Real3 CalculationTools::getCOG(const Molecule::residue& residueGroup, const SimulationState& simstate, \
std::vector<int>& indices_)
{
    auto& simbox = simstate.getSimulationBox();

    // For COM calculation, for each residue, we shift the atoms with respect to the first atom
    // obtain the position of the first atom in the residue group
    int index0 = indices_[0];
    auto& pos1 = residueGroup.atoms_[index0].positions_;

    // Total mass of the atoms of interest
    Real chargeTot = 0;
    Real3 COM_pos = {{0,0,0}};
    

    // iterate over the indices of interest
    for (int j=0;j<indices_.size();j++)
    {
        Real3 distance;
        Real distsq;

        int index = indices_[j];
        Real3 shiftWRTatom1 = simbox.calculateShift(residueGroup.atoms_[index].positions_, pos1);

        Real3 shiftedPos;

        for (int k=0;k<3;k++)
        {
            shiftedPos[k] = residueGroup.atoms_[index].positions_[k] + shiftWRTatom1[k];
        }

        for (int k=0;k<3;k++)
        {
            COM_pos[k] += shiftedPos[k];
        }
    }

    // check if it is outside of the box, if it is, then move it inside the box
    Real3 shift = simbox.calculateShift(COM_pos, simbox.getCenter());

    for (int j=0;j<3;j++)
    {
        COM_pos[j] += shift[j];
    }

    return COM_pos;
}
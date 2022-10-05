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

Calculation::Real3 Calculation::calcCOM(const Molecule::residue& residues){
    if (COM_mode_ == "mass"){return CalculationTools::getCOM(residues, simstate_, COMIndices_);}
    else if (COM_mode_ == "charge"){return CalculationTools::getCOC(residues, simstate_, COMIndices_);}
    else{ASSERT((true == false), "The mode " << COM_mode_ << " is not yet available.");}
}

Calculation::Real3 Calculation::calcCOM(const Molecule::residue& residues, std::vector<int>& COMIndices)
{
    if (COM_mode_ == "mass"){return CalculationTools::getCOM(residues, simstate_, COMIndices);}
    else if (COM_mode_ == "charge"){return CalculationTools::getCOC(residues, simstate_, COMIndices);}
    else if (COM_mode_ == "geometry"){return CalculationTools::getCOG(residues, simstate_, COMIndices);}
    else{ASSERT((true == false), "The mode " << COM_mode_ << " is not yet available.");}
}

void Calculation::initializeNotInProbeVolumes(){
    pack_.ReadVectorString("NotInprobevolumes", ParameterPack::KeyType::Optional, NotInprobevolumeNames_);

    for (int i=0;i<NotInprobevolumeNames_.size();i++){
        auto& pv = simstate_.getProbeVolume(NotInprobevolumeNames_[i]);
        NotInprobevolumes_.push_back(&pv);
    }

}
void Calculation::initializeProbeVolumes(){
    // read in the names of the probevolumes
    pack_.ReadVectorString("probevolumes", ParameterPack::KeyType::Optional, probevolumeNames_);

    for (int i=0;i<probevolumeNames_.size();i++){
        auto& pv = simstate_.getProbeVolume(probevolumeNames_[i]);
        probevolumes_.push_back(&pv);
    }
}


void Calculation::registerOutputFileOutputs(std::string name, OutputValue::ValueFunction func){
    OutputValue val(name, func);

    output_.insert(std::make_pair(name, val));
}

void Calculation::addAtomgroup(std::string name){
    // Check if the AtomGroupName is not in the map
    auto it = MapAtomGroupNameToIndex_.find(name);
    ASSERT((it == MapAtomGroupNameToIndex_.end()), "The AtomGroup with name " << name << " is registered twice.");

    // Get the AtomGroup reference from simulation state
    AtomGroup& ag = simstate_.getAtomGroup(name);

    int index = AtomGroups_.size();
    AtomGroups_.push_back(&ag);

    MapAtomGroupNameToIndex_.insert(std::make_pair(name, index));
}

const AtomGroup& Calculation::getAtomGroup(std::string name) const
{
    auto it  = MapAtomGroupNameToIndex_.find(name);

    ASSERT((it != MapAtomGroupNameToIndex_.end()), "The residue with name " << name << " is not registered.");

    return *AtomGroups_[it->second];
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
    COMIndices_ = COMIndices_ - 1;

    // resize the COM 
    COM_.resize(res.size());
}

void Calculation::initializeResidueGroup(const std::string& residueName, std::string COMName, std::vector<int>& COMIndices, \
std::vector<Real3>& COM)
{
    ASSERT((! residueName.empty()), "The residue name is not provided.");

    // add the residue group
    addResidueGroup(residueName);

    // initialize the COM indice
    auto& res = getResidueGroup(residueName).getResidues();
    COMIndices.resize(res[0].atoms_.size());
    std::iota(COMIndices.begin(), COMIndices.end(), 1);

    // COMIndices are in 1-based counting 
    pack_.ReadVectorNumber(COMName, ParameterPack::KeyType::Optional, COMIndices);
    COMIndices = COMIndices - 1;

    // resize the COM 
    COM.resize(res.size());
}

bool Calculation::isInPV(Real3& pos)
{
    // if it's in the probe volumes that it's not supposed to be in 
    // then we should return false  --> if that makes any sense
    for (auto pv : NotInprobevolumes_){
        auto out = pv -> calculate(pos);

        if (out.hx_ == 1){
            return false;
        }
    }

    for (auto pv : probevolumes_){
        auto out = pv -> calculate(pos);
        if (out.hx_ == 1){
            return true;
        }
    }

    return false;
}

bool Calculation::isInPV(Real3& pos, Real& htildex)
{
    htildex = 1.0;
    for (auto pv : NotInprobevolumes_){
        auto out = pv -> calculate(pos);
        if (out.htilde_x_ > 0){
            return false;
        }
    }

    for (auto pv : probevolumes_){
        auto out = pv -> calculate(pos);
        if (out.htilde_x_ > 0){
            htildex = out.htilde_x_;
            return true;
        }
    }

    return false;

}

void Calculation::printOutput()
{
    for (int i=0;i<vectorOutputs_.size();i++){
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
    for (int i=0;i<ofsVector_.size();i++){
        ofsVector_[i]->close();
    }
}

CalculationTools::Real3 CalculationTools::getCOM(const Molecule::residue& residueGroup, const SimulationState& simstate, \
std::vector<int>& indices_){
    auto& simbox = simstate.getSimulationBox();

    // For COM calculation, for each residue, we shift the atoms with respect to the first atom
    // obtain the position of the first atom in the residue group
    int index0 = indices_[0];
    auto& pos1 = residueGroup.atoms_[index0].positions_;

    // Total mass of the atoms of interest
    Real massTot = 0;
    Real3 COM_pos = {{0,0,0}};
    

    // iterate over the indices of interest
    for (int j=0;j<indices_.size();j++){
        Real3 distance;
        Real distsq;

        int index = indices_[j];
        Real3 shiftWRTatom1 = simbox.calculateShift(residueGroup.atoms_[index].positions_, pos1);

        Real3 shiftedPos = residueGroup.atoms_[index].positions_ + shiftWRTatom1;
        COM_pos = COM_pos + shiftedPos * residueGroup.atoms_[index].mass_;
        massTot += residueGroup.atoms_[index].mass_;
    }

    // calculate the center of mass 
    COM_pos = COM_pos / massTot;

    // check if it is outside of the box, if it is, then move it inside the box
    Real3 shift = simbox.calculateShift(COM_pos, simbox.getCenter());
    COM_pos = COM_pos + shift;

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
    Real3 COM_pos = {{0,0,0}};
    
    // iterate over the indices of interest
    for (int j=0;j<indices_.size();j++){
        Real3 distance;
        Real distsq;

        int index = indices_[j];
        Real3 shiftWRTatom1 = simbox.calculateShift(residueGroup.atoms_[index].positions_, pos1);

        Real3 shiftedPos = residueGroup.atoms_[index].positions_ + shiftWRTatom1;
        COM_pos = COM_pos + shiftedPos * residueGroup.atoms_[index].charge_;
    }

    // check if it is outside of the box, if it is, then move it inside the box
    Real3 shift = simbox.calculateShift(COM_pos, simbox.getCenter());
    COM_pos = COM_pos + shift;

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
            if (isInPV(pos[i]))
            {
                buffer_.push_back(i);
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
            if (isInPV(pos[i]))
            {
                buffer_.push_back(i);
            }
            else
            {
                outsidebuffer.push_back(i);
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

void Calculation::ReadResidueIndices(const std::string& residuename, std::string IndicesName, std::vector<int>& Indices)
{
    // cleat the indices 
    Indices.clear();

    const auto res = getResidueGroup(residuename);
    int atomsize   = res.getResidues()[0].atoms_.size(); 
    Indices.resize(atomsize);
    std::iota(Indices.begin(), Indices.end(), 1.0);

    pack_.ReadVectorNumber(IndicesName, ParameterPack::KeyType::Optional, Indices);
    Indices = Indices - 1;
}


CalculationTools::Real3 CalculationTools::getCOG(const Molecule::residue& residueGroup, const SimulationState& simstate, \
std::vector<int>& indices)
{
    auto& simbox = simstate.getSimulationBox();

    // For COM calculation, for each residue, we shift the atoms with respect to the first atom
    // obtain the position of the first atom in the residue group
    int index0 = indices[0];
    auto& pos1 = residueGroup.atoms_[index0].positions_;

    // Total mass of the atoms of interest
    Real chargeTot = 0;
    Real3 COM_pos = {{0,0,0}};

    // iterate over the indices of interest
    for (int j=0;j<indices.size();j++){
        Real3 distance;
        Real distsq;

        int index = indices[j];
        Real3 shiftWRTatom1 = simbox.calculateShift(residueGroup.atoms_[index].positions_, pos1);

        Real3 shiftedPos = residueGroup.atoms_[index].positions_ + shiftWRTatom1;
        COM_pos = COM_pos + shiftedPos;
    }

    COM_pos = COM_pos / indices.size();

    // check if it is outside of the box, if it is, then move it inside the box
    Real3 shift = simbox.calculateShift(COM_pos, simbox.getCenter());
    COM_pos = COM_pos + shift;

    return COM_pos;
}

void CalculationTools::CalculateUsrBetweenPair(const Molecule::residue& residue1, const Molecule::residue& residue2, \
                                 const SimulationState& simstate, Real r, \
                                 const std::vector<int>& indices1, const std::vector<int>& indices2, Real beta, \
                                 Real& attr, Real& repul, Real& total)
{
    attr=0.0;
    repul=0.0;
    total=0.0;

    Real f = 138.935485;

    for (int i=0;i<indices1.size();i++)
    {
        for (int j=0;j<indices2.size();j++)
        {
            Real qi = residue1.atoms_[i].charge_;
            Real qj = residue2.atoms_[j].charge_;

            Real qiqj = qi * qj;
            Real e = f * std::erfc(beta * r)/r * qiqj;

            if (qiqj < 0){attr += e;}
            if (qiqj > 0){repul += e;}
            total += e;
        }
    }
}

CalculationTools::INT3 CalculationTools::NearestLatticeIndex(const Real3& pos, const Real3& dL){
    Real3 index  = pos / dL;
    INT3 LatticeIndex;
    for (int i=0;i<3;i++){
        LatticeIndex[i]  = std::round(index[i]);
    }

    return LatticeIndex;
}

CalculationTools::INT3 CalculationTools::NearestLatticeIndex(const Real3& pos, const Real3& dL, const INT3& lattice_shape){
    INT3 nonPBCIndex = NearestLatticeIndex(pos, dL);

    return correctPBCLatticeIndex(nonPBCIndex, lattice_shape);
}

CalculationTools::Real CalculationTools::corase_grain_function(Real rsq, Real sigma){
    Real sigma2 = sigma*sigma;
    Real prefactor = std::pow(2*Constants::PI*sigma2, -1.5);

    return prefactor * std::exp(-rsq/(2*sigma2));
}

CalculationTools::INT3 CalculationTools::correctPBCLatticeIndex(const INT3& latticeIndex, const INT3& lattice_shape){
    INT3 pbcIndex;
    for (int i=0;i<3;i++){
        if (latticeIndex[i] < 0){pbcIndex[i] = latticeIndex[i] + lattice_shape[i];}
        else{pbcIndex[i] = latticeIndex[i] % lattice_shape[i];}
    }

    return pbcIndex;
}
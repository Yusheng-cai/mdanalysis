#include "Driver.h"

Driver::Driver(std::string filename, CommandLineArguments& cmd)
{
    // Check if path is provided, if so, then all the relative paths in the input file follows rpath
    apath_ =  FileSystem::getCurrentPath();
    cmd.readString("apath", CommandLineArguments::Keys::Optional,apath_);

    InputParser ip;
    ip.ParseFile(filename, pack_);

    // Read the xdr file inputted, this must be provided
    initializeXdr();

    // read the topology and gro file inputted, these are optional
    initializeGroFile();

    // initialize topology
    initializeTop();

    // initialize the driver, this is optional
    initializeDriverPack();

    // initialize atom groups
    initializeAtomGroups();

    // initialize residue group
    initializeResidueGroups();

    // initialize probe volume
    initializeProbeVolume();

    // initialize calculation
    initializeCalculation();

    // initialize the Order Parameters
    initializeOP();

    RegisterOuputValues();

    initializeOutputFiles();
}

void Driver::initializeResidueGroups(){
    auto resPack    = pack_.findParamPacks("residuegroup", ParameterPack::KeyType::Optional);
    VectorResNames_.clear();

    if (resPack.size() > 0){
        ASSERT((topology_read_), "Topology file is required if residue group were to be provided.");
    }

    for (int i=0;i<resPack.size();i++){
        const auto res = resPack[i];
        std::string resname;
        res -> ReadString("name", ParameterPack::KeyType::Required,resname);
        VectorResNames_.push_back(resname);

        ResidueInput input = {const_cast<ParameterPack&>(*res), grofile_, top_};

        ResidueGroup resgroup(input);
        
        simstate_.registerResidueGroup(resname, resgroup);
    }
}

void Driver::initializeCalculation(){
    auto calcPack  = pack_.findParamPacks("calculation", ParameterPack::KeyType::Optional);
    Calc_.clear();

    for (int i=0;i<calcPack.size();i++)
    {
        std::string type;
        calcPack[i] -> ReadString("type", ParameterPack::KeyType::Required, type);
        CalculationInput input = {const_cast<ParameterPack&>(*calcPack[i]),simstate_};

        Calc_.push_back(calcptr(CalculationRegistry::Factory::instance().create(type, input)));
    }
}

void Driver::initializeTop(){
    auto topPack    = pack_.findParamPack("topology", ParameterPack::KeyType::Optional);

    if (topPack != nullptr){
        topology_read_  = true;
        std::string topPath_;
        topPack->ReadString("path", ParameterPack::KeyType::Required, topPath_);
        
        std::string topAbsPath = FileSystem::joinPath(apath_, topPath_);

        top_.Parse(topAbsPath);
    }
}

bool Driver::isValidStep(int step)
{
    if (step >= startingFrame_){
        int a = (step - startingFrame_)%(skip_+1);
        if (a == 0){
            return true;
        }
        else{
            return false;
        }
    }

    return false;
}

void Driver::initializeDriverPack()
{
    auto driverpack = pack_.findParamPack("driver", ParameterPack::KeyType::Optional);

    if(driverpack != nullptr){
        driverpack -> ReadNumber("startingframe", ParameterPack::KeyType::Optional, startingFrame_);
        driverpack -> ReadNumber("skip", ParameterPack::KeyType::Optional, skip_);
    }

    int nframes = Xdr_ -> getNframes();
    ASSERT((nframes >= startingFrame_), "The starting frame specified " << startingFrame_ << " is more than the total number of frames " << nframes);

    int frames = startingFrame_;
    int numframes = 0;
    while (frames <= nframes){
        numframes ++;
        frames += skip_ + 1;
    }

    startingFrame_ -= 1;
    simstate_.setTotalFrames(numframes);
}

void Driver::initializeOutputFiles(){
    auto output_pack = pack_.findParamPacks("outputfile", ParameterPack::KeyType::Optional);

    for (int i=0;i<output_pack.size();i++){
        auto pack = output_pack[i];
        OutputFiles_.push_back(outputptr(new OutputStream(*pack, *this)));
    }
}

void Driver::initializeAtomGroups()
{
    auto agpack   = pack_.findParamPacks("atomgroup", ParameterPack::KeyType::Optional);

    for (int i=0;i<agpack.size();i++){
        auto pack = agpack[i];

        std::string ag_name;
        pack -> ReadString("name", ParameterPack::KeyType::Required,ag_name);

        AtomGroupInput agInput = { const_cast<ParameterPack&>(*pack), grofile_};
        AtomGroup ag(agInput);

        simstate_.registerAtomGroup(ag_name, ag);
        VectorAgNames_.push_back(ag_name);
    }
}

void Driver::RegisterOuputValues()
{
    for (int i=0; i<OP_.size();i++){
        auto registry = OP_[i]->getOutputRegistry();
        for (auto it = registry.begin(); it != registry.end();it++){
            std::string op_name = OP_[i]->getName();
            std::string output_name = it -> first;
            std::string full_name = op_name + "." + output_name;

            outputValueRegistry_.insert(full_name, it -> second);
        }
    }

    for (int i=0;i<Calc_.size();i++){
        auto registry = Calc_[i] -> getOutputRegistry();
        for (auto it = registry.begin() ;it != registry.end(); it++){
            std::string calc_name = Calc_[i]->getName();
            std::string output_name = it -> first;
            std::string full_name = calc_name + "." + output_name;

            outputValueRegistry_.insert(full_name, it -> second);
        }
    }
}

const OutputValue& Driver::getOutputValue(std::string name) const
{
    return outputValueRegistry_.find(name);
}

void Driver::initializeGroFile()
{
    auto gropack    = pack_.findParamPack("grofile", ParameterPack::KeyType::Optional);

    if (gropack != nullptr)
    {
        std::string gropath_;
        gropack->ReadString("path", ParameterPack::KeyType::Required,gropath_);

        std::string groAbsPath = FileSystem::joinPath(apath_, gropath_);

        grofile_.Open(groAbsPath);
    }
}

void Driver::initializeXdr()
{
    auto xdrpack    = pack_.findParamPack("xdrfile", ParameterPack::KeyType::Required);

    ASSERT((xdrpack != nullptr), "Xdr file such as trr or xtc file is not provided.");

    std::string path;
    std::string mode;

    // Read the path, name, mode of the file
    xdrpack->ReadString("path", ParameterPack::KeyType::Required, path);

    // split the path to obtain the type of xdr file we are working with, i.e test.xtc
    int found = path.find_last_of(".");
    std::string type = path.substr(found+1);

    // Create a pointer to the Xdr file that we are working with
    XdrInput input = { *xdrpack, apath_};
    Xdr_ = XdrPtr(XdrFiles::factory::instance().create(type, input));
    Xdr_->open();

    total_atom_positions_.reserve(Xdr_->getNumAtoms());
}

void Driver::initializeProbeVolume()
{
    auto PVpack = pack_.findParamPacks("probevolume", ParameterPack::KeyType::Optional);

    for (int i=0; i< PVpack.size();i++){
        auto pack = PVpack[i];
        std::string pv_type;
        std::string pv_name;
        pack->ReadString("type", ParameterPack::KeyType::Required, pv_type);
        pack->ReadString("name", ParameterPack::KeyType::Required, pv_name);

        ProbeVolumeNames_.push_back(pv_name);

        // TODO: change the PV input initialization
        ProbeVolumeInput PV_input{*pack, simstate_}; 
        auto pv_ptr = ProbeVolumes::Factory::instance().create(pv_type, PV_input);

        simstate_.registerProbeVolume(pv_name, pv_ptr);
    }
}

void Driver::initializeOP()
{
    auto OPpack = pack_.findParamPacks("orderparameter", ParameterPack::KeyType::Optional);

    for (int i=0;i<OPpack.size();i++){
        auto pack = OPpack[i];
        std::string op_type;

        pack -> ReadString("type", ParameterPack::KeyType::Required, op_type); 
        OrderParametersInput OP_input{*pack, simstate_};

        OPptr op_pointer = OPptr(OrderParametersRegistry::Factory::instance().create(op_type, OP_input));

        OP_.push_back(std::move(op_pointer));
    }
}

bool Driver::readFrame(int FrameNum)
{
    bool read = Xdr_->readFrame(FrameNum);

    // set the frame number in simulation state
    simstate_.setFrameNumber(FrameNum);

    if (read){return true;}
    else{return false;}
}

void Driver::update()
{
    // get the total positions
    const auto& total_atom_positions_ = Xdr_->getPositions();
    simstate_.setTotalNumberAtoms(total_atom_positions_.size());

    // update atom groups
    for (int i=0;i<VectorAgNames_.size();i++){
        auto& ag = simstate_.getAtomGroup(VectorAgNames_[i]);

        ag.update(total_atom_positions_);
    }

    // update residue groups
    for (int i=0;i<VectorResNames_.size();i++){
        auto& res = simstate_.getResidueGroup(VectorResNames_[i]);

        res.update(total_atom_positions_);
    }

    // update the simulation box
    auto& box = simstate_.getSimulationBox();
    simstate_.setSimulationBox(Xdr_->getSimulationBox());
    simstate_.setTime(Xdr_->getTime());
    simstate_.setStep(Xdr_->getStep());
    simstate_.setTotalAtomPos(total_atom_positions_);

    // update the probe volume
    for (int i=0;i<ProbeVolumeNames_.size();i++){
        auto& pv = simstate_.getProbeVolume(ProbeVolumeNames_[i]);

        pv.update();
    }

    // update the Order Parameters
    for (int i=0;i<OP_.size();i++){
        OP_[i] -> update();
    }

    // update the calculation
    for (int i=0;i<Calc_.size();i++){
        Calc_[i] -> update();
    }
}

void Driver::calculate()
{
    // perform OP calculations
    for (int i = 0;i<OP_.size();i++){
        auto& op = OP_[i];
        op ->calculate();
    }

    // perform calculations
    for (int i=0;i<Calc_.size();i++){
        Calc_[i] -> calculate();
    }

    // output calculations
    for (int i=0;i<Calc_.size();i++){
        Calc_[i] -> printOutputOnStep();
    }

    for (int i=0; i< OutputFiles_.size();i++){
        auto& out = OutputFiles_[i];
        out ->printIfOnStep();
    }
}

void Driver::finishCalculate()
{
    for (int i=0;i<Calc_.size();i++){
        Calc_[i] ->finishCalculate();
    }
}

void Driver::printOutput()
{
    for (int i=0;i<Calc_.size();i++){
        Calc_[i] -> printOutput();
    }
}
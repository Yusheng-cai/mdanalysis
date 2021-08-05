#include "Driver.h"

Driver::Driver(std::string filename, CommandLineArguments& cmd)
{
    // Check if path is provided, if so, then all the relative paths in the input file follows rpath
    apath_ =  FileSystem::getCurrentPath();
    cmd.readString("apath", CommandLineArguments::Keys::Optional,apath_);

    InputParser ip;
    ip.ParseFile(filename, pack_);

    auto pv_pack = pack_.findParamPacks("probevolume", ParameterPack::KeyType::Optional);
    auto op_pack = pack_.findParamPacks("orderparameter", ParameterPack::KeyType::Optional);

    // Make sure that there is only one xdr file
    auto xdr_pack= pack_.findParamPack("xdrfile", ParameterPack::KeyType::Required);
    auto ag_pack = pack_.findParamPacks("atomgroup", ParameterPack::KeyType::Optional);
    auto output_pack = pack_.findParamPacks("outputfile", ParameterPack::KeyType::Optional);
    auto gro_pack= pack_.findParamPack("grofile", ParameterPack::KeyType::Optional);

    // Read the xdr file inputted, this must be provided
    initializeXdr(xdr_pack);
    initializeGroFile(gro_pack);

    if (ag_pack.size() != 0)
    {
        initializeAtomGroups(ag_pack);
    }


    if (pv_pack.size() != 0)
    {
        initializeProbeVolume(pv_pack);
    }

    if (op_pack.size() != 0)
    {
        initializeOP(op_pack);
        // only register output values
        RegisterOuputValues();
    }

    if (output_pack.size() != 0)
    {
        initializeOutputFiles(output_pack);
    }
}

void Driver::initializeOutputFiles(const std::vector<const ParameterPack*>& output_pack)
{
    for (int i=0;i<output_pack.size();i++)
    {
        auto pack = output_pack[i];

        OutputFiles_.push_back(outputptr(new OutputStream(*pack, *this)));
    }
}

void Driver::initializeAtomGroups(const std::vector<const ParameterPack*>& agpack)
{
    for (int i=0;i<agpack.size();i++)
    {
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
    for (int i=0; i<OP_.size();i++)
    {
        auto registry = OP_[i]->getOutputRegistry();

        for (auto it = registry.begin(); it != registry.end();it++)
        {
            std::string op_name = OP_[i]->getName();
            std::string output_name = it -> first;
            std::string full_name = op_name + "." + output_name;

            outputValueRegistry_.insert(full_name, it -> second);
        }
    }
}

const OutputValue& Driver::getOutputValue(std::string name) const
{
    return outputValueRegistry_.find(name);
}

void Driver::initializeGroFile(const ParameterPack* gropack)
{
    if (gropack != nullptr)
    {
        std::string gropath_;
        gropack->ReadString("path", ParameterPack::KeyType::Required,gropath_);

        std::string groAbsPath = FileSystem::joinPath(apath_, gropath_);

        grofile_.Open(groAbsPath);
    }
}

void Driver::initializeXdr(const ParameterPack* xdrpack)
{
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

void Driver::initializeProbeVolume(const std::vector<const ParameterPack*>& PVpack)
{
    for (int i=0; i< PVpack.size();i++)
    {
        auto pack = PVpack[i];
        std::string pv_type;
        std::string pv_name;
        pack->ReadString("type", ParameterPack::KeyType::Required, pv_type);
        pack->ReadString("name", ParameterPack::KeyType::Required, pv_name);

        ProbeVolumeInput PV_input{*pack, simstate_}; 
        auto pv_ptr = ProbeVolumes::Factory::instance().create(pv_type, PV_input);

        pv_registry_.registerProbeVolume(pv_name, pv_ptr);
    }
}

void Driver::initializeOP(const std::vector<const ParameterPack*>& OPpack)
{
    for (int i=0;i<OPpack.size();i++)
    {
        auto pack = OPpack[i];
        std::string op_type;

        pack -> ReadString("type", ParameterPack::KeyType::Required, op_type); 
        OrderParametersInput OP_input{*pack, simstate_, pv_registry_};

        OPptr op_pointer = OPptr(OrderParametersRegistry::Factory::instance().create(op_type, OP_input));

        OP_.push_back(std::move(op_pointer));
    }
}

void Driver::update()
{
    bool read = Xdr_->readNextFrame();
    if (read)
    {
        is_Active_ = false;
    }

    const auto& total_atom_positions_ = Xdr_->getPositions();

    for (int i=0;i<VectorAgNames_.size();i++)
    {
        auto start = std::chrono::high_resolution_clock::now();
        auto& ag = simstate_.getAtomGroup(VectorAgNames_[i]);

        ag.update(total_atom_positions_);

        auto stop = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout << "Time it took for ag update is " << duration.count() << " microseconds" << std::endl;
    }


    // update the simulation box
    auto& box = simstate_.getSimulationBox();
    simstate_.setSimulationBox(Xdr_->getSimulationBox());
    simstate_.setTime(Xdr_->getTime());
    simstate_.setStep(Xdr_->getStep());

    // update the Order Parameters
    for (int i=0;i<OP_.size();i++)
    {
        auto& op = OP_[i];
        auto start = std::chrono::high_resolution_clock::now();
        op->update();
        auto stop = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout << "Time it took for op update is " << duration.count() << " microseconds" << std::endl;
    }
}

void Driver::calculate()
{
    for (int i = 0;i<OP_.size();i++)
    {
        auto& op = OP_[i];
        auto start = std::chrono::high_resolution_clock::now();
        op ->calculate();
        auto stop = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout << "Time it took for OP calculation is " << duration.count() << " microseconds" << std::endl;
    }

    auto start = std::chrono::high_resolution_clock::now();
    for (int i=0; i< OutputFiles_.size();i++)
    {
        auto& out = OutputFiles_[i];
        out ->printIfOnStep();
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
    std::cout << "Time it took for writing file is " << duration.count() << " microseconds" << std::endl;
}
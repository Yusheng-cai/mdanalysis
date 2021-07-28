#include "Driver.h"

Driver::Driver(std::string filename)
{
    InputParser ip;
    ip.ParseFile(filename, pack);

    auto pv_pack = pack.findParamPacks("probevolume", ParameterPack::KeyType::Optional);
    auto op_pack = pack.findParamPacks("orderparameter", ParameterPack::KeyType::Optional);
    auto xdr_pack= pack.findParamPacks("xdrfiles", ParameterPack::KeyType::Optional);

    initializeXdr(xdr_pack);
}

void Driver::initializeXdr(std::vector<const ParameterPack*> xdrpack)
{
    std::vector<std::string> names;
    for (int i=0;i<xdrpack.size();i++)
    {
        auto pack = xdrpack[i];

        std::string name;
        pack->ReadString("name", ParameterPack::KeyType::Required, name);
        names.push_back(name);

        int found = name.find_first_of(".");
        std::string type = name.substr(found+1);

        XdrPtr x = XdrPtr(XdrFiles::factory::instance().create(type));
        Xdr.push_back(x);

        // For each xdr file we have a corresponding simulation state
        statePtr sim_ptr(new SimulationState());
        sim_state_.push_back(sim_ptr);
    }

    for (int i=0;i<Xdr.size();i++)
    {
        auto pack = xdrpack[i];

        std::string mode;
        std::string name = names[i];
        bool read = pack->ReadString("mode", ParameterPack::KeyType::Optional, mode);

        if (read == false)
        {
            Xdr[i]->open(name, XdrWrapper::Mode::Read);
        }
        else
        {
            if (mode == "read")
            {
                Xdr[i] ->open(name, XdrWrapper::Mode::Read);
            }
            else if (mode == "append")
            {
                Xdr[i]->open(name, XdrWrapper::Mode::Append);
            }   
            else if (mode == "write")
            {
                Xdr[i]->open(name, XdrWrapper::Mode::Write);
            }
            else 
            {
                ASSERT((true == false), "Provided option " << mode << " is not in write/append/read.");
            }
        }
    }
}

void initializeProbeVolume(std::vector<const ParameterPack*> PVpack)
{
    for (int i=0; i< PVpack.size();i++)
    {
        auto pack = PVpack[i];
        std::string pv_type;
        pack->ReadString("type", ParameterPack::KeyType::Required, pv_type);
    }
}
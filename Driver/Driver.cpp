#include "Driver.h"

Driver::Driver(std::string filename)
{
    InputParser ip;
    ip.ParseFile(filename, pack_);

    auto pv_pack = pack_.findParamPacks("probevolume", ParameterPack::KeyType::Optional);
    auto op_pack = pack_.findParamPacks("orderparameter", ParameterPack::KeyType::Optional);
    auto xdr_pack= pack_.findParamPack("xdrfiles", ParameterPack::KeyType::Optional);

    if (xdr_pack != nullptr)
        initializeXdr(xdr_pack);
}

void Driver::initializeXdr(const ParameterPack* xdrpack)
{
    std::vector<std::string> paths;

    std::string path;
    std::string name;
    std::string mode;
    xdrpack->ReadString("path", ParameterPack::KeyType::Required, path);
    xdrpack->ReadString("name", ParameterPack::KeyType::Required, name);
    bool read = xdrpack->ReadString("mode", ParameterPack::KeyType::Optional, mode);

    // split the path to obtain the type of xdr file we are working with, i.e test.xtc
    int found = path.find_first_of(".");
    std::string type = path.substr(found+1);

    Xdr_ = XdrPtr(XdrFiles::factory::instance().create(type));
    SimState_ = statePtr(new SimulationState());

    // default is reading mode
    if (read == false)
    {
        Xdr_->open(path, XdrWrapper::Mode::Read);
    }
    else
    {
        if (mode == "read")
        {
            Xdr_ ->open(path, XdrWrapper::Mode::Read);
        }
        else if (mode == "append")
        {
            Xdr_->open(path, XdrWrapper::Mode::Append);
        }   
        else if (mode == "write")
        {
            Xdr_->open(path, XdrWrapper::Mode::Write);
        }
        else 
        {
            ASSERT((true == false), "Provided option " << mode << " is not in write/append/read.");
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
#include "tools/InputParser.h"
#include "xdr/XdrWrapper.h"

#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    std::cout << "Welcome to main." << std::endl;
    InputParser ip;
    std::string fname = argv[1];

    ParameterPack pp;
    ip.ParseFile(fname,pp);

    auto xtcpp = pp.findParamPack("xtc_file", ParameterPack::KeyType::Required);

    std::string xtcstr;
    std::string type;
    xtcpp->ReadString("type", ParameterPack::KeyType::Required,type);
    xtcpp->ReadString("name", ParameterPack::KeyType::Required, xtcstr);

    XdrWrapper* xx = XdrFiles::factory::instance().create(type);

    return 0;
};
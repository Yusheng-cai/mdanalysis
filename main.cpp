#include "tools/InputParser.h"
#include "tools/XdrWrapper.h"

#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    std::cout << "type" << std::endl;
    InputParser ip;
    std::string fname = argv[1];

    ParameterPack pp;
    ip.ParseFile(fname,pp);

    auto xtcpp = pp.findParamPack("xtc_file", ParameterPack::KeyType::Required);

    std::string xtcstr;
    std::string type;
    xtcpp->ReadString("type", ParameterPack::KeyType::Required,type);
    xtcpp->ReadString("name", ParameterPack::KeyType::Required, xtcstr);

    return 0;
};
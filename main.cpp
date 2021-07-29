#include "tools/InputParser.h"
#include "xdr/XdrWrapper.h"
#include "OrderParameters/AtomGroup.h"

#include <iostream>
#include <string>
#include <sstream>

int main(int argc, char** argv)
{
    InputParser ip;
    std::string fname = argv[1];

    // Parse the input file
    ParameterPack pp;
    ip.ParseFile(fname,pp);

    // Read the xdr files
    auto xdrpp = pp.findParamPack("xdrfile", ParameterPack::KeyType::Required);
    std::string xdrstr;
    std::string type;
    xdrpp->ReadString("name", ParameterPack::KeyType::Required, xdrstr);
    int found_pos = xdrstr.find_first_of(".");
    type = xdrstr.substr(found_pos + 1);

    XdrWrapper* xx = XdrFiles::factory::instance().create(type);
    xx->open(xdrstr, XdrWrapper::Mode::Read);
    xx->readNextFrame();
    auto& pos = xx->getPositions();

    for (int i=0;i<xx->getNumAtoms();i++)
    {
        for(int j=0;j<3;j++)
        {
            std::cout << pos[i][j] << "\t";
        }
        std::cout << "\n";
    }
    std::cout << "Read sucessful." << std::endl;

    auto agg = pp.findParamPack("atomgroup", ParameterPack::KeyType::Required);
    AtomGroup ag_(*agg);
    return 0;
};
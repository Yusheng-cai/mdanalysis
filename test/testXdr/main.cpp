#include "xdr/TopologyReader.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"

#include <string>
#include <vector>
#include <array>

int main(int argc, char** argv)
{
    std::string fname(argv[1]);

    TopologyReader top;
    top.Parse(fname);
}
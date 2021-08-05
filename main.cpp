#include "OrderParameters/Driver.h"
#include "tools/CommandLineArguments.h"

#include <iostream>
#include <string>
#include <sstream>

int main(int argc, char** argv)
{
    std::string fname="";

    CommandLineArguments cmd(argc, argv);
    if (cmd.has_key("f"))
    {
        cmd.readString("f", CommandLineArguments::Keys::Required,fname);
    }
    else
    {
        ASSERT((argc > 1), "There must be an input file.");
        fname = argv[1];
    }

    ASSERT((!fname.empty()), "Missing input file.");

    Driver d(fname,cmd);

    auto start = std::chrono::high_resolution_clock::now();
    for(int i=0;i<d.getNframes();i++)
    {
        d.update();
        d.calculate();
    }
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Entire operation took " << duration.count() << " microseconds"<< std::endl;
 
    return 0;
};
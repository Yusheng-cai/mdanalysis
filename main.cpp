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
    int step = 0;

    auto start = std::chrono::high_resolution_clock::now();
    while (! d.readNextFrame())
    {
        d.update();
        d.calculate();
        step ++;
        std::cout << "step = " << step << std::endl;
    }
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Entire operation took " << duration.count() << " microseconds"<< std::endl;
 
    return 0;
};
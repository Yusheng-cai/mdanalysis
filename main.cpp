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

    auto starttot = std::chrono::high_resolution_clock::now();
    for (int i=0;i<d.getNframes();i++)
    {
        if (d.isValidStep(i))
        {
            d.readFrame(i);
            d.update(); 
            d.calculate();
        }
    }
    auto endtot = std::chrono::high_resolution_clock::now();

    auto durationtot = std::chrono::duration_cast<std::chrono::microseconds>(endtot- starttot);

    std::cout << "Entire operation took " << durationtot.count() << " microseconds"<< std::endl;

    d.finishCalculate();
 
    return 0;
};
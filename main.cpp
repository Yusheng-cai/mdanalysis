#include "OrderParameters/Driver.h"

#include <iostream>
#include <string>
#include <sstream>

int main(int argc, char** argv)
{
    std::string fname = argv[1];

    Driver d(fname);

    for(int i=0;i<d.getNframes();i++)
    {
        d.update();
        d.calculate();
    }
 
    return 0;
};
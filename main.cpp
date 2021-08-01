#include "OrderParameters/Driver.h"

#include <iostream>
#include <string>
#include <sstream>

int main(int argc, char** argv)
{
    std::string fname = argv[1];

    Driver d(fname);

    while (d.isActive())
    {
        d.update();
        d.calculate();
    }
 
    return 0;
};
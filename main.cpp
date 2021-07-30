#include "OrderParameters/Driver.h"

#include <iostream>
#include <string>
#include <sstream>

int main(int argc, char** argv)
{
    std::string fname = argv[1];

    Driver d(fname);

    d.update();
 
    return 0;
};
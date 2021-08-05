#include "tools/CommandLineArguments.h"
#include "tools/Assert.h"

#include <vector>
#include <string>

int main(int argc, char** argv)
{
    CommandLineArguments cmd(argc, argv);

    cmd.print();

    return 0;
}
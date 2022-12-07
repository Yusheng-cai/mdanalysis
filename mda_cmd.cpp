#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <map>
#include <functional>

#include "tools/CommandLineArguments.h"
#include "tools/CommonTypes.h"
#include "tools/Assert.h"

#include "Actions/mda_actions.hpp"

using ActionFunction = std::function<void(CommandLineArguments& cmd)>;
using mapFunction    = std::map<std::string, ActionFunction>;
using mapUsage       = std::map<std::string, std::string>;

void RegisterAction(std::string name, std::string usage, ActionFunction func, mapFunction& mapF, mapUsage& mapU);
void RegisterAllActions(mapFunction& mapF, mapUsage& mapU);

int main(int argc, char** argv){
    CommandLineArguments cmd(argc, argv);

    std::string help_key_ = "help";
    std::string operation_;
    mapFunction MapNameToAction;
    mapUsage    MapNameToUsage;

    RegisterAllActions(MapNameToAction, MapNameToUsage);

    if (cmd.has_key("help") && ! cmd.has_key("op")){
        std::cout << "Usage: mesh_op -op OPERATION " << "\n";
        std::cout << "operation  Usage" << "\n";
        for (auto it = MapNameToUsage.begin(); it != MapNameToUsage.end(); it ++){
            std::cout << it -> first << "\t" << it -> second << "\n";
        }
        return 0;
    }


    // find the operation
    cmd.readString("op",CommandLineArguments::Keys::Required,operation_);

    if (cmd.has_key("help") && cmd.has_key("op")){
        auto it = MapNameToUsage.find(operation_);
        std::cout << it -> first << "\t" << it -> second << "\n";

        return 0;
    }

    // perform the operation
    auto it = MapNameToAction.find(operation_);
    ASSERT((it != MapNameToAction.end()), "Operation " << operation_ << " is not registered.");

    it -> second(cmd);

    return 0;
}

void RegisterAction(std::string name, std::string usage, ActionFunction func, mapFunction& mapF, mapUsage& mapU)
{
    auto Fit = mapF.find(name);
    ASSERT((Fit == mapF.end()), "The operation " << name << " is registered more than once.");

    auto Uit = mapU.find(name);
    ASSERT((Uit == mapU.end()), "The operation " << name << " is registered more than once.");

    mapF.insert(std::make_pair(name,func));
    mapU.insert(std::make_pair(name, usage));
}

void RegisterAllActions(mapFunction& mapF, mapUsage& mapU)
{
    RegisterAction("NPgen", "-i input.out", \
                            [](CommandLineArguments& cmd) -> void {mda_actions::generateNP(cmd);}, 
                            mapF, mapU);
}

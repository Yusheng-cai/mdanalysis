#include "mda_actions.hpp"

void mda_actions::generateNP(CommandLineArguments& cmd){
    std::string inputfname, outputfname="sulfur.out";

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    std::string ext = StringTools::ReadFileExtension(inputfname);
    std::vector<Real3> pos;

    if (ext == "gro"){
        GroFile gro;
        gro.Open(inputfname);
        pos = gro.getPosition();
    }
    else{
        std::vector<std::vector<Real>> data;
        std::vector<Real3> pos;
        StringTools::ReadTabulatedData(inputfname, data);
        ASSERT((data.size() > 0), "The positions of the nanoparticle must be larger than to 0");

        pos.resize(data.size());
        for (int i=0;i<data.size();i++){
            ASSERT((data[i].size() == 3), "The positions must be 3d.");
            Real3 p = {{data[i][0], data[i][1], data[i][2]}};
            pos[i] = p;
        }
    }

    NanoparticleGeneration NPgen(pos);
    NPgen.Generate();
    NPgen.writeGroFile(outputfname);
}

void mda_actions::FindSurfaceSiO2(CommandLineArguments& cmd){
    std::string inputfname, outputfname="index.out";

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    std::string ext = StringTools::ReadFileExtension(inputfname);
    ASSERT((ext == "gro"), "Input must be a gro file.");
}

void mda_actions::TileCrystals(CommandLineArguments& cmd){
    std::string inputfname;
    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);

    Crystal c(inputfname);
    c.calculate();
    c.writeOutput("test.gro");
}

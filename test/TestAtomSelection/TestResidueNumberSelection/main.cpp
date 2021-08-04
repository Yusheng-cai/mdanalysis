#include "OrderParameters/AtomGroup.h"
#include "tools/Assert.h"
#include "tools/InputParser.h"
#include "tools/CommandLineArguments.h"

#include <string>
#include <vector>

void run(std::string fname, CommandLineArguments& cmd);

int main(int argc, char** argv)
{
    // read in the argc & argv
    CommandLineArguments cmd(argc, argv);

    std::string fname(argv[1]); 

    run(fname, cmd);

    return 0;
}

void run(std::string fname, CommandLineArguments& cmd)
{
    InputParser ip;
    ParameterPack pack;

    ASSERT((cmd.get_num_keys() == 2), "Two keys must be provided, -f for input file name and -t for test folder")
    std::string fpath;
    std::string tpath;

    // read the file name + its path
    cmd.readString("f", CommandLineArguments::Keys::Required, fpath);
    // This tells you where you can find the test folder that contains all the data resources
    cmd.readString("t", CommandLineArguments::Keys::Required, tpath);
    
    // parse the file name using the file path
    ip.ParseFile(fpath, pack);
    auto atomPack = pack.findParamPacks("atomgroup", ParameterPack::KeyType::Required);
    auto groPack  = pack.findParamPack("grofile", ParameterPack::KeyType::Required);

    std::string groPath;
    groPack->ReadString("path", ParameterPack::KeyType::Required, groPath);
    groPath = tpath + "/" + groPath;

    GroFile grofile;
    grofile.Open(groPath);
    
    // construct the correct answer
    std::vector<std::vector<int>> correct_answers;

    // For this one we are doing 1-500
    std::vector<int> ans1;
    int begin = 1 - 1;
    int end = 9500 - 1;
    int skip = 1;
    while (begin <= end)
    {
        ans1.push_back(begin);
        begin += skip;
    }
    correct_answers.push_back(ans1);

    // 1 - 500:2
    std::vector<int> ans2;
    begin = 1 - 1;
    end = 9500 - 1;
    skip = 1;
    while ( begin <= end)
    {
        ans2.push_back(begin);
        begin += skip;
    }
    correct_answers.push_back(ans2);

    // 1 , 2 , 3
    std::vector<int> ans3 = { 0, 1, 2};
    correct_answers.push_back(ans3);

    // 1 - 1000 , 3000 - 4000:2
    std::vector<int> ans4;
    begin = 1 - 1;
    end = 1000 -1;
    skip = 1;

    while(begin <= end)
    {
        ans4.push_back(begin);
        begin+=skip;
    }
    

    begin = 3000 - 1;
    end = 4000 - 1;
    skip = 2;
    while (begin <= end)
    {
        ans4.push_back(begin);
        begin += skip;
    }

    correct_answers.push_back(ans4);

    ASSERT((atomPack.size() == 4), "We are only doing 4 tests in this file.");
    ASSERT((correct_answers.size() == 4), "There aren't enough correct answers to be compared to.");

    for (int i =0;i<atomPack.size();i++)
    {
        auto pack = atomPack[i];
        AtomGroupInput input = {const_cast<ParameterPack&>(*pack), grofile};
        AtomGroup ag(input);

        auto GroupIndices = ag.getAtomGroupIndices();
        for (int i=0; i< GroupIndices.size();i++)
        {
            std::cout << GroupIndices[i] << std::endl;
        }

        auto correct_ans = correct_answers[i];

        ASSERT((GroupIndices == correct_ans), "The answer for " << i << " does not match.");
    }
}
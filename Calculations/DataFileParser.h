#pragma once
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

class DataFileParser
{
    public:
        DataFileParser(){};
        template<typename T>
        void ParseFile(std::string& name, std::vector<std::vector<T>>& input);
    
    private:
};

template<typename T>
void DataFileParser::ParseFile(std::string& name, std::vector<std::vector<T>>& input)
{
    std::ifstream file;
    std::stringstream ss;

    // open the file
    file.open(name);

    ASSERT((file.is_open()), "The file with name " << name << " is not opened.");

    std::string sentence;

    while (std::getline(file, sentence))
    {
        // Find the comment symbol
        int found = sentence.find_first_of("#");

        if (! sentence.empty() && found == std::string::npos)
        {
            T num;
            std::vector<T> vec;

            ss.str(sentence);
            while (ss >> num)
            {
                vec.push_back(num);
            }

            input.push_back(vec);

            ss.clear();
        }
    }

    file.close();
}
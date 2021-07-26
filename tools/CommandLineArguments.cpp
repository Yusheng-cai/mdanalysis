#include "CommandLineArguments.h"

CommandLineArguments::CommandLineArguments(int argc, char** argv)
:argc_(argc)
{
    argv_.clear();
    for (int i=1; i< argc_; i++)
    {
        argv_.push_back(std::string(argv[i]));
    }

    ParseArgv();
}

void CommandLineArguments::ParseArgv()
{
    // Keeps track of the key indices 
    std::vector<int> key_indices;

    // Read the key indices 
    for (int i=0;i<argv_.size();i++)
    {
        bool isNum = StringTools::isNumber(argv_[i]);

        // First check if the argv is empty, then check if the first char is "-", then check if it is a number
        if (argv_[i].size() > 0 && argv_[i][0] == '-' &&  (! isNum))
        {
            key_indices.push_back(i);
        }
    }

    num_keys_ = key_indices.size();

    // Read the actual keys 
    for (int i=0; i< num_keys_;i++)
    {
        // Actual index of the key in the argv_
        int key_idx  = key_indices[i];

        // Find the position of the "-" in the key string
        std::size_t dash_pos = argv_[key_idx].find_first_not_of("-");

        // Find the actual key
        std::string key = argv_[key_idx].substr(dash_pos);

        // Create insert the item into the map
        std::vector<std::string> a_;
        auto args_it = args_map_.find(key);
        ASSERT((args_it == args_map_.end()), "The argument " << key << " is repeated.");

        args_map_.insert(std::pair<std::string, std::vector<std::string>>(key, a_));

        // begin is the next charater after the one with "-", end is either end of string or the next "-"
        int begin = key_idx + 1;
        int end;

        if ( i == num_keys_ - 1)
        {
            end = argc_ - 1;
        }
        else
        {
            end = key_indices[i+1];
        }

        auto it_begin  = argv_.begin() + begin;
        auto it_end    = argv_.begin() + end;
        auto& this_arg_ = args_map_.find(key)->second; 

        // Assign the values of interest to the vector
        this_arg_.assign(it_begin, it_end);
    }
}

void CommandLineArguments::print() const 
{
    std::cout << "Key" << " " << "Values" << std::endl;
    for (auto it: args_map_)
    {
        std::string key = it.first; 
        std::vector<std::string> val = it.second;

        std::cout << key << " ";

        for (int i=0;i<val.size();i++)
        {
            std::cout << val[i] << " ";
        }

        std::cout << "\n";
    }
}

bool CommandLineArguments::has_key(const std::string& key) const
{
    auto it = args_map_.find(key);

    if (it != args_map_.end())
    {
        return true;
    }
    else{
        return false;
    }
}
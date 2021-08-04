#include "InputParser.h"

void StringTools::to_lower(std::string& str)
{
    std::for_each(str.begin(), str.end(), [](char& c){c = std::tolower(c);});
}

bool StringTools::isNumber(std::string str)
{
    std::stringstream ss(str);

    float a;
    ss >> a;

    return (! ss.fail());
}

void StringTools::RemoveBlankInString(std::string& str)
{
    str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
}

                            ////// Parameter packs ///////// 
std::string& ParameterPack::insert(const std::string& key, const std::string& value)
{
    auto new_it = value_.insert(std::make_pair(key,value));

    return new_it->second;
}

std::vector<std::string>& ParameterPack::insert(const std::string& key, const std::vector<std::string>& value)
{
    auto new_it = vectors_.insert(std::make_pair(key,value));

    return new_it->second;
}

ParameterPack& ParameterPack::insert(const std::string& key, const ParameterPack& parampack)
{
    auto new_it = parampacks_.insert(std::make_pair(key, parampack));

    return new_it->second;
}

const std::string* ParameterPack::findValue(const std::string& key, KeyType keytype) const
{
    std::vector<const std::string*> vec_str = findValues(key, keytype);
    
    if(keytype == ParameterPack::KeyType::Required)
    {
        // If the size is equal to 0, it would have been caught in findValues
        ASSERT((vec_str.size() == 1), "There's multiple definition for key: " << key << " ");
    }

    if (vec_str.size() == 0)
    {
        return nullptr;
    }
    else
    {
        return vec_str[0];
    }

}


std::vector<const std::string*> ParameterPack::findValues(const std::string& key, KeyType keytype) const
{
    auto mapit = value_.equal_range(key);

    int num_matches = std::distance(mapit.first, mapit.second);

    if (keytype == ParameterPack::KeyType::Required)
    {
        ASSERT((num_matches), "The required value for key: '" << key << "' is not provided.");
    }

    std::vector<const std::string*> ret;
    for(auto it = mapit.first;it != mapit.second;it++)
    {
        ret.push_back(&(it->second));
    }

    return ret;
}

std::vector<const std::vector<std::string>*> ParameterPack::findVectors(const std::string& key, KeyType keytype) const
{
    auto mapit = this->vectors_.equal_range(key);

    int num_matches = std::distance(mapit.first, mapit.second);

    if (num_matches == 0)
    {
        ASSERT((keytype != ParameterPack::KeyType::Required), "The required key for vector is not provided.");
    }

    std::vector<const std::vector<std::string>*> vector_ptrs;

    for ( auto it = mapit.first; it != mapit.second; ++it ) {
        vector_ptrs.push_back(&(it->second));
    }

    return vector_ptrs;
}

const std::vector<std::string>* ParameterPack::findVector(const std::string& key, const KeyType keytype) const
{
    std::vector<const std::vector<std::string>*> vecptrs = findVectors(key, keytype);

     if(keytype == ParameterPack::KeyType::Required)
    {
        // If the size is equal to 0, it would have been caught in findValues
        ASSERT((vecptrs.size() == 1), "There's multiple definition for key: " << key << " ");
    }

    if (vecptrs.size() == 0)
    {
        return nullptr;
    }
    else
    {
        return vecptrs[0];
    }

}

std::vector<const ParameterPack*> ParameterPack::findParamPacks(const std::string& key, const KeyType keytype) const
{
    auto mapit = this->parampacks_.equal_range(key);

    int num_matches = std::distance(mapit.first, mapit.second);

    if (num_matches == 0)
    {
        ASSERT((keytype != ParameterPack::KeyType::Required), "The required key for ParamPack is not provided.");
    }

    std::vector<const ParameterPack*> vec_param;

    for (auto it = mapit.first; it != mapit.second; ++it)
    {
        vec_param.push_back(&(it -> second));
    }

    return vec_param;
}

const ParameterPack* ParameterPack::findParamPack(const std::string& key, const KeyType keytype) const
{
    auto packs = findParamPacks(key, keytype);

    if (keytype == ParameterPack::KeyType::Required)
    {
        ASSERT((packs.size() == 1), "Multiple instance found for parameter pack with name '"\
        << key << "'");
    }

    if (packs.size() == 0)
    {
        return nullptr;
    }
    else
    {
        return packs[0];
    }
}


bool ParameterPack::ReadString(const std::string& key, const KeyType keytype, std::string& str) const
{
    const std::string* strPtr = findValue(key, keytype);

    // This should only evaluate if it's an optional keytype & not found 
    if(strPtr != nullptr)
    {
        str = *strPtr;
        return true;
    }

    return false;
}

bool ParameterPack::Readbool(const std::string& key, const KeyType keytype, bool& boolean) const
{
    std::string temp_str;

    bool read = ReadString(key, keytype, temp_str);
    // StringTools::to_lower(temp_str);

    if (read == true)
    {
        ASSERT((temp_str.compare("true") == 0 || temp_str.compare("false") == 0), \
        "Failed to read a boolean as the input value is " << temp_str << " for key " << key);

        if (temp_str.compare("true") == 0)
        {
            boolean  = true;
            return true;
        }

        if (temp_str.compare("false") == 0)
        {
            boolean = false;
            return true;
        }
    }

    return false;
}

bool ParameterPack::ReadVectorString(const std::string& key, const KeyType keytype,std::vector<std::string>& vecstr) const
{
    vecstr.clear();
    auto strvec = findVector(key, keytype);

    if (strvec != nullptr)
    {
        for (int i=0;i< strvec->size();i++)
        {
            std::string str = strvec->at(i);
            // StringTools::to_lower(str);
            vecstr.push_back(str);
        }
        return true;
    }

    return false;
}

                                        //// TokenStream ////
TokenStream::Status TokenStream::ReadNextToken(std::string& token)
{
    // clear the string for safety 
    token.clear();
    std::string line;

    // try if we can read in a token from the line_stream_
    if (line_stream_ >> token)
    {
        // If the first token of the line_stream_ does not contain "#", it is a valid token
        if (token.find_first_of(comment_str_) == std::string::npos)
        {
            return TokenStream::Status::Success;
        }
        // else we read another line
        else
        {
            line_stream_.str("");
            line_stream_.clear();

            return this ->ReadNextToken(token);
        }
    }
    // elseif we can read a line
    else if(std::getline(ifstream_,line))
    {
        line_stream_.clear(); // reset the ss state
        line_stream_.str(line);
        return this -> ReadNextToken(token);
    }
    // See if EOF
    else if(ifstream_.eof())
    {
        return Status::EndOfFile;
    }
    // Hopefully it doesn't get here
    else
    {
        return Status::Failure;
    }
}

                    ///////// Input Parser ///////////
void InputParser::ParseFile(const std::string& filename, ParameterPack& parampack)
{
    // instantiate the ifstream as well as the tokenstream 
    std::ifstream ifs(filename);

    // First Assert that the file has to be opened
    ASSERT(ifs.is_open(), "File " << filename << " is not opened.");

    // instantiate the TokenStream
    TokenStream toks(ifs);
    TokenStream::Status status;

    while(1)
    {
        status = ParseNextToken(toks, parampack); 
        if(status == TokenStream::Status::EndOfFile)
        {
            break;
        }
    }
}

TokenStream::Status InputParser::ParseNextToken(TokenStream& toks, ParameterPack& parampack)
{
    TokenStream::Status status;

    std::string key, delimiter, value;

    // Read in the key 
    status = toks.ReadNextToken(key);
    ASSERT((status != TokenStream::Status::Failure), "Reading token failed.");
    if (key == "}")
    {
        status = TokenStream::Status::Close_Brace;
        return status;
    }
    else if (status == TokenStream::Status::EndOfFile)
    {
        status = TokenStream::Status::EndOfFile;
        return status;
    }
    else if (status == TokenStream::Status::Close_bracket)
    {
        status = TokenStream::Status::Close_bracket;
        return status;
    }

    // Read in the delimiter
    status = toks.ReadNextToken(delimiter);
    ASSERT((status != TokenStream::Status::Failure), "Reading delimiter failed.");
    ASSERT((delimiter == "="), "Delimiter missing.");

    // Read in the value
    status = toks.ReadNextToken(value);
    ASSERT((status != TokenStream::Status::Failure), "Reading value failed.");
    if (value == "{")
    {
        // StringTools::to_lower(key);
        auto& new_param = parampack.insert(key, ParameterPack(key));
        status = ParseParamPack(toks, new_param);

        ASSERT((status == TokenStream::Status::Close_Brace), "Missing Ending Brace.");

        return status;
    }
    else if (value == "[")
    {
        std::vector<std::string> vecvals;
        status = ParseVector(toks, vecvals);

        ASSERT((status == TokenStream::Status::Close_bracket), "Missing Ending bracket");

        // make the key lower case too allow for flexibility
        // StringTools::to_lower(key);
        parampack.insert(key, vecvals);
    }
    else
    {
        //StringTools::to_lower(key);
        //StringTools::to_lower(value);
        parampack.insert(key,value);
    }

    return TokenStream::Status::Success;
}

TokenStream::Status InputParser::ParseParamPack(TokenStream& toks, ParameterPack& parampack)
{
    std::string token;
    TokenStream::Status status;

    while(1)
    { 
        // Read the next Token 
        status = ParseNextToken(toks, parampack);

        if (status == TokenStream::Status::Close_Brace || status == TokenStream::EndOfFile)
        {
            break;
        }
    }

    return status;
}

TokenStream::Status InputParser::ParseVector(TokenStream& toks, std::vector<std::string>& vecval)
{
    vecval.clear();

    TokenStream::Status status;

    while(1)
    { 
        std::string token;

        // Read the next Token 
        status = toks.ReadNextToken(token);
        
        if (token == "]")
        {
            status = TokenStream::Close_bracket;
            break;
        }
        else if (status == TokenStream::EndOfFile)
        {
            break;
        }

        vecval.push_back(token);
    }

    return status;
}
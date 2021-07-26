#pragma once
#include <string>
#include <iostream>

#include <exception>
#include <sstream>
#include <string>
#include <stdexcept>
#include <map>

#define ASSERT(test,message) if (! test)  \
    {                                     \
    std::stringstream err_s;              \
    err_s << std::string(__FILE__) << ": "\
    << std::to_string(__LINE__) << " in " \
    << std::string(__PRETTY_FUNCTION__) <<\
    " and " << message << "\n";        \
    throw std::runtime_error(err_s.str());}

template<class Key, class Value>
void assertIfInMap(std::map<Key,Value>& map, std::string str, std::string msg)
{
    typename std::map<Key, Value>::iterator it = map.find(str);
    ASSERT((it != map.end()), msg);
}


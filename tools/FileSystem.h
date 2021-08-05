#pragma once

#include <iostream>
#include <filesystem>
#include <stdio.h>
#include <unistd.h>
#include <string>

namespace FileSystem
{
    std::string getCurrentPath(); 

    std::string joinPath(const std::string& path1, const std::string& path2);
};
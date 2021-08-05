#include "FileSystem.h"

std::string FileSystem::getCurrentPath()
{
    char temp[256];
    std::string path =  getcwd(temp, sizeof(temp));

    return path;
} 

std::string FileSystem::joinPath(const std::string& path1, const std::string& path2)
{
    std::string joinedPath;

    std::string path1_copy(path1);
    std::string path2_copy(path2);

    int found1 = path1_copy.find_last_of("/");
    int found2 = path2_copy.find_first_of("/");

    if (found1 != std::string::npos)
    {
        if (found1 == path1_copy.size() - 1)
        {
            path1_copy.erase(found1,1);
        }
    }

    if (found2 != std::string::npos)
    {
        if (found2 == 0)
        {
            path2_copy.erase(0,1);
        }
    }

    joinedPath = path1_copy + "/" + path2_copy;

    return joinedPath;
}
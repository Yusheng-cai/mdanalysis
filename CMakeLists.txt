cmake_minimum_required(VERSION 3.18.1)

project(mdanalysis)

# Require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
#set(CMAKE_CXX_FLAGS_DEBUG_INIT "-O0 -Wall -g")
#set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O3 -Wall")

add_subdirectory(tools)

add_executable(analysis main.cpp)

target_include_directories(analysis PUBLIC ${CMAKE_SOURCE_DIR})
target_link_libraries(analysis PUBLIC tools)
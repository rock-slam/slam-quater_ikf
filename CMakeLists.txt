# CMakeLists.txt has to be located in the project folder and cmake has to be
# executed from 'project/build' with 'cmake ../'.
cmake_minimum_required(VERSION 2.6)
if(ROCK_CMAKE_MACROS)
else()
    include($ENV{ROCK_CMAKE_MACROS}/Rock.cmake)
endif()
find_package(Rock)
rock_init(quater_ikf 0.1)
rock_standard_layout()

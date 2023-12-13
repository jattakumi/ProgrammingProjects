cmake_minimum_required(VERSION 3.25)
project(Project_04)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)
add_executable(Project_04 hartreefock.cpp molecule.cpp main.cpp mp2.h mp2.cpp)

target_link_libraries(Project_04 Eigen3::Eigen)
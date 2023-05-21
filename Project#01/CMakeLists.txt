cmake_minimum_required(VERSION 3.25)
project(Project_01)

set(CMAKE_CXX_STANDARD 17)

# Set the path to the Eigen library
find_package(Eigen3 REQUIRED)

add_executable(Project_01 main.cpp molecule.cpp)

#Include the Eigen Directory
target_link_libraries(Project_01 Eigen3::Eigen)


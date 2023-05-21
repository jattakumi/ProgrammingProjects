cmake_minimum_required(VERSION 3.25)
project(Project_02)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)
add_executable(Project_02 hessian.cpp molecule.cpp projs.cpp)

target_link_libraries(Project_02 Eigen3::Eigen)
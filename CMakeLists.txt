cmake_minimum_required(VERSION 2.8)
project(crawdad_tutorials CXX)

add_library(molecule diag.cc molecule.hpp molecule.cpp)

add_executable(project1_acetaldehyde project1_acetaldehyde.cpp)
add_executable(project1_benzene project1_benzene.cpp)
add_executable(project1_allene project1_allene.cpp)
add_executable(project2_h2o project2_h2o.cpp)

target_link_libraries(project1_acetaldehyde molecule)
target_link_libraries(project1_benzene molecule)
target_link_libraries(project1_allene molecule)
target_link_libraries(project2_h2o molecule)
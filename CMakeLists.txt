cmake_minimum_required(VERSION 2.8)
project(crawdad_tutorials CXX)

set(CMAKE_BUILD_TYPE debug)

add_library(molecule 
  diag.cpp 
  utils.hpp 
  utils.cpp 
  molecule.hpp 
  molecule.cpp
  )

add_custom_command(
  OUTPUT "acetaldehyde.txt"
  COMMAND "${CMAKE_COMMAND}" -E copy 
  "${crawdad_tutorials_SOURCE_DIR}/acetaldehyde.txt"
  "${crawdad_tutorials_BINARY_DIR}/acetaldehyde.txt"
  DEPENDS "${crawdad_tutorials_SOURCE_DIR}/acetaldehyde.txt"
  )

add_custom_command(
  OUTPUT "benzene.txt"
  COMMAND "${CMAKE_COMMAND}" -E copy 
  "${crawdad_tutorials_SOURCE_DIR}/benzene.txt"
  "${crawdad_tutorials_BINARY_DIR}/benzene.txt"
  DEPENDS "${crawdad_tutorials_SOURCE_DIR}/benzene.txt"
  )

add_custom_command(
  OUTPUT "allene.txt"
  COMMAND "${CMAKE_COMMAND}" -E copy 
  "${crawdad_tutorials_SOURCE_DIR}/allene.txt"
  "${crawdad_tutorials_BINARY_DIR}/allene.txt"
  DEPENDS "${crawdad_tutorials_SOURCE_DIR}/allene.txt"
  )

add_custom_command(
  OUTPUT "h2o_geom.txt"
  COMMAND "${CMAKE_COMMAND}" -E copy 
  "${crawdad_tutorials_SOURCE_DIR}/h2o_geom.txt"
  "${crawdad_tutorials_BINARY_DIR}/h2o_geom.txt"
  DEPENDS "${crawdad_tutorials_SOURCE_DIR}/h2o_geom.txt"
  )

add_custom_command(
  OUTPUT "h2o_hessian.txt"
  COMMAND "${CMAKE_COMMAND}" -E copy 
  "${crawdad_tutorials_SOURCE_DIR}/h2o_hessian.txt"
  "${crawdad_tutorials_BINARY_DIR}/h2o_hessian.txt"
  DEPENDS "${crawdad_tutorials_SOURCE_DIR}/h2o_hessian.txt"
  )

add_custom_command(
  OUTPUT 
  "h2o_sto3g_enuc.dat"
  "h2o_sto3g_s.dat"
  "h2o_sto3g_t.dat"
  "h2o_sto3g_v.dat"
  COMMAND "${CMAKE_COMMAND}" -E copy "${crawdad_tutorials_SOURCE_DIR}/h2o_sto3g_enuc.dat" "${crawdad_tutorials_BINARY_DIR}/h2o_sto3g_enuc.dat"
  COMMAND "${CMAKE_COMMAND}" -E copy "${crawdad_tutorials_SOURCE_DIR}/h2o_sto3g_s.dat" "${crawdad_tutorials_BINARY_DIR}/h2o_sto3g_s.dat"
  COMMAND "${CMAKE_COMMAND}" -E copy "${crawdad_tutorials_SOURCE_DIR}/h2o_sto3g_t.dat" "${crawdad_tutorials_BINARY_DIR}/h2o_sto3g_t.dat"
  COMMAND "${CMAKE_COMMAND}" -E copy "${crawdad_tutorials_SOURCE_DIR}/h2o_sto3g_v.dat" "${crawdad_tutorials_BINARY_DIR}/h2o_sto3g_v.dat"
  DEPENDS
  "${crawdad_tutorials_SOURCE_DIR}/h2o_sto3g_enuc.dat"
  "${crawdad_tutorials_SOURCE_DIR}/h2o_sto3g_s.dat"
  "${crawdad_tutorials_SOURCE_DIR}/h2o_sto3g_t.dat"
  "${crawdad_tutorials_SOURCE_DIR}/h2o_sto3g_v.dat"
  )

add_custom_target(inp_p1_acetaldehyde.txt DEPENDS "${crawdad_tutorials_BINARY_DIR}/acetaldehyde.txt")
add_custom_target(inp_p1_benzene.txt DEPENDS "${crawdad_tutorials_BINARY_DIR}/benzene.txt")
add_custom_target(inp_p1_allene.txt DEPENDS "${crawdad_tutorials_BINARY_DIR}/allene.txt")
add_custom_target(inp_p2_h2o_geom.txt DEPENDS "${crawdad_tutorials_BINARY_DIR}/h2o_geom.txt")
add_custom_target(inp_p2_h2o_hessian.txt DEPENDS "${crawdad_tutorials_BINARY_DIR}/h2o_hessian.txt")
add_custom_target(inp_p3_h2o_sto3g 
  DEPENDS 
  "${crawdad_tutorials_BINARY_DIR}/h2o_sto3g_enuc.dat"
  "${crawdad_tutorials_BINARY_DIR}/h2o_sto3g_s.dat"
  "${crawdad_tutorials_BINARY_DIR}/h2o_sto3g_t.dat"
  "${crawdad_tutorials_BINARY_DIR}/h2o_sto3g_v.dat"
  )

add_executable(project1_acetaldehyde project1_acetaldehyde.cpp)
add_executable(project1_benzene project1_benzene.cpp)
add_executable(project1_allene project1_allene.cpp)
add_executable(project2_h2o project2_h2o.cpp)
add_executable(project3 project3.cpp)

add_dependencies(project1_acetaldehyde inp_p1_acetaldehyde.txt)
add_dependencies(project1_benzene inp_p1_benzene.txt)
add_dependencies(project1_allene inp_p1_allene.txt)
add_dependencies(project2_h2o inp_p2_h2o_geom.txt inp_p2_h2o_hessian.txt)
add_dependencies(project3 inp_p3_h2o_sto3g)

target_link_libraries(project1_acetaldehyde molecule)
target_link_libraries(project1_benzene molecule)
target_link_libraries(project1_allene molecule)
target_link_libraries(project2_h2o molecule)
target_link_libraries(project3 molecule)
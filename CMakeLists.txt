cmake_minimum_required(VERSION 2.8)
project(crawdad_tutorials CXX)

set(CMAKE_BUILD_TYPE debug)

find_package(Armadillo)
if (ARMADILLO_FOUND)
  include_directories(${ARMADILLO_INCLUDE_DIRS})
endif (ARMADILLO_FOUND)

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
  "${crawdad_tutorials_SOURCE_DIR}/project1/acetaldehyde.txt"
  "${crawdad_tutorials_BINARY_DIR}/acetaldehyde.txt"
  DEPENDS "${crawdad_tutorials_SOURCE_DIR}/project1/acetaldehyde.txt"
  )

add_custom_command(
  OUTPUT "benzene.txt"
  COMMAND "${CMAKE_COMMAND}" -E copy 
  "${crawdad_tutorials_SOURCE_DIR}/project1/benzene.txt"
  "${crawdad_tutorials_BINARY_DIR}/benzene.txt"
  DEPENDS "${crawdad_tutorials_SOURCE_DIR}/project1/benzene.txt"
  )

add_custom_command(
  OUTPUT "allene.txt"
  COMMAND "${CMAKE_COMMAND}" -E copy 
  "${crawdad_tutorials_SOURCE_DIR}/project1/allene.txt"
  "${crawdad_tutorials_BINARY_DIR}/allene.txt"
  DEPENDS "${crawdad_tutorials_SOURCE_DIR}/project1/allene.txt"
  )

add_custom_command(
  OUTPUT "h2o_geom.txt"
  COMMAND "${CMAKE_COMMAND}" -E copy 
  "${crawdad_tutorials_SOURCE_DIR}/project2/h2o_geom.txt"
  "${crawdad_tutorials_BINARY_DIR}/h2o_geom.txt"
  DEPENDS "${crawdad_tutorials_SOURCE_DIR}/project2/h2o_geom.txt"
  )

add_custom_command(
  OUTPUT "h2o_hessian.txt"
  COMMAND "${CMAKE_COMMAND}" -E copy 
  "${crawdad_tutorials_SOURCE_DIR}/project2/h2o_hessian.txt"
  "${crawdad_tutorials_BINARY_DIR}/h2o_hessian.txt"
  DEPENDS "${crawdad_tutorials_SOURCE_DIR}/project2/h2o_hessian.txt"
  )

add_custom_command(
  OUTPUT 
  "h2o_sto3g_enuc.dat"
  "h2o_sto3g_s.dat"
  "h2o_sto3g_t.dat"
  "h2o_sto3g_v.dat"
  "h2o_sto3g_eri.dat"
  COMMAND "${CMAKE_COMMAND}" -E copy "${crawdad_tutorials_SOURCE_DIR}/project3/h2o_sto3g_enuc.dat" "${crawdad_tutorials_BINARY_DIR}/h2o_sto3g_enuc.dat"
  COMMAND "${CMAKE_COMMAND}" -E copy "${crawdad_tutorials_SOURCE_DIR}/project3/h2o_sto3g_s.dat" "${crawdad_tutorials_BINARY_DIR}/h2o_sto3g_s.dat"
  COMMAND "${CMAKE_COMMAND}" -E copy "${crawdad_tutorials_SOURCE_DIR}/project3/h2o_sto3g_t.dat" "${crawdad_tutorials_BINARY_DIR}/h2o_sto3g_t.dat"
  COMMAND "${CMAKE_COMMAND}" -E copy "${crawdad_tutorials_SOURCE_DIR}/project3/h2o_sto3g_v.dat" "${crawdad_tutorials_BINARY_DIR}/h2o_sto3g_v.dat"
  COMMAND "${CMAKE_COMMAND}" -E copy "${crawdad_tutorials_SOURCE_DIR}/project3/h2o_sto3g_eri.dat" "${crawdad_tutorials_BINARY_DIR}/h2o_sto3g_eri.dat"
  DEPENDS
  "${crawdad_tutorials_SOURCE_DIR}/project3/h2o_sto3g_enuc.dat"
  "${crawdad_tutorials_SOURCE_DIR}/project3/h2o_sto3g_s.dat"
  "${crawdad_tutorials_SOURCE_DIR}/project3/h2o_sto3g_t.dat"
  "${crawdad_tutorials_SOURCE_DIR}/project3/h2o_sto3g_v.dat"
  "${crawdad_tutorials_SOURCE_DIR}/project3/h2o_sto3g_eri.dat"
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
  "${crawdad_tutorials_BINARY_DIR}/h2o_sto3g_eri.dat"
  )

add_executable(project1_acetaldehyde project1/project1_acetaldehyde.cpp)
add_executable(project1_benzene project1/project1_benzene.cpp)
add_executable(project1_allene project1/project1_allene.cpp)
add_executable(project2_h2o project2/project2_h2o.cpp)
add_executable(project3 project3/project3.cpp)

if (ARMADILLO_FOUND)
  add_executable(project3_arma project3/project3_arma.cpp)
  add_executable(project4_arma project4/project4_arma.cpp)
  add_executable(project3_arma_nocomments project3/project3_arma_nocomments.cpp)
endif (ARMADILLO_FOUND)

add_dependencies(project1_acetaldehyde inp_p1_acetaldehyde.txt)
add_dependencies(project1_benzene inp_p1_benzene.txt)
add_dependencies(project1_allene inp_p1_allene.txt)
add_dependencies(project2_h2o inp_p2_h2o_geom.txt inp_p2_h2o_hessian.txt)
add_dependencies(project3 inp_p3_h2o_sto3g)

if (ARMADILLO_FOUND)
  add_dependencies(project3_arma inp_p3_h2o_sto3g)
  add_dependencies(project3_arma_nocomments inp_p3_h2o_sto3g)
endif (ARMADILLO_FOUND)

target_link_libraries(project1_acetaldehyde molecule)
target_link_libraries(project1_benzene molecule)
target_link_libraries(project1_allene molecule)
target_link_libraries(project2_h2o molecule)
target_link_libraries(project3 molecule)

if (ARMADILLO_FOUND)
  target_link_libraries(project3_arma
    molecule
    ${ARMADILLO_LIBRARIES}
    )

  target_link_libraries(project4_arma
    molecule
    ${ARMADILLO_LIBRARIES}
    )

  target_link_libraries(project3_arma_nocomments
    molecule
    ${ARMADILLO_LIBRARIES}
    )
endif (ARMADILLO_FOUND)

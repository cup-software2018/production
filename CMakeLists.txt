# The WCMD (FADC) and PSMD (IADC) producers.
#
#   source env.sh
#   cmake -S . -B build && cmake --build build
#
# Built here and run on the cluster nodes as they are: ROOT and CUPDataModel
# are found under /home/cupsoft/Products, which is the same path, and for
# CUPDataModel V1.1 the same files, on both sides, and the RPATH points there.
# yaml-cpp is only needed as headers at build time; the nodes carry its runtime.
cmake_minimum_required(VERSION 3.16)
project(production CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)   # -O3 -DNDEBUG, what ACLiC used
endif()
set(CMAKE_CXX_STANDARD 17)        # ROOT 6.32 was built with it
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(ROOT REQUIRED COMPONENTS Core RIO Tree MathCore)
find_package(CUPDataModel REQUIRED
             HINTS /home/cupsoft/Products/CUPDATAMODEL/V1.1/lib64/cmake/CUPDataModel)
find_package(yaml-cpp REQUIRED)

foreach(det fadc iadc)
  add_executable(produce_${det} produce_${det}.cc)
  target_compile_options(produce_${det} PRIVATE -Wall -Wextra)
  target_link_libraries(produce_${det} PRIVATE
    ROOT::Core ROOT::RIO ROOT::Tree ROOT::MathCore CUPDataModel::RawObjs yaml-cpp)
endforeach()

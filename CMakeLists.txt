cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(slvana)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_BUILD_TYPE Debug)

find_package(EDepSim REQUIRED)
#include_directories(/dune/app/users/mylab/NDsofts/c7_edepsim_install_modify/include/EDepSim)
#LINK_DIRECTORIES(/dune/app/users/mylab/NDsofts/c7_edepsim_install_modify/lib)

find_package(ROOT REQUIRED COMPONENTS Geom EG)


#include(${ROOT_USE_FILE})

ROOT_GENERATE_DICTIONARY(G__Event evtinfo.h LINKDEF Linkdef.h)
add_library(evtinfo SHARED G__Event.cxx)
target_link_libraries(evtinfo ${ROOT_LIBRARIES})


add_library(Musep SHARED musep.C)
target_link_libraries(Musep PUBLIC EDepSim::edepsim_io)

add_library(Digi SHARED digitize.C)


add_executable(sandSmearGo smearSAND.C)
target_link_libraries(sandSmearGo PUBLIC Musep Digi evtinfo)

add_executable(ReadAnaEvt read_anaevt.C)
target_link_libraries(ReadAnaEvt PUBLIC evtinfo)

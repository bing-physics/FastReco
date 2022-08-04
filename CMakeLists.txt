cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(SandFastReco LANGUAGES CXX VERSION 1.0.0)

set(CMAKE_BUILD_TYPE Debug)

# Locate the ROOT package and define a number of useful targets and variables.
find_package(ROOT REQUIRED)

# Locate EDep-sim
find_package(EDepSim REQUIRED)


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
ROOT_GENERATE_DICTIONARY(G__Struct struct.h LINKDEF include/Linkdef.h)

add_library(Struct SHARED G__Struct.cxx)
target_link_libraries(Struct PUBLIC)

add_library(Utils SHARED src/utils.cpp)
target_link_libraries(Utils PUBLIC EDepSim::edepsim_io)



# Creates Digitize executable.
add_executable(FastReco src/fastreco.cpp)
target_link_libraries(FastReco Struct Utils)


install(TARGETS Struct FastReco Utils
        EXPORT SandFastRecoTargets
        RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
        
install(
   DIRECTORY "${CMAKE_BINARY_DIR}/."
   DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
   FILES_MATCHING 
   PATTERN "*.pcm"
   PATTERN "*.rootmap"
   PATTERN "lib" EXCLUDE
   PATTERN "bin" EXCLUDE
   PATTERN "CMakeFiles" EXCLUDE)

install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include 
  DESTINATION "${CMAKE_INSTALL_PREFIX}"
  PATTERN "Linkdef.h" EXCLUDE)



# ======================================================================
#  duneana main build file
# ======================================================================


# use cmake 3.14 or later
cmake_minimum_required (VERSION 3.20 FATAL_ERROR) 
# cmake_policy(VERSION 3.19) 

#project(duneana VERSION 09.48.00.d00 LANGUAGES CXX)
find_package(cetmodules REQUIRED) 
project(duneana LANGUAGES CXX)
set(${PROJECT_NAME}_CMAKE_PROJECT_VERSION_STRING 10.03.01d01)

# cetbuildtools contains our cmake modules

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

find_package( nufinder REQUIRED )
find_package( larfinder REQUIRED )
find_package( art REQUIRED ) 
find_package( art_root_io REQUIRED ) 
find_package( canvas_root_io REQUIRED ) 
find_package( Boost REQUIRED ) 
find_package(ROOT REQUIRED) 
find_package( dunepdlegacy REQUIRED EXPORT )
find_package( dunecore REQUIRED EXPORT ) 
find_package( dunereco REQUIRED EXPORT ) 
find_package( duneopdet REQUIRED EXPORT ) 
find_package( duneanaobj REQUIRED EXPORT ) 
find_package( larana REQUIRED EXPORT )
find_package( larcoreobj REQUIRED EXPORT ) 
find_package( larreco REQUIRED EXPORT ) 
find_package( larsim REQUIRED EXPORT ) 
find_package( larevt REQUIRED EXPORT ) 
find_package( larpandora REQUIRED EXPORT ) 
#find_package( srproxy REQUIRED ) 
find_package( systematicstools REQUIRED EXPORT ) 
find_package(GENIE REQUIRED EXPORT)
find_package( Geant4 REQUIRED EXPORT ) 

# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

add_subdirectory(duneana)

cet_cmake_config() 

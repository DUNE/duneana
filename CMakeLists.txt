# ======================================================================
#  duneana main build file
# ======================================================================


# use cmake 3.14 or later
cmake_minimum_required (VERSION 3.19 FATAL_ERROR)

#project(duneana VERSION 09.48.00.d00 LANGUAGES CXX)
project(duneana LANGUAGES CXX)
set(${PROJECT_NAME}_CMAKE_PROJECT_VERSION_STRING 09.79.00d02)

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

find_ups_product( art )
find_ups_product( art_root_io )
find_ups_product( canvas_root_io )
find_ups_boost( )
find_ups_root()
find_ups_product( cetbuildtools )
find_ups_product( dunecore )
find_ups_product( dunereco )
find_ups_product( duneopdet )
find_ups_product( duneanaobj )
find_ups_product( larcoreobj )
find_ups_product( larreco )
find_ups_product( larsim )
find_ups_product( larevt )
find_ups_product( larpandora )
find_ups_product( srproxy )
find_ups_product( systematicstools )
find_package(GENIE REQUIRED EXPORT)
#find_ups_product( dunepdlegacy )
find_ups_geant4( )

# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

add_subdirectory(duneana)
add_subdirectory(ups)

include(UseCPack)

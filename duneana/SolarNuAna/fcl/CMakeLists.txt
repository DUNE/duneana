install_fhicl()
add_subdirectory(dunefd)
add_subdirectory(dunefdvd)
# Also put a copy in the source tree.

FILE(GLOB fcl_files *.fcl)
install_source( EXTRAS ${fcl_files} )

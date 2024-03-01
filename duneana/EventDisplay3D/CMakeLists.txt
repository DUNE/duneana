include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make(BASENAME_ONLY
   LIB_LIBRARIES
    ${ROOT_CORE}
    ${ROOT_EVE}
    cetlib_except::cetlib_except
    cetlib::cetlib 
    larcorealg::Geometry
    lardataob::_RawData
    art::Persistency_Common canvas::canvas
   SERVICE_LIBRARIES
    duneana::EventDisplay3D
    ${ART_ROOT_IO_ROOTINPUT_SOURCE}
   MODULE_LIBRARIES
    nusimdata::SimulationBase
    larreco::RecoAlg
    duneana::EventDisplay3D
    EvtDisplayService_service
)

install_headers()
install_fhicl()
install_source()

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make(BASENAME_ONLY
   LIB_LIBRARIES
    ${ROOT_CORE}
    ${ROOT_EVE}
    cetlib_except
    cetlib::cetlib 
    larcorealg_Geometry
    lardataobj_RawData
    art::Persistency_Common canvas
   SERVICE_LIBRARIES
    duneana_EventDisplay3D
    ${ART_ROOT_IO_ROOTINPUT_SOURCE}
   MODULE_LIBRARIES
    nusimdata_SimulationBase
    larreco_RecoAlg
    duneana_EventDisplay3D
    EvtDisplayService_service
)

install_headers()
install_fhicl()
install_source()

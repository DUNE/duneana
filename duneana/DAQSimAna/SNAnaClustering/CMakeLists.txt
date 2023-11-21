art_make( BASENAME_ONLY MODULE_LIBRARIES
  lardataalg::DetectorInfo
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardataobj::RawData
  lardata::DetectorPropertiesService
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_PhotonBackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  larsim::Simulation
  nug4::ParticleNavigation
  lardataobj::Simulation
  nusimdata::SimulationBase	
  ART_FRAMEWORK_CORE
  ART_FRAMEWORK_PRINCIPAL
  ART_FRAMEWORK_SERVICES_REGISTRY
  ART_ROOT_IO_TFILE_SUPPORT
  ROOT::Core
  ART_ROOT_IO_TFILESERVICE_SERVICE
  art::Persistency_Common canvas
  art::Persistency_Provenance canvas
  art::Utilities canvas
  messagefacility::MF_MessageLogger
  ROOT_GEOM
  cetlib::cetlib cetlib_except
  ROOT_BASIC_LIB_LIST
  )

install_fhicl()
install_headers()
install_source()

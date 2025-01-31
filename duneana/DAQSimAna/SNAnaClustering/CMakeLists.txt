add_subdirectory(fcl_old)

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
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support
  ROOT::Core
  art_root_io::TFileService_service
  art::Persistency_Common 
  art::Persistency_Provenance 
  art::Utilities canvas::canvas
  messagefacility::MF_MessageLogger
  ROOT::Geom
  cetlib::cetlib cetlib_except::cetlib_except
  ROOT::Core ROOT::Hist ROOT::Tree
  )

install_fhicl()
install_headers()
install_source()

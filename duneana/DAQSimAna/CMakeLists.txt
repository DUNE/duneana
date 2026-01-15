include(MakeDuneToolBuilder)
add_subdirectory(SNAnaClustering)
add_subdirectory(RunningSumHitFinder)
add_subdirectory(CalibrationTree)
add_subdirectory(DAQQuickClustering)
add_subdirectory(DataHit)
add_subdirectory(fcl)
add_subdirectory(PackedDump)
add_subdirectory(TriggerPrimitiveFinder)
add_subdirectory(WaveformDump)

art_make( BASENAME_ONLY MODULE_LIBRARIES
  dunepdlegacy::Services_ChannelMap_PdspChannelMapService_service
  larcore::Geometry_Geometry_service
  lardataalg::DetectorInfo
  larcorealg::Geometry
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  dunecore::DuneObj
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  larsim::MCCheater_PhotonBackTrackerService_service
  larsim::Simulation nug4::ParticleNavigation lardataobj::Simulation
  lardata::ArtDataHelper
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
  cetlib::cetlib cetlib_except::cetlib_except
  ROOT::Core ROOT::Hist ROOT::Tree
  )

install_fhicl()
install_headers()
install_source()

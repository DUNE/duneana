make_dune_tool_builder(dunecore::AbsRunningSumTPFinderTool)

art_make( BASENAME_ONLY MODULE_LIBRARIES
  dunepdlegacy_Services_ChannelMap_PdspChannelMapService_service
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
  larsim::Simulation nug4::ParticleNavigation lardataobj_Simulation
  lardata::ArtDataHelper
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
  cetlib::cetlib cetlib_except
  ROOT_BASIC_LIB_LIST
  EXCLUDE read_packed.cxx AbsRunningSumTPFinderPass1_tool.cc
  )

cet_build_plugin(AbsRunningSumTPFinderPass1 art::tool
  fhiclcpp cetlib::cetlib cetlib_except

)

install_fhicl()
install_headers()
install_source()

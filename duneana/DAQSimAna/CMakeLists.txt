include(MakeDuneToolBuilder)
add_subdirectory(SNAnaClustering)
add_subdirectory(SNPointingTools)
add_subdirectory(AbsRunningSumHitFinder)

make_dune_tool_builder(dunecore::TriggerPrimitiveFinderTool)

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
  EXCLUDE read_packed.cxx TriggerPrimitiveFinderPass1_tool.cc TriggerPrimitiveFinderPass2_tool.cc TriggerPrimitiveFinderTemplate_tool.cc
  )

cet_build_plugin(TriggerPrimitiveFinderPass1 art::tool
  fhiclcpp cetlib::cetlib cetlib_except
  
)

cet_build_plugin(TriggerPrimitiveFinderPass2 art::tool
  fhiclcpp cetlib::cetlib cetlib_except
  
  duneana_DAQSimAna_TriggerPrimitiveFinderPass1_tool
)

cet_build_plugin(TriggerPrimitiveFinderTemplate art::tool
  fhiclcpp cetlib::cetlib cetlib_except
  
)

cet_make_exec(read_packed SOURCE read_packed.cxx)

install_fhicl()
install_headers()
install_source()

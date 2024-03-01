include(MakeDuneToolBuilder)
add_subdirectory(SNAnaClustering)
add_subdirectory(RunningSumHitFinder)
add_subdirectory(AbsRunningSumHitFinder)

make_dune_tool_builder(dunecore::TriggerPrimitiveFinderTool)

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
  EXCLUDE read_packed.cxx TriggerPrimitiveFinderPass1_tool.cc TriggerPrimitiveFinderPass2_tool.cc TriggerPrimitiveFinderTemplate_tool.cc
  )

cet_build_plugin(TriggerPrimitiveFinderPass1 art::tool LIBRARIES
  fhiclcpp::fhiclcpp cetlib::cetlib cetlib_except::cetlib_except
  
)

cet_build_plugin(TriggerPrimitiveFinderPass2 art::tool LIBRARIES
  fhiclcpp::fhiclcpp cetlib::cetlib cetlib_except::cetlib_except
  
  duneana_DAQSimAna_TriggerPrimitiveFinderPass1_tool
)

cet_build_plugin(TriggerPrimitiveFinderTemplate art::tool LIBRARIES
  fhiclcpp::fhiclcpp cetlib::cetlib cetlib_except::cetlib_except
  
)

cet_make_exec( NAME read_packed SOURCE read_packed.cxx)

install_fhicl()
install_headers()
install_source()

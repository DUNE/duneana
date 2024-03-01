make_dune_tool_builder(dunecore::RunningSumTPFinderTool)

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
  art_root_io::TFileService_service
  art::Persistency_Common
  art::Persistency_Provenance 
  art::Utilities canvas::canvas
  messagefacility::MF_MessageLogger
  cetlib::cetlib cetlib_except::cetlib_except
  ROOT::Core ROOT::Hist ROOT::Tree
  EXCLUDE read_packed.cxx RunningSumTPFinderPass1_tool.cc RunningSumTPFinderPass2_tool.cc RunningSumTPFinderPass3_tool.cc RunningSumTPFinderPass4_tool.cc
  )

cet_build_plugin(RunningSumTPFinderPass1 art::tool LIBRARIES
  fhiclcpp::fhiclcpp cetlib::cetlib cetlib_except::cetlib_except
  
)
cet_build_plugin(RunningSumTPFinderPass2 art::tool LIBRARIES
  fhiclcpp::fhiclcpp cetlib::cetlib cetlib_except::cetlib_except
  
)
cet_build_plugin(RunningSumTPFinderPass3 art::tool LIBRARIES
  fhiclcpp::fhiclcpp cetlib::cetlib cetlib_except::cetlib_except
  
)
cet_build_plugin(RunningSumTPFinderPass4 art::tool LIBRARIES
  fhiclcpp::fhiclcpp cetlib::cetlib cetlib_except::cetlib_except
  
)

install_fhicl()
install_headers()
install_source()

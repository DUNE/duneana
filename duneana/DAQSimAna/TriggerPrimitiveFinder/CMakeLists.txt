include(MakeDuneToolBuilder)

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
  larsim::Simulation nug4::ParticleNavigation
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
  EXCLUDE TriggerPrimitiveFinderPass1_tool.cc TriggerPrimitiveFinderPass2_tool.cc TriggerPrimitiveFinderTemplate_tool.cc
  )

cet_build_plugin(TriggerPrimitiveFinderPass1 art::tool
  fhiclcpp::fhiclcpp cetlib::cetlib cetlib_except::cetlib_except
  
)

cet_build_plugin(TriggerPrimitiveFinderPass2 art::tool
  fhiclcpp::fhiclcpp cetlib::cetlib cetlib_except::cetlib_except

  duneana_DAQSimAna_TriggerPrimitiveFinder_TriggerPrimitiveFinderPass1_tool

)

cet_build_plugin(TriggerPrimitiveFinderTemplate art::tool
  fhiclcpp::fhiclcpp cetlib::cetlib cetlib_except::cetlib_except
  
)

install_fhicl()
install_headers()
install_source()



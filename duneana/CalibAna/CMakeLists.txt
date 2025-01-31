add_subdirectory(fcl)
build_dictionary()

include_directories("$ENV{LAREVT_DIR}/include")
include_directories("$ENV{ART_ROOT_IO_DIR}/include")

cet_build_plugin(
  CalibAnaTreeSelectAllTracks art::tool
  LIBRARIES
  #canvas
  lardataalg::DetectorInfo
  lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
  larcore::Geometry_Geometry_service
)

cet_build_plugin(
  CalibAnaTreeSelectAnode2CathodeTrack art::tool
  LIBRARIES
  #canvas
  lardataalg::DetectorInfo
  lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
  larcore::Geometry_Geometry_service
)

cet_build_plugin(
  CalibAnaTreeSelectStoppingTrack art::tool
  LIBRARIES
  #canvas
  lardataalg::DetectorInfo
  lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
  larcore::Geometry_Geometry_service
)
#target_include_directories(
#  duneana::CalibAna_CalibAnaTreeSelectStoppingTrack_tool
#  PRIVATE
#  "$ENV{LAREVT_DIR}/include"
#)

cet_build_plugin(
  CalibAnaTree art::module
  LIBRARIES
  lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
  lardataalg::DetectorInfo
  dunecore::DAQTriggerSim_TriggerDataProducts
  larcorealg::Geometry
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  nusimdata::SimulationBase
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_PhotonBackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  nug4::ParticleNavigation
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art_root_io::tfile_support
  ROOT::Core
  art_root_io::TFileService_service
  art::Persistency_Common
  art::Persistency_Provenance
  art::Utilities
  #canvas
  messagefacility::MF_MessageLogger
  cetlib::cetlib
  cetlib_except::cetlib_except
  ${ROOT_BASIC_LIB_LIST}
  ROOT::Tree
  ROOT::Minuit
  ROOT::X3d   	 
  larcorealg::GeoAlgo
)
#target_include_directories(
#  duneana::CalibAna_CalibAnaTree_module
#  PRIVATE
#  "$ENV{ART_ROOT_IO_DIR}/include"
#)


install_fhicl()
install_headers()
install_source()

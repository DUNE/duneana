art_make(BASENAME_ONLY
	 MODULE_LIBRARIES
	 larreco::RecoAlg
	 lardataobj::RecoBase
         lardataobj::RawData
         lardata::ArtDataHelper
	 larcore::Geometry_Geometry_service
         larcorealg::Geometry
	 lardataalg::DetectorInfo
         nusimdata::SimulationBase
	 larsim::MCCheater_BackTrackerService_service
         larsim::MCCheater_ParticleInventoryService_service
	 dunereco::HitFinderDUNE
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
)

install_headers()
install_fhicl()
install_scripts()

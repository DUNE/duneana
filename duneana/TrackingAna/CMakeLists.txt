art_make(BASENAME_ONLY
	 MODULE_LIBRARIES
         lardataalg::DetectorInfo
	 lardataobj::RawData
	 lardataobj::RecoBase
	 larcorealg::Geometry
	 larcore::Geometry_Geometry_service
	 lardata::Utilities
	 larsim::MCCheater_BackTrackerService_service
	 larsim::MCCheater_ParticleInventoryService_service
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
	 ROOT_GEOM
)

install_headers()
install_fhicl()
install_scripts()

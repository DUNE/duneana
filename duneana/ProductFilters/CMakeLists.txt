
art_make( BASENAME_ONLY
         MODULE_LIBRARIES
         lardataobj::RecoBase
         lardata::ArtDataHelper
         larcorealg::Geometry
         larcore::Geometry_Geometry_service
         larsim::MCCheater_BackTrackerService_service
         larsim::MCCheater_ParticleInventoryService_service
         lardata::Utilities
	 nurandom::RandomUtils_NuRandomService_service
	 nusimdata::SimulationBase
 	 art::Framework_Core
	 art::Framework_Principal
	 art::Framework_Services_Registry
	 art_root_io::tfile_support
	 ROOT::Core
	 art_root_io::TFileService_service
	 art::Persistency_Common
	 art::Utilities 
	 canvas::canvas
	 messagefacility::MF_MessageLogger		   
	 fhiclcpp::fhiclcpp
	 cetlib::cetlib 
	 cetlib_except::cetlib_except
         ROOT::Core ROOT::Hist ROOT::Tree
         ROOT::Geom
         )

install_headers()
install_fhicl()
install_source()

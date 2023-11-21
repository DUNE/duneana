
cet_build_plugin(ShowerAnalysis art::module
			     lardataobj::RawData
		      	     lardataobj::RecoBase
		      	     lardata::ArtDataHelper
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
 			     art::Persistency_Common 
			     art::Persistency_Provenance 
			     art::Utilities 
			     canvas
			     messagefacility::MF_MessageLogger
			     cetlib::cetlib 
			     cetlib_except
			     ROOT_BASIC_LIB_LIST
)


install_headers()
install_fhicl()
install_source()
install_scripts()


# art_make(LIBRARY_NAME  ShowerAna
# 	 LIB_LIBRARIES lardataobj_RawData
# 	 	       lardataobj_RecoBase	  
# 	 	       larcorealg_Geometry  
# 	 	       larcore_Geometry_Geometry_service
# 	 	       lardata_Utilities
# 		       larsim_MCCheater_BackTrackerService_service
#           larsim_MCCheater_ParticleInventoryService_service
# 	 	       nusimdata::SimulationBase
#          	       art::Framework_Core
# 	 	       art::Framework_Principal
# 	 	       art::Framework_Services_Registry
#          	       art_root_io::tfile_support ROOT::Core
#          	       art_root_io::TFileService_service
#  	 	       art::Persistency_Common canvas
# 	 	       art::Persistency_Provenance canvas
# 	 	       art::Utilities canvas
# 	 	       ${messagefacility::MF_MessageLogger}
# 	 	       
# 	 	       cetlib::cetlib cetlib_except
# 	 	       ROOT_BASIC_LIB_LIST
# 		       ShowerAna
# )

# install_headers()
# install_fhicl()
# install_source()
# install_scripts()

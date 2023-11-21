art_make(BASENAME_ONLY
	       MODULE_LIBRARIES
         dunecore_DAQTriggerSim_TriggerDataProducts
	       larreco_RecoAlg
	       lardataobj_RecoBase
         lardataobj_RawData
         lardata_ArtDataHelper
	       larcore_Geometry_Geometry_service
         larcorealg_Geometry
	       lardataalg_DetectorInfo
         nusimdata::SimulationBase
         larsim_MCCheater_BackTrackerService_service
         larsim_MCCheater_ParticleInventoryService_service
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
         ART_PERSISTENCY_COMMON
         ART_PERSISTENCY_PROVENANCE
         ART_UTILITIES

         cetlib::cetlib cetlib_except
         ROOT_BASIC_LIB_LIST
)

install_headers()
install_fhicl()
install_scripts()

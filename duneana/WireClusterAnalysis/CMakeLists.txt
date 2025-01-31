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
         art::Persistency_Common
         art::Persistency_Provenance
         cetlib::cetlib cetlib_except::cetlib_except
         ROOT::Core ROOT::Hist ROOT::Tree
)

install_headers()
install_fhicl()
install_scripts()

add_subdirectory(fcl)
art_make(BASENAME_ONLY
	       MODULE_LIBRARIES
         dunecore::DAQTriggerSim_TriggerDataProducts
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
         art::Framework_Core
         art::Framework_Principal
         art::Framework_Services_Registry
         art_root_io::tfile_support
         ROOT::Core
         art_root_io::TFileService_service
         art::Persistency_Common
         art::Persistency_Provenance
         art::Utilities 
	 canvas::canvas
         art::Persistency_Common
         art::Persistency_Provenance
         messagefacility::MF_MessageLogger
         cetlib::cetlib 
	 cetlib_except::cetlib_except
         ROOT::Core ROOT::Hist ROOT::Tree
         #hep_hpc_hdf5
         #HDF5
)

install_headers()
install_fhicl()
install_scripts()

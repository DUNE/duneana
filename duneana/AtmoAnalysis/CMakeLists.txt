art_make( MODULE_LIBRARIES 
	  larcorealg::Geometry
	  larcore::Geometry_Geometry_service
	  larsim::Simulation nug4::ParticleNavigation 
	  lardataobj::Simulation
	  larsim::MCCheater_BackTrackerService_service
	  larsim::MCCheater_ParticleInventoryService_service
          larsim::Utils
	  lardata::Utilities
	  larevt::Filters
	  lardataobj::RawData
	  lardataobj::RecoBase
	  lardataobj::AnalysisBase
	  lardata::RecoObjects
    	  lardata::ArtDataHelper
	  larreco::RecoAlg
	  larpandora::LArPandoraInterface
	  dunecore::SignalShapingServiceDUNE10kt_service
	  dunecore::SignalShapingServiceDUNE35t_service
          nusimdata::SimulationBase
	  dunereco::AnaUtils
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
	  fhiclcpp::fhiclcpp
	  cetlib::cetlib cetlib_except::cetlib_except
          ROOT::Geom
          ROOT::XMLIO
          ROOT::Gdml
	  ROOT::Core ROOT::Hist ROOT::Tree
        )

install_headers()
install_fhicl()
install_source()
install_scripts()

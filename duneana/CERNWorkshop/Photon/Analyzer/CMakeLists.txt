art_make( MODULE_LIBRARIES
 
        duneopdet::OpticalDetector
	nusimdata::SimulationBase
	art::Framework_Core
	dunereco::FDSensOpt_FDSensOptData
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

  # lardataobj_Simulation

  # duneopdet_OpticalDetector
)

install_headers()
install_fhicl()
install_source()

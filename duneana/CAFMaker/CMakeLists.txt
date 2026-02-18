if( DEFINED ENV{TENSORFLOW_DIR} )
include_directories($ENV{SRPROXY_INC})
cet_build_plugin(CAFMaker   art::module LIBRARIES
                        duneanaobj::StandardRecord
                        duneanaobj::StandardRecordFlat
                        nugen::NuReweight
                        nugen::NuReweight_art
			nugen::EventGeneratorBase_GENIE
			nusimdata::SimulationBase
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
                        cetlib::cetlib cetlib_except::cetlib_except
                        ROOT::Core ROOT::Hist ROOT::Tree
			${GENIE_LIB_LIST}
                        dunereco::CVN_func
			larpandora::LArPandoraInterface
                        larreco::Calorimetry
                        dunereco::AnaUtils
                        larsim::Utils
              BASENAME_ONLY
)
endif()

install_headers()
install_fhicl()
install_source()
install_scripts()

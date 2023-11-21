include_directories ( $ENV{GENIE_INC}/GENIE )

# nusystematics
##cet_find_library( SYSTTOOLS NAMES systematicstools_interface systematicstools_interpreters systematicstools_utility PATHS ENV SYSTEMATICSTOOLS_LIB NO_DEFAULT_PATH )

cet_build_plugin(CAFMaker   art::module
                        duneanaobj_StandardRecord
                        duneanaobj_StandardRecordFlat
                        nugen::NuReweight
                        nugen::NuReweight_art
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
			${GENIE_LIB_LIST}
                        systematicstools::interface
			systematicstools::interpreters
			systematicstools::utility
                        dunereco::CVN_func
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
install_scripts()


cet_build_plugin(DuneEventFilter art::module
			dunecore::ArtSupport
			lardataobj::RawData
                        lardataobj::RecoBase
	                larcorealg::Geometry
	                larcore::Geometry_Geometry_service
                        lardata::Utilities    
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
	                fhiclcpp::fhiclcpp
			cetlib::cetlib cetlib_except
                        ROOT_GEOM
                        ROOT_XMLIO
                        ROOT_GDML
			ROOT_BASIC_LIB_LIST
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()

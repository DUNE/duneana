
cet_build_plugin(dunezsanalysis art::module
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardataobj::RawData
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
  BASENAME_ONLY
  )

cet_build_plugin(dEdx art::module
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardataobj::RawData
  lardataobj::RecoBase
  lardata::Utilities
  larsim::Simulation 
  nug4::ParticleNavigation 
  lardataobj::Simulation
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
  BASENAME_ONLY
  )

cet_build_plugin(EMEnergyCalib art::module
  lardataalg::DetectorInfo
  lardataobj::RawData
  lardataobj::RecoBase
  larcorealg::Geometry
  larcore::Geometry_Geometry_service
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  larsim::Simulation 
  nug4::ParticleNavigation 
  lardataobj::Simulation
  lardata::Utilities
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
  BASENAME_ONLY
  )

cet_build_plugin(EMPi0Energy art::module
  lardataobj::RawData
  lardataobj::RecoBase
  larcorealg::Geometry
  larcore::Geometry_Geometry_service
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  larsim::Simulation nug4::ParticleNavigation lardataobj_Simulation
  lardata::Utilities
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
  BASENAME_ONLY
  )

install_headers()
install_fhicl()
install_source()
install_scripts()

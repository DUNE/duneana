#add_subdirectory(OpHitFinder)

cet_make_library(LIBRARY_NAME OpDetResponseService INTERFACE
  SOURCE OpDetResponseInterface.h
  LIBRARIES INTERFACE
  larcore::Geometry_Geometry_service
  art_plugin_types::serviceDeclaration
  art::Framework_Services_Registry
)

cet_write_plugin_builder(lar::OpDetResponseService art::service Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  larana::OpDetResponseService
)

cet_make_library(SOURCE
  FlashHypothesis.cxx
  FlashHypothesisAnaAlg.cxx
  FlashHypothesisCalculator.cxx
  FlashHypothesisComparison.cxx
  FlashHypothesisCreator.cxx
  FlashUtilities.cxx
  OpFlashAlg.cxx
  OpFlashAnaAlg.cxx
  SimPhotonCounter.cxx
  SimPhotonCounterAlg.cxx
  LIBRARIES
  PUBLIC
  larcorealg::headers
  lardataobj::MCBase
  lardataobj::RecoBase
  lardataobj::Simulation
  larsim::MCSTReco
  larsim::PhotonVisibilityTypes
  ROOT::Hist
  ROOT::Physics
  PRIVATE
  larsim::PhotonPropagation_PhotonVisibilityService_service
  larcore::Geometry_Geometry_service
  lardataalg::DetectorInfo
  lardataalg::UtilitiesHeaders
  lardataobj::MCBase
  lardataobj::RecoBase
  larcorealg::Geometry
  larcoreobj::SimpleTypesAndConstants
  fhiclcpp::fhiclcpp
  ROOT::RIO
  ROOT::Tree
)

include(lar::OpDetResponseService)

cet_build_plugin(DefaultOpDetResponse lar::OpDetResponseService
  LIBRARIES PRIVATE
  lardata::LArPropertiesService
  larcore::ServiceUtil
  messagefacility::MF_MessageLogger
)

cet_build_plugin(OpDigiProperties art::service
  LIBRARIES
  PUBLIC
  larcore::Geometry_Geometry_service
  lardataobj::headers
  PRIVATE
  messagefacility::MF_MessageLogger
  CLHEP::Random
)

cet_build_plugin(MicrobooneOpDetResponse lar::OpDetResponseService
  LIBRARIES PRIVATE
  lardataobj::Simulation
)

cet_build_plugin(FIFOHistogramAna art::EDAnalyzer
  LIBRARIES PRIVATE
  lardataobj::headers
  art_root_io::TFileService_service
  art_root_io::tfile_support
  art::Framework_Principal
  art::Framework_Services_Registry
  fhiclcpp::fhiclcpp
  ROOT::Hist
)

cet_build_plugin(FlashHypothesisAna art::EDAnalyzer
  LIBRARIES PRIVATE
  larana::OpticalDetector
  lardata::LArPropertiesService
  lardataobj::headers
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  fhiclcpp::fhiclcpp
  ROOT::Hist
  ROOT::Tree
)


cet_build_plugin(OpDigiAna art::EDAnalyzer
  LIBRARIES PRIVATE
  larana::OpticalDetector_OpDigiProperties_service
  lardataobj::RawData
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  fhiclcpp::fhiclcpp
  ROOT::Hist
)

cet_build_plugin(OpFlashAna art::EDAnalyzer
  LIBRARIES PRIVATE
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
  ROOT::Hist
  ROOT::Tree
)

cet_build_plugin(OpHitAnaNew art::EDAnalyzer
  LIBRARIES PRIVATE
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
  ROOT::Hist
  ROOT::Tree
)

cet_build_plugin(OpFlashFinder art::EDProducer
  LIBRARIES PRIVATE
  larana::OpticalDetector
  lardata::DetectorClocksService
  lardata::AssociationUtil
  larcore::ServiceUtil
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
)

cet_build_plugin(OpFlashMCTruthAna art::EDAnalyzer
  LIBRARIES PRIVATE
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  fhiclcpp::fhiclcpp
  ROOT::Tree
)

cet_build_plugin(OpFlashSimpleAna art::EDAnalyzer
  LIBRARIES PRIVATE
  larana::OpticalDetector
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  fhiclcpp::fhiclcpp
  ROOT::Tree
)

cet_build_plugin(OpHitAna art::EDAnalyzer
  LIBRARIES PRIVATE
  larana::OpticalDetector_OpDigiProperties_service
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Hist
  ROOT::Physics
  ROOT::Tree
)



cet_build_plugin(OpMCDigi art::EDProducer
  LIBRARIES PRIVATE
  larana::OpDetResponseService
  larana::OpticalDetector_OpDigiProperties_service
  larsim::Simulation
  lardataobj::RawData
  lardataobj::Simulation
  nurandom::RandomUtils_NuRandomService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  fhiclcpp::fhiclcpp
  CLHEP::Random
)

cet_build_plugin(OptDetDigitizer art::EDProducer
  LIBRARIES PRIVATE
  larana::OpDetResponseService
  larana::OpticalDetector_OpDigiProperties_service
  larsim::Simulation
  lardataobj::headers
  lardataobj::Simulation
  nurandom::RandomUtils_NuRandomService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  fhiclcpp::fhiclcpp
  CLHEP::Random
)

cet_build_plugin(OpticalRawDigitReformatter art::EDProducer
  LIBRARIES PRIVATE
  lardata::DetectorClocksService
  lardataobj::headers
  lardataobj::RawData
  lardataobj::Simulation
  art::Framework_Principal
  fhiclcpp::fhiclcpp
)

cet_build_plugin(OpticalRecoAna art::EDAnalyzer
  LIBRARIES PRIVATE
  larana::OpticalDetector_OpDigiProperties_service
  larsim::MCCheater_ParticleInventoryService_service
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Hist
  ROOT::Tree
)


cet_build_plugin(SimPhotonCounter art::EDAnalyzer
  LIBRARIES PRIVATE
  larana::OpDetResponseService
  larsim::MCCheater_ParticleInventoryService_service
  larsim::PhotonPropagation_PhotonVisibilityService_service
  larsim::Simulation
  larcore::Geometry_Geometry_service
  lardataobj::MCBase
  lardataobj::Simulation
  larcorealg::Geometry
  nug4::ParticleNavigation
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Hist
  ROOT::Physics
  ROOT::Tree
)

cet_build_plugin(PhotonInf art::EDAnalyzer
  LIBRARIES PRIVATE
  larana::OpDetResponseService
  larsim::MCCheater_ParticleInventoryService_service
  larsim::PhotonPropagation_PhotonVisibilityService_service
  larsim::Simulation
  larcore::Geometry_Geometry_service
  lardataobj::MCBase
  lardataobj::Simulation
  larcorealg::Geometry
  nug4::ParticleNavigation
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Hist
  ROOT::Physics
  ROOT::Tree
)

# --- BEGIN -- hit finding algorithm maker tools -------------------------------

# --- END ---- hit finding algorithm maker tools -------------------------------

# --- BEGIN -- pedestal estimator algorithm maker tools ------------------------

# --- END ---- pedestal estimator algorithm maker tools ------------------------

install_headers()
install_fhicl()
install_source()
#install_fw(LIST toyWaveform.txt SUBDIRNAME OpticalDetector)

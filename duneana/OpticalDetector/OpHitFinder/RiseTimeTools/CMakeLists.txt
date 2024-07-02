cet_enable_asserts()

cet_make_library(LIBRARY_NAME RiseTimeCalculatorTool INTERFACE
  SOURCE 
     RiseTimeCalculatorBase.h
  LIBRARIES
  ROOT::Hist
)

cet_write_plugin_builder(lar::RiseTimeCalculatorTool art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  larana::RiseTimeCalculatorTool
)

include(lar::RiseTimeCalculatorTool)

cet_build_plugin(RiseTimeThreshold lar::RiseTimeCalculatorTool
  LIBRARIES PRIVATE
  fhiclcpp::fhiclcpp
)

cet_build_plugin(RiseTimeGaussFit lar::RiseTimeCalculatorTool
  LIBRARIES PRIVATE
  fhiclcpp::fhiclcpp
  messagefacility::MF_MessageLogger
  ROOT::Hist
)


install_headers()
install_source()
install_fhicl()

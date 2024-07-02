/**
 * @file   larana/OpticalDetector/PedAlgoMakerToolBase.h
 * @brief  Base class wrapping hit finding algorithms into _art_ tools.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   February 12, 2023
 *
 * This library is header-only.
 */

#ifndef LARANA_OPTICALDETECTOR_PEDALGOMAKERTOOLBASE_H
#define LARANA_OPTICALDETECTOR_PEDALGOMAKERTOOLBASE_H

// LArSoft libraries
#include "larana/OpticalDetector/IPedAlgoMakerTool.h"
#include "larana/OpticalDetector/OpHitFinder/PMTPedestalBase.h"
#include "larana/OpticalDetector/OpHitFinder/RiseTimeTools/RiseTimeCalculatorBase.h"

// framework libraries
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/OptionalDelegatedParameter.h"

namespace opdet {
  template <class PedAlgoClass>
  class PedAlgoMakerToolBase;
}

/**
 * @brief Base _art_ tool class wrapping a pedestal estimator algorithm.
 * @tparam PedAlgoClass the pedestal estimator algorithm class being created
 *
 * Algorithms of the pedestal estimation mini-framework in `larana` follow a
 * factory pattern which is not the one native in _art_.
 *
 * This base class provides the backbone to a _art_ tool wrapping one of the
 * pedestal estimator algorithms. The only function of these tools is to create
 * and configure a pedestal estimator algorithm object: the tools do not offer
 * any pedestal estimation functionality by themselves.
 *
 * Note that these tools all implement a single _art_ tool interface, which
 * is defined as `opdet::IPedAlgoMakerTool`.
 *
 * The current _art_ plugin system will be content with a simple declaration
 * like:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * DEFINE_ART_CLASS_TOOL(opdet::HitAlgoMakerToolBase<pmtana::PedAlgoRollingMean>)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for defining the tool that creates a `pmtana::PedAlgoRollingMean` algorithm.
 * This definition needs to live in a source file with an appropriate name
 * (e.g. `PedAlgoRollingMeanMaker_tool.cc`), as the source name will be used to
 * identify which library to load at run time.
 *
 *
 * Pedestal estimator class (`PedAlgoClass`) requirements
 * -----------------------------------------------
 *
 * The constructor of the pedestal estimator class, `PedAlgoClass`, must support
 * one parameter:
 *
 *     PedAlgoClass::PedAlgoClass(fhicl::ParameterSet const& config);
 *
 *
 */
template <class PedAlgoClass>
class opdet::PedAlgoMakerToolBase : public opdet::IPedAlgoMakerTool {

public:
  struct Config {

    fhicl::DelegatedParameter PedAlgoPset{
      fhicl::Name{"PedAlgoPset"},
      fhicl::Comment{"configuration of the pedestal estimator  algorithm"}};

  }; // struct Config

  using Parameters = art::ToolConfigTable<Config>;

  /// Constructor: copies and stores the configuration for the algorithm.
  PedAlgoMakerToolBase(Parameters const& params);

  /// Creates and returns the algorithm from the configuration on construction.
  virtual std::unique_ptr<pmtana::PMTPedestalBase> makeAlgo() override;

protected:
  Config fConfig; ///< Tool configuration cache.

}; // opdet::PedAlgoMakerToolBase

// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <class PedAlgoClass>
opdet::PedAlgoMakerToolBase<PedAlgoClass>::PedAlgoMakerToolBase(Parameters const& params)
  : fConfig{params()}
{}

// -----------------------------------------------------------------------------
template <class PedAlgoClass>
std::unique_ptr<pmtana::PMTPedestalBase> opdet::PedAlgoMakerToolBase<PedAlgoClass>::makeAlgo()
{

  return std::make_unique<PedAlgoClass>(fConfig.PedAlgoPset.template get<fhicl::ParameterSet>());

} // opdet::PedAlgoMakerToolBase::makeAlgo()

// -----------------------------------------------------------------------------

#endif // LARANA_OPTICALDETECTOR_PEDALGOMAKERTOOLBASE_H

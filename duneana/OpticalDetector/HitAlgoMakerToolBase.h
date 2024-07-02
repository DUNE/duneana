/**
 * @file   larana/OpticalDetector/HitAlgoMakerToolBase.h
 * @brief  Base class wrapping hit finding algorithms into _art_ tools.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   February 12, 2023
 *
 * This library is header-only.
 */

#ifndef LARANA_OPTICALDETECTOR_HITALGOMAKERTOOLBASE_H
#define LARANA_OPTICALDETECTOR_HITALGOMAKERTOOLBASE_H

// LArSoft libraries
#include "larana/OpticalDetector/IHitAlgoMakerTool.h"
#include "larana/OpticalDetector/OpHitFinder/PMTPulseRecoBase.h"
#include "larana/OpticalDetector/OpHitFinder/RiseTimeTools/RiseTimeCalculatorBase.h"

// framework libraries
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/OptionalDelegatedParameter.h"

namespace opdet {
  template <class HitAlgoClass>
  class HitAlgoMakerToolBase;
}

/**
 * @brief Base _art_ tool class wrapping a hit algorithm.
 * @tparam HitAlgoClass the hit finder algorithm class being created
 *
 * Algorithms of the hit finding mini-framework in `larana` follow a factory
 * pattern which is not the one native in _art_.
 *
 * This base class provides the backbone to a _art_ tool wrapping one of the hit
 * finder algorithms. The only function of these tools is to create and
 * configure a hit finder algorithm object: the tools do not offer any hit
 * finding functionality by themselves.
 *
 * Note that these tools all implement a single _art_ tool interface, which
 * is defined as `opdet::IHitAlgoMakerTool`.
 *
 * The current _art_ plugin system will be content with a simple declaration
 * like:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * DEFINE_ART_CLASS_TOOL(opdet::HitAlgoMakerToolBase<pmtana::AlgoSlidingWindow>)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * for defining the tool that creates a `pmtana::AlgoSlidingWindow` algorithm.
 * This definition needs to live in a source file with an appropriate name
 * (e.g. `AlgoSlidingWindowMaker_tool.cc`), as the source name will be used to
 * identify which library to load at run time.
 *
 * Hit finder class (`HitAlgoClass`) requirements
 * -----------------------------------------------
 *
 * The constructor of the hit finder class, `HitAlgoClass`, must support two
 * parameters:
 *
 *     HitAlgoClass::HitAlgoClass(
 *       fhicl::ParameterSet const& config,
 *       std::unique_ptr<pmtana::RiseTimeCalculatorBase>&& riseCalcAlgo
 *       );
 *
 *
 */
template <class HitAlgoClass>
class opdet::HitAlgoMakerToolBase : public opdet::IHitAlgoMakerTool {

public:
  struct Config {

    fhicl::DelegatedParameter HitAlgoPset{
      fhicl::Name{"HitAlgoPset"},
      fhicl::Comment{"configuration of the hit finder algorithm"}};

    fhicl::OptionalDelegatedParameter RiseTimeCalculator{
      fhicl::Name{"RiseTimeCalculator"},
      fhicl::Comment{"configuration of the rise time calculator algorithm"}};

  }; // struct Config

  using Parameters = art::ToolConfigTable<Config>;

  /// Constructor: copies and stores the configuration for the algorithm.
  HitAlgoMakerToolBase(Parameters const& params);

  /// Creates and returns the algorithm from the configuration on construction.
  virtual std::unique_ptr<pmtana::PMTPulseRecoBase> makeAlgo() override;

protected:
  Config fConfig; ///< Tool configuration cache.

}; // opdet::HitAlgoMakerToolBase

// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <class HitAlgoClass>
opdet::HitAlgoMakerToolBase<HitAlgoClass>::HitAlgoMakerToolBase(Parameters const& params)
  : fConfig{params()}
{}

// -----------------------------------------------------------------------------
template <class HitAlgoClass>
std::unique_ptr<pmtana::PMTPulseRecoBase> opdet::HitAlgoMakerToolBase<HitAlgoClass>::makeAlgo()
{

  std::unique_ptr<pmtana::RiseTimeCalculatorBase> riseCalcAlgo =
    fConfig.RiseTimeCalculator.hasValue() ?
      art::make_tool<pmtana::RiseTimeCalculatorBase>(
        fConfig.RiseTimeCalculator.template get_if_present<fhicl::ParameterSet>().value()) :
      nullptr;

  return std::make_unique<HitAlgoClass>(fConfig.HitAlgoPset.template get<fhicl::ParameterSet>(),
                                        std::move(riseCalcAlgo));

} // opdet::HitAlgoMakerToolBase::makeAlgo()

// -----------------------------------------------------------------------------

#endif // LARANA_OPTICALDETECTOR_HITALGOMAKERTOOLBASE_H

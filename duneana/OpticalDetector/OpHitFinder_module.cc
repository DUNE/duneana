// -*- mode: c++; c-basic-offset: 2; -*-
// Gleb Sinev, Duke, 2016
//
// This module finds periods of time-localized activity
// on each channel of the optical system and creates OpHits.
//
// Split from OpFlashFinder_module.cc
// by Ben Jones, MIT, 2013
//

// LArSoft includes
#include "larana/OpticalDetector/IHitAlgoMakerTool.h"
#include "larana/OpticalDetector/IPedAlgoMakerTool.h"
#include "larana/OpticalDetector/OpHitFinder/OpHitAlg.h"
#include "larana/OpticalDetector/OpHitFinder/PMTPedestalBase.h"
#include "larana/OpticalDetector/OpHitFinder/PMTPulseRecoBase.h"
#include "larana/OpticalDetector/OpHitFinder/PulseRecoManager.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "larreco/Calibrator/IPhotonCalibrator.h"
#include "larreco/Calibrator/IPhotonCalibratorService.h"
#include "larreco/Calibrator/PhotonCalibratorStandard.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes

// C++ Includes
#include <map>
#include <memory>
#include <string>

namespace {

  /**
   * @brief Adapts module FHiCL configuration to a algorithm maker tool config.
   * @param baseConfig the configuration of the module
   * @param configKey name of algorithm configuration table within `baseConfig`
   * @param algoClassPrefix prefix of the algorithm class name
   * @param algoNameKey (default: `"Name"`) name of the algorithm name config
   * @return the adapted configuration
   *
   * Input configuration is expected to include also:
   *  * `configKey` table: the configuration table passed to the actual
   *     algorithm; its content must include at least one of:
   *      * `Name` (string): the name of the algorithm; this is usually the name
   *        of the class trimmed of a prefix, e.g. `"SlidingWindow"` for
   *        `pmtana::AlgoSlidingWindow`.
   *      * `tool_type` (string): the name of the algorithm maker _art_ tool;
   *        this is usually the name of the algorithm class, with the `Maker`
   *        suffix, e.g. `"AlgoSlidingWindowMaker"` for the _art_ tool creating
   *        a `pmtana::AlgoSlidingWindow`.
   *
   * The output configuration is structured as:
   *  * `configKey`: it is the original table from the input, except that
   *    `tool_type` key is removed if it was present.
   *  * `tool_type`: a new string that contains a copy of the original
   *    `<configKey>.tool_type` if it was present, or a string made with
   *    the value of the original `<configKey>.Name`, plus the prefix in
   *    `algoClassPrefix` and the suffix `"Maker"` (e.g., if `algoClassPrefix`
   *    is `"Algo"`, name `"SlidingWindow"` will become the _art_ tool name
   *    `"AlgoSlidingWindowMaker"`).
   *    If `<configKey>.Name` is also missing, `tool_type` atom is not added,
   *    and this will result in an error when using this configuration to create
   *    an _art_ tool.
   *
   */
  fhicl::ParameterSet makeAlgoToolConfig(fhicl::ParameterSet const& baseConfig,
                                         std::string const& configKey,
                                         std::string const& algoClassPrefix = "",
                                         std::string const& algoNameKey = "Name")
  {

    fhicl::ParameterSet toolConfig;

    // add algo configuration
    fhicl::ParameterSet algoConfig = baseConfig.get<fhicl::ParameterSet>(configKey);

    if (auto toolType = algoConfig.get_if_present<std::string>("tool_type")) {
      toolConfig.put("tool_type", *toolType);
      algoConfig.erase("tool_type");
    }
    else if (auto algoName = algoConfig.get_if_present<std::string>(algoNameKey)) {
      toolConfig.put("tool_type", algoClassPrefix + *algoName + "Maker");
      // we leave the algorithm Name in its configuration for compatibility
    }
    toolConfig.put(configKey, std::move(algoConfig));

    return toolConfig;
  } // makeAlgoToolConfig()

  /**
   * @brief Adapts module FHiCL configuration to a pedestal maker tool config.
   * @param baseConfig the configuration of the module
   * @return the adapted configuration
   * @see `makeAlgoToolConfig()`
   *
   * This adapter operates simply like `makeAlgoToolConfig()`, using as
   * algorithm configuration table key `"PedAlgoPset"`.
   *
   */
  fhicl::ParameterSet makePedAlgoToolConfig(fhicl::ParameterSet const& baseConfig)
  {
    return makeAlgoToolConfig(baseConfig, "PedAlgoPset", "PedAlgo");
  }

  /**
   * @brief Adapts module FHiCL configuration to a hit finder maker tool config.
   * @param baseConfig the configuration of the module
   * @return the adapted configuration
   * @see `makeAlgoToolConfig()`
   *
   * Input configuration is expected to include also:
   *  * `HitAlgoPset` table: the configuration table passed to the actual hit
   *     finder algorithm; see `makeAlgoToolConfig()` for the requirements.
   *  * `RiseTimeCalculator` table: if present, it's the complete _art_ tool
   *     configuration for the rise time calculator algorithm.
   *
   * The output configuration is structured as:
   *  * `HitAlgoPset`: it is the original table from the input, except that
   *    `tool_type` key is removed.
   *  * `RiseTimeCalculator`: a verbatim copy of the original table; if the
   *    original table was missing, this one is not specified at all.
   *  * `tool_type`: a new string with the algorithm maker tool name;
   *    see `makeAlgoToolConfig()` for the details.
   *
   */
  fhicl::ParameterSet makeHitAlgoToolConfig(fhicl::ParameterSet const& baseConfig)
  {
    fhicl::ParameterSet toolConfig = makeAlgoToolConfig(baseConfig, "HitAlgoPset", "Algo");

    // add rise time calculator configuration ("no configuration" is ok)
    if (auto riseCalcCfg = baseConfig.get_if_present<fhicl::ParameterSet>("RiseTimeCalculator")) {
      toolConfig.put("RiseTimeCalculator", std::move(*riseCalcCfg));
    }

    return toolConfig;
  } // makeHitAlgoToolConfig()

}

namespace opdet {

  class OpHitFinder : public art::EDProducer {
  public:
    // Standard constructor and destructor for an ART module.
    explicit OpHitFinder(const fhicl::ParameterSet&);

    // The producer routine, called once per event.
    void produce(art::Event&);

  private:
    std::map<int, int> GetChannelMap();
    std::vector<double> GetSPEScales();
    std::vector<double> GetSPEShifts();

    // The parameters we'll read from the .fcl file.
    std::string fInputModule; // Input tag for OpDetWaveform collection
    std::string fGenModule;
    std::vector<std::string> fInputLabels;
    std::set<unsigned int> fChannelMasks;

    pmtana::PulseRecoManager fPulseRecoMgr;
    std::unique_ptr<pmtana::PMTPulseRecoBase> const fThreshAlg;
    std::unique_ptr<pmtana::PMTPedestalBase> const fPedAlg;

    Float_t fHitThreshold;
    unsigned int fMaxOpChannel;
    bool fUseStartTime;

    calib::IPhotonCalibrator const* fCalib = nullptr;
  };

}

namespace opdet {
  DEFINE_ART_MODULE(OpHitFinder)
}

namespace opdet {

  //----------------------------------------------------------------------------
  // Constructor
  OpHitFinder::OpHitFinder(const fhicl::ParameterSet& pset)
    : EDProducer{pset}
    , fPulseRecoMgr()
    , fThreshAlg{art::make_tool<opdet::IHitAlgoMakerTool>(makeHitAlgoToolConfig(pset))->makeAlgo()}
    , fPedAlg{art::make_tool<opdet::IPedAlgoMakerTool>(makePedAlgoToolConfig(pset))->makeAlgo()}
  {
    // Indicate that the Input Module comes from .fcl
    fInputModule = pset.get<std::string>("InputModule");
    fGenModule = pset.get<std::string>("GenModule");
    fInputLabels = pset.get<std::vector<std::string>>("InputLabels");
    fUseStartTime = pset.get<bool>("UseStartTime", false);

    for (auto const& ch :
         pset.get<std::vector<unsigned int>>("ChannelMasks", std::vector<unsigned int>()))
      fChannelMasks.insert(ch);

    fHitThreshold = pset.get<float>("HitThreshold");
    bool useCalibrator = pset.get<bool>("UseCalibrator", false);

    auto const& geometry(*lar::providerFrom<geo::Geometry>());
    fMaxOpChannel = geometry.MaxOpChannel();

    if (useCalibrator) {
      // If useCalibrator, get it from ART
      fCalib = lar::providerFrom<calib::IPhotonCalibratorService>();
    }
    else {
      // If not useCalibrator, make an internal one based
      // on fhicl settings to hit finder.
      bool areaToPE = pset.get<bool>("AreaToPE");
      float SPEArea = pset.get<float>("SPEArea");
      float SPEShift = pset.get<float>("SPEShift", 0.);

      // Reproduce behavior from GetSPEScales()
      if (!areaToPE) SPEArea = 20;

      // Delete and replace if we are reconfiguring
      if (fCalib) { delete fCalib; }

      fCalib = new calib::PhotonCalibratorStandard(SPEArea, SPEShift, areaToPE);
    }

    produces<std::vector<recob::OpHit>>();

    fPulseRecoMgr.AddRecoAlgo(fThreshAlg.get());
    fPulseRecoMgr.SetDefaultPedAlgo(fPedAlg.get());

    // show the algorithm selection on screen
    mf::LogInfo{"OpHitFinder"} << "Pulse finder algorithm: '" << fThreshAlg->Name() << "'"
                               << "\nPedestal algorithm:     '" << fPedAlg->Name() << "'";
  }

  //----------------------------------------------------------------------------
  void OpHitFinder::produce(art::Event& evt)
  {

    // These is the storage pointer we will put in the event
    std::unique_ptr<std::vector<recob::OpHit>> HitPtr(new std::vector<recob::OpHit>);

    std::vector<const sim::BeamGateInfo*> beamGateArray;
    try {
      evt.getView(fGenModule, beamGateArray);
    }
    catch (art::Exception const& err) {
      if (err.categoryCode() != art::errors::ProductNotFound) throw;
    }

    auto const& geometry(*lar::providerFrom<geo::Geometry>());
    auto const clock_data =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const& calibrator(*fCalib);
    //
    // Get the pulses from the event
    //

    // Load pulses into WaveformVector
    if (fChannelMasks.empty() && fInputLabels.size() < 2) {
      art::Handle<std::vector<raw::OpDetWaveform>> wfHandle;
      if (fInputLabels.empty())
        evt.getByLabel(fInputModule, wfHandle);
      else
        evt.getByLabel(fInputModule, fInputLabels.front(), wfHandle);
      assert(wfHandle.isValid());
      RunHitFinder(*wfHandle,
                   *HitPtr,
                   fPulseRecoMgr,
                   *fThreshAlg,
                   geometry,
                   fHitThreshold,
                   clock_data,
                   calibrator,
                   fUseStartTime);
    }
    else {

      // Reserve a large enough array
      int totalsize = 0;
      for (auto label : fInputLabels) {
        art::Handle<std::vector<raw::OpDetWaveform>> wfHandle;
        evt.getByLabel(fInputModule, label, wfHandle);
        if (!wfHandle.isValid()) continue; // Skip non-existent collections
        totalsize += wfHandle->size();
      }

      std::vector<raw::OpDetWaveform> WaveformVector;
      WaveformVector.reserve(totalsize);

      for (auto label : fInputLabels) {
        art::Handle<std::vector<raw::OpDetWaveform>> wfHandle;
        evt.getByLabel(fInputModule, label, wfHandle);
        if (!wfHandle.isValid()) continue; // Skip non-existent collections

        //WaveformVector.insert(WaveformVector.end(),
        //                      wfHandle->begin(), wfHandle->end());
        for (auto const& wf : *wfHandle) {
          if (fChannelMasks.find(wf.ChannelNumber()) != fChannelMasks.end()) continue;
          WaveformVector.push_back(wf);
        }
      }

      RunHitFinder(WaveformVector,
                   *HitPtr,
                   fPulseRecoMgr,
                   *fThreshAlg,
                   geometry,
                   fHitThreshold,
                   clock_data,
                   calibrator,
                   fUseStartTime);
    }
    // Store results into the event
    evt.put(std::move(HitPtr));
  }

} // namespace opdet

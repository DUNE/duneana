////////////////////////////////////////////////////////////////////////
//
//  PedAlgoUB source
//
////////////////////////////////////////////////////////////////////////

#include "PedAlgoUB.h"

#include "fhiclcpp/ParameterSet.h"

namespace pmtana {

  //************************************************
  PedAlgoUB::PedAlgoUB(const std::string name) : PMTPedestalBase(name), _beamgatealgo(name)
  //************************************************
  {}

  //*************************************************************
  PedAlgoUB::PedAlgoUB(const fhicl::ParameterSet& pset,
                       //PedAlgoUB::PedAlgoUB(const ::fcllite::PSet &pset,
                       const std::string name)
    : PMTPedestalBase(name)
    //, _beamgatealgo(pset.get_pset("BeamGateAlgo"),"BeamGateAlgo")
    , _beamgatealgo(pset, "BeamGateAlgo")
  //*************************************************************
  {
    _beam_gate_samples = pset.get<unsigned int>("BeamGateSamples");
  }

  //*********************************************************************
  bool PedAlgoUB::ComputePedestal(const pmtana::Waveform_t& wf,
                                  pmtana::PedestalMean_t& mean_v,
                                  pmtana::PedestalSigma_t& sigma_v)
  //*********************************************************************
  {

    if (wf.size() < _beam_gate_samples) {

      double ped_mean = wf.front(); //first sample
      double ped_sigma = 0;

      for (auto& v : mean_v)
        v = ped_mean;
      for (auto& v : sigma_v)
        v = ped_sigma;

      return true;
    }

    else {

      _beamgatealgo.Evaluate(wf);
      mean_v = _beamgatealgo.Mean();
      sigma_v = _beamgatealgo.Sigma();

      return true;
    }
  }

}

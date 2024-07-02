////////////////////////////////////////////////////////////////////////
//
//  PulseRecoManager source
//
////////////////////////////////////////////////////////////////////////

#include "PulseRecoManager.h"
#include "OpticalRecoException.h"
#include "larana/OpticalDetector/OpHitFinder/PMTPedestalBase.h"
#include "larana/OpticalDetector/OpHitFinder/PMTPulseRecoBase.h"

#include <sstream>

namespace pmtana {

  //*******************************************************
  PulseRecoManager::PulseRecoManager() : _ped_algo(nullptr)
  //*******************************************************
  {
    _reco_algo_v.clear();
  }

  //************************************************************************************
  void PulseRecoManager::AddRecoAlgo(PMTPulseRecoBase* algo, PMTPedestalBase* ped_algo)
  //************************************************************************************
  {
    if (!algo) throw OpticalRecoException("Invalid PulseReco algorithm!");

    _reco_algo_v.push_back(std::make_pair(algo, ped_algo));
  }

  //**************************************************************
  void PulseRecoManager::SetDefaultPedAlgo(PMTPedestalBase* algo)
  //**************************************************************
  {
    if (!algo) throw OpticalRecoException("Invalid Pedestal algorithm!");
    _ped_algo = algo;
  }

  //**********************************************************************
  bool PulseRecoManager::Reconstruct(const pmtana::Waveform_t& wf) const
  //**********************************************************************
  {
    if (_reco_algo_v.empty() && !_ped_algo)

      throw OpticalRecoException("No Pulse/Pedestal reconstruction to run!");

    bool ped_status = true;

    if (_ped_algo) ped_status = _ped_algo->Evaluate(wf);

    bool pulse_reco_status = ped_status;

    for (auto& algo_pair : _reco_algo_v) {

      auto& pulse_algo = algo_pair.first;
      auto& ped_algo = algo_pair.second;

      if (ped_algo) {

        ped_status = ped_status && ped_algo->Evaluate(wf);

        pulse_reco_status = (ped_status && pulse_reco_status &&
                             pulse_algo->Reconstruct(wf, ped_algo->Mean(), ped_algo->Sigma()));
      }
      else {

        if (!_ped_algo) {
          std::stringstream ss;
          ss << "No pedestal algorithm available for pulse algo " << pulse_algo->Name();
          throw OpticalRecoException(ss.str());
        }

        pulse_reco_status =
          (pulse_reco_status && pulse_algo->Reconstruct(wf, _ped_algo->Mean(), _ped_algo->Sigma()));
      }
    }

    return pulse_reco_status;
  }

}

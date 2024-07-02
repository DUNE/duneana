////////////////////////////////////////////////////////////////////////
//
//  PedAlgoEdges source
//
////////////////////////////////////////////////////////////////////////

#include "PedAlgoEdges.h"
#include "OpticalRecoException.h"
#include "UtilFunc.h"

#include "fhiclcpp/ParameterSet.h"

namespace pmtana {

  //************************************************
  PedAlgoEdges::PedAlgoEdges(const std::string name) : PMTPedestalBase(name)
  //************************************************
  {
    _nsample_front = 3;
    _nsample_tail = 5;
    _method = kHEAD;
  }

  //*************************************************************
  PedAlgoEdges::PedAlgoEdges(const fhicl::ParameterSet& pset,
                             //PedAlgoEdges::PedAlgoEdges(const ::fcllite::PSet &pset,
                             const std::string name)
    : PMTPedestalBase(name)
  //*************************************************************
  {

    _nsample_front = pset.get<size_t>("NumSampleFront");
    _nsample_tail = pset.get<size_t>("NumSampleTail");
    int method = pset.get<int>("Method");

    if (method < 0 || method > kBOTH)
      throw OpticalRecoException("PedAlgoEdges received invalid \"Method\" parameter value!");

    _method = (PED_METHOD)method;
  }

  //*********************************************************************
  bool PedAlgoEdges::ComputePedestal(const pmtana::Waveform_t& wf,
                                     pmtana::PedestalMean_t& mean_v,
                                     pmtana::PedestalSigma_t& sigma_v)
  //*********************************************************************
  {

    double ped_mean = 0;
    double ped_sigma = 0;
    switch (_method) {
    case kHEAD:
      ped_mean = mean(wf, 0, _nsample_front);
      ped_sigma = std(wf, ped_mean, 0, _nsample_front);
      for (auto& v : mean_v)
        v = ped_mean;
      for (auto& v : sigma_v)
        v = ped_sigma;
      break;
    case kTAIL:
      ped_mean = mean(wf, (wf.size() - _nsample_tail), _nsample_tail);
      ped_sigma = std(wf, ped_mean, (wf.size() - _nsample_tail), _nsample_tail);
      for (auto& v : mean_v)
        v = ped_mean;
      for (auto& v : sigma_v)
        v = ped_sigma;
      break;
    case kBOTH:
      double ped_mean_head = mean(wf, 0, _nsample_front);
      double ped_sigma_head = std(wf, ped_mean_head, 0, _nsample_front);
      double ped_mean_tail = mean(wf, (wf.size() - _nsample_tail), _nsample_tail);
      double ped_sigma_tail = std(wf, ped_mean_tail, (wf.size() - _nsample_tail), _nsample_tail);

      ped_mean = ped_mean_head;
      ped_sigma = ped_sigma_head;
      if (ped_sigma_tail < ped_sigma) {
        ped_mean = ped_mean_tail;
        ped_sigma = ped_sigma_tail;
      }
      for (auto& v : mean_v)
        v = ped_mean;
      for (auto& v : sigma_v)
        v = ped_sigma;
      break;
    }
    return true;
  }

}

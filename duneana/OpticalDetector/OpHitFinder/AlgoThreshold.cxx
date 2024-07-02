////////////////////////////////////////////////////////////////////////
//
//  AlgoThreshold source
//
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"

#include "AlgoThreshold.h"

namespace pmtana {

  //***************************************************************************
  AlgoThreshold::AlgoThreshold(const std::string name) : PMTPulseRecoBase(name)
  //***************************************************************************
  {
    //_adc_thres = 3;
    //_nsigma = 5;
    Reset();
  }

  //************************************************************
  AlgoThreshold::AlgoThreshold(const fhicl::ParameterSet& pset,
                               std::unique_ptr<pmtana::RiseTimeCalculatorBase> risetimecalculator,
                               //AlgoThreshold::AlgoThreshold(const ::fcllite::PSet &pset,
                               const std::string name)
    : PMTPulseRecoBase(name)
  //*******************************************************
  {

    _start_adc_thres = pset.get<double>("StartADCThreshold");
    _end_adc_thres = pset.get<double>("EndADCThreshold");

    //_nsigma = pset.get<double>("NSigmaThreshold");

    _nsigma_start = pset.get<double>("NSigmaThresholdStart");
    _nsigma_end = pset.get<double>("NSigmaThresholdEnd");

    _risetime_calc_ptr = std::move(risetimecalculator);

    Reset();
  }

  //***************************************************************
  void AlgoThreshold::Reset()
  //***************************************************************
  {
    PMTPulseRecoBase::Reset();
  }

  //***************************************************************
  bool AlgoThreshold::RecoPulse(const Waveform_t& wf,
                                const PedestalMean_t& mean_v,
                                const PedestalSigma_t& sigma_v)
  //***************************************************************
  {
    bool fire = false;

    double counter = 0;

    double ped_mean = mean_v.front();
    double ped_rms = sigma_v.front();

    //double threshold = ( _adc_thres > (_nsigma * ped_rms) ? _adc_thres : (_nsigma * ped_rms) );
    auto start_threshold =
      (_start_adc_thres > (_nsigma_start * ped_rms) ? _start_adc_thres : (_nsigma_start * ped_rms));
    auto end_threshold =
      (_end_adc_thres > (_nsigma_end * ped_rms) ? _end_adc_thres : (_nsigma_end * ped_rms));

    //    threshold += ped_mean

    start_threshold += ped_mean;
    end_threshold += ped_mean;

    Reset();

    for (auto const& value : wf) {

      if (!fire && ((double)value) >= start_threshold) {

        // Found a new pulse

        fire = true;

        _pulse.ped_mean = ped_mean;
        _pulse.ped_sigma = ped_rms;

        //vic: i move t_start back one, this helps with porch

        _pulse.t_start = counter - 1 > 0 ? counter - 1 : counter;
        //std::cout << "counter: " << counter << " tstart : " << _pulse.t_start << "\n";
      }

      if (fire && ((double)value) < end_threshold) {

        // Found the end of a pulse

        fire = false;

        //vic: i move t_start forward one, this helps with tail
        _pulse.t_end = counter < wf.size() ? counter : counter - 1;

        if (_risetime_calc_ptr)
          _pulse.t_rise = _risetime_calc_ptr->RiseTime(
            {wf.begin() + _pulse.t_start, wf.begin() + _pulse.t_end},
            {mean_v.begin() + _pulse.t_start, mean_v.begin() + _pulse.t_end},
            true);

        _pulse_v.push_back(_pulse);

        _pulse.reset_param();
      }

      //std::cout << "\tFire=" << fire << std::endl;

      if (fire) {

        // Add this adc count to the integral

        _pulse.area += ((double)value - (double)ped_mean);

        if (_pulse.peak < ((double)value - (double)ped_mean)) {

          // Found a new maximum

          _pulse.peak = ((double)value - (double)ped_mean);

          _pulse.t_max = counter;
        }
      }

      counter++;
    }

    if (fire) {

      // Take care of a pulse that did not finish within the readout window.

      fire = false;

      _pulse.t_end = counter - 1;

      if (_risetime_calc_ptr)
        _pulse.t_rise = _risetime_calc_ptr->RiseTime(
          {wf.begin() + _pulse.t_start, wf.begin() + _pulse.t_end},
          {mean_v.begin() + _pulse.t_start, mean_v.begin() + _pulse.t_end},
          true);

      _pulse_v.push_back(_pulse);

      _pulse.reset_param();
    }

    return true;
  }

}

////////////////////////////////////////////////////////////////////////
//
//  AlgoCFD source
//
////////////////////////////////////////////////////////////////////////

#include "AlgoCFD.h"
#include "UtilFunc.h"

#include "fhiclcpp/ParameterSet.h"

#include <unordered_map>

namespace pmtana {

  //*********************************************************************
  AlgoCFD::AlgoCFD(const std::string name) : PMTPulseRecoBase(name)
  //*********************************************************************
  {}

  //*********************************************************************
  AlgoCFD::AlgoCFD(const fhicl::ParameterSet& pset,
                   std::unique_ptr<pmtana::RiseTimeCalculatorBase> risetimecalculator,
                   //AlgoCFD::AlgoCFD(const ::fcllite::PSet &pset,
                   const std::string name)
    : PMTPulseRecoBase(name)
  //*********************************************************************
  {

    _F = pset.get<float>("Fraction");
    _D = pset.get<int>("Delay");

    //_number_presample = pset.get<int>   ("BaselinePreSample");
    _peak_thresh = pset.get<double>("PeakThresh");
    _start_thresh = pset.get<double>("StartThresh");
    _end_thresh = pset.get<double>("EndThresh");

    _risetime_calc_ptr = std::move(risetimecalculator);

    Reset();
  }

  //***************************************************************
  void AlgoCFD::Reset()
  //***************************************************************
  {
    PMTPulseRecoBase::Reset();
  }

  //***************************************************************
  bool AlgoCFD::RecoPulse(const pmtana::Waveform_t& wf,
                          const pmtana::PedestalMean_t& mean_v,
                          const pmtana::PedestalSigma_t& sigma_v)
  //***************************************************************
  {

    Reset();

    std::vector<double> cfd;
    cfd.reserve(wf.size());

    // follow cfd procedure: invert waveform, multiply by constant fraction
    // add to delayed waveform.
    for (unsigned int k = 0; k < wf.size(); ++k) {

      auto delayed = -1.0 * _F * ((float)wf.at(k) - mean_v.at(k));

      if ((int)k < _D)

        cfd.push_back(delayed);

      else

        cfd.push_back(delayed + ((float)wf.at(k - _D) - mean_v.at(k)));
    }

    // Get the zero point crossings, how can I tell which are meaningful?
    // go to each crossing, see if waveform is above pedestal (high above pedestal)

    auto crossings = LinearZeroPointX(cfd);

    // lambda criteria to determine if inside pulse

    auto in_peak = [&wf, &sigma_v, &mean_v](int i, float thresh) -> bool {
      return wf.at(i) > sigma_v.at(i) * thresh + mean_v.at(i);
    };

    // loop over CFD crossings
    for (const auto& cross : crossings) {

      if (in_peak(cross.first, _peak_thresh)) {
        _pulse.reset_param();

        int i = cross.first;

        //backwards
        while (in_peak(i, _start_thresh)) {
          i--;
          if (i < 0) {
            i = 0;
            break;
          }
        }
        _pulse.t_start = i;

        //walk a little further backwards to see if we can get 5 low RMS
        // while ( !in_peak(i,_start_thresh) ) {
        //   if (i == ( _pulse.t_start - _number_presample ) ) break;
        //   i--;
        //   if ( i < 0 ) { i = 0; break; }
        // }

        // auto before_mean = double{0.0};

        // if ( _pulse.t_start - i > 0 )
        //   before_mean = std::accumulate(std::begin(mean_v) + i,
        // 				std::begin(mean_v) + _pulse.t_start, 0.0) / ((double) (_pulse.t_start - i));

        i = _pulse.t_start + 1;

        //forwards
        while (in_peak(i, _end_thresh)) {
          i++;
          if (i > (int)(wf.size()) - 1) {
            i = (int)(wf.size()) - 1;
            break;
          }
        }

        _pulse.t_end = i;

        // //walk a little further forwards to see if we can get 5 low RMS
        // while ( !in_peak(i,_end_thresh) ) {
        //   if (i == ( _pulse.t_end + _number_presample ) ) break;
        //   i++;
        //   if ( i > wf.size() - 1 ) { i = wf.size() - 1; break; }
        // }

        // auto after_mean = double{0.0};

        // if( i - _pulse.t_end > 0)
        //   after_mean = std::accumulate(std::begin(mean_v) + _pulse.t_end + 1,
        // 			       std::begin(mean_v) + i + 1, 0.0) / ((double) (i - _pulse.t_end));

        //how to decide before or after? set before for now
        //if ( wf.size() < 1500 ) //it's cosmic discriminator
        //before_mean = mean_v.front();

        // if( after_mean <= 0 and before_mean <= 0 ) {
        //   std::cerr << "\033[93m<<" << __FUNCTION__ << ">>\033[00m Could not find good pedestal for CDF"
        // 	    << " both before_mean and after_mean are zero or less? Ignoring this crossing." << std::endl;
        //   continue;
        // }

        //x

        auto start_ped = mean_v.at(_pulse.t_start);
        auto end_ped = mean_v.at(_pulse.t_end);

        //just take the "smaller one"
        _pulse.ped_mean = start_ped <= end_ped ? start_ped : end_ped;

        if (wf.size() < 50) _pulse.ped_mean = mean_v.front(); //is COSMIC DISCRIMINATOR

        auto it = std::max_element(std::begin(wf) + _pulse.t_start, std::begin(wf) + _pulse.t_end);

        _pulse.t_max = it - std::begin(wf);
        _pulse.peak = *it - _pulse.ped_mean;
        _pulse.t_cfdcross = cross.second;

        for (auto k = _pulse.t_start; k <= _pulse.t_end; ++k) {
          auto a = wf.at(k) - _pulse.ped_mean;
          if (a > 0) _pulse.area += a;
        }

        if (_risetime_calc_ptr)
          _pulse.t_rise = _risetime_calc_ptr->RiseTime(
            {wf.begin() + _pulse.t_start, wf.begin() + _pulse.t_end},
            {mean_v.begin() + _pulse.t_start, mean_v.begin() + _pulse.t_end},
            true);

        _pulse_v.push_back(_pulse);
      }
    }

    // Vic:
    // Very close in time pulses have multiple CFD
    // crossing points. Should we check that pulses now have
    // some multiplicity? No lets just delete them.

    auto pulses_copy = _pulse_v;
    _pulse_v.clear();

    std::unordered_map<unsigned, pulse_param> delta;

    //unsigned width = 0;
    for (const auto& p : pulses_copy) {

      if (delta.count(p.t_start)) {
        if ((p.t_end - p.t_start) > (delta[p.t_start].t_end - delta[p.t_start].t_start))
          delta[p.t_start] = p;
        else
          continue;
      }
      else {
        delta[p.t_start] = p;
      }
    }

    for (const auto& p : delta)
      _pulse_v.push_back(p.second);

    //do the same now ensure t_final's are all unique
    //width = 0;

    pulses_copy.clear();
    pulses_copy = _pulse_v;

    _pulse_v.clear();
    delta.clear();

    for (const auto& p : pulses_copy) {

      if (delta.count(p.t_end)) {
        if ((p.t_end - p.t_start) > (delta[p.t_end].t_end - delta[p.t_end].t_start))
          delta[p.t_end] = p;
        else
          continue;
      }
      else {
        delta[p.t_end] = p;
      }
    }

    for (const auto& p : delta)
      _pulse_v.push_back(p.second);

    //there should be no overlapping pulses now...

    return true;
  }

  // currently returns ALL zero point crossings, we really just want ones associated with peak...
  const std::map<unsigned, double> AlgoCFD::LinearZeroPointX(const std::vector<double>& trace)
  {

    std::map<unsigned, double> crossing;

    //step through the trace and find where slope is POSITIVE across zero
    for (unsigned i = 0; i < trace.size() - 1; ++i) {

      auto si = ::pmtana::sign(trace.at(i));
      auto sf = ::pmtana::sign(trace.at(i + 1));

      if (si == sf) //no sign flip, no zero cross
        continue;

      if (sf < si) //this is a negative slope, continue
        continue;

      //calculate the crossing X based on linear interpolation bt two pts

      crossing[i] = (double)i - trace.at(i) * (1.0 / (trace.at(i + 1) - trace.at(i)));
    }

    return crossing;
  }

}

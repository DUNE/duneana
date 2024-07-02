////////////////////////////////////////////////////////////////////////
//
//  AlgoSlidingWindow source
//
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"

#include "AlgoSlidingWindow.h"

namespace pmtana {

  //*********************************************************************
  AlgoSlidingWindow::AlgoSlidingWindow(const std::string name) : PMTPulseRecoBase(name)
  //*********************************************************************
  {}

  //*********************************************************************
  AlgoSlidingWindow::AlgoSlidingWindow(
    const fhicl::ParameterSet& pset,
    std::unique_ptr<pmtana::RiseTimeCalculatorBase> risetimecalculator,
    //AlgoSlidingWindow::AlgoSlidingWindow(const ::fcllite::PSet &pset,
    const std::string name)
    : PMTPulseRecoBase(name)
  //*********************************************************************
  {
    _positive = pset.get<bool>("PositivePolarity", true);

    _adc_thres = pset.get<float>("ADCThreshold");

    _tail_adc_thres = pset.get<float>("TailADCThreshold", _adc_thres);

    _end_adc_thres = pset.get<float>("EndADCThreshold");

    _nsigma = pset.get<float>("NSigmaThreshold");

    _tail_nsigma = pset.get<float>("TailNSigma", _nsigma);

    _end_nsigma = pset.get<float>("EndNSigmaThreshold");

    _verbose = pset.get<bool>("Verbosity");

    _num_presample = pset.get<size_t>("NumPreSample");

    _num_postsample = pset.get<size_t>("NumPostSample", 0);

    _min_width = pset.get<size_t>("MinPulseWidth", 0);

    _risetime_calc_ptr = std::move(risetimecalculator);

    Reset();
  }

  //***************************************************************
  void AlgoSlidingWindow::Reset()
  //***************************************************************
  {
    PMTPulseRecoBase::Reset();
  }

  //***************************************************************
  bool AlgoSlidingWindow::RecoPulse(const pmtana::Waveform_t& wf,
                                    const pmtana::PedestalMean_t& mean_v,
                                    const pmtana::PedestalSigma_t& sigma_v)
  //***************************************************************
  {

    bool fire = false;

    bool in_tail = false;

    bool in_post = false;

    double pulse_tail_threshold = 0;

    double pulse_end_threshold = 0;

    double pulse_start_baseline = 0;

    int post_integration = 0;

    assert(wf.size() == mean_v.size() && wf.size() == sigma_v.size());

    //double threshold = ( _adc_thres > (_nsigma * _ped_rms) ? _adc_thres : (_nsigma * _ped_rms) );

    //threshold += _ped_mean;

    Reset();

    for (size_t i = 0; i < wf.size(); ++i) {

      double value = 0.;
      if (_positive)
        value = ((double)(wf[i])) - mean_v[i];
      else
        value = mean_v[i] - ((double)(wf[i]));

      float start_threshold = 0.;
      float tail_threshold = 0.;
      if (sigma_v[i] * _nsigma < _adc_thres)
        start_threshold = _adc_thres;
      else
        start_threshold = sigma_v[i] * _nsigma;

      if (sigma_v[i] * _tail_nsigma < _tail_adc_thres)
        tail_threshold = _tail_adc_thres;
      else
        tail_threshold = sigma_v[i] * _tail_nsigma;

      // End pulse if significantly high peak found (new pulse)
      if ((!fire || in_tail || in_post) && ((double)value > start_threshold)) {

        // If there's a pulse, end it
        if (in_tail) {
          _pulse.t_end = i - 1;

          // Register if width is acceptable
          if ((_pulse.t_end - _pulse.t_start) >= _min_width) {
            if (_risetime_calc_ptr)
              _pulse.t_rise = _risetime_calc_ptr->RiseTime(
                {wf.begin() + _pulse.t_start, wf.begin() + _pulse.t_end},
                {mean_v.begin() + _pulse.t_start, mean_v.begin() + _pulse.t_end},
                _positive);

            _pulse_v.push_back(_pulse);
          }

          _pulse.reset_param();

          if (_verbose)
            std::cout << "\033[93mPulse End\033[00m: "
                      << "baseline: " << mean_v[i] << " ... "
                      << " ... adc above: " << value << " T=" << i << std::endl;
        }

        //
        // Found a new pulse ... try to get a few samples prior to this
        //

        pulse_tail_threshold = tail_threshold;
        pulse_start_baseline = mean_v[i];

        pulse_end_threshold = 0.;
        if (sigma_v[i] * _end_nsigma < _end_adc_thres)
          pulse_end_threshold = _end_adc_thres;
        else
          pulse_end_threshold = sigma_v[i] * _end_nsigma;

        int buffer_num_index = 0;
        if (_pulse_v.size())
          buffer_num_index = (int)i - _pulse_v.back().t_end - 1;
        else
          buffer_num_index = std::min(_num_presample, i);

        if (buffer_num_index > (int)_num_presample) buffer_num_index = _num_presample;

        if (buffer_num_index < 0) {
          std::cerr << "\033[95m[ERROR]\033[00m Logic error! Negative buffer_num_index..."
                    << std::endl;
          throw std::exception();
        }

        // If there's a pulse, end we where in in_post, end the previous pulse first
        if (in_post) {
          // Find were
          _pulse.t_end = static_cast<int>(i) - buffer_num_index;
          if (_pulse.t_end > 0) --_pulse.t_end; // leave a gap, if we can

          // Register if width is acceptable
          if ((_pulse.t_end - _pulse.t_start) >= _min_width) {
            if (_risetime_calc_ptr)
              _pulse.t_rise = _risetime_calc_ptr->RiseTime(
                {wf.begin() + _pulse.t_start, wf.begin() + _pulse.t_end},
                {mean_v.begin() + _pulse.t_start, mean_v.begin() + _pulse.t_end},
                _positive);

            _pulse_v.push_back(_pulse);
          }

          _pulse.reset_param();

          if (_verbose)
            std::cout << "\033[93mPulse End\033[00m: new pulse starts during in_post: "
                      << "baseline: " << mean_v[i] << " ... "
                      << " ... adc above: " << value << " T=" << i << std::endl;
        }

        _pulse.t_start = i - buffer_num_index;
        _pulse.ped_mean = pulse_start_baseline;
        _pulse.ped_sigma = sigma_v[i];

        for (size_t pre_index = _pulse.t_start; pre_index < i; ++pre_index) {

          double pre_adc = wf[pre_index];
          if (_positive)
            pre_adc -= pulse_start_baseline;
          else
            pre_adc = pulse_start_baseline - pre_adc;

          if (pre_adc > 0.) _pulse.area += pre_adc;
        }

        if (_verbose)
          std::cout << "\033[93mPulse Start\033[00m: "
                    << "baseline: " << mean_v[i] << " ... threshold: " << start_threshold
                    << " ... adc above baseline: " << value << " ... pre-adc sum: " << _pulse.area
                    << " T=" << i << std::endl;

        fire = true;
        in_tail = false;
        in_post = false;
      }

      if (fire && value < pulse_tail_threshold) {
        fire = false;
        in_tail = true;
        in_post = false;
      }

      if ((fire || in_tail || in_post) && _verbose) {
        std::cout << (fire ? "\033[93mPulsing\033[00m: " : "\033[93mIn-tail\033[00m: ")
                  << "baseline: " << mean_v[i] << " std: " << sigma_v[i]
                  << " ... adc above baseline " << value << " T=" << i << std::endl;
      }

      if ((fire || in_tail) && value < pulse_end_threshold) {
        in_post = true;
        fire = in_tail = false;
        post_integration = _num_postsample;
      }

      if (in_post && post_integration < 1) {
        // Found the end of a pulse
        _pulse.t_end = i - 1;

        // Register if width is acceptable
        if ((_pulse.t_end - _pulse.t_start) >= _min_width) {
          if (_risetime_calc_ptr)
            _pulse.t_rise = _risetime_calc_ptr->RiseTime(
              {wf.begin() + _pulse.t_start, wf.begin() + _pulse.t_end},
              {mean_v.begin() + _pulse.t_start, mean_v.begin() + _pulse.t_end},
              _positive);

          _pulse_v.push_back(_pulse);
        }

        if (_verbose)
          std::cout << "\033[93mPulse End\033[00m: "
                    << "baseline: " << mean_v[i] << " ... adc: " << value << " T=" << i
                    << " ... area sum " << _pulse.area << std::endl;

        _pulse.reset_param();

        fire = false;
        in_tail = false;
        in_post = false;
      }

      if (fire || in_tail || in_post) {

        //_pulse.area += ((double)value - (double)mean_v[i]);
        _pulse.area += value;

        if (_pulse.peak < value) {

          // Found a new maximum
          _pulse.peak = value;

          _pulse.t_max = i;
        }

        if (in_post) --post_integration;
      }
    }

    if (fire || in_tail || in_post) {

      // Take care of a pulse that did not finish within the readout window.

      fire = false;
      in_tail = false;

      _pulse.t_end = wf.size() - 1;

      // Register if width is acceptable
      if ((_pulse.t_end - _pulse.t_start) >= _min_width) {
        if (_risetime_calc_ptr)
          _pulse.t_rise = _risetime_calc_ptr->RiseTime(
            {wf.begin() + _pulse.t_start, wf.begin() + _pulse.t_end},
            {mean_v.begin() + _pulse.t_start, mean_v.begin() + _pulse.t_end},
            _positive);
        _pulse_v.push_back(_pulse);
      }

      _pulse.reset_param();
    }

    return true;
  }

}

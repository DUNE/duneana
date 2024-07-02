////////////////////////////////////////////////////////////////////////
//
//  PedAlgoRmsSlider source
//
////////////////////////////////////////////////////////////////////////

#include "PedAlgoRmsSlider.h"
#include "OpticalRecoException.h"
#include "UtilFunc.h"
#include "fhiclcpp/ParameterSet.h"

#include <fstream>
#include <iostream>
#include <numeric>

namespace pmtana {

  //*****************************************************************
  PedAlgoRmsSlider::PedAlgoRmsSlider(const std::string name) : PMTPedestalBase(name)
  //*****************************************************************
  {
    srand(static_cast<unsigned int>(time(0)));
  }

  //**************************************************************************
  PedAlgoRmsSlider::PedAlgoRmsSlider(const fhicl::ParameterSet& pset, const std::string name)
    : PMTPedestalBase(name)
  //############################################################
  {

    _sample_size = pset.get<size_t>("SampleSize", 7);
    _threshold = pset.get<double>("Threshold", 0.6);
    _max_sigma = pset.get<float>("MaxSigma", 0.5);
    _ped_range_max = pset.get<float>("PedRangeMax", 2150);
    _ped_range_min = pset.get<float>("PedRangeMin", 100);
    _num_presample = pset.get<int>("NumPreSample", 0);
    _num_postsample = pset.get<int>("NumPostSample", 0);
    _verbose = pset.get<bool>("Verbose", true);
    _n_wf_to_csvfile = pset.get<int>("NWaveformsToFile", 12);

    if (_n_wf_to_csvfile > 0) {
      _csvfile.open("wf_pedalgormsslider.csv", std::ofstream::out | std::ofstream::trunc);
      _csvfile << "n,time,wf,wf_ped_mean,wf_ped_rms" << std::endl;
    }
  }

  //*******************************************
  void PedAlgoRmsSlider::PrintInfo()
  //*******************************************
  {
    std::cout << "PedAlgoRmsSlider setting:"
              << "\n\t SampleSize:       " << _sample_size
              << "\n\t Threshold:        " << _threshold << "\n\t Verbose:          " << _verbose
              << "\n\t NWaveformsToFile: " << _n_wf_to_csvfile << std::endl;
  }

  //****************************************************************************
  double PedAlgoRmsSlider::CalcMean(const std::vector<double>& wf, size_t start, size_t nsample)
  //****************************************************************************
  {
    if (!nsample) nsample = wf.size();
    if (start > wf.size() || (start + nsample) > wf.size())
      throw OpticalRecoException("Invalid start/end index!");

    double sum = std::accumulate(wf.begin() + start, wf.begin() + start + nsample, 0.0);

    sum /= ((double)nsample);

    return sum;
  }

  //****************************************************************************
  double PedAlgoRmsSlider::CalcStd(const std::vector<double>& wf,
                                   const double ped_mean,
                                   size_t start,
                                   size_t nsample)
  //****************************************************************************
  {
    if (!nsample) nsample = wf.size();
    if (start > wf.size() || (start + nsample) > wf.size())
      throw OpticalRecoException("Invalid start/end index!");

    double sigma = 0;

    for (size_t index = start; index < (start + nsample); ++index) {
      sigma += pow((wf[index] - ped_mean), 2);
    }

    sigma = sqrt(sigma / ((double)(nsample)));

    return sigma;
  }

  //****************************************************************************
  bool PedAlgoRmsSlider::ComputePedestal(const pmtana::Waveform_t& wf,
                                         pmtana::PedestalMean_t& mean_v,
                                         pmtana::PedestalSigma_t& sigma_v)
  //****************************************************************************
  {

    if (_verbose) this->PrintInfo();

    if (wf.size() <= (_sample_size * 2)) return false;

    // Prepare output
    mean_v.resize(wf.size(), 0);
    sigma_v.resize(wf.size(), 0);

    // **********
    // To start, set the pedestal equal to
    // the wf itself
    // **********

    pmtana::PedestalMean_t mean_temp_v;
    mean_temp_v.resize(wf.size(), 0);

    for (size_t i = 0; i < wf.size(); ++i) {
      mean_temp_v[i] = wf[i];
      sigma_v[i] = 0;
    }

    // **********
    // Now look for rms variations
    // and change the mean and rms accordingly
    // **********
    int last_good_index = -1;
    double local_mean, local_rms;
    std::vector<double> local_mean_v(wf.size(), -1.);
    std::vector<double> local_sigma_v(wf.size(), -1.);

    for (size_t i = 0; i < wf.size() - _sample_size; i++) {

      local_mean = mean(wf, i, _sample_size);
      local_rms = std(wf, local_mean, i, _sample_size);

      if (_verbose)
        std::cout << "\033[93mPedAlgoRmsSlider\033[00m: i " << i << "  local_mean: " << local_mean
                  << "  local_rms: " << local_rms << std::endl;

      if (local_rms < _threshold) {

        local_mean_v[i] = local_mean;
        local_sigma_v[i] = local_rms;

        if (_verbose)
          std::cout << "\033[93mBelow threshold\033[00m: "
                    << "at i " << i << " last good index was: " << last_good_index << std::endl;
      }
    }

    // find the gaps (regions to be interpolated
    last_good_index = -1;
    std::vector<bool> ped_interapolated(wf.size(), false);
    for (size_t i = 0; i < wf.size() - _sample_size; i++) {

      if (local_mean_v[i] > -0.1) {
        // good pedestal!

        if ((last_good_index + 1) < (int)i) {
          // finished the gap. try interpolation
          // 0) find where to start/end interpolation
          int start_tick = last_good_index;
          int end_tick = i;
          int start_bound = std::max(last_good_index - _num_presample, 0);
          int end_bound = std::min(i + _num_postsample, (int)(wf.size()) - _sample_size);
          for (int j = start_tick; j >= start_bound; --j) {
            if (local_mean_v[j] < 0) continue;
            start_tick = j;
          }
          for (int j = end_tick; j <= end_bound; ++j) {
            if (local_mean_v[j] < 0) continue;
            end_tick = j;
          }

          //this should become generic interpolation function, for now lets leave.
          float slope =
            (local_mean_v[end_tick] - local_mean_v[start_tick]) / (float(end_tick - start_tick));

          for (int j = start_tick + 1; j < end_tick; ++j) {
            mean_temp_v[j] = slope * (float(j - start_tick)) + local_mean_v[start_tick];
            // for sigma, put either the sigma in the region before the pulse or
            // after the pulse, depending on which one if != 0. If both are !=0 put the one after
            // the pulse (just default), if both are zero then put zero
            sigma_v[j] = (local_sigma_v[end_tick] != 0 ?
                            local_sigma_v[end_tick] :
                            local_sigma_v[start_tick]); // todo: fluctuate baseline
            ped_interapolated[j] = true;
          }
        }

        last_good_index = i;
      }
    }

    /*

    for (size_t i = 0; i < wf.size() - _sample_size; i++) {

      local_mean = mean(wf, i, _sample_size);
      local_rms  = std(wf, local_mean, i, _sample_size);

      if(_verbose) std::cout << "\033[93mPedAlgoRmsSlider\033[00m: i " << i << "  local_mean: " << local_mean << "  local_rms: " << local_rms << std::endl;

      if (local_rms < _threshold) {

	local_mean_v[i] = local_mean;
	local_sigma_v[i] = local_sigma;

        if(_verbose)
          std::cout << "\033[93mBelow threshold\033[00m: "
                    << "at i " << i
                    << " last good index was: " << last_good_index
                    << std::endl;

        if(last_good_index<0) {
          last_good_index = (int)i;
          last_local_mean = local_mean;
          last_local_rms  = local_rms;
          continue;
        }


        if( ( last_good_index + 1 ) < (int)i ) {

          //this should become generic interpolation function, for now lets leave.
          float slope = (local_mean - last_local_mean) / (float(i - last_good_index));

          for(size_t j = last_good_index + 1; j < i && j < wf.size(); ++j) {
            mean_temp_v.at(j)  = slope * ( float(j - last_good_index) ) + mean_temp_v.at(last_good_index);
            // for sigma, put either the sigma in the region before the pulse or
            // after the pulse, depending on which one if != 0. If both are !=0 put the one after
            // the pulse (just default), if both are zero then put zero
            sigma_v.at(j) = (local_rms != 0 ? local_rms : last_local_rms); // todo: fluctuate baseline
            ped_interapolated.at(j) = true;
          }
        }

	// record this mean & rms as good mean value
        last_good_index = i;
        last_local_mean = local_mean;
        last_local_rms  = local_rms;
	// if _num_postsample is specified, go back in time to look for it
	if(_num_postsample >0 && (i>_num_postsample)) {
	  int loop_start = std::max(((int)i) - _num_postsample, 0);
	  for(int j=loop_start; j>=0; --j) {
	    if(local_mean_v[j] <0) continue;
	    last_good_index = j;
	    last_local_mean = local_mean_v[j];
	    last_local_rms  = local_sigma_v[j];
	    break;
	  }
	}
      }


    }
     */

    // **********
    // Now look at special cases, if wf starts or
    // ends with a pulse
    // **********

    // At start

    bool end_found = false;

    local_mean = mean(wf, 0, _sample_size);
    local_rms = std(wf, local_mean, 0, _sample_size);

    if (local_rms >= _threshold) {

      for (size_t i = 1; i < wf.size() - _sample_size; i++) {

        local_mean = mean(wf, i, _sample_size);
        local_rms = std(wf, local_mean, i, _sample_size);

        if (local_rms < _threshold) {

          end_found = true;

          for (size_t j = 0; j < i; j++) {
            mean_temp_v[j] = local_mean;
            sigma_v[j] = local_rms;
            ped_interapolated[j] = true;
          }
          break;
        }
      }

      if (!end_found) {
        std::cerr << "\033[93m<<" << __FUNCTION__
                  << ">>\033[00m Could not find good pedestal for CDF"
                  << "There is pulse on first sample and baseline never went back down. Returning "
                     "false here.";
        return false;
      }
    }

    // At end

    bool start_found = false;

    local_mean = mean(wf, wf.size() - 1 - _sample_size, _sample_size);
    local_rms = std(wf, local_mean, wf.size() - 1 - _sample_size, _sample_size);

    if (local_rms >= _threshold) {

      size_t i = wf.size() - 1 - _sample_size;
      while (i-- > 0) {
        local_mean = mean(wf, i, _sample_size);
        local_rms = std(wf, local_mean, i, _sample_size);

        if (local_rms < _threshold) {

          start_found = true;

          for (size_t j = wf.size() - 1; j > i; j--) {
            mean_temp_v[j] = local_mean;
            sigma_v[j] = local_rms;
            ped_interapolated[j] = true;
          }
          break;
        }
      }

      if (!start_found) {
        std::cerr << "\033[93m<<" << __FUNCTION__
                  << ">>\033[00m Could not find good pedestal for CDF"
                  << "There is pulse on last sample and baseline never went back down. Returning "
                     "false here.";
        return false;
      }
    }

    // **********
    // Now smooth it to estimate the final pedestal
    // **********

    const size_t window_size = _sample_size * 2;

    // middle mean
    for (size_t i = 0; i < mean_temp_v.size(); ++i) {

      if (i < _sample_size || i >= (wf.size() - _sample_size)) continue;

      mean_v[i] = this->CalcMean(mean_temp_v, i - _sample_size, window_size);
      if (!ped_interapolated[i]) {
        sigma_v[i] = this->CalcStd(mean_temp_v, mean_v[i], i - _sample_size, window_size);
      }
    }

    // front mean
    for (size_t i = 0; i < _sample_size; ++i) {

      mean_v[i] = mean_v[_sample_size];
      if (!ped_interapolated[i]) { sigma_v[i] = sigma_v[_sample_size]; }
    }

    // tail mean
    for (size_t i = (mean_temp_v.size() - _sample_size); i < mean_temp_v.size(); ++i) {

      mean_v[i] = mean_v[wf.size() - _sample_size - 1];
      if (!ped_interapolated[i]) { sigma_v[i] = sigma_v[wf.size() - _sample_size - 1]; }
    }

    // Save to file
    if (_wf_saved + 1 <= _n_wf_to_csvfile) {
      _wf_saved++;
      for (size_t i = 0; i < wf.size(); i++) {
        _csvfile << _wf_saved - 1 << "," << i << "," << wf[i] << "," << mean_v[i] << ","
                 << sigma_v[i] << std::endl;
      }
    }

    bool is_sane = this->CheckSanity(mean_v, sigma_v);

    return is_sane;
  }

  //*******************************************
  bool PedAlgoRmsSlider::CheckSanity(pmtana::PedestalMean_t& mean_v,
                                     pmtana::PedestalSigma_t& sigma_v)
  //*******************************************
  {

    float best_sigma = 1.1e9;
    size_t best_sigma_index = 0;
    size_t num_good_adc = 0;

    for (size_t i = 0; i < sigma_v.size(); ++i) {
      // Only consider adcs which mean is in the allowed range
      auto const& mean = mean_v[i];

      if (mean < _ped_range_min || mean > _ped_range_max) continue;

      auto const& sigma = sigma_v[i];
      if (sigma < best_sigma) {
        best_sigma = sigma;
        best_sigma_index = i;
      }

      if (sigma < _max_sigma) num_good_adc += 1;
    }

    if (num_good_adc < 1) {
      std::cerr << "\033[93m<<" << __FUNCTION__
                << ">>\033[00m Could not find good pedestal at all..." << std::endl;
      return false;
    }

    // If not enough # of good mean indices, use the best guess within this waveform
    if (best_sigma > _max_sigma || num_good_adc < 3) {

      if (_verbose) {
        std::cout << "\033[93mPedAlgoRmsSlider\033[00m: Not enough number of good mean indices."
                  << "Using the best guess within this waveform." << std::endl;
      }

      for (size_t i = 0; i < mean_v.size(); ++i) {
        mean_v[i] = mean_v[best_sigma_index];
        sigma_v[i] = sigma_v[best_sigma_index];
      }
    }

    return true;
  }
}

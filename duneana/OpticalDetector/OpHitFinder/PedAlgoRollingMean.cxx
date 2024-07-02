////////////////////////////////////////////////////////////////////////
//
//  PedAlgoRollingMean source
//
////////////////////////////////////////////////////////////////////////

#include "PedAlgoRollingMean.h"
#include "UtilFunc.h"
#include "fhiclcpp/ParameterSet.h"

#include <iostream>

namespace pmtana {

  //*****************************************************************
  PedAlgoRollingMean::PedAlgoRollingMean(const std::string name) : PMTPedestalBase(name)
  //*****************************************************************
  {
    srand(static_cast<unsigned int>(time(0)));
  }

  //**************************************************************************
  PedAlgoRollingMean::PedAlgoRollingMean(
    const fhicl::ParameterSet& pset,
    //PedAlgoRollingMean::PedAlgoRollingMean(const ::fcllite::PSet &pset,
    const std::string name)
    : PMTPedestalBase(name)
  //############################################################
  {

    _sample_size = pset.get<size_t>("SampleSize");
    _max_sigma = pset.get<float>("MaxSigma");
    _ped_range_max = pset.get<float>("PedRangeMax");
    _ped_range_min = pset.get<float>("PedRangeMin");

    // _range          = pset.get<int>   ("RandomRange");
    // _divisions      = pset.get<double>("RandomRangeDivisions");
    _threshold = pset.get<double>("Threshold");
    _diff_threshold = pset.get<double>("DiffBetweenGapsThreshold");
    _diff_adc_count = pset.get<double>("DiffADCCounts");

    _n_presamples = pset.get<int>("NPrePostSamples");

    //_random_shift   = pset.get<double>("RandomRangeShift");
    // Random seed number generator
    //srand(static_cast<unsigned int>(time(0)));
  }

  //****************************************************************************
  bool PedAlgoRollingMean::ComputePedestal(const pmtana::Waveform_t& wf,
                                           pmtana::PedestalMean_t& mean_v,
                                           pmtana::PedestalSigma_t& sigma_v)
  //****************************************************************************
  {

    // parameters
    if (wf.size() <= (_sample_size * 2)) return false;

    const size_t window_size = _sample_size * 2;

    // middle mean
    for (size_t i = 0; i < wf.size(); ++i) {

      mean_v[i] = 0;
      sigma_v[i] = 0;

      if (i < _sample_size || i >= (wf.size() - _sample_size)) continue;

      mean_v[i] = mean(wf, i - _sample_size, window_size);
      sigma_v[i] = std(wf, mean_v[i], i - _sample_size, window_size);
    }

    // front mean
    for (size_t i = 0; i < _sample_size; ++i) {

      mean_v[i] = mean_v[_sample_size];
      sigma_v[i] = sigma_v[_sample_size];
    }
    // tail mean
    for (size_t i = (wf.size() - _sample_size); i < wf.size(); ++i) {

      mean_v[i] = mean_v[wf.size() - _sample_size - 1];
      sigma_v[i] = sigma_v[wf.size() - _sample_size - 1];
    }

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
      for (size_t i = 0; i < mean_v.size(); ++i) {
        mean_v[i] = mean_v.at(best_sigma_index);
        sigma_v[i] = sigma_v.at(best_sigma_index);
      }

      return true;
    }

    // Else do extrapolation, or random seed depending on what we find...

    unsigned nbins = 1000;

    //////////////////seg faulting...
    const auto mode_mean = BinnedMaxOccurrence(mean_v, nbins);
    const auto mode_sigma = BinnedMaxOccurrence(sigma_v, nbins);

    //auto mode_mean  = BinnedMaxTH1D(mean_v ,nbins);
    //auto mode_sigma = BinnedMaxTH1D(sigma_v,nbins);

    //std::cout<<mode_mean<<" +/- "<<mode_sigma<<std::endl;

    double const diff_cutoff = std::max(_diff_threshold * mode_sigma, _diff_adc_count);

    int last_good_index = -1;

    for (size_t i = 0; i < wf.size(); ++i) {

      auto const mean = mean_v[i];
      auto const sigma = sigma_v[i];

      // if(sigma <= _max_sigma && mean < _ped_range_max && mean > _ped_range_min) {
      // not sure if this works well for basline that is "linear" seen by David K

      if (sigma <= _threshold * mode_sigma && fabs(mean - mode_mean) <= _threshold * mode_sigma) {

        if (last_good_index < 0) {
          last_good_index = (int)i;
          continue;
        }

        if ((last_good_index + 1) < (int)i) {

          auto diff = fabs(mean_v.at(last_good_index) - mean);

          if (diff > diff_cutoff) {
            //this should become generic interpolation function, for now lets leave.
            float slope = (mean - mean_v.at(last_good_index)) / (float(i - last_good_index));

            for (size_t j = last_good_index + 1; j < i; ++j) {
              mean_v.at(j) = slope * (float(j - last_good_index)) + mean_v.at(last_good_index);
              sigma_v.at(j) = mode_sigma;
            }
          }
          else {
            //difference roughly the same lets fill with constant baseline

            auto presample_mean =
              edge_aware_mean(wf, last_good_index - _n_presamples, last_good_index);
            auto postsample_mean = edge_aware_mean(wf, i, i + _n_presamples);

            auto diff_pre = fabs(presample_mean - mode_mean);
            auto diff_post = fabs(postsample_mean - mode_mean);

            auto constant_mean = diff_pre <= diff_post ? presample_mean : postsample_mean;

            for (size_t j = last_good_index + 1; j < i; ++j) {
              //mean_v.at(j)  = floor( mean_v.at(last_good_index) ) + _random_shift + (double) ( rand() % _range) / _divisions;
              mean_v.at(j) = constant_mean;
              sigma_v.at(j) = mode_sigma;
            }
          }
        }
        last_good_index = (int)i;
      }
    }

    // Next do extrapolation to the first and end (if needed)
    // vic: for now we leave this i'm not sure if this really needs
    //      to be tuned until I can make unit test
    // update: yes this needs work...

    if (sigma_v.front() > mode_sigma) {

      int first_index = -1;
      int second_index = -1;

      for (size_t i = 0; i < wf.size(); ++i) {
        if (sigma_v.at(i) < mode_sigma) {
          if (first_index < 0)
            first_index = (int)i;
          else if (second_index < 0) {
            second_index = (int)i;
            break;
          }
        }
      }

      if (first_index < 0 || second_index < 0) {
        std::cerr << "\033[93m<<" << __FUNCTION__
                  << ">>\033[00m Could not find good pedestal for CDF"
                  << "\n"
                  << "first_index:  " << first_index << "\n"
                  << "second_index: " << second_index << "\n"
                  << "If one of these is less than zero, this means there is pulse\n"
                  << "on first sample and baseline never went back down. Returning false here.";
        return false;
      }

      auto diff = fabs(mean_v.at(first_index) - mean_v.at(second_index));

      if (diff > diff_cutoff) {

        float slope =
          (mean_v.at(second_index) - mean_v.at(first_index)) / (float(second_index - first_index));

        for (int j = 0; j < first_index; ++j) {
          mean_v.at(j) = mean_v.at(first_index) - slope * (first_index - j);
          sigma_v.at(j) = _max_sigma;
        }
      }
      else {

        auto postsample_mean = edge_aware_mean(wf, first_index, first_index + _n_presamples);

        for (int j = 0; j < first_index; ++j) {
          //mean_v.at(j)  = floor( mean_v.at(second_index) ) + _random_shift + (double) ( rand() % _range) / _divisions;
          mean_v.at(j) = postsample_mean;
          sigma_v.at(j) = mode_sigma;
        }
      }
    }

    if (sigma_v.back() > mode_sigma) {

      int first_index = -1;
      int second_index = -1;

      for (int i = wf.size() - 1; i >= 0; --i) {
        if (sigma_v.at(i) < mode_sigma) {
          if (second_index < 0)
            second_index = (int)i;
          else if (first_index < 0) {
            first_index = (int)i;
            break;
          }
        }
      }

      if (first_index < 0 || second_index < 0) {
        std::cerr << "\033[93m<<" << __FUNCTION__
                  << ">>\033[00m Could not find good pedestal for CDF"
                  << "\n"
                  << "first_index:  " << first_index << "\n"
                  << "second_index: " << second_index << "\n"
                  << "If one of these is less than zero, this means there is pulse\n"
                  << "on the last sample and baseline never went back down. Returning false here.";
        return false;
      }

      auto diff = fabs(mean_v.at(second_index) - mean_v.at(first_index));

      if (diff > diff_cutoff) {

        float slope =
          (mean_v.at(second_index) - mean_v.at(first_index)) / (float(second_index - first_index));
        for (int j = second_index + 1; j < int(wf.size()); ++j) {
          mean_v.at(j) = mean_v.at(second_index) + slope * (j - second_index);
          sigma_v.at(j) = _max_sigma;
        }
      }
      else {

        auto presample_mean = edge_aware_mean(wf, second_index - _n_presamples, second_index);

        for (int j = second_index + 1; j < int(wf.size()); ++j) {
          //mean_v.at(j)  = floor( mean_v.at(first_index) ) + _random_shift + (double) ( rand() % _range) / _divisions;
          mean_v.at(j) = presample_mean;
          sigma_v.at(j) = mode_sigma;
        }
      }
    }

    return true;
  }

}

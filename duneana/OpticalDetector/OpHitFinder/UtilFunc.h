#ifndef larana_OPTICALDETECTOR_UTILFUNC_H
#define larana_OPTICALDETECTOR_UTILFUNC_H

#include "OpticalRecoTypes.h"

#include <cstddef>
#include <vector>

namespace pmtana {

  double mean(const std::vector<short>& wf, size_t start = 0, size_t nsample = 0);

  double edge_aware_mean(const std::vector<short>& wf, int start, int end);

  double std(const std::vector<short>& wf,
             const double ped_mean,
             size_t start = 0,
             size_t nsample = 0);

  double BinnedMaxOccurrence(const PedestalMean_t& mean_v, const size_t nbins);

  double BinnedMaxTH1D(const std::vector<double>& v, int bins);

  int sign(double val);

}

#endif

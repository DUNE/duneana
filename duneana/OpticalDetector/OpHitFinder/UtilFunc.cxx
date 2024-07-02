#include "UtilFunc.h"
#include "OpticalRecoException.h"
#include "larana/OpticalDetector/OpHitFinder/OpticalRecoTypes.h"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "TH1D.h"

namespace pmtana {

  double mean(const std::vector<short>& wf, size_t start, size_t nsample)
  {
    if (!nsample) nsample = wf.size();
    if (start > wf.size() || (start + nsample) > wf.size())
      throw OpticalRecoException("Invalid start/end index!");

    double sum =
      std::accumulate(wf.begin() + start, wf.begin() + start + nsample, 0.0) / ((double)nsample);

    return sum;
  }

  double edge_aware_mean(const std::vector<short>& wf, int start, int end)
  {

    auto m = double{0.0};
    auto n_t = unsigned{0};

    for (int k = start; k < end; ++k) {
      if (k < 0 or k > (int)(wf.size()) - 1) continue;
      m += wf.at(k);
      ++n_t;
    }

    if (n_t > 0) m /= n_t;
    n_t = 0;

    return m;
  }

  double std(const std::vector<short>& wf, const double ped_mean, size_t start, size_t nsample)
  {
    if (!nsample) nsample = wf.size();
    if (start > wf.size() || (start + nsample) > wf.size())
      throw OpticalRecoException("Invalid start/end index!");

    double sigma = 0;

    for (size_t index = start; index < (start + nsample); ++index)

      sigma += pow((wf[index] - ped_mean), 2);

    sigma = sqrt(sigma / ((double)(nsample)));

    return sigma;
  }

  double BinnedMaxOccurrence(const PedestalMean_t& mean_v, const size_t nbins)
  {
    if (nbins < 1) throw OpticalRecoException("Cannot have 0 binning");

    auto res = std::minmax_element(std::begin(mean_v), std::end(mean_v));

    double bin_width = ((*res.second) - (*res.first)) / ((double)nbins);

    if (nbins == 1 || bin_width == 0) return ((*res.first) + bin_width / 2.);

    //std::cout<<"Min: "<<(*res.first)<<" Max: "<<(*res.second)<<" Width: "<<bin_width<<std::endl;

    // Construct array of nbins
    static std::vector<size_t> ctr_v(nbins, 0);
    for (auto& v : ctr_v)
      v = 0;
    for (auto const& v : mean_v) {

      size_t index = int((v - (*res.first)) / bin_width);
      //std::cout<<"adc = "<<v<<" width = "<<bin_width<< " ... "
      //<<index<<" / "<<ctr_v.size()<<std::endl;

      ctr_v[index]++;
    }

    // Find max occurrence
    auto max_it = std::max_element(std::begin(ctr_v), std::end(ctr_v));

    // Get the mean of max-occurrence bins
    double mean_max_occurrence = 0;
    double num_occurrence = 0;
    for (size_t bin = 0; bin < ctr_v.size(); ++bin) {

      if (ctr_v[bin] != (*max_it)) continue;

      mean_max_occurrence += ((*res.first) + bin_width / 2. + bin_width * bin);

      num_occurrence += 1.0;
    }

    return (mean_max_occurrence / num_occurrence);
  }

  // template<typename W>
  int sign(double val)
  {

    if (val > 0) return 1;
    if (val < 0) return -1;
    return 0;
  }

  double BinnedMaxTH1D(const std::vector<double>& v, int bins)
  {

    auto max_it = std::max_element(std::begin(v), std::end(v));
    auto min_it = std::min_element(std::begin(v), std::end(v));

    TH1D th("th", ";;", bins, *min_it, *max_it);

    for (const auto& m : v)
      th.Fill(m);

    return th.GetXaxis()->GetBinCenter(th.GetMaximumBin());
  }

}

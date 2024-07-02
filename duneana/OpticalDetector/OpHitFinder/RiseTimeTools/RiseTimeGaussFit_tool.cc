/**
 * \file RiseTimeGaussFit_tool.cc
 *
 * \brief Gaussian fit method implemented to  compute Rise time as
 * the value of the center of the first local maximum.
 * Fixed min threshold required set in the fhicl file
 *
 * @author Rodrigo Alvarez, Sept 2022
 */

#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "RiseTimeCalculatorBase.h"

// ROOT includes
#include "TF1.h"
#include "TH1F.h"

namespace pmtana {

  class RiseTimeGaussFit : RiseTimeCalculatorBase {

  public:
    //Configuration parameters
    struct Config {

      fhicl::Atom<double> MinAmp{fhicl::Name("MinAmp")};
      fhicl::Atom<double> InitSigma{fhicl::Name("InitSigma")};
      fhicl::Atom<int> Nbins{fhicl::Name("Nbins")};
      fhicl::Atom<double> Tolerance{fhicl::Name("Tolerance")};
    };

    // Default constructor
    explicit RiseTimeGaussFit(art::ToolConfigTable<Config> const& config);

    // Method to calculate the OpFlash t0
    double RiseTime(const pmtana::Waveform_t& wf_pulse,
                    const pmtana::PedestalMean_t& ped_pulse,
                    bool _positive) const override;
    // Method to fit the first local max of the wvf above fixed threshold
    std::size_t findFirstMax(const std::vector<double>& arr, double threshold) const;

  private:
    double fMinAmp;
    double fInitSigma;
    double fTolerance;
    int fNbins;
  };

  RiseTimeGaussFit::RiseTimeGaussFit(art::ToolConfigTable<Config> const& config)
    : fMinAmp{config().MinAmp()}
    , fInitSigma{config().InitSigma()}
    , fTolerance{config().Tolerance()}
    , fNbins{config().Nbins()}
  {}

  double RiseTimeGaussFit::RiseTime(const pmtana::Waveform_t& wf_pulse,
                                    const pmtana::PedestalMean_t& ped_pulse,
                                    bool _positive) const
  {

    // Pedestal-subtracted pulse
    std::vector<double> wf_aux(ped_pulse);
    if (_positive) {
      for (size_t ix = 0; ix < wf_aux.size(); ix++) {
        wf_aux[ix] = ((double)wf_pulse[ix]) - wf_aux[ix];
      }
    }
    else {
      for (size_t ix = 0; ix < wf_aux.size(); ix++) {
        wf_aux[ix] = wf_aux[ix] - ((double)wf_pulse[ix]);
      }
    }

    // Find first local maximum
    size_t first_max = findFirstMax(wf_aux, fMinAmp);

    // Create & fill th1
    TH1F* h_aux = new TH1F("aux", "aux", wf_aux.size(), -0.5, wf_aux.size() - 0.5);
    for (long unsigned int j = 0; j < wf_aux.size(); j++)
      h_aux->SetBinContent(j + 1, wf_aux[j]);

    // Function & initial values for the fit, ROOT Gaussian:  [p0]*e**(-0.5 (x-[p1])**2 / [p2]**2)
    TF1* f = new TF1("f", "gaus", double(first_max) - fNbins, double(first_max) + fNbins);
    f->SetParameters(wf_aux[first_max], first_max, fInitSigma);
    // Fit
    h_aux->Fit(f, "q", "SAME", first_max - fNbins, first_max + fNbins);
    double t_fit = f->GetParameter(1);
    double peak_time;
    if (
      std::abs(t_fit - first_max) <
      fTolerance) { //check fit is close in distance to the original max, use max peak as time otherwise
      peak_time = t_fit;
    }
    else {
      peak_time = first_max;
      mf::LogInfo("RiseTimeGaussFit") << "No good fit found, keeping 1st max bin instead";
    }

    return peak_time;
  }

  std::size_t RiseTimeGaussFit::findFirstMax(const std::vector<double>& arr, double threshold) const
  { /**
    * Linear search, O(N), 1st peak should be close to the start of the vector
    * for scintillation LAr Signals. Returns the position of the first local max
    **/
    double max = arr[0];
    for (std::size_t i = 1, n = arr.size(); i < n; ++i) {
      if (arr[i] >= max)
        max = arr[i];

      else if (max < threshold) //didn't pass the threshold, keep searching
        max = arr[i];

      else
        return i - 1;
    }

    //No local maxima found.
    mf::LogInfo("RiseTimeGaussFit") << "No local max found above fixed threshold: " << threshold;
    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(pmtana::RiseTimeGaussFit)

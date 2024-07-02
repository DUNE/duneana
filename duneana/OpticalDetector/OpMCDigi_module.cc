// \file OpMCDigi.h
// \author Ben Jones and Christie Chiu, MIT, Sept 2012
//   bjpjones@mit.edu, cschiu@mit.edu
//
// This module starts from MC truth sim::OnePhoton objects
// and produces a digitized waveform.
//
// It is assumed that the electronics response is linear,
// and the 1PE waveform can be described by a discreate
// response shape.  The many PE response is then the linear
// superposition of the relevant function at the appropriate
// arrival times.
//

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes
#include "larana/OpticalDetector/OpDetResponseInterface.h"
#include "larana/OpticalDetector/OpDigiProperties.h"
#include "lardataobj/RawData/OpDetPulse.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larsim/Simulation/SimListUtils.h"

// CLHEP includes
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"

// nurandom
#include "nurandom/RandomUtils/NuRandomService.h"

// C++ language includes
#include <cstring>

namespace opdet {

  class OpMCDigi : public art::EDProducer {
  public:
    explicit OpMCDigi(const fhicl::ParameterSet&);

  private:
    void produce(art::Event&) override;

    // The parameters we'll read from the .fcl file.
    std::string fInputModule; // Input tag for OpDet collection
    float fSampleFreq;        // in MHz
    float fTimeBegin;         // in us
    float fTimeEnd;           // in us
    //float fQE;                             // quantum efficiency of opdet
    float fSaturationScale; // adc count w/ saturation occurs

    float fDarkRate; // Noise rate in Hz

    std::vector<double> fSinglePEWaveform;

    CLHEP::HepRandomEngine& fEngine;
    CLHEP::RandFlat fFlatRandom;
    CLHEP::RandPoisson fPoissonRandom;

    void AddTimedWaveform(int time, std::vector<double>& OldPulse, std::vector<double>& NewPulse);
  };
}

// Debug flag; only used during code development.
// const bool debug = true;

namespace opdet {

  OpMCDigi::OpMCDigi(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , fInputModule{pset.get<std::string>("InputModule")} //, fQE{pset.get<double>("QE")}
    , fSaturationScale{pset.get<float>("SaturationScale")}
    , fDarkRate{pset.get<float>("DarkRate")}
    // create a default random engine; obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed"
    , fEngine(art::ServiceHandle<rndm::NuRandomService>()->registerAndSeedEngine(createEngine(0),
                                                                                 pset,
                                                                                 "Seed"))
    , fFlatRandom{fEngine}
    , fPoissonRandom{fEngine}
  {
    produces<std::vector<raw::OpDetPulse>>();

    art::ServiceHandle<OpDigiProperties> odp;
    fSampleFreq = odp->SampleFreq();
    fTimeBegin = odp->TimeBegin();
    fTimeEnd = odp->TimeEnd();
    fSinglePEWaveform = odp->SinglePEWaveform();
  }

  //-------------------------------------------------

  void OpMCDigi::AddTimedWaveform(int binTime,
                                  std::vector<double>& OldPulse,
                                  std::vector<double>& NewPulse)
  {

    if ((binTime + NewPulse.size()) > OldPulse.size()) {
      OldPulse.resize(binTime + NewPulse.size());
    }

    // Add shifted NewWaveform to Waveform at pointer
    for (size_t i = 0; i != NewPulse.size(); ++i) {
      OldPulse.at(binTime + i) += NewPulse.at(i);
    }
  }

  //-------------------------------------------------

  void OpMCDigi::produce(art::Event& evt)
  {
    auto StoragePtr = std::make_unique<std::vector<raw::OpDetPulse>>();

    bool const fUseLitePhotons = art::ServiceHandle<sim::LArG4Parameters const>
    {
      } -> UseLitePhotons();

    // Service for determining opdet responses
    art::ServiceHandle<opdet::OpDetResponseInterface const> odresponse;

    double const TimeBegin_ns = fTimeBegin * 1000;
    double const TimeEnd_ns = fTimeEnd * 1000;
    double const SampleFreq_ns = fSampleFreq / 1000;

    int const nSamples = (TimeEnd_ns - TimeBegin_ns) * SampleFreq_ns;
    int const NOpChannels = odresponse->NOpChannels();

    // This vector will store all the waveforms we will make
    std::vector<std::vector<double>> PulsesFromDetPhotons(NOpChannels,
                                                          std::vector<double>(nSamples, 0.0));

    if (!fUseLitePhotons) {
      // Read in the Sim Photons
      sim::SimPhotonsCollection ThePhotCollection =
        sim::SimListUtils::GetSimPhotonsCollection(evt, fInputModule);
      // For every OpDet:
      for (auto const& pr : ThePhotCollection) {
        const sim::SimPhotons& ThePhot = pr.second;

        int const Ch = ThePhot.OpChannel();
        int readoutCh;

        // For every photon in the hit:
        for (const sim::OnePhoton& Phot : ThePhot) {
          // Sample a random subset according to QE
          if (!odresponse->detected(Ch, Phot, readoutCh)) { continue; }

          // Convert photon arrival time to the appropriate bin,
          // dictated by fSampleFreq. Photon arrival time is in ns,
          // beginning time in us, and sample frequency in MHz. Notice
          // that we have to accommodate for the beginning time
          if ((Phot.Time > TimeBegin_ns) && (Phot.Time < TimeEnd_ns)) {
            auto const binTime = static_cast<int>((Phot.Time - TimeBegin_ns) * SampleFreq_ns);
            AddTimedWaveform(binTime, PulsesFromDetPhotons[readoutCh], fSinglePEWaveform);
          }
        } // for each Photon in SimPhotons
      }
    }
    else {
      auto const photons = *evt.getValidHandle<std::vector<sim::SimPhotonsLite>>("largeant");
      // For every OpDet:
      for (auto const& photon : photons) {
        int const Ch = photon.OpChannel;
        int readoutCh;

        std::map<int, int> PhotonsMap = photon.DetectedPhotons;

        // For every photon in the hit:
        for (auto const& pr : photon.DetectedPhotons) {
          for (int i = 0; i < pr.second; i++) {
            // Sample a random subset according to QE
            if (odresponse->detectedLite(Ch, readoutCh)) {
              // Convert photon arrival time to the appropriate bin, dictated by fSampleFreq.
              // Photon arrival time is in ns, beginning time in us, and sample frequency in MHz.
              // Notice that we have to accommodate for the beginning time
              if ((pr.first > TimeBegin_ns) && (pr.first < TimeEnd_ns)) {
                auto const binTime = static_cast<int>((pr.first - TimeBegin_ns) * SampleFreq_ns);
                AddTimedWaveform(binTime, PulsesFromDetPhotons[readoutCh], fSinglePEWaveform);
              }
            } // random QE cut
          }
        } // for each Photon in SimPhotons
      }
    }

    // Create vector of output objects, add dark noise and apply
    //  saturation

    std::vector<raw::OpDetPulse*> ThePulses(NOpChannels);
    for (int iCh = 0; iCh != NOpChannels; ++iCh) {
      PulsesFromDetPhotons[iCh].resize((TimeEnd_ns - TimeBegin_ns) * SampleFreq_ns);

      // Add dark noise
      double const MeanDarkPulses = fDarkRate * (fTimeEnd - fTimeBegin) / 1000000;
      unsigned const int NumberOfPulses = fPoissonRandom.fire(MeanDarkPulses);

      for (size_t i = 0; i != NumberOfPulses; ++i) {
        double const PulseTime = (fTimeEnd - fTimeBegin) * fFlatRandom.fire(1.0);
        int const binTime = static_cast<int>(PulseTime * fSampleFreq);

        AddTimedWaveform(binTime, PulsesFromDetPhotons[iCh], fSinglePEWaveform);
      }

      // Apply saturation for large signals
      for (size_t i = 0; i != PulsesFromDetPhotons[iCh].size(); ++i) {
        if (PulsesFromDetPhotons[iCh].at(i) > fSaturationScale)
          PulsesFromDetPhotons[iCh].at(i) = fSaturationScale;
      }

      // Produce ADC pulse of integers rather than doubles

      std::vector<short> shortvec;

      for (size_t i = 0; i != PulsesFromDetPhotons[iCh].size(); ++i) {
        // Throw randoms to fairly sample +ve and -ve side of doubles
        int ThisSample = PulsesFromDetPhotons[iCh].at(i);
        if (ThisSample > 0) {
          if (fFlatRandom.fire(1.0) > (ThisSample - int(ThisSample)))
            shortvec.push_back(int(ThisSample));
          else
            shortvec.push_back(int(ThisSample) + 1);
        }
        else {
          if (fFlatRandom.fire(1.0) > (int(ThisSample) - ThisSample))
            shortvec.push_back(int(ThisSample));
          else
            shortvec.push_back(int(ThisSample) - 1);
        }
      }

      StoragePtr->emplace_back(iCh, shortvec, 0, fTimeBegin);

    } // for each OpDet in SimPhotonsCollection

    evt.put(std::move(StoragePtr));
  }
}

DEFINE_ART_MODULE(opdet::OpMCDigi)

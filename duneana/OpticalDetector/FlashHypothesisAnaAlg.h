#ifndef FLASHHYPOTHESISANAALG_H
#define FLASHHYPOTHESISANAALG_H

/*!
 * Title:   FlashHypothesisAnaAlg Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description:
 * Alg that compares the flash hypotheses with truth photons and stores the
 * results in a TTree.
 *
*/

namespace fhicl {
  class ParameterSet;
}

namespace phot {
  class PhotonVisibilityServce;
}

namespace opdet {
  class OpDigiProperties;
}

#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "larsim/MCSTReco/MCTrackCollectionAnaAlg.h"

#include "FlashHypothesisComparison.h"
#include "FlashHypothesisCreator.h"
#include "SimPhotonCounterAlg.h"

class TH1F;
class TTree;

namespace geo {
  class Geometry;
}

namespace opdet {

  class FlashHypothesisAnaAlg {

  public:
    using Providers_t = FlashHypothesisCreator::Providers_t;

    FlashHypothesisAnaAlg(fhicl::ParameterSet const& p)
      : fCounterIndex(p.get<unsigned int>("SimPhotonCounterIndex", 0))
      , fdEdx(p.get<float>("dEdx", 2.1))
      , fXOffset(p.get<float>("HypothesisXOffset", 0.0))
      , fSPCAlg(p.get<fhicl::ParameterSet>("SimPhotonCounterAlgParams"))
    {}

    void SetOutputObjects(TTree*,
                          TH1F*,
                          TH1F*,
                          TH1F*,
                          TH1F*,
                          TH1F*,
                          TH1F*,
                          TH1F*,
                          TH1F*,
                          TH1F*,
                          geo::Geometry const&);

    void FillOpDetPositions(geo::Geometry const&);

    void RunComparison(const unsigned int run,
                       const unsigned int event,
                       std::vector<sim::MCTrack> const&,
                       std::vector<sim::SimPhotons> const&,
                       Providers_t providers,
                       opdet::OpDigiProperties const& opdigip,
                       phot::PhotonVisibilityService const& pvs);

  private:
    unsigned int fCounterIndex;
    float fdEdx;
    float fXOffset;

    TTree* fTree;

    FlashHypothesisCreator fFHCreator;
    SimPhotonCounterAlg fSPCAlg;
    FlashHypothesisComparison fFHCompare;
    sim::MCTrackCollectionAnaAlg fMCTAlg;

    std::vector<float> fOpDetPositions_Y;
    std::vector<float> fOpDetPositions_Z;
  };

}

#endif

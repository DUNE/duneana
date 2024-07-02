#ifndef SIMPHOTONCOUNTERALG_H
#define SIMPHOTONCOUNTERALG_H

/*!
 * Title:   SimPhotonCounterALG Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: Alg class that counts up sim photons, leading towards
 *              comparisons with flashes and flash hypotheses.
*/

#include "SimPhotonCounter.h"

#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/Simulation/SimPhotons.h"

namespace geo {
  class GeometryCore;
}
namespace opdet {
  class OpDigiProperties;
}

#include <vector>

namespace opdet {

  class SimPhotonCounterAlg {

  public:
    SimPhotonCounterAlg(fhicl::ParameterSet const&);

    void InitializeCounters(geo::GeometryCore const&, opdet::OpDigiProperties const&);

    void AddSimPhotonCollection(sim::SimPhotonsCollection const&);
    void AddSimPhotonsVector(std::vector<sim::SimPhotons> const&);

    void ClearCounters();

    SimPhotonCounter const& GetSimPhotonCounter(size_t);
    std::vector<float> const& PromptPhotonVector(size_t);
    std::vector<float> const& LatePhotonVector(size_t);

  private:
    std::vector<std::vector<float>> fWavelengthRanges;
    std::vector<std::vector<float>> fTimeRanges;
    std::vector<SimPhotonCounter> fCounters;

    void FillAllRanges(std::vector<fhicl::ParameterSet> const&);
    void FillRanges(fhicl::ParameterSet const&);
  };

}

#endif

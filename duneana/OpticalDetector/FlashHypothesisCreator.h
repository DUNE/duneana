#ifndef FLASHHYPOTHESISCREATOR_H
#define FLASHHYPOTHESISCREATOR_H

/*!
 * Title:   FlashHypothesis Creator Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: class that produces a flash hypothesis for a trajectory.
 * Input:       Trajectory (std::vector<TVector3> objects)
 * Output:      FlashHypotheses
*/

#include <vector>

namespace detinfo {
  class LArProperties;
}
namespace recob {
  class Track;
}
namespace sim {
  class MCTrack;
}
namespace geo {
  class GeometryCore;
}
namespace opdet {
  class OpDigiProperties;
}
namespace phot {
  class PhotonVisibilityService;
}

#include "larcorealg/CoreUtils/ProviderPack.h"

#include "TVector3.h"

#include "FlashHypothesis.h"
#include "FlashHypothesisCalculator.h"

namespace opdet {

  class FlashHypothesisCreator {

  public:
    /// Set of service providers used in the common(est) interface
    using Providers_t = lar::ProviderPack<geo::GeometryCore, detinfo::LArProperties>;

    FlashHypothesisCreator() {}

    FlashHypothesisCollection GetFlashHypothesisCollection(recob::Track const& track,
                                                           std::vector<float> const& dEdxVector,
                                                           Providers_t providers,
                                                           phot::PhotonVisibilityService const& pvs,
                                                           opdet::OpDigiProperties const& opdigip,
                                                           float XOffset = 0);

    FlashHypothesisCollection GetFlashHypothesisCollection(sim::MCTrack const& mctrack,
                                                           std::vector<float> const& dEdxVector,
                                                           Providers_t providers,
                                                           phot::PhotonVisibilityService const& pvs,
                                                           opdet::OpDigiProperties const& opdigip,
                                                           float XOffset = 0);

    FlashHypothesisCollection GetFlashHypothesisCollection(std::vector<TVector3> const& trajVector,
                                                           std::vector<float> const& dEdxVector,
                                                           Providers_t providers,
                                                           phot::PhotonVisibilityService const& pvs,
                                                           opdet::OpDigiProperties const& opdigip,
                                                           float XOffset = 0);

    FlashHypothesisCollection GetFlashHypothesisCollection(TVector3 const& pt1,
                                                           TVector3 const& pt2,
                                                           float const& dEdx,
                                                           Providers_t providers,
                                                           phot::PhotonVisibilityService const& pvs,
                                                           opdet::OpDigiProperties const& opdigip,
                                                           float XOffset = 0);

  private:
    FlashHypothesisCollection CreateFlashHypothesesFromSegment(
      TVector3 const& pt1,
      TVector3 const& pt2,
      float const& dEdx,
      Providers_t providers,
      phot::PhotonVisibilityService const& pvs,
      opdet::OpDigiProperties const& opdigip,
      float XOffset);

    FlashHypothesisCalculator _calc;
  };

}

#endif

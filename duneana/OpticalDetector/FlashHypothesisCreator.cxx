/*!
 * Title:   FlashHypothesis Creator Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: Class that produces a flash hypothesis for a trajectory.
 * Input:       Trajectory (std::vector<TVector3> objects)
 * Output:      FlashHypotheses
*/

#include "FlashHypothesisCreator.h"
#include "OpDigiProperties.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataalg/DetectorInfo/LArProperties.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/Track.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"

opdet::FlashHypothesisCollection opdet::FlashHypothesisCreator::GetFlashHypothesisCollection(
  recob::Track const& track,
  std::vector<float> const& dEdxVector,
  Providers_t providers,
  phot::PhotonVisibilityService const& pvs,
  opdet::OpDigiProperties const& opdigip,
  float XOffset)
{
  bool interpolate_dEdx = false;
  if (track.NumberTrajectoryPoints() == dEdxVector.size())
    interpolate_dEdx = true;
  else if (track.NumberTrajectoryPoints() == dEdxVector.size() + 1)
    interpolate_dEdx = false;
  else
    throw "ERROR in FlashHypothesisCreator: dEdx vector size not compatible with track size.";

  auto const* geom = providers.get<geo::GeometryCore>();
  FlashHypothesisCollection fhc(geom->NOpDets());
  for (size_t pt = 1; pt < track.NumberTrajectoryPoints(); pt++) {
    if (interpolate_dEdx)
      fhc = fhc + CreateFlashHypothesesFromSegment(track.LocationAtPoint<TVector3>(pt - 1),
                                                   track.LocationAtPoint<TVector3>(pt),
                                                   0.5 * (dEdxVector[pt] + dEdxVector[pt - 1]),
                                                   providers,
                                                   pvs,
                                                   opdigip,
                                                   XOffset);
    else
      fhc = fhc + CreateFlashHypothesesFromSegment(track.LocationAtPoint<TVector3>(pt - 1),
                                                   track.LocationAtPoint<TVector3>(pt),
                                                   dEdxVector[pt - 1],
                                                   providers,
                                                   pvs,
                                                   opdigip,
                                                   XOffset);
  }
  return fhc;
}

opdet::FlashHypothesisCollection opdet::FlashHypothesisCreator::GetFlashHypothesisCollection(
  sim::MCTrack const& mctrack,
  std::vector<float> const& dEdxVector,
  Providers_t providers,
  phot::PhotonVisibilityService const& pvs,
  opdet::OpDigiProperties const& opdigip,
  float XOffset)
{
  bool interpolate_dEdx = false;
  if (mctrack.size() == dEdxVector.size())
    interpolate_dEdx = true;
  else if (mctrack.size() == dEdxVector.size() + 1)
    interpolate_dEdx = false;
  else
    throw "ERROR in FlashHypothesisCreator: dEdx vector size not compatible with mctrack size.";

  auto const* geom = providers.get<geo::GeometryCore>();
  FlashHypothesisCollection fhc(geom->NOpDets());
  for (size_t pt = 1; pt < mctrack.size(); pt++) {
    if (interpolate_dEdx)
      fhc = fhc + CreateFlashHypothesesFromSegment(mctrack[pt - 1].Position().Vect(),
                                                   mctrack[pt].Position().Vect(),
                                                   0.5 * (dEdxVector[pt] + dEdxVector[pt - 1]),
                                                   providers,
                                                   pvs,
                                                   opdigip,
                                                   XOffset);
    else
      fhc = fhc + CreateFlashHypothesesFromSegment(mctrack[pt - 1].Position().Vect(),
                                                   mctrack[pt].Position().Vect(),
                                                   dEdxVector[pt - 1],
                                                   providers,
                                                   pvs,
                                                   opdigip,
                                                   XOffset);
  }
  return fhc;
}

opdet::FlashHypothesisCollection opdet::FlashHypothesisCreator::GetFlashHypothesisCollection(
  std::vector<TVector3> const& trajVector,
  std::vector<float> const& dEdxVector,
  Providers_t providers,
  phot::PhotonVisibilityService const& pvs,
  opdet::OpDigiProperties const& opdigip,
  float XOffset)
{
  bool interpolate_dEdx = false;
  if (trajVector.size() == dEdxVector.size())
    interpolate_dEdx = true;
  else if (trajVector.size() == dEdxVector.size() + 1)
    interpolate_dEdx = false;
  else
    throw "ERROR in FlashHypothesisCreator: dEdx vector size not compatible with trajVector size.";

  auto const* geom = providers.get<geo::GeometryCore>();
  FlashHypothesisCollection fhc(geom->NOpDets());
  for (size_t pt = 1; pt < trajVector.size(); pt++) {
    if (interpolate_dEdx)
      fhc = fhc + CreateFlashHypothesesFromSegment(trajVector[pt - 1],
                                                   trajVector[pt],
                                                   0.5 * (dEdxVector[pt] + dEdxVector[pt - 1]),
                                                   providers,
                                                   pvs,
                                                   opdigip,
                                                   XOffset);
    else
      fhc =
        fhc +
        CreateFlashHypothesesFromSegment(
          trajVector[pt - 1], trajVector[pt], dEdxVector[pt - 1], providers, pvs, opdigip, XOffset);
  }
  return fhc;
}

opdet::FlashHypothesisCollection opdet::FlashHypothesisCreator::GetFlashHypothesisCollection(
  TVector3 const& pt1,
  TVector3 const& pt2,
  float const& dEdx,
  Providers_t providers,
  phot::PhotonVisibilityService const& pvs,
  opdet::OpDigiProperties const& opdigip,
  float XOffset)
{
  return CreateFlashHypothesesFromSegment(pt1, pt2, dEdx, providers, pvs, opdigip, XOffset);
}

opdet::FlashHypothesisCollection opdet::FlashHypothesisCreator::CreateFlashHypothesesFromSegment(
  TVector3 const& pt1,
  TVector3 const& pt2,
  float const& dEdx,
  Providers_t providers,
  phot::PhotonVisibilityService const& pvs,
  opdet::OpDigiProperties const& opdigip,
  float XOffset)
{
  auto const* geom = providers.get<geo::GeometryCore>();
  auto const* larp = providers.get<detinfo::LArProperties>();
  auto const nOpDets = geom->NOpDets();
  FlashHypothesisCollection fhc(nOpDets);

  FlashHypothesis prompt_hyp = FlashHypothesis(nOpDets);

  std::vector<double> xyz_segment(_calc.SegmentMidpoint(pt1, pt2, XOffset));

  //get the visibility vector
  auto const& PointVisibility = pvs.GetAllVisibilities(&xyz_segment[0]);

  //check visibility pointer, as it may be null if given a y/z outside some range
  if (!PointVisibility) return fhc;

  //klugey ... right now, set a qe_vector that gives constant qe across all opdets
  std::vector<float> qe_vector(nOpDets, opdigip.QE());
  _calc.FillFlashHypothesis(larp->ScintYield() * larp->ScintYieldRatio(),
                            dEdx,
                            pt1,
                            pt2,
                            qe_vector,
                            PointVisibility,
                            prompt_hyp);

  fhc.SetPromptHypAndPromptFraction(prompt_hyp, larp->ScintYieldRatio());
  return fhc;
}

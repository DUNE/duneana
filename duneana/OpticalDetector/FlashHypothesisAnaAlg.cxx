/*!
 * Title:   FlashHypothesisAnaAlg Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description:
 * Alg that compares the flash hypotheses with truth photons and stores the
 * results in a TTree.
 *
*/

#include "FlashHypothesisAnaAlg.h"
#include "FlashHypothesis.h"
#include "FlashHypothesisComparison.h"
#include "FlashHypothesisCreator.h"
#include "SimPhotonCounterAlg.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "lardataalg/DetectorInfo/LArProperties.h"
#include "lardataobj/MCBase/MCTrack.h"

#include "TTree.h"

void opdet::FlashHypothesisAnaAlg::SetOutputObjects(TTree* tree,
                                                    TH1F* h_h_p,
                                                    TH1F* h_s_p,
                                                    TH1F* h_c_p,
                                                    TH1F* h_h_l,
                                                    TH1F* h_s_l,
                                                    TH1F* h_c_l,
                                                    TH1F* h_h_t,
                                                    TH1F* h_s_t,
                                                    TH1F* h_c_t,
                                                    geo::Geometry const& geom)
{
  fTree = tree;
  fFHCompare.SetOutputObjects(
    tree, h_h_p, h_s_p, h_c_p, h_h_l, h_s_l, h_c_l, h_h_t, h_s_t, h_c_t, geom.NOpDets(), false);

  fMCTAlg.SetOutputTree(tree, false);
}

void opdet::FlashHypothesisAnaAlg::FillOpDetPositions(geo::Geometry const& geom)
{

  fOpDetPositions_Y.resize(geom.NOpDets());
  fOpDetPositions_Z.resize(geom.NOpDets());

  for (size_t i_opdet = 0; i_opdet < geom.NOpDets(); i_opdet++) {
    auto const xyz = geom.Cryostat().OpDet(i_opdet).GetCenter();
    fOpDetPositions_Y[i_opdet] = (float)xyz.Y();
    fOpDetPositions_Z[i_opdet] = (float)xyz.Z();
  }
}

void opdet::FlashHypothesisAnaAlg::RunComparison(const unsigned int run,
                                                 const unsigned int event,
                                                 std::vector<sim::MCTrack> const& mctrackVec,
                                                 std::vector<sim::SimPhotons> const& simPhotonsVec,
                                                 Providers_t providers,
                                                 opdet::OpDigiProperties const& opdigip,
                                                 phot::PhotonVisibilityService const& pvs)
{
  auto const* geom = providers.get<geo::GeometryCore>();

  FlashHypothesisCollection fhc(geom->NOpDets());
  for (auto const& mctrack : mctrackVec) {
    if (mctrack.size() == 0) continue;
    std::vector<float> dEdxVector(mctrack.size() - 1, fdEdx);
    fhc = fhc + fFHCreator.GetFlashHypothesisCollection(
                  mctrack, dEdxVector, providers, pvs, opdigip, fXOffset);
  }

  fSPCAlg.InitializeCounters(*geom, opdigip);
  fSPCAlg.AddSimPhotonsVector(simPhotonsVec);

  fFHCompare.RunComparison(run,
                           event,
                           fhc,
                           fSPCAlg.GetSimPhotonCounter(fCounterIndex),
                           fOpDetPositions_Y,
                           fOpDetPositions_Z);

  fMCTAlg.FillTree(run, event, mctrackVec);

  fTree->Fill();
}

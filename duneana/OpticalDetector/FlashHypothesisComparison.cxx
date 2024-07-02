/*!
 * Title:   FlashHypothesisComparison Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description:
 * Class for comparing a flash hypothesis to MC truth (via SimPhotonCounter).
 * Needs a flash hypothesis and a SimPhotonCounter object as input.
 * Outputs a Tree with relevent info.
 */

#include "FlashHypothesisComparison.h"

#include "larana/OpticalDetector/FlashHypothesis.h"
#include "larana/OpticalDetector/FlashUtilities.h"
#include "larana/OpticalDetector/SimPhotonCounter.h"

#include "TH1F.h"
#include "TTree.h"

#include <iostream>

void opdet::FlashHypothesisComparison::SetOutputObjects(TTree* tree,
                                                        TH1F* h_h_p,
                                                        TH1F* h_s_p,
                                                        TH1F* h_c_p,
                                                        TH1F* h_h_l,
                                                        TH1F* h_s_l,
                                                        TH1F* h_c_l,
                                                        TH1F* h_h_t,
                                                        TH1F* h_s_t,
                                                        TH1F* h_c_t,
                                                        const unsigned int n_opdet,
                                                        bool fill)
{
  fTree = tree;
  fFillTree = fill;
  fTree->SetName("ctree");

  fTree->Branch("run", &fRun, "run/i");
  fTree->Branch("event", &fEvent, "event/i");

  fHypHist_p = h_h_p;
  fSimHist_p = h_s_p;
  fCompareHist_p = h_c_p;

  fHypHist_p->SetBins(n_opdet, -0.5, (float)n_opdet - 0.5);
  fSimHist_p->SetBins(n_opdet, -0.5, (float)n_opdet - 0.5);
  fCompareHist_p->SetBins(n_opdet, -0.5, (float)n_opdet - 0.5);

  fHypHist_p->SetNameTitle("hHypHist_p", "Hypothesis (Prompt);Opdet;PEs");
  fSimHist_p->SetNameTitle("hSimHist_p", "SimPhoton (Prompt);Opdet;PEs");
  fCompareHist_p->SetNameTitle("hCompareHist_p", "Comparison (Hyp - Sim) (Prompt);Opdet;PEs");

  fTree->Branch("hyp_PEs_p", &fHypPEs_p, "hyp_PEs_p/F");
  fTree->Branch("hyp_PEsError_p", &fHypPEsError_p, "hyp_PEsError_p/F");
  fTree->Branch("sim_PEs_p", &fSimPEs_p, "sim_PEs_p/F");

  fTree->Branch("hyp_Y_p", &fHypY_p, "hyp_Y_p/F");
  fTree->Branch("sim_Y_p", &fSimY_p, "sim_Y_p/F");

  fTree->Branch("hyp_RMSY_p", &fHypRMSY_p, "hyp_RMSY_p/F");
  fTree->Branch("sim_RMSY_p", &fSimRMSY_p, "sim_RMSY_p/F");

  fTree->Branch("hyp_Z_p", &fHypZ_p, "hyp_Z_p/F");
  fTree->Branch("sim_Z_p", &fSimZ_p, "sim_Z_p/F");

  fTree->Branch("hyp_RMSZ_p", &fHypRMSZ_p, "hyp_RMSZ_p/F");
  fTree->Branch("sim_RMSZ_p", &fSimRMSZ_p, "sim_RMSZ_p/F");

  fTree->Branch("comp_total_p", &fCompare_p, "comp_total_p/F");

  fTree->Branch("hHypHist_p", &fHypHist_p);
  fTree->Branch("hSimHist_p", &fSimHist_p);
  fTree->Branch("hCompareHist_p", &fCompareHist_p);

  fHypHist_l = h_h_l;
  fSimHist_l = h_s_l;
  fCompareHist_l = h_c_l;

  fHypHist_l->SetBins(n_opdet, -0.5, (float)n_opdet - 0.5);
  fSimHist_l->SetBins(n_opdet, -0.5, (float)n_opdet - 0.5);
  fCompareHist_l->SetBins(n_opdet, -0.5, (float)n_opdet - 0.5);

  fHypHist_l->SetNameTitle("hHypHist_l", "Hypothesis (Late);Opdet;PEs");
  fSimHist_l->SetNameTitle("hSimHist_l", "SimPhoton (Late);Opdet;PEs");
  fCompareHist_l->SetNameTitle("hCompareHist_l", "Comparison (Hyp - Sim) (Late);Opdet;PEs");

  fTree->Branch("hyp_PEs_l", &fHypPEs_l, "hyp_PEs_l/F");
  fTree->Branch("hyp_PEsError_l", &fHypPEsError_l, "hyp_PEsError_l/F");
  fTree->Branch("sim_PEs_l", &fSimPEs_l, "sim_PEs_l/F");

  fTree->Branch("hyp_Y_l", &fHypY_l, "hyp_Y_l/F");
  fTree->Branch("sim_Y_l", &fSimY_l, "sim_Y_l/F");

  fTree->Branch("hyp_RMSY_l", &fHypRMSY_l, "hyp_RMSY_l/F");
  fTree->Branch("sim_RMSY_l", &fSimRMSY_l, "sim_RMSY_l/F");

  fTree->Branch("hyp_Z_l", &fHypZ_l, "hyp_Z_l/F");
  fTree->Branch("sim_Z_l", &fSimZ_l, "sim_Z_l/F");

  fTree->Branch("hyp_RMSZ_l", &fHypRMSZ_l, "hyp_RMSZ_l/F");
  fTree->Branch("sim_RMSZ_l", &fSimRMSZ_l, "sim_RMSZ_l/F");

  fTree->Branch("comp_total_l", &fCompare_l, "comp_total_l/F");

  fTree->Branch("hHypHist_l", &fHypHist_l);
  fTree->Branch("hSimHist_l", &fSimHist_l);
  fTree->Branch("hCompareHist_l", &fCompareHist_l);

  fHypHist_t = h_h_t;
  fSimHist_t = h_s_t;
  fCompareHist_t = h_c_t;

  fHypHist_t->SetBins(n_opdet, -0.5, (float)n_opdet - 0.5);
  fSimHist_t->SetBins(n_opdet, -0.5, (float)n_opdet - 0.5);
  fCompareHist_t->SetBins(n_opdet, -0.5, (float)n_opdet - 0.5);

  fHypHist_t->SetNameTitle("hHypHist_t", "Hypothesis (Total);Opdet;PEs");
  fSimHist_t->SetNameTitle("hSimHist_t", "SimPhoton (Total);Opdet;PEs");
  fCompareHist_t->SetNameTitle("hCompareHist_t", "Comparison (Hyp - Sim) (Total);Opdet;PEs");

  fTree->Branch("hyp_PEs_t", &fHypPEs_t, "hyp_PEs_t/F");
  fTree->Branch("hyp_PEsError_t", &fHypPEsError_t, "hyp_PEsError_t/F");
  fTree->Branch("sim_PEs_t", &fSimPEs_t, "sim_PEs_t/F");

  fTree->Branch("hyp_Y_t", &fHypY_t, "hyp_Y_t/F");
  fTree->Branch("sim_Y_t", &fSimY_t, "sim_Y_t/F");

  fTree->Branch("hyp_RMSY_t", &fHypRMSY_t, "hyp_RMSY_t/F");
  fTree->Branch("sim_RMSY_t", &fSimRMSY_t, "sim_RMSY_t/F");

  fTree->Branch("hyp_Z_t", &fHypZ_t, "hyp_Z_t/F");
  fTree->Branch("sim_Z_t", &fSimZ_t, "sim_Z_t/F");

  fTree->Branch("hyp_RMSZ_t", &fHypRMSZ_t, "hyp_RMSZ_t/F");
  fTree->Branch("sim_RMSZ_t", &fSimRMSZ_t, "sim_RMSZ_t/F");

  fTree->Branch("comp_total_t", &fCompare_t, "comp_total_t/F");

  fTree->Branch("hHypHist_t", &fHypHist_t);
  fTree->Branch("hSimHist_t", &fSimHist_t);
  fTree->Branch("hCompareHist_t", &fCompareHist_t);
}

void opdet::FlashHypothesisComparison::RunComparison(const unsigned int run,
                                                     const unsigned int event,
                                                     const FlashHypothesisCollection& fhc,
                                                     const SimPhotonCounter& spc,
                                                     const std::vector<float>& posY,
                                                     const std::vector<float>& posZ)
{
  if (fhc.GetVectorSize() != (unsigned int)fHypHist_p->GetNbinsX() ||
      fhc.GetVectorSize() != spc.PromptPhotonVector().size() ||
      fhc.GetVectorSize() != posY.size() || fhc.GetVectorSize() != posZ.size()) {
    std::cout << (unsigned int)fHypHist_p->GetNbinsX() << " " << spc.PromptPhotonVector().size()
              << " " << posY.size() << " " << posZ.size() << std::endl;
    throw std::runtime_error("ERROR in FlashHypothesisComparison: Mismatch in vector sizes.");
  }
  fRun = run;
  fEvent = event;

  FillFlashHypothesisInfo(fhc, posY, posZ);
  FillSimPhotonCounterInfo(spc, posY, posZ);
  FillComparisonInfo(fhc, spc);

  if (fFillTree) fTree->Fill();
}

void opdet::FlashHypothesisComparison::FillFlashHypothesisInfo(const FlashHypothesisCollection& fhc,
                                                               const std::vector<float>& posY,
                                                               const std::vector<float>& posZ)
{
  fHypPEs_p = fhc.GetPromptHypothesis().GetTotalPEs();
  fHypPEsError_p = fhc.GetPromptHypothesis().GetTotalPEsError();
  fUtil.GetPosition(fhc.GetPromptHypothesis().GetHypothesisVector(), posY, fHypY_p, fHypRMSY_p);
  fUtil.GetPosition(fhc.GetPromptHypothesis().GetHypothesisVector(), posZ, fHypZ_p, fHypRMSZ_p);

  for (size_t i = 0; i < fhc.GetVectorSize(); i++)
    fHypHist_p->SetBinContent(i + 1, fhc.GetPromptHypothesis().GetHypothesis(i));

  fHypPEs_l = fhc.GetLateHypothesis().GetTotalPEs();
  fHypPEsError_l = fhc.GetLateHypothesis().GetTotalPEsError();
  fUtil.GetPosition(fhc.GetLateHypothesis().GetHypothesisVector(), posY, fHypY_l, fHypRMSY_l);
  fUtil.GetPosition(fhc.GetLateHypothesis().GetHypothesisVector(), posZ, fHypZ_l, fHypRMSZ_l);

  for (size_t i = 0; i < fhc.GetVectorSize(); i++)
    fHypHist_l->SetBinContent(i + 1, fhc.GetLateHypothesis().GetHypothesis(i));

  fHypPEs_t = fhc.GetTotalHypothesis().GetTotalPEs();
  fHypPEsError_t = fhc.GetTotalHypothesis().GetTotalPEsError();
  fUtil.GetPosition(fhc.GetTotalHypothesis().GetHypothesisVector(), posY, fHypY_t, fHypRMSY_t);
  fUtil.GetPosition(fhc.GetTotalHypothesis().GetHypothesisVector(), posZ, fHypZ_t, fHypRMSZ_t);

  for (size_t i = 0; i < fhc.GetVectorSize(); i++)
    fHypHist_t->SetBinContent(i + 1, fhc.GetLateHypothesis().GetHypothesis(i));
}

void opdet::FlashHypothesisComparison::FillSimPhotonCounterInfo(const SimPhotonCounter& spc,
                                                                const std::vector<float>& posY,
                                                                const std::vector<float>& posZ)
{
  fSimPEs_p = spc.PromptPhotonTotal();
  fUtil.GetPosition(spc.PromptPhotonVector(), posY, fSimY_p, fSimRMSY_p);
  fUtil.GetPosition(spc.PromptPhotonVector(), posZ, fSimZ_p, fSimRMSZ_p);

  for (size_t i = 0; i < spc.PromptPhotonVector().size(); i++)
    fSimHist_p->SetBinContent(i + 1, spc.PromptPhotonVector()[i]);

  fSimPEs_l = spc.LatePhotonTotal();
  fUtil.GetPosition(spc.LatePhotonVector(), posY, fSimY_l, fSimRMSY_l);
  fUtil.GetPosition(spc.LatePhotonVector(), posZ, fSimZ_l, fSimRMSZ_l);

  for (size_t i = 0; i < spc.LatePhotonVector().size(); i++)
    fSimHist_l->SetBinContent(i + 1, spc.LatePhotonVector()[i]);

  fSimPEs_t = fSimPEs_p + fSimPEs_l;
  fUtil.GetPosition(spc.TotalPhotonVector(), posY, fSimY_t, fSimRMSY_t);
  fUtil.GetPosition(spc.TotalPhotonVector(), posZ, fSimZ_t, fSimRMSZ_t);

  for (size_t i = 0; i < spc.GetVectorSize(); i++)
    fSimHist_t->SetBinContent(i + 1, spc.TotalPhotonVector(i));
}

void opdet::FlashHypothesisComparison::FillComparisonInfo(const FlashHypothesisCollection& fhc,
                                                          const SimPhotonCounter& spc)
{
  std::vector<float> result_p, result_l, result_t;
  fCompare_p = fUtil.CompareByError(fhc.GetPromptHypothesis(), spc.PromptPhotonVector(), result_p);
  fCompare_l = fUtil.CompareByError(fhc.GetLateHypothesis(), spc.LatePhotonVector(), result_l);
  fCompare_t = fUtil.CompareByError(fhc.GetTotalHypothesis(), spc.TotalPhotonVector(), result_t);

  for (size_t i = 0; i < result_p.size(); i++) {
    fCompareHist_p->SetBinContent(i + 1, result_p[i]);
    fCompareHist_l->SetBinContent(i + 1, result_l[i]);
    fCompareHist_t->SetBinContent(i + 1, result_t[i]);
  }
}

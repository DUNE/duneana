#ifndef FLASHHYPOTHESISCOMPARISON_H
#define FLASHHYPOTHESISCOMPARISON_H

/*!
 * Title:   FlashHypothesisComparison Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description:
 * Class for comparing a flash hypothesis to MC truth (via SimPhotonCounter).
 * Needs a flash hypothesis and a SimPhotonCounter object as input.
 * Outputs a Tree with relevent info.
 */

#include "FlashHypothesis.h"
#include "FlashUtilities.h"
#include "SimPhotonCounter.h"

class TTree;
class TH1F;

#include <vector>

namespace opdet {

  class SimPhotonCounter;

  class FlashHypothesisComparison {

  public:
    FlashHypothesisComparison() {}

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
                          const unsigned int,
                          bool fill = true);

    void RunComparison(const unsigned int,
                       const unsigned int,
                       const FlashHypothesisCollection&,
                       const SimPhotonCounter&,
                       const std::vector<float>&,
                       const std::vector<float>&);

  private:
    FlashUtilities fUtil;

    void FillFlashHypothesisInfo(const FlashHypothesisCollection&,
                                 const std::vector<float>&,
                                 const std::vector<float>&);

    void FillSimPhotonCounterInfo(const SimPhotonCounter&,
                                  const std::vector<float>&,
                                  const std::vector<float>&);

    void FillComparisonInfo(const FlashHypothesisCollection&, const SimPhotonCounter&);

    bool fFillTree;
    TTree* fTree;

    TH1F* fHypHist_p;
    TH1F* fSimHist_p;
    TH1F* fCompareHist_p;
    TH1F* fHypHist_l;
    TH1F* fSimHist_l;
    TH1F* fCompareHist_l;
    TH1F* fHypHist_t;
    TH1F* fSimHist_t;
    TH1F* fCompareHist_t;

    unsigned int fRun;
    unsigned int fEvent;

    float fHypPEs_p;
    float fHypPEsError_p;
    float fSimPEs_p;
    float fHypY_p;
    float fSimY_p;
    float fHypRMSY_p;
    float fSimRMSY_p;
    float fHypZ_p;
    float fSimZ_p;
    float fHypRMSZ_p;
    float fSimRMSZ_p;
    float fCompare_p;

    float fHypPEs_l;
    float fHypPEsError_l;
    float fSimPEs_l;
    float fHypY_l;
    float fSimY_l;
    float fHypRMSY_l;
    float fSimRMSY_l;
    float fHypZ_l;
    float fSimZ_l;
    float fHypRMSZ_l;
    float fSimRMSZ_l;
    float fCompare_l;

    float fHypPEs_t;
    float fHypPEsError_t;
    float fSimPEs_t;
    float fHypY_t;
    float fSimY_t;
    float fHypRMSY_t;
    float fSimRMSY_t;
    float fHypZ_t;
    float fSimZ_t;
    float fHypRMSZ_t;
    float fSimRMSZ_t;
    float fCompare_t;
  };

}

#endif

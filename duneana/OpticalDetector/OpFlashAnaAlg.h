#ifndef OPFLASHANAALG_H
#define OPFLASHANAALG_H

/*!
 * Title:   OpFlashAnaAlg Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description:
 * Alg that stores flash information in a handy tree.
 *
*/

#include <string>
#include <vector>

#include "TH1.h"
class TTree;

#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"

namespace opdet {

  class OpFlashAnaAlg {

  public:
    OpFlashAnaAlg() { fMakeOpDetPEHist = false; }
    void SetOpFlashTree(TTree*, bool makeOpDetPEHist = true);
    void SetOpHitTree(TTree*);

    void FillOpFlashes(const std::vector<recob::OpFlash>&);
    void FillOpHits(const std::vector<recob::OpHit>&);

  private:
    struct FlashAnaStruct {
      double FlashTime;
      double FlashTimeWidth;
      double FlashAbsTime;
      double FlashY;
      double FlashYWidth;
      double FlashZ;
      double FlashZWidth;
      double FlashFastToTotal;
      double FlashTotalPE;
      unsigned int FlashFrame;
      int FlashOnBeamTime;
      bool FlashInBeamFrame;

      TH1D* FlashOpDetPEHist;

      std::string LeafList;
      FlashAnaStruct()
        : LeafList("time/D:timewidth/D:abstime/D:y/D:ywidth/D:z/D:zwidth/D:fasttototal/D:totalpe/"
                   "D:onbeamtime/I:frame/i:inbeamframe/O")
      {
        FlashOpDetPEHist = new TH1D();
      }
      ~FlashAnaStruct() { delete FlashOpDetPEHist; }
    };

    FlashAnaStruct fOpFlashAnaStruct;
    bool fMakeOpDetPEHist;

    TTree* fOpFlashTree;
    void FillOpFlashTree(const std::vector<recob::OpFlash>&);

    struct HitAnaStruct {
      double HitPeakTime;
      double HitPeakTimeAbs;
      double HitWidth;
      double HitArea;
      double HitAmplitude;
      double HitFastToTotal;
      double HitPE;
      unsigned int HitFrame;
      int HitOpChannel;

      std::string LeafList;
      HitAnaStruct()
        : LeafList(
            "time/D:abstime/D:width/D:area/D:amplitude/D:fasttototal/D:pe/D:frame/I:opchannel/i")
      {}
    };
    HitAnaStruct fOpHitAnaStruct;

    TTree* fOpHitTree;
    void FillOpHitTree(const std::vector<recob::OpHit>&);
  };

}

#endif

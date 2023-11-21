//
// - Serves as the link between ROOT "events" (e.g. mouse-clicks) and the ART
//   event display service by providing a receiver slot for signals generated
//   by the ROOT events.  A ROOT dictionary needs to be generated for this.
//

#ifndef EventDisplay3D_EvtDisplayUtils_hh
#define EventDisplay3D_EvtDisplayUtils_hh

// ... libCore
#include <TObject.h>
#include <TApplication.h>
#include <TStyle.h>

// ... libGui
#include <TGTextBuffer.h>

// ... libRGL
#include <TGLViewer.h>
#include <TGLCameraOverlay.h>

// ... libEve
#include <TEveManager.h>
#include <TEveGedEditor.h>
#include <TEveViewer.h>
#include <TEveScene.h>
#include <TEvePointSet.h>
#include <TEveQuadSet.h>
#include <TEveRGBAPalette.h>
#include <TEveRGBAPaletteOverlay.h>
#include <TEveFrameBox.h>
#include <TEveTrans.h>

#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "canvas/Persistency/Common/Ptr.h"
#include <iostream>

namespace dune
{
  class EvtDisplayUtils
  {
    public:
      explicit EvtDisplayUtils();
      void PrevEvent();
      void NextEvent();
      void GotoEvent();
      TGTextBuffer *fTbRun;
      TGTextBuffer *fTbEvt;

      void TrjPtSelAction(Int_t n);
      const geo::GeometryCore *geomcore;
      double fSampleRate;
      double fRecipDriftVel[3];
      int fNSamplesReadout;
      int fWinWid;
      int fWinHgt;
      std::map<raw::ChannelID_t, art::Ptr<raw::RawDigit>> rawdigitMap;

      TEveScene *fUplaneScene,*fVplaneScene,*fZplaneScene;
      TEveViewer *fUplaneViewer,*fVplaneViewer,*fZplaneViewer;
      TEveQuadSet* fQuadSetU,*fQuadSetV,*fQuadSetZ;
      TEveRGBAPalette *fRGBAPaletteU,*fRGBAPaletteV,*fRGBAPaletteZ;
   };
}
#endif /* EventDisplay3D_EvtDisplayUtils_hh */

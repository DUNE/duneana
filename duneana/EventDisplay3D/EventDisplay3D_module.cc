#include <random>

// ROOT includes

// ... libCore
#include <TApplication.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TStyle.h>

// ... libGui
#include <TGLabel.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TGTab.h>
#include <TG3DLine.h>

// ... libGeom
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>

// ... libRGL
#include <TGLViewer.h>
#include <TGLCameraOverlay.h>

// ... libEve
#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveBrowser.h>
#include <TEveGeoNode.h>
#include <TEveViewer.h>
#include <TEveGedEditor.h>
#include <TEvePointSet.h>
#include <TEveLine.h>
#include <TEveTrans.h>
#include <TEveFrameBox.h>
#include <TEveRGBAPalette.h>
#include <TEveScene.h>
#include <TEveRGBAPaletteOverlay.h>

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/View.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft libraries
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// DUNETPC specific includes
#include "dunecore/DuneInterface/Data/AdcTypes.h"
#include "dunecore/DuneInterface/Service/SimChannelExtractService.h"

#include "Style.h"
#include "EvtDisplayUtils.h"

// ... Anonymous namespace for helpers.
namespace {

  // ... Helper for setting color and transparency of detector objects
  void setRecursiveColorTransp(TGeoVolume *vol, Int_t color, Int_t transp)
  {
     if(color>=0)vol->SetLineColor(color);
     if(transp>=0)vol->SetTransparency(transp);
     Int_t nd = vol->GetNdaughters();
     for (Int_t i=0; i<nd; i++) {
        setRecursiveColorTransp(vol->GetNode(i)->GetVolume(), color, transp);
     }
  }

  // ... Helper for setting color and transparency of detector objects
  void setRecursiveTransp(TGeoVolume *vol, Int_t transp)
  {
     if(transp>=0)vol->SetTransparency(transp);
     Int_t nd = vol->GetNdaughters();
     for (Int_t i=0; i<nd; i++) {
        setRecursiveTransp(vol->GetNode(i)->GetVolume(), transp);
     }
  }
}

namespace dune
{
  class EventDisplay3D : public art::EDAnalyzer {

  public:

    explicit EventDisplay3D(fhicl::ParameterSet const& p);

    EventDisplay3D(EventDisplay3D const&) = delete;
    EventDisplay3D(EventDisplay3D&&) = delete;
    EventDisplay3D& operator=(EventDisplay3D const&) = delete;
    EventDisplay3D& operator=(EventDisplay3D&&) = delete;

    void reconfigure(fhicl::ParameterSet const & p);
    void beginJob() override;
    void endJob() override;
    void beginRun( const art::Run& run ) override;
    void analyze(art::Event const& e) override;

  private:

    // Set by parameter set variables.
    std::string	g4ModuleLabel_;
    std::string	digitModuleLabel_;
    bool          drawGenTrajPoints_;
    bool          drawGenPolyLines_;
    Double_t	camRotateCenterH_;
    Double_t	camRotateCenterV_;
    Double_t	camDollyDelta_;

    art::ServiceHandle<geo::Geometry> fgeom;
    art::ServiceHandle<detinfo::DetectorClocksService> fclock;
    art::ServiceHandle<detinfo::DetectorPropertiesService> fdetprop;

    std::unique_ptr<dune::EvtDisplayUtils>visutil_;

    TApplication* fTApplication;
    TGTextEntry *fTeRun,*fTeEvt;
    TGLabel     *fTlRun,*fTlEvt;

    TEveElementList *fPartTrjPtsList;
    TEveElementList *fPartPlinesList;

    void makeNavPanel();
    int GetParticle(const art::Event& evt, std::vector<const simb::MCParticle*>& plist);
  };

//-----------------------------------------------------------------------
EventDisplay3D::EventDisplay3D(fhicl::ParameterSet const& p):
  EDAnalyzer{p},
  g4ModuleLabel_    ( p.get<std::string>("g4ModuleLabel") ),
  digitModuleLabel_ ( p.get<std::string>("DigitModuleLabel", "simWire")),
  drawGenTrajPoints_( p.get<bool>       ("drawGenTrajPoints",true) ),
  drawGenPolyLines_ ( p.get<bool>       ("drawGenPolyLines",true)  ),
  camRotateCenterH_ ( p.get<Double_t>	("camRotateCenterH",-0.26) ),
  camRotateCenterV_ ( p.get<Double_t>	("camRotateCenterV",-2.  ) ),
  camDollyDelta_    ( p.get<Double_t>	("camDollyDelta",500.) ),
  visutil_(new dune::EvtDisplayUtils()),
  fTeRun(nullptr),fTeEvt(nullptr),
  fTlRun(nullptr),fTlEvt(nullptr),
  fPartTrjPtsList(0),fPartPlinesList(0)
{
  this->reconfigure(p);
}

//-----------------------------------------------------------------------
void EventDisplay3D::reconfigure(fhicl::ParameterSet const & p)
{
  return;
}

//-----------------------------------------------------------------------
void EventDisplay3D::makeNavPanel()
{
  // Create control panel for event navigation
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TEveBrowser* browser = gEve->GetBrowser();
  browser->StartEmbedding(TRootBrowser::kLeft); // insert nav frame as new tab in left pane

  TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
  //frmMain->SetWindowName("EVT NAV");
  frmMain->SetCleanup(kDeepCleanup);

  TGHorizontalFrame* navFrame = new TGHorizontalFrame(frmMain);
  TGVerticalFrame* evtidFrame = new TGVerticalFrame(frmMain);
  {
    TString icondir(TString::Format("%s/icons/", gSystem->Getenv("ROOTSYS")) );
    TGPictureButton* b = 0;

    // ... Create back button and connect to "PrevEvent" rcvr in visutils
    b = new TGPictureButton(navFrame, gClient->GetPicture(icondir + "GoBack.gif"));
    navFrame->AddFrame(b);
    b->Connect("Clicked()", "dune::EvtDisplayUtils", visutil_.get(), "PrevEvent()");

    // ... Create forward button and connect to "NextEvent" rcvr in visutils
    b = new TGPictureButton(navFrame, gClient->GetPicture(icondir + "GoForward.gif"));
    navFrame->AddFrame(b);
    b->Connect("Clicked()", "dune::EvtDisplayUtils", visutil_.get(), "NextEvent()");

    // ... Create run num text entry widget and connect to "GotoEvent" rcvr in visutils
    TGHorizontalFrame* runoFrame = new TGHorizontalFrame(evtidFrame);
    fTlRun = new TGLabel(runoFrame,"Run Number");
    fTlRun->SetTextJustify(kTextLeft);
    fTlRun->SetMargins(5,5,5,0);
    runoFrame->AddFrame(fTlRun);
    
    fTeRun = new TGTextEntry(runoFrame, visutil_->fTbRun = new TGTextBuffer(5));
    visutil_->fTbRun->AddText(0, "1");
    fTeRun->Connect("ReturnPressed()","dune::EvtDisplayUtils", visutil_.get(),"GotoEvent()");
    runoFrame->AddFrame(fTeRun,new TGLayoutHints(kLHintsExpandX));

    // ... Create evt num text entry widget and connect to "GotoEvent" rcvr in visutils
    TGHorizontalFrame* evnoFrame = new TGHorizontalFrame(evtidFrame);
    fTlEvt = new TGLabel(evnoFrame,"Evt Number");
    fTlEvt->SetTextJustify(kTextLeft);
    fTlEvt->SetMargins(5,5,5,0);
    evnoFrame->AddFrame(fTlEvt);

    fTeEvt = new TGTextEntry(evnoFrame, visutil_->fTbEvt = new TGTextBuffer(5));
    visutil_->fTbEvt->AddText(0, "1");
    fTeEvt->Connect("ReturnPressed()","dune::EvtDisplayUtils", visutil_.get(),"GotoEvent()");
    evnoFrame->AddFrame(fTeEvt,new TGLayoutHints(kLHintsExpandX));

    // ... Add horizontal run & event number subframes to vertical evtidFrame
    evtidFrame->AddFrame(runoFrame,new TGLayoutHints(kLHintsExpandX));
    evtidFrame->AddFrame(evnoFrame,new TGLayoutHints(kLHintsExpandX));

    // ... Add navFrame and evtidFrame to MainFrame
    frmMain->AddFrame(navFrame);
    TGHorizontal3DLine *separator = new TGHorizontal3DLine(frmMain);
    frmMain->AddFrame(separator, new TGLayoutHints(kLHintsExpandX));
    frmMain->AddFrame(evtidFrame);

    frmMain->MapSubwindows();
    frmMain->Resize();
    frmMain->MapWindow();

    browser->StopEmbedding();
    browser->SetTabTitle("Event Nav", 0);
  }
}

//-----------------------------------------------------------------------
int EventDisplay3D::GetParticle(const art::Event& evt, std::vector<const simb::MCParticle*>& plist)
{
  plist.clear();

  if (evt.isRealData()) return 0;

  std::vector<const simb::MCParticle*> temp;

  art::View<simb::MCParticle> plcol;
  // use get by Type because there should only be one collection of these in the event
  try {
    evt.getView(g4ModuleLabel_, plcol);
    for (unsigned int i = 0; i < plcol.vals().size(); ++i) {
      temp.push_back(plcol.vals().at(i));
    }
    temp.swap(plist);
  }
  catch (cet::exception& err) {
    throw cet::exception("EventDisplay3D::GetParticle")
          << "error: unable to get MCParticle collection: " << err << std::endl;
  }

  return plist.size();
}

//-----------------------------------------------------------------------
void EventDisplay3D::beginJob()
{

  // ... Set up the physical constants
  
  //  .. calculate the drift velocity for each region:
  //     0: main drift region
  //     1: between first two planes
  //     2: between last two planes
  auto const detProp =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  for (int i = 0; i < 3; ++i) {
    double driftVelocity = detProp.DriftVelocity(detProp.Efield(i),
  			                         detProp.Temperature())
			   *1.e-3; // cm/us -> cm/ns (standard LArSoft units)
    visutil_->fRecipDriftVel[i] = 1. / driftVelocity;
  }

  //  .. get the sampling rate (time between ticks)
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  auto const detProp2 =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);
  visutil_->fSampleRate = sampling_rate(clockData);
  visutil_->fNSamplesReadout = detProp.NumberTimeSamples();
  //visutil_->fNSamplesReadout = detProp.ReadOutWindowSize();

  //gDebug=3;
  
  // ... No need to start up TApplication here since EvtDisplayService will do that
  //     via EnsureTApplication.  All we need to do here is initialize the graphics
  //     and switch to non-batch mode
  TApplication::NeedGraphicsLibs();
  gApplication->InitializeGraphics();
  gROOT->SetBatch(kFALSE);

  // Initialize global Eve application manager (return gEve)
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TEveManager::Create(kTRUE);
  //TEveManager::Create(kTRUE,"FIVV");

  // Get the detector geometry
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ... Get top volume and find main detector enclosure underneath it
  TGeoVolume* topvol = gGeoManager->GetTopVolume();
  TEveGeoTopNode* detnode = new TEveGeoTopNode(gGeoManager,topvol->FindNode("volDetEnclosure_0"));

  detnode->SetVisLevel(0);
  setRecursiveTransp(detnode->GetNode()->GetVolume(), 60);
  TEveTrans& tr0 = detnode->RefMainTrans();
  tr0.SetPos(0., 120.5, 571.); // 695.

  // ... Add static detector geometry to global scene
  gEve->AddGlobalElement(detnode);

  // Create multiview tab
  // ~~~~~~~~~~~~~~~~~~~~
  TEveWindowSlot *slot = 0;
  TEveViewer *v = 0;

  slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
  TEveWindowPack* pack1 = slot->MakePack();
  pack1->SetShowTitleBar(kFALSE);
  pack1->SetHorizontal();
  pack1->SetElementName("Plane Views");

  // Embedded viewer.
  slot = pack1->NewSlot();
  v = new TEveViewer("TPCViewer");
  v->SpawnGLEmbeddedViewer(gEve->GetEditor());
  slot->ReplaceWindow(v);
  v->SetElementName("Embedded TPCViewer");

  gEve->GetViewers()->AddElement(v);
  v->AddScene(gEve->GetEventScene());
  v->AddScene(gEve->GetGlobalScene());

  TGLViewer *glv1 = v->GetGLViewer();
  glv1->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);  // other option: kAxesEdge
  Double_t camcent[3] = { 0., 0., 697. };
  glv1->SetPerspectiveCamera(TGLViewer::kCameraPerspXOZ, 37.8, 1000., camcent, -0.34, 0.52);
  glv1->SetDrawCameraCenter(kTRUE);

  slot = pack1->NewSlot();
  TEveWindowPack* pack2 = slot->MakePack();
  pack2->SetShowTitleBar(kFALSE);

  // .. set up the wire plane viewers
  //    =============================

  // .. U-induction plane viewer
  //    ------------------------
  slot = pack2->NewSlot();
  visutil_->fUplaneViewer = new TEveViewer("U-Plane Viewer");
  visutil_->fUplaneViewer->SpawnGLViewer(gEve->GetEditor());
  slot->ReplaceWindow(visutil_->fUplaneViewer);
  gEve->GetViewers()->AddElement(visutil_->fUplaneViewer);
  visutil_->fUplaneScene = gEve->SpawnNewScene("Scene - Uplane");
  visutil_->fUplaneViewer->AddScene(visutil_->fUplaneScene);

  TGLViewer* glvwr = visutil_->fUplaneViewer->GetGLViewer();
  glvwr->SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
  TGLCameraOverlay* co = glvwr->GetCameraOverlay();
  co->SetShowOrthographic(kTRUE);
  co->SetOrthographicMode(TGLCameraOverlay::kGridFront);

  visutil_->fRGBAPaletteU = new TEveRGBAPalette(-300, 300, kTRUE, kTRUE, kFALSE);
  visutil_->fRGBAPaletteU->ClearColorArray();
  visutil_->fRGBAPaletteU->SetupColorArray();

  TEveRGBAPaletteOverlay *po = new TEveRGBAPaletteOverlay(visutil_->fRGBAPaletteU, 0.3, 0.15, 0.4, 0.05);
  glvwr->AddOverlayElement(po);

  // .. V-induction plane viewer
  //    ------------------------
  slot = pack2->NewSlot();
  visutil_->fVplaneViewer = new TEveViewer("V-Plane Viewer");
  visutil_->fVplaneViewer->SpawnGLViewer(gEve->GetEditor());
  slot->ReplaceWindow(visutil_->fVplaneViewer);
  gEve->GetViewers()->AddElement(visutil_->fVplaneViewer);
  visutil_->fVplaneScene = gEve->SpawnNewScene("Scene - Vplane");
  visutil_->fVplaneViewer->AddScene(visutil_->fVplaneScene);

  glvwr = visutil_->fVplaneViewer->GetGLViewer();
  glvwr->SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
  co = glvwr->GetCameraOverlay();
  co->SetShowOrthographic(kTRUE);
  co->SetOrthographicMode(TGLCameraOverlay::kGridFront);

  visutil_->fRGBAPaletteV = new TEveRGBAPalette(-300, 300, kTRUE, kTRUE, kFALSE);
  visutil_->fRGBAPaletteV->ClearColorArray();
  visutil_->fRGBAPaletteV->SetupColorArray();

  po = new TEveRGBAPaletteOverlay(visutil_->fRGBAPaletteV, 0.3, 0.15, 0.4, 0.05);
  glvwr->AddOverlayElement(po);

  // .. Z-collection plane viewer
  //    ------------------------
  slot = pack2->NewSlot();
  visutil_->fZplaneViewer = new TEveViewer("Z-Plane Viewer");
  visutil_->fZplaneViewer->SpawnGLViewer(gEve->GetEditor());
  slot->ReplaceWindow(visutil_->fZplaneViewer);
  gEve->GetViewers()->AddElement(visutil_->fZplaneViewer);
  visutil_->fZplaneScene = gEve->SpawnNewScene("Scene - Zplane");
  visutil_->fZplaneViewer->AddScene(visutil_->fZplaneScene);

  glvwr = visutil_->fZplaneViewer->GetGLViewer();
  glvwr->SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
  co = glvwr->GetCameraOverlay();
  co->SetShowOrthographic(kTRUE);
  co->SetOrthographicMode(TGLCameraOverlay::kGridFront);

  visutil_->fRGBAPaletteZ = new TEveRGBAPalette(-300, 300, kTRUE, kTRUE, kFALSE);
  visutil_->fRGBAPaletteZ->ClearColorArray();
  visutil_->fRGBAPaletteZ->SetupColorArray();

  po = new TEveRGBAPaletteOverlay(visutil_->fRGBAPaletteZ, 0.3, 0.15, 0.4, 0.05);
  glvwr->AddOverlayElement(po);

  gEve->GetBrowser()->GetTabRight()->SetTab(0);

  // Create navigation panel
  // ~~~~~~~~~~~~~~~~~~~~~~~~
  makeNavPanel();

  // Add new Eve event into the "Event" scene and make it the current event
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // (Subsequent elements added using "AddElements" will be added to this event)
  gEve->AddEvent(new TEveEventManager("Event", "DUNE Detector Event"));

  // ... Set up initial camera orientation in main 3d view
  TGLViewer *glv = gEve->GetDefaultGLViewer();
  glv->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);  // other option: kAxesEdge
  //Double_t camcent[3] = { 0., 0., 697. };
  glv->SetPerspectiveCamera(TGLViewer::kCameraPerspXOZ, 37.8, 1000., camcent, -0.34, 0.52);
  glv->SetDrawCameraCenter(kTRUE);

}

//-----------------------------------------------------------------------
void EventDisplay3D::endJob()
{
}

//-----------------------------------------------------------------------
void EventDisplay3D::beginRun( const art::Run& )
{
  //if(gGeoManager){
  //  gGeoManager->GetListOfNodes()->Delete();
  //  gGeoManager->GetListOfVolumes()->Delete();
  //  gGeoManager->GetListOfShapes()->Delete();
  //}
  //gEve->GetGlobalScene()->DestroyElements();
}

//-----------------------------------------------------------------------
void EventDisplay3D::analyze(art::Event const& evt)
{
  std::cout << "Event "
            << " " << evt.id().run()
            << " " << evt.id().subRun()
            << " " << evt.id().event() << std::endl;

  // ... Read in the digit List object(s).
  art::Handle<std::vector<raw::RawDigit>> digitVecHandle;
  std::vector<art::Ptr<raw::RawDigit>> rawdigitlist;
  if (evt.getByLabel(digitModuleLabel_, digitVecHandle)) {
    art::fill_ptr_vector(rawdigitlist, digitVecHandle);
  }

  const auto& digitVec0 = rawdigitlist[0];
  if (int(digitVec0->Samples()) != visutil_->fNSamplesReadout) {
    std::cout << "!!!!! Bad dataSize: " << digitVec0->Samples() << std::endl;
    return;
  }
  std::cout << "dataSize: " << digitVec0->Samples() << std::endl;
  std::cout << "rawdigitlist size: " << rawdigitlist.size() << std::endl;

  // ... Build a map from channel number -> rawdigitVec
  std::map<raw::ChannelID_t, art::Ptr<raw::RawDigit>> rawdigitMap;
  raw::ChannelID_t chnum = raw::InvalidChannelID;
  for (unsigned int i = 0; i < rawdigitlist.size(); ++i) {
    const auto& digitVec = rawdigitlist[i];
    chnum = digitVec->Channel();
    if (chnum == raw::InvalidChannelID) continue;
    visutil_->rawdigitMap[chnum] = digitVec;
  }

  visutil_->geomcore = fgeom->provider();

  // ... Update the run and event numbers in the TGTextEntry widgets in the Navigation panel
  std::ostringstream sstr;
  sstr << evt.id().run();
  visutil_->fTbRun->Clear();
  visutil_->fTbRun->AddText(0,sstr.str().c_str());
  gClient->NeedRedraw(fTeRun);

  sstr.str("");
  sstr << evt.id().event();
  visutil_->fTbEvt->Clear();
  visutil_->fTbEvt->AddText(0,sstr.str().c_str());
  gClient->NeedRedraw(fTeEvt);

  // ... Delete visualization structures associated with previous event
  gEve->GetViewers()->DeleteAnnotations();
  gEve->GetCurrentEvent()->DestroyElements();

  std::vector<const simb::MCParticle*> plist;
  if ( drawGenTrajPoints_ || drawGenPolyLines_ ){
    this->GetParticle(evt, plist);
  }

  // Draw individual points of generated particle trajectories
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (drawGenTrajPoints_) {

    if (fPartTrjPtsList == 0) {
      fPartTrjPtsList = new TEveElementList("Particle trajectory points"); 
      fPartTrjPtsList->IncDenyDestroy();	      // protect element against destruction
    }
    else {
      fPartTrjPtsList->DestroyElements();	      // destroy children of the element
    }

    for (size_t p = 0; p < plist.size(); ++p) {

      const simb::MCParticle* mcPart = plist[p];
      const simb::MCTrajectory& mcTraj = mcPart->Trajectory();

      int colorIdx(dune::Style::ColorFromPDG(mcPart->PdgCode()));

      if (!mcTraj.empty()) {

        TEvePointSet* trjpt = new TEvePointSet(Form("Particle %d traj pts",int(p+1))); 
        trjpt->SetTitle(Form("TrjPts for Particle %d",int(p+1)));

        int numTrajPoints = mcTraj.size();
        for (int hitIdx = 0; hitIdx < numTrajPoints; hitIdx++) {
          double xPos = mcTraj.X(hitIdx);
          double yPos = mcTraj.Y(hitIdx);
          double zPos = mcTraj.Z(hitIdx);

          trjpt->SetNextPoint(xPos,yPos,zPos);
          //trjpt->SetTitle(Form("E = %.3f MeV",mcTraj.E(hitIdx)*1000.));
	  trjpt->SetPointId(new TNamed(Form("TrjPt %d, E=%.3f MeV", hitIdx,mcTraj.E(hitIdx)*1000.), ""));
	  trjpt->SetMarkerSize(1);
          if ( hitIdx==0 ) {
            trjpt->SetMarkerColor(kRed);
          } else {
            trjpt->SetMarkerColor(colorIdx);
          }
          trjpt->SetPickable(1);
	  trjpt->Connect(trjpt,"PointSelected(Int_t)", "dune::EvtDisplayUtils", visutil_.get(), "TrjPtSelAction(Int_t)");
        }
        fPartTrjPtsList->AddElement(trjpt);
      }
    }
    gEve->AddElement(fPartTrjPtsList);
  }

  // Draw polylines representing generated particle trajectories
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (drawGenPolyLines_) {

    if (fPartPlinesList == 0) {
      fPartPlinesList = new TEveElementList("Particle trajectory plines"); 
      fPartPlinesList->IncDenyDestroy();	      // protect element against destruction
    }
    else {
      fPartPlinesList->DestroyElements();	      // destroy children of the element
    }

    for (size_t p = 0; p < plist.size(); ++p) {

      const simb::MCParticle* mcPart = plist[p];
      const simb::MCTrajectory& mcTraj = mcPart->Trajectory();

      int colorIdx(dune::Style::ColorFromPDG(mcPart->PdgCode()));
      int linestyleIdx(dune::Style::LineStyleFromPDG(mcPart->PdgCode()));
      double partEnergy = mcPart->E();

      if (!mcTraj.empty()) {

	TEveLine* gentrkseg = new TEveLine(Form("Particle %d traj pline",int(p+1)));
	gentrkseg->SetTitle(Form("PartID: %s\nE_init = %.3f MeV\nProcess: %s\nEndProcess: %s",
	                         dune::Style::LatexName(mcPart->PdgCode()),partEnergy*1000.,
                                 mcPart->Process().c_str(),mcPart->EndProcess().c_str() ));

        int numTrajPoints = mcTraj.size();
        for (int hitIdx = 0; hitIdx < numTrajPoints; hitIdx++) {
          double xPos = mcTraj.X(hitIdx);
          double yPos = mcTraj.Y(hitIdx);
          double zPos = mcTraj.Z(hitIdx);

	  gentrkseg->SetMainColor(colorIdx);
	  gentrkseg->SetLineStyle(linestyleIdx);
	  gentrkseg->SetLineWidth(3);
	  gentrkseg->SetNextPoint(xPos,yPos,zPos);
        }
        fPartPlinesList->AddElement(gentrkseg);
      }
    }
    gEve->AddElement(fPartPlinesList);

  }
}
}
DEFINE_ART_MODULE(dune::EventDisplay3D)

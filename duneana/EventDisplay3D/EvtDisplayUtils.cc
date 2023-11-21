//
// - Serves as the link between ROOT "events" (e.g. mouse-clicks) and the ART
//   event display service by providing a receiver slot for signals generated
//   by the ROOT events.  A ROOT dictionary needs to be generated for this.
//

#include "NavState.h"
#include "EvtDisplayUtils.h"
#include <string>

namespace dune
{
  EvtDisplayUtils::EvtDisplayUtils():
    fTbRun(0),fTbEvt(0),fQuadSetU(0),fQuadSetV(0),fQuadSetZ(0){
  }
  void EvtDisplayUtils::PrevEvent(){
    NavState::Set(kPREV_EVENT);
  }
  void EvtDisplayUtils::NextEvent(){
    NavState::Set(kNEXT_EVENT);
  }
  void EvtDisplayUtils::GotoEvent(){
    int run = std::stoi(fTbRun->GetString());
    int event = std::stoi(fTbEvt->GetString());
    NavState::SetTarget(run, event);
    NavState::Set(kGOTO_EVENT);
  }
  void EvtDisplayUtils::TrjPtSelAction(Int_t n){

    Float_t x, y, z;

    // .. get the TEvePointSet that was selected and get its coordinates
    auto ps = dynamic_cast<TEvePointSet*> (reinterpret_cast<TQObject*>(gTQSender)); // see eve/paramlist.C
    ps->GetPoint(n,x,y,z);

    geo::Point_t xu(x,y,z);	// Point_t is of type ROOT::Math::PositionVector3D
    geo::Point_t xv(x,y,z);
    geo::Point_t xz(x,y,z);
    geo::TPCID tID = geomcore->FindTPCAtPosition(xu);    // which TPC is this point in?

    if (tID.TPC!=geo::TPCID::InvalidID) {
      geo::PlaneID const pIDu(tID, 0);  //larcoreobj/SimpleTypesAndConstants/geo_types.h
      geo::PlaneID const pIDv(tID, 1);
      geo::PlaneID const pIDz(tID, 2);

      geo::PlaneGeo const& planeU = geomcore->Plane(pIDu); //larcorealg/Geometry/Plane.h
      geo::PlaneGeo const& planeV = geomcore->Plane(pIDv);
      geo::PlaneGeo const& planeZ = geomcore->Plane(pIDz);

      double du = planeU.DistanceFromPlane(xu);
      double dv = planeV.DistanceFromPlane(xv);
      double dz = planeZ.DistanceFromPlane(xz);
    
      planeU.DriftPoint(xu, du);                         // drift points to planes
      planeU.DriftPoint(xv, dv);
      planeU.DriftPoint(xz, dz);

      geo::WireID const wIDu = planeU.NearestWireID(xu); // find wire closest to drifted pt
      geo::WireID const wIDv = planeV.NearestWireID(xv);
      geo::WireID const wIDz = planeZ.NearestWireID(xz);

      double tu = du * fRecipDriftVel[0];
      double tv = tu + (dv - du) * fRecipDriftVel[1];
      double tz = tv + (dz - dv) * fRecipDriftVel[2];
      int ticku  = tu/fSampleRate;
      int tickv  = tv/fSampleRate;
      int tickz  = tz/fSampleRate;

      std::cout << "Number of samples in readout window: " << fNSamplesReadout << std::endl;
      fWinWid = 256;
      fWinHgt = 64;

      int tick1u = ticku - fWinWid/2;
      int tick1v = tickv - fWinWid/2;
      int tick1z = tickz - fWinWid/2;

      if (tick1u < 0) tick1u = 0;
      if (tick1v < 0) tick1v = 0;
      if (tick1z < 0) tick1z = 0;
      if ((tick1u+fWinWid) > fNSamplesReadout) tick1u = fNSamplesReadout - fWinWid;
      if ((tick1v+fWinWid) > fNSamplesReadout) tick1v = fNSamplesReadout - fWinWid;
      if ((tick1z+fWinWid) > fNSamplesReadout) tick1z = fNSamplesReadout - fWinWid;
      
      int wire1u = wIDu.Wire - fWinHgt/2;
      int wire1v = wIDv.Wire - fWinHgt/2;
      int wire1z = wIDz.Wire - fWinHgt/2;
      if (wire1u < 0) wire1u = 0;
      if (wire1v < 0) wire1v = 0;
      if (wire1z < 0) wire1z = 0;
      if ((wire1u+fWinHgt) > int(planeU.Nwires())) wire1u = planeU.Nwires() - fWinHgt;
      if ((wire1v+fWinHgt) > int(planeV.Nwires())) wire1v = planeV.Nwires() - fWinHgt;
      if ((wire1z+fWinHgt) > int(planeZ.Nwires())) wire1z = planeZ.Nwires() - fWinHgt;

      std::cout << "U-plane ticks range from " << tick1u << " to " << tick1u+fWinWid-1 <<std::endl;
      std::cout << "U-plane wires range from " << wire1u << " to " << wire1u+fWinHgt-1 <<std::endl;
      std::cout << "V-plane ticks range from " << tick1v << " to " << tick1v+fWinWid-1 <<std::endl;
      std::cout << "V-plane wires range from " << wire1v << " to " << wire1v+fWinHgt-1 <<std::endl;
      std::cout << "Z-plane ticks range from " << tick1z << " to " << tick1z+fWinWid-1 <<std::endl;
      std::cout << "Z-plane wires range from " << wire1z << " to " << wire1z+fWinHgt-1 <<std::endl;
      
      // -----------------------------------------------------------------------

      gStyle->SetPalette(kLake); // 51,57,97

      TEveFrameBox *box = new TEveFrameBox();
      box->SetAAQuadXY(0, 0, 0, fWinWid, fWinHgt);
      box->SetFrameColor(kGray);
      short minadc, maxadc;
      size_t k;

      // .... U-Plane
      //      -------
      if (fQuadSetU!=0) fQuadSetU->DestroyElements(); 	    // destroy children of the element
      fQuadSetU = new TEveQuadSet("RectangleXY"); 
      fQuadSetU->IncDenyDestroy();  	    // protect element against destruction
      fQuadSetU->SetOwnIds(kTRUE);
      fQuadSetU->SetFrame(box);
      fQuadSetU->Reset(TEveQuadSet::kQT_RectangleXY, kFALSE, 64);

      minadc = 32767;
      maxadc = -32768;
      k = 0;
      for (int iw=wire1u; iw<wire1u+fWinHgt; ++iw){
        geo::WireID const wID (pIDu, iw);
        raw::ChannelID_t chnum = geomcore->PlaneWireToChannel(wID);

        auto search = rawdigitMap.find(chnum);
        if (search == rawdigitMap.end()) {
          continue;
        }      

        art::Ptr<raw::RawDigit> digitVec = (*search).second;
        std::vector<short> rawadc(fNSamplesReadout);
        raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->GetPedestal(), digitVec->Compression());
        size_t j = 0;
        for (size_t i = size_t(tick1u); i < size_t(tick1u+fWinWid); ++i) {
          short adcvec = rawadc[i] - digitVec->GetPedestal();
	  if (adcvec<minadc) minadc=adcvec;
	  if (adcvec>maxadc) maxadc=adcvec;
          fQuadSetU->AddQuad(j, k, 0, 1, 1);
          fQuadSetU->QuadValue(adcvec);
          fQuadSetU->QuadId(new TNamed(Form("QuadIdx %d", adcvec),
  			  "ADC value."));
	  j++;
        }
	k++;
      }
      fQuadSetU->RefitPlex();
      std::cout << " !!!! Uplane Min ADC: " << minadc << std::endl;
      std::cout << " !!!! Uplane Max ADC: " << maxadc << std::endl;

      fRGBAPaletteU->SetMinMax(minadc,maxadc);
      fRGBAPaletteU->SetLimits(minadc,maxadc);
      fQuadSetU->SetPalette(fRGBAPaletteU);
    
      TEveTrans& xfu = fQuadSetU->RefMainTrans();
      xfu.RotateLF(1, 3, 0.5*TMath::Pi());
      xfu.SetPos(0, 0, 0);

      fQuadSetU->SetPickable(1);
      fQuadSetU->SetAlwaysSecSelect(1);

      fUplaneScene->DestroyElements();
      fUplaneScene->AddElement(fQuadSetU);
      fUplaneViewer->Redraw(kTRUE);


      // .... V-Plane
      //      -------
      if (fQuadSetV!=0) fQuadSetV->DestroyElements(); 	    // destroy children of the element
      fQuadSetV = new TEveQuadSet("RectangleXY"); 
      fQuadSetV->IncDenyDestroy();  	    // protect element against destruction
      fQuadSetV->SetOwnIds(kTRUE);
      fQuadSetV->SetFrame(box);
      fQuadSetV->Reset(TEveQuadSet::kQT_RectangleXY, kFALSE, 64);

      minadc = 32767;
      maxadc = -32768;
      k = 0;
      for (int iw=wire1v; iw<wire1v+fWinHgt; ++iw){
        geo::WireID const wID (pIDv, iw);
        raw::ChannelID_t chnum = geomcore->PlaneWireToChannel(wID);

        auto search = rawdigitMap.find(chnum);
        if (search == rawdigitMap.end()) {
          continue;
        }      

        art::Ptr<raw::RawDigit> digitVec = (*search).second;
        std::vector<short> rawadc(fNSamplesReadout);
        raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->GetPedestal(), digitVec->Compression());
        size_t j = 0;
        for (size_t i = size_t(tick1v); i < size_t(tick1v+fWinWid); ++i) {
          short adcvec = rawadc[i] - digitVec->GetPedestal();
	  if (adcvec<minadc) minadc=adcvec;
	  if (adcvec>maxadc) maxadc=adcvec;
          fQuadSetV->AddQuad(j, k, 0, 1, 1);
          fQuadSetV->QuadValue(adcvec);
          fQuadSetV->QuadId(new TNamed(Form("ADC = %d", adcvec),
  			  "ADC value."));
	  j++;
        }
	k++;
      }
      fQuadSetV->RefitPlex();
      std::cout << " !!!! Vplane Min ADC: " << minadc << std::endl;
      std::cout << " !!!! Vplane Max ADC: " << maxadc << std::endl;

      fRGBAPaletteV->SetMinMax(minadc,maxadc);
      fRGBAPaletteV->SetLimits(minadc,maxadc);
      fQuadSetV->SetPalette(fRGBAPaletteV);
    
      TEveTrans& xfv = fQuadSetV->RefMainTrans();
      xfv.RotateLF(1, 3, 0.5*TMath::Pi());
      xfv.SetPos(0, 0, 0);

      fQuadSetV->SetPickable(1);
      fQuadSetV->SetAlwaysSecSelect(1);

      fVplaneScene->DestroyElements();
      fVplaneScene->AddElement(fQuadSetV);
      fVplaneViewer->Redraw(kTRUE);

      // .... Z-Plane
      //      -------
      if (fQuadSetZ!=0) fQuadSetZ->DestroyElements(); 	    // destroy children of the element
      fQuadSetZ = new TEveQuadSet("RectangleXY"); 
      fQuadSetZ->IncDenyDestroy();  	    // protect element against destruction
      fQuadSetZ->SetOwnIds(kTRUE);
      fQuadSetZ->SetFrame(box);
      fQuadSetZ->Reset(TEveQuadSet::kQT_RectangleXY, kFALSE, 64);

      minadc = 32767;
      maxadc = -32768;
      k = 0;
      for (int iw=wire1z; iw<wire1z+fWinHgt; ++iw){
        geo::WireID const wID (pIDz, iw);
        raw::ChannelID_t chnum = geomcore->PlaneWireToChannel(wID);

        auto search = rawdigitMap.find(chnum);
        if (search == rawdigitMap.end()) {
          continue;
        }      

        art::Ptr<raw::RawDigit> digitVec = (*search).second;
        std::vector<short> rawadc(fNSamplesReadout);
        raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->GetPedestal(), digitVec->Compression());
        size_t j = 0;
        for (size_t i = size_t(tick1z); i < size_t(tick1z+fWinWid); ++i) {
          short adcvec = rawadc[i] - digitVec->GetPedestal();
	  if (adcvec<minadc) minadc=adcvec;
	  if (adcvec>maxadc) maxadc=adcvec;
          fQuadSetZ->AddQuad(j, k, 0, 1, 1);
          fQuadSetZ->QuadValue(adcvec);
          fQuadSetZ->QuadId(new TNamed(Form("ADC = %d", adcvec),
  			  "ADC value."));
	  j++;
        }
	k++;
      }
      fQuadSetZ->RefitPlex();
      std::cout << " !!!! Zplane Min ADC: " << minadc << std::endl;
      std::cout << " !!!! Zplane Max ADC: " << maxadc << std::endl;

      fRGBAPaletteZ->SetMinMax(minadc,maxadc);
      fRGBAPaletteZ->SetLimits(minadc,maxadc);
      fQuadSetZ->SetPalette(fRGBAPaletteZ);
    
      TEveTrans& xfz = fQuadSetZ->RefMainTrans();
      xfz.RotateLF(1, 3, 0.5*TMath::Pi());
      xfz.SetPos(0, 0, 0);

      fQuadSetZ->SetPickable(1);
      fQuadSetZ->SetAlwaysSecSelect(1);

      fZplaneScene->DestroyElements();
      fZplaneScene->AddElement(fQuadSetZ);
      fZplaneViewer->Redraw(kTRUE);



    } else {
      std::cout << " !!!! Invalid TPCID, point is not in a TPC" << std::endl;
    }

  }
}

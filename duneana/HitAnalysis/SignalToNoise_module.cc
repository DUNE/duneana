////////////////////////////////////////////////////////////////////////
// Class:       SignalToNoise
// Module Type: analyzer
// File:        SignalToNoise_module.cc
//
// Generated at Sat Sep 10 22:01:34 2016 by Tingjun Yang using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcorealg/CoreUtils/NumericUtils.h"
#include "larcorealg/Geometry/Intersections.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/ArtDataHelper/TrackUtils.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "larreco/RecoAlg/LinFitAlg.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TTree.h"

#include <vector>

namespace dune {
  class SignalToNoise;
}

class dune::SignalToNoise : public art::EDAnalyzer {
public:
  explicit SignalToNoise(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SignalToNoise(SignalToNoise const &) = delete;
  SignalToNoise(SignalToNoise &&) = delete;
  SignalToNoise & operator = (SignalToNoise const &) = delete;
  SignalToNoise & operator = (SignalToNoise &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // Declare member data here.
  art::InputTag fTrackModuleLabel;
  art::InputTag fExternalCounterModuleLabel;
  art::InputTag fRawDigitModuleLabel;
  art::InputTag fCalDataModuleLabel;

  bool fIsData;

  trkf::LinFitAlg fLinFitAlg;

  double cx[93], cy[93], cz[93];
  
  TH1D *dcos;
  TH2D *hyz, *hxz;
  TH1D *hsig[8][3][4];
  TH1D *hbkg[8][3][4];
  TH1D *hstb[8][3][4];
  TH1D *hdx;
  TH2D *hphx;
  TH1D *hdt;

  TTree *ftree;
  int run;
  int event;
  int tpc;
  int wire;
  int trackid;
  double dqdx;
  double x;
};


dune::SignalToNoise::SignalToNoise(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fTrackModuleLabel(p.get< art::InputTag >("TrackModuleLabel")),
  fExternalCounterModuleLabel(p.get< art::InputTag >("ExternalCounterModuleLabel")),
  fRawDigitModuleLabel(p.get< art::InputTag >("RawDigitModuleLabel")),
  fCalDataModuleLabel(p.get< art::InputTag >("CalDataModuleLabel")),
  fIsData(p.get<bool>("IsData",true))
{
}

void dune::SignalToNoise::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  run = e.run();
  event = e.id().event();

  // * tracks
  std::vector<art::Ptr<recob::Track> > tracklist;
  auto trackListHandle = e.getHandle< std::vector<recob::Track> >(fTrackModuleLabel);
  if (trackListHandle)
    art::fill_ptr_vector(tracklist, trackListHandle);

  // * External Counters
  std::vector<art::Ptr<raw::ExternalTrigger> > countlist;
  auto countListHandle = e.getHandle< std::vector<raw::ExternalTrigger> >(fExternalCounterModuleLabel);
  if (countListHandle)
    art::fill_ptr_vector(countlist, countListHandle);
 
  // * Raw Digits
  std::vector<art::Ptr<raw::RawDigit> > rawlist;
  auto rawListHandle = e.getHandle<std::vector<raw::RawDigit> >(fRawDigitModuleLabel);
  if (rawListHandle)
    art::fill_ptr_vector(rawlist, rawListHandle);

  // * Noise removed
  std::vector<art::Ptr<recob::Wire> > wirelist;
  auto wireListHandle = e.getHandle<std::vector<recob::Wire> >(fCalDataModuleLabel);
  if (wireListHandle)
    art::fill_ptr_vector(wirelist, wireListHandle);

  // * associations
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, e, fTrackModuleLabel);

  geo::WireReadoutGeom const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  const lariov::ChannelStatusProvider* fCSP = &art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();

  std::vector<int> c1;  //counter index
  std::vector<int> c2;
  std::vector<int> ci1; //trigger id
  std::vector<int> ci2;
 
  for (size_t i = 0; i<countlist.size(); ++i){
    for (size_t j = 0; j<countlist.size(); ++j){
      if (i==j) continue;
      if (std::abs(countlist[i]->GetTrigTime()-countlist[j]->GetTrigTime()) < 2){
        // for c2: GetTrigID() is an unsigned int and always >= 0
        //if ( (    countlist[i]->GetTrigID() >= 6  && countlist[i]->GetTrigID() <= 15 && countlist[j]->GetTrigID() >= 28 && countlist[j]->GetTrigID() <= 37 ) // East Lower, West Upper 
         //    || ( countlist[i]->GetTrigID() >= 0  && countlist[i]->GetTrigID() <= 5  && countlist[j]->GetTrigID() >= 22 && countlist[j]->GetTrigID() <= 27 ) // South Lower, North Upper
         //    || ( countlist[i]->GetTrigID() >= 16 && countlist[i]->GetTrigID() <= 21 && countlist[j]->GetTrigID() >= 38 && countlist[j]->GetTrigID() <= 43 ) // North Lower, South Upper
        if ( (    countlist[i]->GetTrigID() >= 6  && countlist[i]->GetTrigID() <= 15 && countlist[j]->GetTrigID() >= 28 && countlist[j]->GetTrigID() <= 37 ) // East Lower, West Upper 
             || ( countlist[i]->GetTrigID() <= 5  && countlist[j]->GetTrigID() >= 22 && countlist[j]->GetTrigID() <= 27 ) // South Lower, North Upper
             || ( countlist[i]->GetTrigID() >= 16 && countlist[i]->GetTrigID() <= 21 && countlist[j]->GetTrigID() >= 38 && countlist[j]->GetTrigID() <= 43 ) // North Lower, South Upper
             ) {
          //              std::cout << "I have a match..."
          //                        << "i " << i << ", ID " << countlist[i]->GetTrigID() << ", time " << countlist[i]->GetTrigTime() << "..."
          //                        << "j " << j << ", ID " << countlist[j]->GetTrigID() << ", Time " <<  countlist[j]->GetTrigTime() 
          //                        << std::endl;
          c1.push_back(countlist[i]->GetTrigID());
          c2.push_back(countlist[j]->GetTrigID());
          ci1.push_back(i);
          ci2.push_back(j);
        }
      }
    }
  }
  std::cout<<"counter array size = "<<c1.size()<<" tracklist size = "<<tracklist.size()<<std::endl;
  if (c1.size()&&tracklist.size()){
    double ccosx = (cx[c1[0]]-cx[c2[0]])/sqrt(pow(cx[c1[0]]-cx[c2[0]],2)+pow(cy[c1[0]]-cy[c2[0]],2)+pow(cz[c1[0]]-cz[c2[0]],2));
    double ccosy = (cy[c1[0]]-cy[c2[0]])/sqrt(pow(cx[c1[0]]-cx[c2[0]],2)+pow(cy[c1[0]]-cy[c2[0]],2)+pow(cz[c1[0]]-cz[c2[0]],2));
    double ccosz = (cz[c1[0]]-cz[c2[0]])/sqrt(pow(cx[c1[0]]-cx[c2[0]],2)+pow(cy[c1[0]]-cy[c2[0]],2)+pow(cz[c1[0]]-cz[c2[0]],2));   

    recob::Track::Vector_t larStart;
    recob::Track::Vector_t larEnd;
 
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);
    for (size_t i = 0; i<tracklist.size(); ++i){
      recob::Track::Point_t trackStart, trackEnd;
      std::tie(trackStart, trackEnd) = tracklist[i]->Extent(); 
      larStart = tracklist[i]->VertexDirection();
      larEnd = tracklist[i]->EndDirection();

      double dc = std::abs(ccosx*larStart.X()+ccosy*larStart.Y()+ccosz*larStart.Z());
      dcos->Fill(dc);
      /*
      if (i==4){
        for (size_t j = 0; j<tracklist[i]->NPoints(); ++j){
          if (tracklist[i]->HasValidPoint(j)){
            std::cout<<j<<" "<<tracklist[i]->LocationAtPoint(j).X()
                     <<" "<<tracklist[i]->LocationAtPoint(j).Y()
                     <<" "<<tracklist[i]->LocationAtPoint(j).Z()
                     <<" "<<tracklist[i]->DirectionAtPoint(j).X()
                     <<" "<<tracklist[i]->DirectionAtPoint(j).Y()
                     <<" "<<tracklist[i]->DirectionAtPoint(j).Z()
                     <<std::endl;
          }
        }
      }
      */
      if (dc>0.98){//found a match
        double ctime = (countlist[ci1[0]]->GetTrigTime() + countlist[ci2[0]]->GetTrigTime())/(2*32.); //ticks
        double x0 = trackStart.X() + (trackStart.X()>0?-1:1)*ctime*0.5*0.1085;
        //double y0 = trackStart.Y();
        double z0 = trackStart.Z();
        double x1 = trackEnd.X() + (trackStart.X()>0?-1:1)*ctime*0.5*0.1085;
        //double y1 = trackEnd.Y();
        double z1 = trackEnd.Z();

        double dx1 = 1e10;
        double dx2 = 1e10;
        if (c1[0]>=6&&c1[0]<=15){
          dx1 = std::abs(x0+(cz[c1[0]]-z0)*(x0-x1)/(z0-z1)-cx[c1[0]]);
        }
        if (c2[0]>=28&&c2[0]<=37){
          dx2 = std::abs(x0+(cz[c2[0]]-z0)*(x0-x1)/(z0-z1)-cx[c2[0]]);
        }
        //std::cout<<dx1<<" "<<dx2<<std::endl;
        hdx->Fill(x0+(cz[c1[0]]-z0)*(x0-x1)/(z0-z1)-cx[c1[0]]);
        hdx->Fill(x0+(cz[c2[0]]-z0)*(x0-x1)/(z0-z1)-cx[c2[0]]);
        if (dx1>40||dx2>40) continue;
        for (size_t j = 0; j<tracklist[i]->NumberTrajectoryPoints(); ++j){
          auto loc = tracklist[i]->LocationAtPoint(j);
          hyz->Fill(loc.Z(),loc.Y());
          if (loc.X()>0){
            hxz->Fill(loc.Z(),loc.X()-ctime*0.5*0.1085);
          }
          else{
            hxz->Fill(loc.Z(),loc.X()+ctime*0.5*0.1085);
          }
        }
        //map wireid with hit pointer, index on track and x position
        std::map<size_t,art::Ptr<recob::Hit> > hitmap;
        std::map<size_t,int> indexmap;
        std::map<size_t,double> xposmap;
        if (fmthm.isValid()){//fmthm.isValid()
          auto vhit = fmthm.at(i);
          auto vmeta = fmthm.data(i);

          std::map<geo::PlaneID, std::vector<std::pair<int, double>>> planehits;
        
          for (size_t h = 0; h < vhit.size(); ++h){
            auto loc = tracklist[i]->LocationAtPoint(vmeta[h]->Index());
            if (vhit[h]->WireID().TPC==5&&vhit[h]->WireID().Plane==2){
              hitmap[vhit[h]->WireID().Wire] = vhit[h];
              indexmap[vhit[h]->WireID().Wire] = vmeta[h]->Index();
              xposmap[vhit[h]->WireID().Wire] = loc.X()-ctime*0.5*0.1085;
            }
            planehits[vhit[h]->WireID().planeID()].push_back(std::make_pair(vhit[h]->WireID().Wire, vhit[h]->PeakTime()));
            //now study the difference between hit time and prediction from counters
            double xyzStart[3], xyzEnd[3]; //wire ends
            wireReadout.WireEndPoints(vhit[h]->WireID(), xyzStart, xyzEnd);
            double xhit, yhit, zhit; //intersection between wire and counter line
            //std::cout<<"about to calculate intersection"<<std::endl;
            //std::cout<<xyzStart[1]<<" "<<xyzStart[2]<<" "<<xyzEnd[1]<<" "<<xyzEnd[2]<<" "<<cy[c1[0]]<<" "<<cz[c1[0]]<<" "<<cy[c2[0]]<<" "<<cz[c2[0]]<<std::endl;
            if (geo::IntersectLines(xyzStart[1], xyzStart[2],
                                     xyzEnd[1], xyzEnd[2],
                                     cy[c1[0]], cz[c1[0]],
                                     cy[c2[0]], cz[c2[0]],
                                     yhit, zhit)){
              xhit = cx[c1[0]] + (cx[c2[0]]-cx[c1[0]])/(cy[c2[0]]-cy[c1[0]])*(yhit-cy[c1[0]]);
              hdt->Fill(vhit[h]->PeakTime() - detProp.ConvertXToTicks(xhit+(vhit[h]->WireID().TPC%2==1?1:-1)*ctime*0.5*0.1085, vhit[h]->WireID()));
            }
          }//loop over all associated hits

          //make sure there are enough hits on each plane before proceeding
          for (auto &plhits : planehits){//loop over planes
            int w0 = INT_MAX;
            int w1 = INT_MIN;
            std::vector<float> x, y, ey;
            for (auto &hit : plhits.second){
              if (hit.first < w0) w0 = hit.first;
              if (hit.first > w1) w1 = hit.first;
              x.push_back(hit.first);
              y.push_back(hit.second);
              ey.push_back(10); //arbitrary error on time
            }
            if (w0!=INT_MAX&&w1!=INT_MIN&&(w1-w0>5)&&plhits.second.size()>2){//there are at least 3 hits and the track is at least 5 wires long
              float intcpt, slope, intcpterr, slopeerr, chidof;
              fLinFitAlg.LinFit(x, y, ey, intcpt, slope, intcpterr, slopeerr, chidof);
              //Loop over all wires in the range identified
              std::map<raw::ChannelID_t, int> chused;
              for (int w = w0; w<=w1; ++w){
                geo::WireID wid(plhits.first, w);
                auto channel = wireReadout.PlaneWireToChannel(wid);
                if (fCSP->IsBad(channel)) continue;
                if (chused.find(channel)==chused.end()){ //channel not looked at before
                  double tick = intcpt + slope*w;
                  double x = detProp.ConvertTicksToX(tick, wid);//raw x
                  x += (wid.TPC%2==1?-1:1)*ctime*0.5*0.1085; //correct for t0
                  double xyzStart[3], xyzEnd[3]; //wire ends
                  wireReadout.WireEndPoints(wid, xyzStart, xyzEnd);
                  x = std::abs(x-xyzStart[0]); //distance to wire plane
                  if (x<50){//hit is within 50 cm of wire plane
                    int besttime = -1;
                    for (size_t j = 0; j<rawlist.size(); ++j){
                      if (rawlist[j]->Channel() == channel){
                        double pedestal = rawlist[j]->GetPedestal();
                        int maxph = -1;
                        for (size_t k = 0; k<rawlist[j]->NADC(); ++k){
                          if (float(k)>=tick&&float(k)<=tick+20){
                            if (int(rawlist[j]->ADC(k)-pedestal)>maxph){
                              maxph = int(rawlist[j]->ADC(k)-pedestal);
                              besttime = k;
                            }
                          }
                        }
                        double mean=0;
                        double mean2=0;
                        int npts = 0;
                        for (size_t k = 0; k<rawlist[j]->NADC(); ++k){
                          if ((int(k)>besttime-200&&int(k)<besttime-100)||
                              (int(k)>besttime+100&&int(k)<besttime+200)){
                            mean += rawlist[j]->ADC(k)-pedestal;
                            mean2 += pow(rawlist[j]->ADC(k)-pedestal,2);
                            ++npts;
                          }
                        }
                        mean/=npts;
                        mean2/=npts;
                        double angleToVert = wireReadout.WireAngleToVertical(wireReadout.Plane(wid).View(), wid.asPlaneID().asTPCID()) - 0.5*::util::pi<>();
                        //std::cout<<vhit[h]->View()<<" "<<vhit[h]->WireID().TPC<<" "<<vhit[h]->WireID().Cryostat<<" "<<angleToVert<<std::endl;
                        const auto& dir = tracklist[i]->DirectionAtPoint(0);
                        double cosgamma = std::abs(std::sin(angleToVert)*dir.Y() + std::cos(angleToVert)*dir.Z());
                        //std::cout<<maxph*cosgamma<<" "<<sqrt(mean2-mean*mean)<<std::endl;
                        //std::cout<<vhit[h]->PeakTime()<<" "<<vhit[h]->PeakAmplitude()<<" "<<besttime<<" "<<maxph<<std::endl;
                        hsig[wid.TPC][wid.Plane][0]->Fill(maxph*cosgamma);
                        hbkg[wid.TPC][wid.Plane][0]->Fill(sqrt(mean2-mean*mean));
                      }
                    }//loop over rawlist
                    if (besttime<0||besttime>=15000) continue;
                    for (size_t j = 0; j<wirelist.size(); ++j){
                      if (wirelist[j]->Channel() == channel){
                        double maxph1 = wirelist[j]->Signal()[besttime];
                        double mean=0;
                        double mean2=0;
                        int npts = 0;
                        for (size_t k = 0; k<wirelist[j]->NSignal(); ++k){
                          if ((int(k)>besttime-200&&int(k)<besttime-100)||
                              (int(k)>besttime+100&&int(k)<besttime+200)){
                            mean += wirelist[j]->Signal()[k];
                            mean2 += pow(wirelist[j]->Signal()[k],2);
                            ++npts;
                          }
                        }
                        mean/=npts;
                        mean2/=npts;
                        double angleToVert = wireReadout.WireAngleToVertical(wireReadout.Plane(wid).View(), wid.asPlaneID().asTPCID()) - 0.5*::util::pi<>();
                        //std::cout<<vhit[h]->View()<<" "<<vhit[h]->WireID().TPC<<" "<<vhit[h]->WireID().Cryostat<<" "<<angleToVert<<std::endl;
                        const auto& dir = tracklist[i]->DirectionAtPoint(0);
                        double cosgamma = std::abs(std::sin(angleToVert)*dir.Y() + std::cos(angleToVert)*dir.Z());
                        //std::cout<<maxph*cosgamma<<" "<<sqrt(mean2-mean*mean)<<std::endl;
                        //std::cout<<vhit[h]->PeakTime()<<" "<<vhit[h]->PeakAmplitude()<<" "<<besttime<<" "<<maxph<<std::endl;
                        hsig[wid.TPC][wid.Plane][1]->Fill(maxph1*cosgamma);
                        hbkg[wid.TPC][wid.Plane][1]->Fill(sqrt(mean2-mean*mean));
                      }
                    }//loop over wirelist
                  }//within 50 cm of anode
                  chused[channel] = 1;
                }//channel not used 
              }//loop over wires in range
            }//find enough wires
          }//loop over planes
        }//fmthm is valid
        if (hitmap.size()>=10){//found at least 10 hits in TPC5
          for (auto const& wireID : wireReadout.Iterate<geo::WireID>(geo::PlaneID{0, 5, 2})) {//plane 2, tpc 5, cstat 0
            std::size_t const h = wireID.Wire;
            if (fCSP->IsBad(wireReadout.PlaneWireToChannel(wireID))) continue;
            double pt = -1;
            double xpos = -1;
            if (hitmap.find(h)!=hitmap.end()){//found hit on track
              pt = hitmap[h]->PeakTime();
              xpos = xposmap[h];
            }
            else{//found hit time through extrapolation
              std::vector<double> vw;
              std::vector<double> vt;
              std::vector<double> vx;
              for (auto& hv : hitmap){
                // for c2: fix ambiguous call to abs
                //if (std::abs(hv.first-h)<=5){
                if (lar::util::absDiff(hv.first,h)<=5){
                  vw.push_back((hv.second)->WireID().Wire);
                  vt.push_back((hv.second)->PeakTime());
                  vx.push_back(xposmap[hv.first]);
                }
              }
              if (vw.size()>=3){
                TGraph *gr = new TGraph(vw.size(), &vw[0], &vt[0]);
                pt = gr->Eval(h);
                delete gr;
                gr = new TGraph(vw.size(), &vw[0], &vx[0]);
                xpos = gr->Eval(h);
                delete gr;
              }
            }
            if (pt>=0){
              double maxph = -1e10;
              for (size_t j = 0; j<rawlist.size(); ++j){
                if (rawlist[j]->Channel() == wireReadout.PlaneWireToChannel(wireID)){
                  double pedestal = rawlist[j]->GetPedestal();
                  for (size_t k = 0; k<rawlist[j]->NADC(); ++k){
                    if (float(k)>=pt&&float(k)<=pt+20){
                      if (rawlist[j]->ADC(k)-pedestal>maxph){
                        maxph = rawlist[j]->ADC(k)-pedestal;
                      }
                    }
                  }
                }
              }
              double angleToVert = wireReadout.WireAngleToVertical(hitmap.begin()->second->View(), hitmap.begin()->second->WireID().asPlaneID().asTPCID()) - 0.5*::util::pi<>();
              //std::cout<<vhit[h]->View()<<" "<<vhit[h]->WireID().TPC<<" "<<vhit[h]->WireID().Cryostat<<" "<<angleToVert<<std::endl;
              //const TVector3& dir = tracklist[i]->DirectionAtPoint(indexmap.begin()->second);
              const auto& dir = tracklist[i]->VertexDirection();
              double cosgamma = std::abs(std::sin(angleToVert)*dir.Y() + std::cos(angleToVert)*dir.Z());
              hphx->Fill(xpos,maxph*cosgamma/wireReadout.Plane({0, 0, hitmap.begin()->second->View()}).WirePitch());
              tpc = 5;
              wire = h;
              trackid = i;
              x = xpos;
              dqdx = maxph*cosgamma/wireReadout.Plane({0, 0, hitmap.begin()->second->View()}).WirePitch();
              ftree->Fill();
              //std::cout<<h<<" "<<pt<<" "<<xpos<<" "<<maxph*cosgamma/wireReadout.WirePitch(hitmap.begin()->second->View())<<std::endl;
            }
          }
        }
      }//matched track
    }
  }
}

void dune::SignalToNoise::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  
  dcos = tfs->make<TH1D>("dcos",";cos#theta;",100,0,1.1);
  dcos->Sumw2();

  hxz = tfs->make<TH2D>("hxz",";z (cm);x (cm)",1000,0,160,1000,-50,250);
  hyz = tfs->make<TH2D>("hyz",";z (cm);y (cm)",1000,0,160,1000,-90,120);

  for (int i = 0; i<8; ++i){
    for (int j = 0; j<3; ++j){
      for (int k = 0; k<4; ++k){
        hsig[i][j][k]= tfs->make<TH1D>(Form("hsig_%d_%d_%d",i,j,k),Form("TPC=%d, Plane=%d, %d",i,j,k), 200,0,200);
        hbkg[i][j][k]= tfs->make<TH1D>(Form("hbkg_%d_%d_%d",i,j,k),Form("TPC=%d, Plane=%d, %d",i,j,k), 100,0,100);
        hstb[i][j][k]= tfs->make<TH1D>(Form("hstb_%d_%d_%d",i,j,k),Form("TPC=%d, Plane=%d, %d",i,j,k), 100,0,100);
      }
    }
  }
        
  hdx = tfs->make<TH1D>("hdx",";#Delta x (cm)",100,-100,100);
  hdx->Sumw2();
  hphx = tfs->make<TH2D>("hphx",";x (cm); dQ/dx (ADC/cm)",50,0,250,100,0,600);

  hdt = tfs->make<TH1D>("hdt",";#Delta t (ticks);",1000,-1000,1000);

  ftree = tfs->make<TTree>("tree","a tree for sn analysis");
  ftree->Branch("run", &run, "run/I");
  ftree->Branch("event", &event, "event/I");
  ftree->Branch("tpc", &tpc, "tpc/I");
  ftree->Branch("wire", &wire, "wire/I");
  ftree->Branch("trackid", &trackid, "trackid/I");
  ftree->Branch("dqdx", &dqdx, "dqdx/D");
  ftree->Branch("x", &x, "x/D");

  if (fIsData){
    std::ifstream in;
    //in.open("/dune/app/users/mthiesse/olddev/CounterZOffset/work/counterInformation.txt");
    in.open("counterInformation.txt");
    char line[1024];
    while(1){
      in.getline(line,1024);
      if (!in.good()) break;
      int i;
      float x,y,z;
      sscanf(line,"%d %f %f %f",&i,&x,&y,&z);
      //std::cout<<i<<" "<<x<<" "<<y<<" "<<z<<std::endl;
      cx[i] = x;
      cy[i] = y;
      cz[i] = z;
    }
    in.close();
  }
  else{
    auto const& auxDetGeom = art::ServiceHandle<geo::AuxDetGeometry>()->GetProvider();
    for (size_t i = 0; i<auxDetGeom.NAuxDets(); ++i){
      auto& auxdet = auxDetGeom.AuxDet(i);
      //std::cout<<i<<" "<<auxdet.GetCenter().X()<<" "<<auxdet.GetCenter().Y()<<" "<<auxdet.GetCenter().Z()<<std::endl;
      cx[i] = auxdet.GetCenter().X();
      cy[i] = auxdet.GetCenter().Y();
      cz[i] = auxdet.GetCenter().Z();
    }
  }
}

DEFINE_ART_MODULE(dune::SignalToNoise)

/*
 * \file: CalibrationTreeBuilder.h
 * \author: Jason Stock (jason.stock@mines.sdsmt.edu
 * \brief: This is a small Tree made for use in calibration analysis.
 *
 */

#ifndef CALIBRATIONTREEBUILDER_H
#define CALIBRATIONTREEBUILDER_H

//Includes
#include "dunecore/DuneObj/CalibTreeRecord.h"
#include "dunecore/DuneObj/OpDetDivRec.h"

//LArSoft Includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"
namespace detinfo {
  class DetectorClocksData;
}

//FrameworkIncludes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//ROOT includes
#include "TTree.h"
#include "TFile.h"

//CPP includes
#include <vector>
#include <map>


namespace{}//

namespace CalibrationTreeBuilder {

  struct flatfiller{ //ordered for descending size.
    Double_t eve_x ;
    Double_t eve_y ;
    Double_t eve_z ;
    Double_t eve_t ;
    Double_t part_x ;
    Double_t part_y ;
    Double_t part_z ;
    Double_t part_t ;
    Double_t hit_charge;
    Double_t hit_energy;
    Double_t hit_time ;
    Double_t hit_width ;
    Double_t hit_split ;
    Double_t ophit_pes ;
    Double_t ophit_energy;
    Double_t ophit_time ;
    Double_t ophit_width ;
    Double_t ophit_split ;
    int64_t hit_index ;
    int64_t ophit_index ;
    UInt_t run ;
    UInt_t subrun ;
    UInt_t event_n ;
    UInt_t hit_wire ;
    UInt_t ophit_opchan ;
    UInt_t eve_index;
    UInt_t part_index;
    Int_t eve_trackid;
    Int_t eve_pdgid ;
    Int_t part_trackid ;
    Int_t part_pdgid;
    bool part_iseve ;
  };


  class CalibrationTreeBuilder : public art::EDAnalyzer
  {
    public:

      //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      CalibrationTreeBuilder(fhicl::ParameterSet const& pSet);

      virtual void beginJob() override;
      virtual void analyze(const art::Event& evt) override;
      virtual void endJob() override;


    private:
      //  fhicl::Atom<art::InputTag> WavLabel{fhicl::Name("WavLabel"), fhicl::Comment("The default label for the module to use when DivRecs "), "opdigi"};

      bool AddHit(detinfo::DetectorClocksData const& clockData, const art::Ptr<recob::Hit> hit, unsigned int& counter);
      bool AddHit(const art::Ptr<recob::OpHit> hit, unsigned int& counter);
      void PrepDivRec(const art::Event& evt);
      std::pair<std::vector<CalibTreeRecord::EveRecord>::iterator, bool> EmplaceEve(const simb::MCParticle* new_eve);
      std::pair<std::vector<CalibTreeRecord::ParticleRecord>::iterator, bool> EmplaceParticle(const simb::MCParticle* new_part);

      art::ServiceHandle<cheat::ParticleInventoryService> PIS;
      art::ServiceHandle<cheat::BackTrackerService> BTS;
      art::ServiceHandle<cheat::PhotonBackTrackerService> PBS;
      art::ServiceHandle<geo::Geometry> GS;
      geo::WireReadoutGeom const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();

      art::InputTag private_HitLabel;
      art::InputTag private_OpHitLabel;
      art::InputTag fWavLabel;

      art::ServiceHandle<art::TFileService> private_service_tfs;

      TTree* private_CalibrationTree;
      TTree* private_FlatCalibrationTree;
      TBranch* private_CalibrationRecord;
      TBranch* private_FlatCalibrationRecord;

      const art::Ptr< sim::OpDetDivRec > FindDivRec(int const& opDetNum) const;
      const std::vector< sim::SDP > OpHitToChannelWeightedSimSDPs(art::Ptr<recob::OpHit> const& opHit_P) const;

      CalibTreeRecord::CalibTreeRecord private_eventBuffer;
      CalibTreeRecord::CalibTreeRecord private_eventPrep;
      flatfiller fl;

      mutable std::vector<art::Ptr<sim::OpDetDivRec>> priv_DivRecs;
      //std::map<UInt_t, sim::OpDetDivRec> priv_od_to_chanDiv;


      //Also need buffers

  };// end class CalibrationTreeBuilder

}

#endif//endif CALIBRATIONTREEBUILDER_H

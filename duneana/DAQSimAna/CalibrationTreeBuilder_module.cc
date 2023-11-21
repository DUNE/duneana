/*
 * \file: CalibrationTreeBuilder_module.cc
 * \author: JStock (jason.stock@mines.sdsmt.edu)
 * \brief: This is a small analysis Tree made for use in the Calibration group to investigate radiological backgrounds and calibration sources.
 *
 * 2018_May: JStock
 * The Overload of AddHit for ophits is now diverging from the AddHit(Hit) function, as it will also fill the PerOpDet det records as well.
 *
 */

#include "CalibrationTreeBuilder.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include <vector>
#include <TObject.h>

namespace{}


namespace CalibrationTreeBuilder {

  struct genFinder{
    private:
      typedef std::pair<int, std::string> track_id_to_string;
      std::vector<track_id_to_string> track_id_map;
      std::set<std::string> generator_names;
      bool isSorted=false;
//      static bool pairsort(const std::pair<int, std::string>& a, const std::pair<int, std::string>& b){
//        return (a.first<b.first);
//      }

    public:
      void sort_now(){
        std::sort(this->track_id_map.begin(), this->track_id_map.end(), [](const auto &a, const auto &b){return (a.first < b.first) ; } );
//        std::sort(this->track_id_map.begin(), this->track_id_map.end(), pairsort);
        isSorted=true;
      }
      void add(const int& track_id, const std::string& gname){
        this->track_id_map.push_back(std::make_pair(track_id, gname));
        generator_names.emplace(gname);
        isSorted=false;
      }
      bool has_gen(std::string gname){
        return static_cast<bool>(generator_names.count(gname));
      };
      std::string get_gen(int tid){
        if( !isSorted ){
          this->sort_now();
        }
        return std::lower_bound(track_id_map.begin(), track_id_map.end(), tid,[](const auto &a, const auto &b){return (a.first < b) ; } )->second;
      };

  };
  genFinder* gf = new genFinder();

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  CalibrationTreeBuilder::CalibrationTreeBuilder(fhicl::ParameterSet const& pSet)
    :EDAnalyzer(pSet),
    private_HitLabel(pSet.get<art::InputTag>("HitLabel","gaushit")),
    private_OpHitLabel(pSet.get<art::InputTag>("OpHitLabel","ophit")),
    fWavLabel(pSet.get<art::InputTag>("WavLabel", "opdigi"))
  {  }

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  /*
     CalibrationTreeBuilder::CalibrationTree(fhiclConfig const& config)
     :EDAnalyzer(config),
     private_HitLabel(config.HitLabel()),
     private_OpHitLabel(config.OpHitLabel())
     {  }
     */ //Commented out until I have a subtable for EDAnalyzer

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  void CalibrationTreeBuilder::beginJob(){
    private_CalibrationTree = private_service_tfs->make<TTree>("CalibrationRecordTree"," A TTree for storing hit and ophit values under their particles for each event in the simulation");
    private_FlatCalibrationTree = private_service_tfs->make<TTree>("FlatCalibrationRecordTree"," A TTree for storing hit and ophit values with their particles for each event in the simulation");
    private_CalibrationRecord = private_CalibrationTree->Branch("event_records", &private_eventBuffer );
    private_FlatCalibrationRecord = private_FlatCalibrationTree->Branch("flat_event_records", &fl, "eve_x/D:eve_y:eve_z:eve_t:part_x:part_y:part_z:part_t:hit_charge:hit_energy:hit_time:hit_width:hit_split:ophit_pes:ophit_energy:ophit_time:ophit_width:ophit_split:hit_index/L:ophit_index:run/i:subrun:event_n:hit_wire:ophit_opchan:eve_index:part_index:eve_trackid/I:eve_pdgid:part_trackid:part_pdgid:part_iseve/O");

  }

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  void CalibrationTreeBuilder::analyze(const art::Event& evt) {
    //eprV.clear();
    private_eventPrep.Clear();
    private_eventBuffer.Clear();
    priv_DivRecs.clear();
    this->PrepDivRec(evt);


    //Set up BackTracking //Now in header
    /*    art::ServiceHandle<cheat::ParticleInventoryService> PIS;
          art::ServiceHandle<cheat::BackTrackerService> BTS;
          art::ServiceHandle<cheat::PhotonBackTrackerService> PBS;*/
    //Get a list of generator names.
    auto const mcHandles = evt.getMany<std::vector<simb::MCTruth>>();
    std::vector< std::pair<int, std::string>> track_id_to_label;

    for( auto const& mcHandle : mcHandles ){
      const std::string& sModuleLabel = mcHandle.provenance()->moduleLabel();
      art::FindManyP<simb::MCParticle> findMCParts(mcHandle, evt, "largeant");
      std::vector<art::Ptr<simb::MCParticle> > mcParts = findMCParts.at(0);
      for( const art::Ptr<simb::MCParticle> ptr : mcParts){
        int track_id = ptr->TrackId();
        gf->add(track_id, sModuleLabel);
      }
    }

    std::vector<art::Ptr<recob::Hit>> hitList;
    auto hitHandle = evt.getHandle<std::vector<recob::Hit>>(private_HitLabel);
    if ( hitHandle )
      art::fill_ptr_vector(hitList, hitHandle);
    std::vector<art::Ptr<recob::OpHit>> opHitList;
    auto opHitHandle = evt.getHandle<std::vector<recob::OpHit>>(private_OpHitLabel);
    if ( opHitHandle )
      art::fill_ptr_vector(opHitList, opHitHandle);

    //private_eventBuffer is the event record
    private_eventBuffer.bunch_hits=true;
    private_eventBuffer.run=evt.id().run();
    private_eventBuffer.subrun=evt.id().subRun();
    private_eventBuffer.event=evt.id().event();
    //private_eventPrep is the event record
    private_eventPrep.bunch_hits=true;
    private_eventPrep.run=evt.id().run();
    private_eventPrep.subrun=evt.id().subRun();
    private_eventPrep.event=evt.id().event();

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    unsigned int hitCounter=0;
    for ( auto const& hit : hitList){
      hitCounter++;
      if(!AddHit(clockData, hit, hitCounter)){
        mf::LogDebug(__FUNCTION__)<<"Hit "<<hitCounter<<" failed to add to the record. Is this hit caused by noise?";
      }
    }//end for hit in hitlist

    hitCounter=0;
    for ( auto const& hit : opHitList ){
      ++hitCounter;
      if(!AddHit(hit, hitCounter)){
        mf::LogDebug(__FUNCTION__)<<"OpHit "<<hitCounter<<" failed to add to the record. Is this hit caused by noise?";
      }
    } //Add ophits to hitmap.

    //Buffer normalized per hit as it is built.
    //Make the event record
    /*
       if(!private_eventPrep.eves.empty()){
       mf::LogDebug(__FUNCTION__)<<"EventBuffer for CalibrationTree contains "<<private_eventPrep.eves.size()<<" eves.\n";
       if(!private_eventPrep.eves.at(0).particles.empty()){
       mf::LogDebug(__FUNCTION__)<<"The first particle has "<<private_eventPrep.eves.at(0).particles.size()<<" particles.\n";
       if(!private_eventPrep.eves.at(0).particles.at(0).partial_hits.empty()){
       mf::LogDebug(__FUNCTION__)<<"The first of these particles has "<<private_eventPrep.eves.at(0).particles.at(0).partial_hits.size()<<" partial hits.\n";
       }else{
       mf::LogDebug(__FUNCTION__)<<"The particle contains no hits.\n";
       }
       }else{
       mf::LogDebug(__FUNCTION__)<<"Eve contains no particles.\n";
       }
       }else{
       mf::LogDebug(__FUNCTION__)<<"EventBuffer for CalibrationTree contains no eves.\n";
       }
       */

    {//Forced local scope for a few variables. This is a good indication that the code I am writing doesn't belong here, but in another function.
      int eve_pos=0;
      for(auto eve : private_eventPrep.eves){
        if(!private_eventPrep.eves.empty()){ //This is not needed and is redundant.
          int part_pos=0;
          for(auto part : eve.particles){
            { //keep partial_hit_pos local scope for the loop
              int partial_hit_pos=0;
              for(auto partial_hit : part.partial_hits){
                int ind = partial_hit.index;
                std::vector<CalibTreeRecord::HitContributor>::iterator hc_pos = std::find_if(private_eventPrep.hits.begin(), private_eventPrep.hits.end(), [&ind](CalibTreeRecord::HitContributor& hc){ return ind==hc.index;});
                if(hc_pos==private_eventPrep.hits.end()){
                  private_eventPrep.hits.emplace_back(ind);
                  if(hc_pos!=(private_eventPrep.hits.end()-1)){
                    hc_pos=private_eventPrep.hits.end()-1;
                  }
                }
                hc_pos->locations.emplace_back(eve_pos,part_pos,partial_hit_pos);
                //                private_eventPrep.hits.emplace_back(eve_pos,part_pos,partial_hit_pos);
                ++partial_hit_pos;
              }
            }
            {
              int partial_ophit_pos=0;
              for(auto partial_ophit : part.partial_ophits){
                int ind = partial_ophit.index;
                std::vector<CalibTreeRecord::HitContributor>::iterator hc_pos = std::find_if(private_eventPrep.ophits.begin(), private_eventPrep.ophits.end(), [&ind](CalibTreeRecord::HitContributor& hc){ return ind==hc.index;});
                if(hc_pos==private_eventPrep.ophits.end()){
                  private_eventPrep.ophits.emplace_back(ind);
                  if(hc_pos!=(private_eventPrep.ophits.end()-1)){
                    hc_pos=private_eventPrep.ophits.end()-1;
                  }
                }
                hc_pos->locations.emplace_back(eve_pos,part_pos,partial_ophit_pos);
                ++partial_ophit_pos;
              }
            }
            ++part_pos;
          }
        }
        ++eve_pos;
      }
    }

    //This is where the partial hit bunching should happen.
    //I can use the hits records to bunch partials on the same particle. This however would mean I would then need to recalculate the positions for all hits post bunching.
    //Now we loop the private_eventPrep and build private_eventBuffer.

    if( private_eventPrep.bunch_hits==true){
      for(size_t eve_num=0; eve_num<private_eventPrep.eves.size(); ++eve_num ){
        auto eve = private_eventPrep.eves.at(eve_num);
        private_eventBuffer.eves.push_back(eve); //We are gettign them in order, we are filling them in order. This should not break anything.
        //private_eventBuffer.eves.at(eve_num).particles.clear();//clear the particles list from the copeid eve. We will now refill it.
        for(size_t part_num=0; part_num<eve.particles.size(); ++part_num ){
          auto part = eve.particles.at(part_num);
          auto& npart = private_eventBuffer.eves.at(eve_num).particles.at(part_num);
          npart.partial_hits.clear();
          npart.partial_ophits.clear(); //I did not clear the particles from the new list because I don't need to. They will be the same particles. The only thing that chages are the partial records, as we are bunching them.
          //add part to particle buffer.
          std::map<int64_t, CalibTreeRecord::PartialHit> bunched_hits;
          std::map<int64_t, std::list<Double_t>> remembered_energies;
          for(auto partialhit : part.partial_hits){ //No value in using an index loop here because I don't need the hit number (it will be internal, and won't match the counter);
            auto emp_h = bunched_hits.emplace(std::make_pair(partialhit.index,partialhit));
            std::list<Double_t> dummy;
            auto emp_e = remembered_energies.emplace(std::make_pair(partialhit.index, dummy));//very memory inefficient, but I am not in a position to care about that right now. One could easily improve on my code here.
            if(emp_h.second==true){ //emp_e.second and emp_h.second are always the same.
              emp_e.first->second.emplace_back(partialhit.energy);
              emp_h.first->second=partialhit;
            }else{
              if(emp_h.first->second.time != partialhit.time){throw cet::exception("CalibrationTreeBuilder")<<"Buncher failed. Trying to merge two of the same hit at different times?\n";}
              emp_e.first->second.emplace_back(partialhit.energy);
              emp_h.first->second.charge+=partialhit.charge;
              emp_h.first->second.split+=partialhit.split;
              emp_h.first->second.num_electrons+=partialhit.num_electrons;//numelectrons came from the ides that contributed to the given hit. They should stack because the IDEs stack(which is why we are adding them in the first place).
            }
            emp_e.first->second.unique();//Throw out duplicates.
            Double_t energy=0.0;
            for( auto en : emp_e.first->second){
              energy += en;
            }
            emp_h.first->second.energy=energy;
          }
          remembered_energies.clear();
          std::map<int64_t, CalibTreeRecord::PartialOpHit> bunched_ophits;
          for(auto partialop : part.partial_ophits){
            auto emp_o = bunched_ophits.emplace(std::make_pair(partialop.index,partialop));
            std::list<Double_t> dummy;
            auto emp_e = remembered_energies.emplace(std::make_pair(partialop.index, dummy));//very memory inefficient, but I am not in a position to care about that right now. One could easily improve on my code here.
            if(emp_o.second==true){
              emp_e.first->second.emplace_back(partialop.energy);
              emp_o.first->second=partialop;
            }else{
              if(emp_o.first->second.time != partialop.time){throw cet::exception("CalibrationTreeBuilder")<<"Buncher failed. Trying to merge two of the same ophit at different times?\n";}
              emp_e.first->second.emplace_back(partialop.energy);
              emp_o.first->second.pes+=partialop.pes;
              emp_o.first->second.split+=partialop.split;
              emp_o.first->second.num_photons+=partialop.num_photons;
              //numphotons came from the sdps that contributed to the given hit. They should stack because the SDPs stack(which is why we are adding them in the first place). NOTE!! We will not do this when we combine records together for the detectors.
            }
            emp_e.first->second.unique();//Throw out duplicates.
            Double_t energy=0.0;
            for( auto en : emp_e.first->second ){
              energy += en;
            }
            emp_o.first->second.energy=energy;
          }
          for(auto entry : bunched_hits){
            npart.partial_hits.push_back(entry.second);
          }
          for(auto entry : bunched_ophits){
            //This is where I should make the combined per detector records.
            //This new record should find any existing record of the same OpDetNum and Time, and add the PEs, and Photons to that record, if it doesn't exist, make a new one.
            npart.partial_ophits.push_back(entry.second);
          }
        }
      }

      //Reding the bunched record is a good time to make the flattened record
      /*
         mf::LogDebug(__FUNCTION__)<<"EventBuffer for CalibrationTree contains "<<private_eventBuffer.eves.size()<<" eves, the first of which has "<<private_eventBuffer.eves.at(0).particles.size()<<" particles. The first of these particles has "<<private_eventBuffer.eves.at(0).particles.at(0).partial_hits.size()<<" partial hits.\n";
         */
      {//Forced local scope for a few variables. This is a good indication that the code I am writing doesn't belong here, but in another function.
        int eve_pos=0;
        fl.run           = private_eventBuffer.run;
        fl.subrun        = private_eventBuffer.subrun;
        fl.event_n       = private_eventBuffer.event;
        for(auto eve : private_eventBuffer.eves){
          if(!private_eventBuffer.eves.empty()){
            fl.eve_x         = eve.x_pos;
            fl.eve_y         = eve.y_pos;
            fl.eve_z         = eve.z_pos;
            fl.eve_t         = eve.t_pos;
            fl.eve_trackid   = eve.trackId;
            fl.eve_pdgid     = eve.pdgid;
            fl.eve_index     = eve_pos;
            int part_pos=0;
            for(auto part : eve.particles){
              fl.part_x        = part.x_pos;
              fl.part_y        = part.y_pos;
              fl.part_z        = part.z_pos;
              fl.part_t        = part.t_pos;
              fl.part_trackid  = part.trackId;
              fl.part_pdgid    = part.pdgid;
              fl.part_iseve    = part.isEve;
              fl.part_index    = part_pos;
              { //keep partial_hit_pos local scope for the loop
                fl.ophit_pes     = 0.0;
                fl.ophit_energy  = 0.0;
                fl.ophit_time    = 0.0;
                fl.ophit_width   = 0.0;
                fl.ophit_split   = 0.0;
                fl.ophit_index   = 0.0;
                fl.ophit_opchan  = 0.0;
                int partial_hit_pos=0;
                for(auto partial_hit : part.partial_hits){
                  //This section is to fill the flat tree
                  fl.hit_charge    = partial_hit.charge;
                  fl.hit_energy    = partial_hit.energy;
                  fl.hit_time      = partial_hit.time;
                  fl.hit_width     = partial_hit.width;
                  fl.hit_split     = partial_hit.split;
                  fl.hit_index     = partial_hit.index;
                  fl.hit_wire      = partial_hit.wire;
                  private_FlatCalibrationTree->Fill();
                  //Flat tree finished.
                  int ind = partial_hit.index;
                  std::vector<CalibTreeRecord::HitContributor>::iterator hc_pos = std::find_if(private_eventBuffer.hits.begin(), private_eventBuffer.hits.end(), [&ind](CalibTreeRecord::HitContributor& hc){ return ind==hc.index;});
                  if(hc_pos==private_eventBuffer.hits.end()){
                    private_eventBuffer.hits.emplace_back(ind);
                    if(hc_pos!=(private_eventBuffer.hits.end()-1)){
                      hc_pos=private_eventBuffer.hits.end()-1;
                    }
                  }
                  hc_pos->locations.emplace_back(eve_pos,part_pos,partial_hit_pos);
                  //                private_eventBuffer.hits.emplace_back(eve_pos,part_pos,partial_hit_pos);
                  ++partial_hit_pos;
                }
              }
              fl.hit_charge    = 0.0;
              fl.hit_energy    = 0.0;
              fl.hit_time      = 0.0;
              fl.hit_width     = 0.0;
              fl.hit_split     = 0.0;
              fl.hit_index     = 0.0;
              fl.hit_wire      = 0.0;
              {
                int partial_ophit_pos=0;
                for(auto partial_ophit : part.partial_ophits){
                  int ind = partial_ophit.index;
                  fl.ophit_pes     = partial_ophit.pes;
                  fl.ophit_energy  = partial_ophit.energy;
                  fl.ophit_time    = partial_ophit.time;
                  fl.ophit_width   = partial_ophit.width;
                  fl.ophit_split   = partial_ophit.split;
                  fl.ophit_index   = partial_ophit.index;
                  fl.ophit_opchan  = partial_ophit.opchan;
                  private_FlatCalibrationTree->Fill();

                  //Seg fault here. I suspect it is because of an empty private_eventBuffer.ophits
                  std::vector<CalibTreeRecord::HitContributor>::iterator hc_pos = std::find_if(private_eventBuffer.ophits.begin(), private_eventBuffer.ophits.end(), [&ind](CalibTreeRecord::HitContributor& hc){ return ind==hc.index;});
                  if(hc_pos==private_eventBuffer.ophits.end()){
                    private_eventBuffer.ophits.emplace_back(ind);
                    if(hc_pos!=(private_eventBuffer.ophits.end()-1)){
                      hc_pos=private_eventBuffer.ophits.end()-1;
                    }
                  }
                  hc_pos->locations.emplace_back(eve_pos,part_pos,partial_ophit_pos);
                  ++partial_ophit_pos;
                }
              }
              ++part_pos;
            }
          }
          ++eve_pos;
        }
      }
      //Private event buffer eves are now done for this particle.
      //Loop private event buffer to make the hit and ophit partial positions list.
    }else{ private_eventBuffer=private_eventPrep;}


    private_CalibrationTree->Fill();



  } //end analyze

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  void CalibrationTreeBuilder::endJob(){
  }

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  //This function will break up the hit into partial hits
  //using the backtracker, and then will assigne the partial
  //hits to the appropriate particles in the record.
  //We will also add the appropriate references to the hit
  //list to be able to rebuild the "full" hit.
  bool CalibrationTreeBuilder::AddHit(detinfo::DetectorClocksData const& clockData,
                                      const art::Ptr<recob::Hit> hit, unsigned int& counter){
    art::ServiceHandle<geo::Geometry> geom;
    std::vector<const sim::IDE*> ides_ps = BTS->HitToSimIDEs_Ps(clockData, hit);
    std::vector<std::pair<int, CalibTreeRecord::PartialHit>> track_partials;
    Double_t total_charge = 0.0;
    if(ides_ps.empty()){return false;}
    for(auto ide_p : ides_ps){
      int track_id = ide_p->trackID;
      CalibTreeRecord::PartialHit tmp;
      total_charge += ide_p->numElectrons;
      //        tmp.charge = ide_p->numElectrons;
      tmp.num_electrons = ide_p->numElectrons;
      tmp.time   = hit->PeakTime();
      tmp.width  = hit->SigmaPeakTime();
      tmp.split  = 0.0;
      tmp.wire   = hit->Channel();
      tmp.is_collection_wire = (geom->SignalType(tmp.wire)==geo::kCollection)?(true):(false);
      tmp.energy = (tmp.is_collection_wire)?(ide_p->energy):(0.0);
      tmp.index  = counter;
      track_partials.push_back(std::make_pair(track_id, tmp));
    }
    for(auto track_partial : track_partials){
      track_partial.second.split  = track_partial.second.num_electrons / total_charge; //What percentage of the hit is in this partial.
      track_partial.second.charge = hit->Integral()*track_partial.second.split; //Charge now normalized by deposition;
      //Add partial to CalibTreeRecord;
      const simb::MCParticle* part = PIS->TrackIdToParticle_P(track_partial.first);
      //Emplace Particle. (Emplace Eve will be handeled in emplace particle
      std::pair<std::vector<CalibTreeRecord::ParticleRecord>::iterator, bool> part_marker = EmplaceParticle(part);
      part_marker.first->partial_hits.push_back(track_partial.second); //Add post normalized partial hit.
      //Don't forget to add the indicies to the hit vector. Actually. To keep it easy, we can just loop the eve records and add the indicies then. This is inefficient, but easier to code since I don't immediately have access to eve_pos here.
      mf::LogDebug(__FUNCTION__)<<"This partialhit is being added to the record with charge "<<track_partial.second.charge<<"\n";
    }//end for track partial.

    return true;
  }

  bool CalibrationTreeBuilder::AddHit(const art::Ptr<recob::OpHit> hit, unsigned int& counter){
    {//make each loop local so the vectors don't persist
//      try{
        const std::vector< sim::SDP> sdps = this->OpHitToChannelWeightedSimSDPs(hit);
//      }catch(const std::exception &e){
//        return false;
//      }
        std::vector<std::pair<int, CalibTreeRecord::PartialOpHit>> track_partials;
        Double_t total_charge = 0.0;
        if(sdps.empty()){return false;}
        for(auto sdp : sdps){
          int track_id = sdp.trackID;
          CalibTreeRecord::PartialOpHit tmp;
          total_charge += sdp.numPhotons;
          //        tmp.pes = //Not used yet.
          tmp.num_photons = sdp.numPhotons; //This is not PEs. This is incident photons. I am just storing it here for later logic (see final assignment in the next for loop.
          tmp.energy = sdp.energy;
          //tmp.energy = 0.0;// sdp.energy;
          tmp.time   = hit->PeakTime();
          tmp.width  = hit->Width();
          tmp.split  = 0.0;
          tmp.opchan = hit->OpChannel();
          tmp.opdet  = GS->OpDetFromOpChannel(hit->OpChannel());
          tmp.index  = counter;
          track_partials.push_back(std::make_pair(track_id, tmp));
        }
        for(auto track_partial : track_partials){
          track_partial.second.split  = track_partial.second.num_photons / total_charge; //What percentage of the hit is in this partial.
          track_partial.second.pes = hit->PE()*track_partial.second.split; //PEs now normalized by deposition;
          //Add partial to CalibTreeRecord;
          const simb::MCParticle* part = PIS->TrackIdToParticle_P(track_partial.first);
          //Emplace Particle. (Emplace Eve will be handeled in emplace particle
          std::pair<std::vector<CalibTreeRecord::ParticleRecord>::iterator, bool> part_marker = EmplaceParticle(part);
          part_marker.first->partial_ophits.push_back(track_partial.second); //Add post normalized partial hit.
          //Don't forget to add the indicies to the hit vector. Actually. To keep it easy, we can just loop the eve records and add the indicies then. This is inefficient, but easier to code since I don't immediately have access to eve_pos here.
        }//end for track partial.

      }
      return true;
//    }
//    catch(...){
//      return false;
//    }
  }



  std::pair<std::vector<CalibTreeRecord::ParticleRecord>::iterator, bool> CalibrationTreeBuilder::EmplaceParticle(const simb::MCParticle* part){
    int track_id = part->TrackId();
    int eve_id   = PIS->TrackIdToEveTrackId(part->TrackId());
    const simb::MCParticle* eve = PIS->TrackIdToParticle_P(eve_id);
    //This search can be put into EmplaceEve
    //std::vector<EvenRecord::ParticleRecord>::iterator = std::find_if(private_eventPrep.eves.begin(), private_eventPrep.eves.end(), [&eve_id](CalibTreeRecord::EveRecord& eve_rec){ return eve_id==eve_rec.trackId;});
    std::pair<std::vector<CalibTreeRecord::EveRecord>::iterator, bool> eve_marker = EmplaceEve(eve);
    std::vector<CalibTreeRecord::ParticleRecord>::iterator part_ptr = std::find_if(eve_marker.first->particles.begin(), eve_marker.first->particles.end(), [&track_id](CalibTreeRecord::ParticleRecord& part_rec){return track_id == part_rec.trackId;});
    if(part_ptr!=eve_marker.first->particles.end()){
      return std::make_pair(part_ptr, false);
    }else{
      CalibTreeRecord::ParticleRecord tmpRec;
      tmpRec.isEve   = (track_id==eve_id)?true:false; //The conditional operator here is unneccessary, but I think it is more clear to include it.
      tmpRec.trackId = track_id;
      tmpRec.pdgid   = part->PdgCode();
      tmpRec.dP      = ( part->EndMomentum().P() - part->Momentum(0).P() );
      tmpRec.dE      = ( part->EndMomentum().E() - part->Momentum(0).E() );
      mf::LogDebug(__FUNCTION__)<<"The particle now enetered into the particle record has PDGID "<<part->PdgCode()<<"\n";
      tmpRec.x_pos   = part->Position(0).X();
      tmpRec.y_pos   = part->Position(0).Y();
      tmpRec.z_pos   = part->Position(0).Z();
      tmpRec.t_pos   = part->Position(0).T();
      eve_marker.first->particles.push_back(tmpRec);
      if(part_ptr == (eve_marker.first->particles.end()-1) ){
        return std::make_pair(part_ptr, true); //If the vector hasn't been moved (can happen with expanding), just return what we have made).
      }else{//The vector has moved. Return the right answer (this could honestly be the only return, the previous being unnecessary).
        return std::make_pair(eve_marker.first->particles.end()-1, true);
      }
    }//end else (part_ptr!=eve_marker.first->particles.end());
  }//end EmplaceParticle

  std::pair<std::vector<CalibTreeRecord::EveRecord>::iterator, bool> CalibrationTreeBuilder::EmplaceEve(const simb::MCParticle* eve){
    std::vector<CalibTreeRecord::EveRecord>& eves = private_eventPrep.eves;
    int eve_id = eve->TrackId();
    std::vector<CalibTreeRecord::EveRecord>::iterator eve_ptr = std::find_if(eves.begin(), eves.end(), [&eve_id](CalibTreeRecord::EveRecord& eve_rec){return eve_id == eve_rec.trackId;});
    if( eve_ptr != eves.end() ){ return std::make_pair(eve_ptr, false); }else{
      CalibTreeRecord::EveRecord tmpRec;
      tmpRec.trackId = eve_id;
      tmpRec.generator = gf->get_gen(eve_id);
      tmpRec.pdgid = eve->PdgCode();
      mf::LogDebug(__FUNCTION__)<<"The eve now enetered into the particle record has PDGID "<<eve->PdgCode()<<"\n";
      tmpRec.x_pos   = eve->Position(0).X();
      tmpRec.y_pos   = eve->Position(0).Y();
      tmpRec.z_pos   = eve->Position(0).Z();      tmpRec.t_pos   = eve->Position(0).T();
      tmpRec.t_pos   = eve->Position(0).T();
      eves.push_back(tmpRec);
      if(eve_ptr == (eves.end()-1) ){
        return std::make_pair(eve_ptr, true);
      }else{
        return std::make_pair(eves.end()-1, true);
      }
    }
  }


  //DUNE CalibTree specific tools intended for PhotonBackTracker (DivRec must be a larsfot wide product for that to happen).

  //----------------------------------------------------------------
  const art::Ptr< sim::OpDetDivRec > CalibrationTreeBuilder::FindDivRec(int const& opDetNum) const
  {
    art::Ptr< sim::OpDetDivRec > opDet;
    for(size_t detnum = 0; detnum < priv_DivRecs.size(); ++detnum){
      if(priv_DivRecs.at(detnum)->OpDetNum() == opDetNum)
        opDet = priv_DivRecs.at(detnum);
    }
    if(!opDet)
    {
      throw cet::exception("CalibTreeBuilder") << "No sim:: OpDetDivRec corresponding "
        << "to opDetNum: " << opDetNum << "\n";
    }
    return opDet;
  }


  //----------------------------------------------------------------
  //This function weights each of the returned sdps by the correct fractional number of detected photons. This make nPEs and Energy in the SDP make sense. Because these are constructed in place as weighted copies of the actual SDPs, they do not exist in the principle, cannot be returned as pointers, and are not comparable between calls.
  //This function is largely copied from OpHitToSimSDPs_Ps
  const std::vector< sim::SDP > CalibrationTreeBuilder::OpHitToChannelWeightedSimSDPs(art::Ptr<recob::OpHit> const& opHit_P) const
  {
    const double fDelay = PBS->GetDelay();
    std::vector< sim::SDP > retVec;
    double fPeakTime = opHit_P->PeakTime();
    double fWidth = opHit_P->Width();
    UInt_t fChan = opHit_P->OpChannel();
    art::ServiceHandle<geo::Geometry> geom;
    int fDet  = geom->OpDetFromOpChannel(fChan);
    //I should use the timing service for these time conversions.
    sim::OpDetBacktrackerRecord::timePDclock_t start_time = ((fPeakTime- fWidth)*1000.0)-fDelay;
    sim::OpDetBacktrackerRecord::timePDclock_t end_time = ((fPeakTime+ fWidth)*1000.0)-fDelay;
    if(start_time > end_time){throw;}//This is bad. Give a reasonable error message here, and use cet::except

    //BUG!!!fGeom->OpDetFromOpChannel(channel)
    art::Ptr<sim::OpDetBacktrackerRecord> fBTR;
    try{
      fBTR = PBS->FindOpDetBTR(fDet);
    }catch(cet::exception &e){
      if(true){
        std::vector<sim::SDP> ret;
        return ret;
      }else{
        throw e;
      }
    }
    const std::vector<std::pair<double, std::vector<sim::SDP>> >& timeSDPMap
      = fBTR->timePDclockSDPsMap(); //Not guranteed to be sorted.
    art::Ptr<sim::OpDetDivRec> div_rec = this->FindDivRec(fDet);//This is an OpDetDivRec collected from this BTR.


    std::vector<const std::pair<double, std::vector<sim::SDP>>*> timePDclockSDPMap_SortedPointers;
    std::vector<const sim::OpDet_Time_Chans*> div_rec_SortedPointers;
    for ( auto& pair : timeSDPMap ){ timePDclockSDPMap_SortedPointers.push_back(&pair); }
    auto pairSort = [](auto& a, auto& b) { return a->first < b->first ; } ;
    if( !std::is_sorted( timePDclockSDPMap_SortedPointers.begin(), timePDclockSDPMap_SortedPointers.end(), pairSort)){
      std::sort(timePDclockSDPMap_SortedPointers.begin(), timePDclockSDPMap_SortedPointers.end(), pairSort);
    }

    //This section is a hack to make comparisons work right.
    std::vector<sim::SDP> dummyVec;
    std::pair<double, std::vector<sim::SDP>> start_timePair = std::make_pair(start_time, dummyVec);
    std::pair<double, std::vector<sim::SDP>> end_timePair = std::make_pair(end_time, dummyVec);
    auto start_timePair_P = &start_timePair;
    auto end_timePair_P = &end_timePair;

    //First interesting iterator.
    auto map_pdsdp_itr = std::lower_bound(timePDclockSDPMap_SortedPointers.begin(), timePDclockSDPMap_SortedPointers.end(), start_timePair_P, pairSort);
    //Last interesting iterator.
    auto map_pdsdp_last = std::upper_bound(map_pdsdp_itr, timePDclockSDPMap_SortedPointers.end(), end_timePair_P, pairSort);

    //retvec.push_back(map_pdsdp_first.second[0]);
    //This code screams fragile. Really. If you read this, feel free to clean up the methods.
    for(auto& sdp_time_pair =  map_pdsdp_itr; sdp_time_pair != map_pdsdp_last; ++sdp_time_pair){
      //cut div_rec by time
      auto time=(*sdp_time_pair)->first;
      for(auto& sdp : (*sdp_time_pair)->second){
        //      auto sdp = (*sdp_time_pair)->second; //This is a vector
        //auto slice = div_rec->GetSlice(start_time, end_time);
        auto time_divr_pair = div_rec->FindClosestTimeChan(time);
        auto time_divr=time_divr_pair.first;
        auto time_divr_found=time_divr_pair.second;
        sim::SDP tmp = sdp;
        if( time_divr_found && time_divr->time == time){ //This will break if time_divr is the end value of the vector
          tmp.energy = time_divr->GetFrac(fChan) * tmp.energy;
          // ((*map_divrec_itr)->DivChans.eScaleFrac(fChan))*tmp.energy;
        }else{
          //          std::cout<<"I am having trouble reconciling "<<time<<" and "<<time_divr->time<<"\n";
          tmp.energy = 0; //This does not corespond to an energy detected. A photon may or may not have been incident on  the detector, but it didn't get picked up.
        }
        retVec.emplace_back(tmp);

      }
    }

    //while( map_pdsdp_itr != map_pdsdp_last /*&& map_divrec_itr != map_divrec_last*/){ //Stepping through.
    /*for(auto const& sdp :  (*map_pdsdp_itr)->second){
      sim::SDP tmp = sdp;
      tmp.energy = ((*map_divrec_itr)->DivChans.eScaleFrac(fChan))*tmp.energy;
    //              second[fChan].tick_photons_frac)*tmp.energy;
    retVec.emplace_back(tmp);
    }
    ++map_pdsdp_itr;
    ++map_divrec_itr;
    }*/

    return retVec;
  }

void CalibrationTreeBuilder::PrepDivRec(const art::Event& evt)
{
  if( 0 ){ return;} //Insert check for DivRecs here, or don't use validHandle below.
  auto const& divrecHandle = evt.getValidHandle <std::vector<sim::OpDetDivRec>>(fWavLabel);
  if(divrecHandle.failedToGet()){
    return;
  }
  art::fill_ptr_vector(priv_DivRecs, divrecHandle);
  auto compareDivReclambda = [](art::Ptr<sim::OpDetDivRec> a, art::Ptr<sim::OpDetDivRec> b) {return(a->OpDetNum() < b->OpDetNum());};
  if (!std::is_sorted(priv_DivRecs.begin(), priv_DivRecs.end(), compareDivReclambda))
    std::sort(priv_DivRecs.begin(), priv_DivRecs.end(), compareDivReclambda);

}


}//end namespace

DEFINE_ART_MODULE(CalibrationTreeBuilder::CalibrationTreeBuilder)

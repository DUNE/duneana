#include "WireAna_module.h"
//Constructor for cnn struct


wireana::WireAna::WireAna(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset}  ,
  fWireProducerLabel(pset.get< art::InputTag >("InputWireProducerLabel", "caldata")),
  fSimChannelLabel(pset.get< art::InputTag >("SimChannelLabel", "elecDrift")),
  fSimulationProducerLabel(pset.get< art::InputTag >("SimulationProducerLabel", "largeant"))
{
  fLogLevel           = pset.get<int>("LogLevel", 10);
  fDoAssns            = pset.get<bool>("DoAssns", false);
  fNChanPerApa        = pset.get<int>("ChannelPerApa", 2560);
  fNTicksPerWire      = pset.get<int>("TickesPerWire", 6000);
  fMakeCluster        = pset.get<bool>("MakeCluster", true );
  auto const* geo = lar::providerFrom<geo::Geometry>();
  fNPlanes = geo->Nplanes();

  fMinPts  =  pset.get<unsigned int>("DBSCAN_MinPts", 2);
  fEps     =  pset.get<float>("DBSCAN_Eps", 4.5);
  fDrift   =  pset.get<float>("DBSCAN_Drift", 1.6);
  fPitch   =  pset.get<float>("DBSCAN_Pitch", 3.0);

  fDeltaMetric = pset.get<float>("CleanClusterDeltaScore", 0.01);

  image_channel_width  = pset.get<int>("IMAGE_CHANNEL_WIDTH",13); 
  image_tick_width     = pset.get<int>("IMAGE_TICK_WIDTH",400); 
  image_rebin_tick     = pset.get<int>("IMAGE_REBIN_TICK",4); 
  image_size           = image_channel_width*image_tick_width/image_rebin_tick;

  fDumpFileName      = pset.get<std::string>("DUMPFILENAME", "out.npy");
  fDumpMaxRow        = pset.get<int>("DUMPMAXROW", 50000);
  fDumpNClusters     = pset.get<int>("DUMPNCLUSTERS",-1);


  fHistEnergyMax         = pset.get<float>("Histo_EMax", 10);
  fHistChargeMax         = pset.get<float>("Histo_CMax", 10000);

  fTreeName          = pset.get<std::string>("TREENAME", "wireana");

  fViewMap[0]   =   "U";
  fViewMap[1]   =   "V";
  fViewMap[2]   =   "Z";

  //selTypes=std::vector<std::string>({"","_neutrino","_rad"});
  //partTypes=std::vector<std::string>({"","_electron","_proton","_neutron","_photon","_other"});
  if (fLogLevel >= 3 )
  {
    std::cout<<"Histogram Selection Types: "<<std::endl;
    for( auto st: selTypes ) std::cout<<selTypes.size()<<std::endl;
  }
}

void wireana::WireAna::analyze(art::Event const & evt) {
  //reset containers
  reset();
  //get detector property
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

  //get event data
  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();
  std::cout<<"########## EvtNo."<<event<<std::endl;

  MC = !evt.isRealData();
  /// FillTruthInfo to internal data objects
  /// i.e. trkid_to_label_map
  if (MC) FillTrackIDtoLabelMap( evt );

  ////////////////////////////////////////////////////////////
  //Build Wire SimChanel List
  //parse Wire/SimChannel Info
  //Sort wires/channels
  art::Handle<std::vector<recob::Wire>> wireListHandle;
  std::vector<art::Ptr<recob::Wire>> wirelist;
  if (evt.getByLabel(fWireProducerLabel, wireListHandle)) 
    art::fill_ptr_vector(wirelist, wireListHandle);
  if ( fLogLevel >= 3 )
    std::cout<<"Size of Wirelist: "<<wirelist.size()<<std::endl;
  SortWirePtrByChannel( wirelist, true );

  art::Handle<std::vector<sim::SimChannel>> simChannelListHandle;
  std::vector<art::Ptr<sim::SimChannel>> channellist;
  if (evt.getByLabel(fSimChannelLabel, simChannelListHandle)) 
    art::fill_ptr_vector(channellist, simChannelListHandle);
  SortWirePtrByChannel( channellist, true );


  // //Get channel-> wire,simchannel map
  if( fLogLevel >= 3 ) std::cout<<"Fill ch_w_sc"<<std::endl;
  for( auto w: wirelist ) 
    ch_w_sc[ w->Channel() ].first = w;
  for( auto w: channellist )
  {

    //accumulate energy and charge for all channels
    for( auto &tdcide: w->TDCIDEMap() )
    {
      std::vector<float> energies(partTypes.size(),fECMin );
      std::vector<float> charges(partTypes.size(),fECMin );
      std::vector<float> energiesNeut(partTypes.size(),fECMin );
      std::vector<float> chargesNeut(partTypes.size(),fECMin );
      std::vector<float> energiesRad(partTypes.size(),fECMin );
      std::vector<float> chargesRad(partTypes.size(),fECMin );

      for( auto &ide: tdcide.second )
      {
        
        if( fLogLevel >= 3 ) std::cout<<"ide.trackID: "<<ide.trackID<<std::endl;
        bool isSignal = trkid_to_label_map[ ide.trackID ] == "NuEScatter" || trkid_to_label_map[ ide.trackID ] == "marley";
        int pdg = PIS->TrackIdToParticle_P( ide.trackID )->PdgCode();
        float energy = ide.energy;
        float numElectrons = ide.numElectrons;
        energies[kAll]+=energy;
        charges[kAll]+=numElectrons;
        int partType = -1;
        if( abs(pdg) == 11 || abs(pdg) == 13 || abs(pdg) == 15 )
        {
          partType = kElectron;
        } else if( abs(pdg) == 2212)
        {
          partType = kProton;
        } else if(abs(pdg) == 2112)
        {
          partType = kNeutron;
        } else if(abs(pdg) == 22)
        {
          partType = kPhoton;
        } else
        {
          partType = kNuc;
        }
        //parsed particle, accumulate energy
        energies[partType]+=energy; charges[partType]+=numElectrons;
        if( isSignal )
        {
          energiesNeut[partType]+=energy; chargesNeut[partType]+=numElectrons;
        }
        else
        {
          energiesRad[partType]+=energy; chargesRad[partType]+=numElectrons;
        }
      }

      if( fLogLevel >= 3 ) std::cout<<"FillHistogram: begin"<<std::endl;
      FillHistogram( energies, charges, kSAll, false );
      FillHistogram( energiesNeut, chargesNeut, kSNeutrino, false );
      FillHistogram( energiesRad, chargesRad, kSRad, false );
      if( fLogLevel >= 3 ) std::cout<<"FillHistogram: end"<<std::endl;

      //TH2FMap["TrueEnergyChargeDeposited"]->Fill(energy,charge);
      //TH2FMap["TrueEnergyChargeDeposited_electron"]->Fill(energy_e,charge_e);
      //TH2FMap["TrueEnergyChargeDeposited_proton"]->Fill(energy_p,charge_p);
      //TH2FMap["TrueEnergyChargeDeposited_photon"]->Fill(energy_gamma,charge_gamma);
      //TH2FMap["TrueEnergyChargeDeposited_neutron"]->Fill(energy_n,charge_n);
      //TH2FMap["TrueEnergyChargeDeposited_other"]->Fill(energy_o,charge_o);
      TrueEnergyDeposited->Fill(energies[0]);
      TrueChargeDeposited->Fill(charges[0]);
    }

    if ( ch_w_sc.find( w->Channel() ) != ch_w_sc.end() ) 
    {
      ch_w_sc[ w->Channel() ].second = w;
      if( fLogLevel>=10 ) 
      {
        std::cout<<"Filled Channel "<<w->Channel()<<std::endl;
      }
    }
  }

  //We can now print truth particles

  //! Since ROI Tag Truth will supply the leading trkid, now need a method to convert trkid to MCPart and the truthlabel.
  //! The method to use seems to be PIS->TrackIdToMotherParticle(trkid
  //const simb::MCParticle *  TrackIdToParticle_P (int id) const
  //simb::MCParticle   TrackIdToParticle (int const id) const
  //const simb::MCParticle *  TrackIdToMotherParticle_P (int id) const
  //simb::MCParticle   TrackIdToMotherParticle (int const id) const
  //const art::Ptr< simb::MCTruth > &   TrackIdToMCTruth_P (int id) const
  //simb::MCTruth  TrackIdToMCTruth )
  

  /////////////////////////////////////////////////////////////
  // Build ROICluster
  //
  if( fLogLevel>=10 ) std::cout<<"BuildPlaneViewROIMap(): Start"<<std::endl;
  BuildPlaneViewROIMap( wirelist );
  if( fLogLevel>=10 ) std::cout<<"BuildPlaneViewROIMap(): End"<<std::endl;
  /////////////////////////////////////////////////////////////
  // Set Truth Info on roi
  if(MC)
  {
    truth_data.truth_intType = 1;
    TagAllROITruth( clockData );
    FillTruthBranches( evt, fTree, truth_data );
  }
  fTree->Fill();
  std::unique_ptr<std::vector<recob::Wire>> outwires(new std::vector<recob::Wire>);


  if(fMakeCluster)
  {
    if( fLogLevel>=10 ) std::cout<<"BuildInitialROIClusters(): Start"<<std::endl;
    BuildInitialROIClusters();
    if( fLogLevel>=10 ) std::cout<<"BuildInitialROIClusters(): Done"<<std::endl;
    if ( MC ) 
    {
      TagAllROIClusterTruth( clockData );
      if( fLogLevel>=10 ) std::cout<<"TagAllROIClusterTruth(): Done"<<std::endl;
    }



    //Match ROI across views
    ROIMatcher matcher;
    matcher.SetData(plane_view_roicluster_map);
    matcher.MatchROICluster();
    if( fLogLevel>=10 ) std::cout<<"matcher.MatchROICluster(): Done"<<std::endl;
    matcher.CleanDuplicates(fDeltaMetric);
    bool hasCluster = (matcher.GetMatchedClusters().size() != 0 );
    if( hasCluster )
    {
      matcher.SortROIClustersBySize();
      //Logs
      if( fLogLevel >= 2 )
      {
        std::cout<<"  List Matched Clusters: "<<std::endl;
        int i = 0;
        for( auto mc : matcher.GetMatchedClusters() )
        {
          if (i == fDumpNClusters) break;
          std::cout<<"    Item "<<i<<", PlaneID: "<<mc.planeid<<std::endl
                   <<"        metric: "<<mc.metric<<std::endl;
          for( int v=0;v<3;v++ )
          {
           std::cout<<Form("        View: %d, tick(%d, %d, %d), ch(%d, %d, %d)",
                        mc.clusters[v].view,
                        mc.clusters[v].begin_index,
                        mc.clusters[v].end_index,
                        mc.clusters[v].end_index-mc.clusters[v].begin_index,
                        mc.clusters[v].channel_min,
                        mc.clusters[v].channel_max, 
                        mc.clusters[v].channel_max-mc.clusters[v].channel_min
                        )<<std::endl;
          }
          if (MC)
          {
            std::cout<<"    Matched Truth::(GenCode, LabelMatch, GenTrkMatch): "<<std::endl;
            std::cout<<Form("                   (%d, %d, %d) ", mc.gencode(), mc.labelmatch(),mc.trkmatch())<<std::endl;
            std::cout<<Form("                   first label: %s",mc.clusters[0].label.c_str() )<<std::endl;
          }
          ++i;
        }
      }//end Logs

      int nDumpedClusters = 0;
      for( auto  mCluster : matcher.GetMatchedClusters() )
      {
        if (nDumpedClusters == fDumpNClusters ) break;

        //matchedroicluster mCluster = matcher.GetMatchedClusters().front();
        int ch_width=image_channel_width,tick_width= image_tick_width, nticks=image_rebin_tick;
        std::vector<double> u_vec = GetArrayFromWire( wirelist, mCluster.clusters[0], ch_width,tick_width);
        std::vector<double> v_vec = GetArrayFromWire( wirelist, mCluster.clusters[1], ch_width,tick_width);
        std::vector<double> z_vec = GetArrayFromWire( wirelist, mCluster.clusters[2], ch_width,tick_width);

        if( fLogLevel >= 3 ) std::cout<<"  Got all ArrayFromWire: "<<std::endl;
        std::vector<double> u_vecc = CombineTicks( u_vec, ch_width, nticks );
        std::vector<double> v_vecc = CombineTicks( v_vec, ch_width, nticks );
        std::vector<double> z_vecc = CombineTicks( z_vec, ch_width, nticks );
        if( fLogLevel >= 3 ) std::cout<<"  Combined all ArrayFromWire: "<<std::endl;
        std::vector<double> u_veccs = ScaleArray( u_vecc, 1, 225 );
        std::vector<double> v_veccs = ScaleArray( v_vecc, 1, 225 );
        std::vector<double> z_veccs = ScaleArray( z_vecc, 1, 225 );
        if( fLogLevel >= 3 ) std::cout<<"  Scaled all ArrayFromWire: "<<std::endl;
        //Logs
        if( fLogLevel >= 10 )
        {
          std::cout<<"  Print Array Values: "<<std::endl;
          //int t=0;
          for( unsigned int i = 0; i<z_vecc.size(); i++ )
          {
            //if( i % ch_width == 0 ) std::cout<<std::endl<<"line "<<t++<<": ";
            if( i % ch_width == 0 ) std::cout<<std::endl<<"[ ";
            std::cout<<Form("(%d,%d,%d), ",(int)u_veccs[i],(int)v_veccs[i],(int)z_veccs[i]);
            if( i % ch_width == (unsigned int) (ch_width-1) ) std::cout<<"],";
            //if (u_vec[i]==0) continue;
            //auto ct = CalculateCT(i,ch_width,tick_width);
            //std::cout<<Form("    Pts: (%d, %d, %.2f)", ct.first,ct.second,u_vec[i])<<std::endl;
          }
          std::cout<<std::endl<<"  == End Print Array Values: "<<std::endl<<
                                 Form("  are wires the same? uv, uz, vz: %d,%d,%d",
                                     (u_vec==v_vec),
                                     (u_vec==z_vec),
                                     (v_vec==z_vec) )<<std::endl;
        }//end Logs

        WriteNumPy( mCluster, wirelist );
        ++nDumpedClusters;
      }

    }
  }//end fMakeCluster


}

void wireana::WireAna::beginJob() {

  gROOT->SetBatch(1);

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>(fTreeName.c_str() ,fTreeName.c_str() );

  fTree->Branch("run", &run);
  fTree->Branch("subrun", &subrun);
  fTree->Branch("event", &event);
  fTree->Branch("MC", &MC);
  DeclareTruthBranches(fTree, truth_data);


  //setup histograms
  fECMin=-1e-8; //settting energy and charge accumulation to start at this value --> ensures true background, i.e 0 energy/0 charge falls into the underflow bin
  string name1 = Form("TrueEnergyChargeDeposited_%s",fTreeName.c_str() );
  string name2 = Form("TrueEnergyChargeDepositedInROI_%s",fTreeName.c_str() );
  name1 = "TrueEnergyChargeDeposited";
  name2 = "TrueEnergyChargeDepositedInROI";

  for( auto selType : selTypes )
  {
    for( auto partType : partTypes )
    {
      string n1 = "TrueEnergyChargeDeposited"+selType+partType, n2 = "TrueEnergyChargeDepositedInROI"+selType+partType;
      if( fLogLevel >= 3 )
      {
        std::cout<<n1<<", "<<n2<<std::endl;
      }
      TH2FMap[ n1 ] = tfs->make<TH2F>( n1.c_str(), "Energy vs Charge; E (MeV); N",100,0,fHistEnergyMax, 100,0,fHistChargeMax );
      TH2FMap[ n2 ] = tfs->make<TH2F>( n2.c_str(), "Energy vs Charge; E (MeV); N",100,0,fHistEnergyMax, 100,0,fHistChargeMax );
    }
  }

  name1 = Form("TrueEnergyDeposited_%s",fTreeName.c_str() );
  name2 = Form("TrueEnergyDepositedInROI_%s",fTreeName.c_str() );
  TrueEnergyDeposited = tfs->make<TH1F>( name1.c_str(), "Energy vs Charge; E (MeV)",100,0,fHistEnergyMax);
  TrueEnergyDepositedInROI =  tfs->make<TH1F>( name2.c_str(), "Energy vs Charge; E (MeV)",100,0,fHistEnergyMax);

  name1 = Form("TrueChargeDeposited_%s",fTreeName.c_str() );
  name2 = Form("TrueChargeDepositedInROI_%s",fTreeName.c_str() );
  TrueChargeDeposited = tfs->make<TH1F>( name1.c_str(), "Charge vs Charge; E (MeV)",100,0,fHistChargeMax);
  TrueChargeDepositedInROI =  tfs->make<TH1F>( name2.c_str(), "Charge vs Charge; E (MeV)",100,0,fHistChargeMax);


  if(fMakeCluster)
  {
    //setup npywriter
    c2numpy_init(&npywriter, fDumpFileName, fDumpMaxRow);
     

    c2numpy_addcolumn(&npywriter, "RUN", C2NUMPY_UINT16 );
    c2numpy_addcolumn(&npywriter, "SUBRUN", C2NUMPY_UINT16 );
    c2numpy_addcolumn(&npywriter, "EVENT", C2NUMPY_UINT16 );
    c2numpy_addcolumn(&npywriter, "MC", C2NUMPY_UINT16 );
    //setup image format
    c2numpy_addcolumn(&npywriter, "NChannels", C2NUMPY_UINT16 );
    c2numpy_addcolumn(&npywriter, "NTicks", C2NUMPY_UINT16 );
    c2numpy_addcolumn(&npywriter, "RTicks", C2NUMPY_UINT16 );
    // generator code
    c2numpy_addcolumn(&npywriter, "GeneratorCode", C2NUMPY_UINT16 );//0 NueScatter, 1 CC, 2 Radiological
    //check match code
    c2numpy_addcolumn(&npywriter, "LabelMatchCode", C2NUMPY_UINT16 );//views: 0 did not match, 1 matched
    c2numpy_addcolumn(&npywriter, "TrkMatchCode", C2NUMPY_UINT16 );//trkid: 0 did not match, 1 matched

    //setup truth information
    c2numpy_addcolumn(&npywriter, "NPhoton", C2NUMPY_UINT16 );
    c2numpy_addcolumn(&npywriter, "NProton", C2NUMPY_UINT16 );
    c2numpy_addcolumn(&npywriter, "NNeutron", C2NUMPY_UINT16 );
    c2numpy_addcolumn(&npywriter, "NMeson", C2NUMPY_UINT16 );
    c2numpy_addcolumn(&npywriter, "NNucleus", C2NUMPY_UINT16 );

    c2numpy_addcolumn(&npywriter, "Nu_Px", C2NUMPY_FLOAT );
    c2numpy_addcolumn(&npywriter, "Nu_Py", C2NUMPY_FLOAT );
    c2numpy_addcolumn(&npywriter, "Nu_Pz", C2NUMPY_FLOAT );
    c2numpy_addcolumn(&npywriter, "Part_Px", C2NUMPY_FLOAT );
    c2numpy_addcolumn(&npywriter, "Part_Py", C2NUMPY_FLOAT );
    c2numpy_addcolumn(&npywriter, "Part_Pz", C2NUMPY_FLOAT );


    //setup index data
    for( int i = 0; i < image_size; i++ )
    {
      for( int v = 0; v < 3; v++ )
      {
        std::string name=Form("%s_%08d", fViewMap[v].c_str(),i);
        c2numpy_addcolumn(&npywriter, name.c_str(), C2NUMPY_FLOAT);
      }
    }
  }// end fMakeCluster


  //GetTrackInNeutrinoTruth if MC

}

void wireana::WireAna::endJob()
{
  if(fMakeCluster) c2numpy_close(&npywriter);
}



void wireana::WireAna::reset()
{
  ResetTruthBranches(this->truth_data);
  ch_w_sc.clear();
  trkid_to_label_map.clear();
  plane_view_roi_map.clear();
  plane_view_roicluster_map.clear();
}



//======================================================================================
//======================================================================================
void wireana::WireAna::DeclareTruthBranches(TTree*t, DataBlock_Truth &blk)
{
  t->Branch("truth_intType",        &blk.truth_intType);         
  t->Branch("truth_nu_momentum_x",  &blk.truth_nu_momentum_x); 
  t->Branch("truth_nu_momentum_y",  &blk.truth_nu_momentum_y); 
  t->Branch("truth_nu_momentum_z",  &blk.truth_nu_momentum_z); 
  t->Branch("truth_nu_momentum_e",  &blk.truth_nu_momentum_e); 
  t->Branch("truth_nu_vtx_x",       &blk.truth_nu_vtx_x); 
  t->Branch("truth_nu_vtx_y",       &blk.truth_nu_vtx_y); 
  t->Branch("truth_nu_vtx_z",       &blk.truth_nu_vtx_z); 
  t->Branch("truth_nu_vtx_t",       &blk.truth_nu_vtx_t); 
  t->Branch("truth_nu_PDG",         &blk.truth_nu_PDG); 
  t->Branch("truth_lep_momentum_x", &blk.truth_lep_momentum_x); 
  t->Branch("truth_lep_momentum_y", &blk.truth_lep_momentum_y); 
  t->Branch("truth_lep_momentum_z", &blk.truth_lep_momentum_z); 
  t->Branch("truth_lep_momentum_e", &blk.truth_lep_momentum_e); 
  t->Branch("truth_lep_vtx_x",      &blk.truth_lep_vtx_x); 
  t->Branch("truth_lep_vtx_y",      &blk.truth_lep_vtx_y); 
  t->Branch("truth_lep_vtx_z",      &blk.truth_lep_vtx_z); 
  t->Branch("truth_lep_vtx_t",      &blk.truth_lep_vtx_t); 
  t->Branch("truth_lep_PDG",        &blk.truth_lep_PDG); 
  t->Branch("truth_had_momentum_x", &blk.truth_had_momentum_x); 
  t->Branch("truth_had_momentum_y", &blk.truth_had_momentum_y); 
  t->Branch("truth_had_momentum_z", &blk.truth_had_momentum_z); 
  t->Branch("truth_had_momentum_e", &blk.truth_had_momentum_e); 
  t->Branch("truth_had_vtx_x",      &blk.truth_had_vtx_x); 
  t->Branch("truth_had_vtx_y",      &blk.truth_had_vtx_y); 
  t->Branch("truth_had_vtx_z",      &blk.truth_had_vtx_z); 
  t->Branch("truth_had_vtx_t",      &blk.truth_had_vtx_t); 
  t->Branch("truth_had_PDG",        &blk.truth_had_PDG); 

  t->Branch("nROIs", &blk.nROIs);
  t->Branch("roi_has_truth", &blk.roi_has_truth);
  t->Branch("lead_pdg", &blk.lead_pdg);
  t->Branch("lead_pdg_total_energy", &blk.lead_pdg_total_energy);
  t->Branch("lead_pdg_total_charge", &blk.lead_pdg_total_charge);
  t->Branch("total_energy", &blk.total_energy);
  t->Branch("total_charge", &blk.total_charge);
  t->Branch("hit_fraction", &blk.hit_fraction);
}

void wireana::WireAna::ResetTruthBranches(DataBlock_Truth &blk)
{
  blk.truth_intType = DEFAULT_VALUE;         
  blk.truth_nu_momentum_x = DEFAULT_VALUE; 
  blk.truth_nu_momentum_y = DEFAULT_VALUE; 
  blk.truth_nu_momentum_z = DEFAULT_VALUE; 
  blk.truth_nu_momentum_e = DEFAULT_VALUE; 
  blk.truth_nu_vtx_x = DEFAULT_VALUE; 
  blk.truth_nu_vtx_y = DEFAULT_VALUE; 
  blk.truth_nu_vtx_z = DEFAULT_VALUE; 
  blk.truth_nu_vtx_t = DEFAULT_VALUE; 
  blk.truth_nu_PDG = DEFAULT_VALUE; 
  blk.truth_lep_momentum_x = DEFAULT_VALUE; 
  blk.truth_lep_momentum_y = DEFAULT_VALUE; 
  blk.truth_lep_momentum_z = DEFAULT_VALUE; 
  blk.truth_lep_momentum_e = DEFAULT_VALUE; 
  blk.truth_lep_vtx_x = DEFAULT_VALUE; 
  blk.truth_lep_vtx_y = DEFAULT_VALUE; 
  blk.truth_lep_vtx_z = DEFAULT_VALUE; 
  blk.truth_lep_vtx_t = DEFAULT_VALUE; 
  blk.truth_lep_PDG = DEFAULT_VALUE; 
  blk.truth_had_momentum_x = DEFAULT_VALUE; 
  blk.truth_had_momentum_y = DEFAULT_VALUE; 
  blk.truth_had_momentum_z = DEFAULT_VALUE; 
  blk.truth_had_momentum_e = DEFAULT_VALUE; 
  blk.truth_had_vtx_x = DEFAULT_VALUE; 
  blk.truth_had_vtx_y = DEFAULT_VALUE; 
  blk.truth_had_vtx_z = DEFAULT_VALUE; 
  blk.truth_had_vtx_t = DEFAULT_VALUE; 
  blk.truth_had_PDG = DEFAULT_VALUE; 
  blk.nROIs=0;
  blk.lead_pdg.clear();
  blk.lead_pdg_total_energy.clear();
  blk.lead_pdg_total_charge.clear();
  blk.total_energy.clear();
  blk.total_charge.clear();
  blk.hit_fraction.clear();

}


void wireana::WireAna::FillTruthBranches(art::Event const& evt, TTree*t, DataBlock_Truth &blk)
{
  //FillROI Truth
  for( auto &pvv : plane_view_roi_map )
  {
    for ( auto & vv : pvv.second )
    {
      for (auto & roi : vv.second )
      {
        ++blk.nROIs;
        blk.roi_has_truth.push_back( roi.hasTrueSignal );
        blk.lead_pdg.push_back( roi.true_leading_pdg );
        blk.lead_pdg_total_energy.push_back( roi.true_leading_energy_deposit );
        blk.lead_pdg_total_charge.push_back( roi.true_leading_electron_deposit );
        blk.total_energy.push_back( roi.true_energy_deposit );
        blk.total_charge.push_back( roi.true_electron_deposit );
        blk.hit_fraction.push_back( roi.tdc_hit_fraction );
      }
    }
  }


}

template<class T>
void wireana::WireAna::SortWirePtrByChannel( std::vector<art::Ptr<T>> &vec, bool increasing )
{
  if( fLogLevel >= 10 ) 
  {
    std::cout<<"Entering SortWirePtrByChannel, sorting "<<vec.size()<<"channels."<<std::endl;
  }
  if (increasing)
  {
    std::sort(vec.begin(), vec.end(), [](art::Ptr<T> &a, art::Ptr<T> &b) { return a->Channel() < b->Channel(); });
  }
  else
  {
    std::sort(vec.begin(), vec.end(), [](art::Ptr<T> &a, art::Ptr<T> &b) { return a->Channel() > b->Channel(); });
  }
}

template void wireana::WireAna::SortWirePtrByChannel<>( std::vector<art::Ptr<recob::Wire>> &vec, bool increasing );
template void wireana::WireAna::SortWirePtrByChannel<>( std::vector<art::Ptr<sim::SimChannel>> &vec, bool increasing );

std::vector<wireana::wirecluster> 
wireana::WireAna::BuildInitialClusters( std::vector<art::Ptr<recob::Wire>> &vec, int dC, int dTick )
{
  std::vector<wireana::wirecluster> ret;

  bool newcluster = true;
  for( auto it = vec.begin(); it!=vec.end(); )
  {
    if (newcluster) 
    {
      wirecluster clus;
      ret.push_back( clus );
      ret.back().push_back( *it );
      newcluster = false;
      ++it;
      continue;
    }

    if( abs( int( (*it)->Channel() - ret.back().back()->Channel()) ) == 1 )
    {
      ret.back().push_back( *it );
      ++it;
    } else 
    {
      newcluster = true;
    }
  }

  //now ret is a vector of wirecluster object
  //each wirecluster contains set of wires adjacent to each other
  return ret;
}



bool 
wireana::WireAna::HasHit( const art::Ptr<recob::Wire> &wire, int minTick )
{
  if( wire->SignalROI().n_ranges() == 0 ) return false;
  for( auto itr = wire->SignalROI().begin_range(); itr!=wire->SignalROI().end_range(); ++itr )
  {
    int dTick = itr->end_index() - itr->begin_index();
    if (dTick > minTick) return true;
  }
  return false;

}

void wireana::WireAna::PrintClusters( std::vector<wirecluster> &clusters )
{
  for( auto it = clusters.begin(); it != clusters.end(); it++ )
  {
    std::cout<<"Printing Cluster "<<it-clusters.begin()<<std::endl;
    std::cout<<Form("Nwires: %d, Min Channel: %d, Max Channel: %d", it->nWires, it->wire_pointers.front()->Channel(), it->wire_pointers.back()->Channel() )<<std::endl;
  }

}


void
wireana::WireAna::BuildPlaneViewROIMap(  std::vector<art::Ptr<recob::Wire>> &wires )
{
  //First get all rois
  //I assume the wires are already grouped by APA and view, because we don't want to
  //search over the entire DUNE volume all at once.
  for( const auto &wire: wires )
  {
    raw::ChannelID_t chan = wire->Channel();
    int planeID = chan/fNChanPerApa;
    recob::Wire::RegionsOfInterest_t signalROI = wire->SignalROI();
    auto dataranges = signalROI.get_ranges();
    if( fLogLevel>=10 ) std::cout<<"Channel "<<chan<<" has "<<dataranges.size()<<" ranges"<<std::endl;
    for ( auto datarange: dataranges ) {
      if( fLogLevel>=10 ) std::cout<<"BuildPlaneViewROIMap(): create roi"<<std::endl;
      wireana::roi thisroi( wire, datarange );
      if( fLogLevel>=10 ) std::cout<<"BuildPlaneViewROIMap(): created roi"<<std::endl;
      plane_view_roi_map[planeID][wire->View()].push_back( thisroi );
      if( fLogLevel>=10 ) std::cout<<"BuildPlaneViewROIMap(): push_back -- done"<<std::endl;
    }
  }
}

void wireana::WireAna::BuildInitialROIClusters()
{
  wireana::WireAnaDBSCAN scanner;
  for( auto &p: plane_view_roi_map ) //loop over planeid
  {
    //Internal clustering of ROIs
    for( auto &v: p.second ) //loop over view
    {
      scanner.SetParameters(v.second, fMinPts, fEps, fDrift, fPitch );
      scanner.run();
      std::map<int, wireana::roicluster> idclusmap;
      for( unsigned int i = 0; i<v.second.size(); i++ )
      {
        idclusmap[ v.second[i].clusterID ].AddROI( v.second[i], p.first );
      }
      for( auto idclus: idclusmap )
      {
        plane_view_roicluster_map[p.first][v.first].push_back(idclus.second);
      }
      std::sort(plane_view_roicluster_map[p.first][v.first].begin(), plane_view_roicluster_map[p.first][v.first].end(),[]( auto &a, auto &b ){ return a.nWires > b.nWires; } );
  //std::sort( trkID_sum.begin(), trkID_sum.end(), [](auto &a, auto &b){ return a.second.first > b.second.first;} );
    }
  }

  //if( fLogLevel >= 3 )
  //{
  //  for( const auto &p: plane_view_roi_map )
  //  {
  //    std::cout<<"==============================="<<std::endl;
  //    std::cout<<"Print Plane: "<< p.first <<std::endl;
  //    for( const auto &v: p.second )
  //    {
  //      std::cout<<"    View: "<<v.first<<std::endl;
  //      PrintROIs( v.second );
  //    }
  //  }
  //}

  if( fLogLevel >= 10 )
  {
    for( const auto &p: plane_view_roicluster_map)
    {
      std::cout<<"==============================="<<std::endl;
      std::cout<<"Print Plane: "<< p.first <<std::endl;
      for( const auto &v: p.second )
      {
        std::cout<<"    View: "<<v.first<<std::endl;
        for ( const auto &cluster: v.second ) PrintROIs( cluster.ROIs );
      }
    }
  }


  return;
}


void 
wireana::WireAna::PrintROIs( const std::vector<wireana::roi> &ROIs)
{
  std::cout<<"Printing ROIs"<<std::endl;
  for( auto roi: ROIs )
  {
    if (roi.clusterID == NOISE ) continue;
    std::cout<<Form(
        "\t\tClusID: %d \n\t\t\t\tChannel: %d, Index: (%d,%d,%d), |Sum|: %f", 
        roi.clusterID, 
        roi.channel, roi.begin_index, roi.end_index, roi.end_index-roi.begin_index,roi.abs_sum
        )
      <<std::endl;
  }
}



void
wireana::WireAna::FillTrackIDtoLabelMap( art::Event const& evt )
{
  ////////////////////////////////////////////////////////////
  //Build trackid -> truth particle map
  //parse MCParticles
  art::Handle<std::vector<simb::MCParticle>> particleHandle;
  if (!evt.getByLabel(fSimulationProducerLabel, particleHandle)) {
    throw cet::exception("AnalysisExample")
      << " No simb::MCParticle objects in this event - "
      << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  } 

  /////////////////////////////////////////////////////////////
  // Get MCTruth handle
  // Get Trackid->generator info
  auto mcHandles = evt.getMany<std::vector<simb::MCTruth>>();

  //! Get a map of Art::Ptr<MCParticle> to the truth label
  for (auto const& mcHandle : mcHandles) {
    const std::string& sModuleLabel = mcHandle.provenance()->moduleLabel();
    art::FindManyP<simb::MCParticle> findMCParts(mcHandle, evt, fSimulationProducerLabel);
    std::vector<art::Ptr<simb::MCParticle>> mcParts = findMCParts.at(0);
    for (const art::Ptr<simb::MCParticle> ptr : mcParts) {
      trkid_to_label_map[ptr->TrackId()] = sModuleLabel;
    }
  }
}





void 
wireana::WireAna::TagAllROIClusterTruth(const detinfo::DetectorClocksData &clock )
{
  for( auto &pvv : plane_view_roicluster_map )
  {
    for ( auto & vv : pvv.second )
    {
      for (auto & cluster : vv.second )
      {
        TagROIClusterTruth( cluster, clock );
      }
    }
  }
  if (fLogLevel >= 3 )
  {
    std::cout<<"========================================="<<std::endl;
    std::cout<<"Printing Truth:"<<std::endl;
    for( auto &pvv : plane_view_roicluster_map )
    {
      for ( auto & vv : pvv.second )
      {
        std::cout<<Form("Plane %d, View %d: ", pvv.first,vv.first)<<std::endl;
        for (auto & cluster : vv.second )
        {
          std::cout<<"-----------------"<<std::endl;
          std::cout<<
            Form("  Cluster ID: %d, (wMin, wMax):(%d, %d), (tMin, tMax):(%d, %d)",
              cluster.clusterID,
              cluster.channel_min, cluster.channel_max, 
              cluster.begin_index, cluster.end_index )
            <<std::endl;
          if (cluster.label_code != 0 )
          {
            std::cout<<
              Form("   Truth Label: %s, Leading PDG: %d, Leading #electron: %f",
                 cluster.label.c_str(), cluster.pdg_energy_list.front().first,
                  cluster.pdg_energy_list.front().second.first)
              <<std::endl;
          }
          else std::cout<<"   Has no Truth"<<std::endl;
          std::cout<<
            Form("   centroid (Channel,Ticks): (%f, %f). width(c,t):(%i,%i)", 
               cluster.abs_centroidChannel, cluster.abs_centroidIndex,
               cluster.GetWidthChannel(),  cluster.GetWidthTick() )
            <<std::endl;
          if( cluster.truthFromNeutrino )
          {
          std::cout<<
            Form("   Has Neutrino Truth: Pnu(%f,%f,%f,%f), Ppart(%f,%f,%f,%f)",
                cluster.momentum_neutrino.Px(),
                cluster.momentum_neutrino.Py(),
                cluster.momentum_neutrino.Pz(),
                cluster.momentum_neutrino.E(),
                cluster.momentum_part.Px(),
                cluster.momentum_part.Py(),
                cluster.momentum_part.Pz(),
                cluster.momentum_part.E())
            <<std::endl;
          }
        }
      }
    }
  }

}

//! Tagging ROI Truth --
//! Look at each tick and search for true particle deposit and the
//! generator tag
void
wireana::WireAna::TagROIClusterTruth( wireana::roicluster &cluster, const detinfo::DetectorClocksData &clock )
{
  if (fLogLevel>=10) std::cout<<"Entering TagROIClusterTruth"<<std::endl;

  std::vector< std::pair< int, std::pair<float,float> > > trkID_sum; 
  for ( auto& roi: cluster.ROIs ) // loop 1
  {
    int roi_channel = roi.wire->Channel();
    if (fLogLevel>=10) 
    {
      std::cout<<"TagROIClusterTruth:look at channel "<<roi_channel<<std::endl;
      std::cout<<"            Map Size: "<<ch_w_sc.size()<<std::endl;
      std::cout<<"            Channel ID Exists? "<< (ch_w_sc.find(roi_channel) != ch_w_sc.end())<<std::endl;
      std::cout<<"            Wire Exist? "<<ch_w_sc[ roi_channel ].first->Channel()<<std::endl;
      std::cout<<"            SimChannel Exist? "<<ch_w_sc[ roi_channel ].second->Channel()<<std::endl;
    }
    art::Ptr<sim::SimChannel> &simchannel = ch_w_sc[ roi_channel ].second;

    //! Get the TPC time from the ticks in the ROI datarange_t
    double startTDC = clock.TPCTick2TDC( roi.begin_index );
    double endTDC   = clock.TPCTick2TDC( roi.end_index );
    if (fLogLevel>=10) 
    {
      std::cout<<"            SimChannel exist? "<<simchannel->Channel()<<std::endl;
      std::cout<<"            startTDC, endTDC: "<<startTDC<<", "<<endTDC<<std::endl;
    }

    //! Get all the IDE (Ionization at a point of the TPC sensitive volume. ) that deposited energy in the datarange_t
    std::vector< sim::IDE >  ide_vec = simchannel->TrackIDsAndEnergies(startTDC, endTDC);
    if (fLogLevel>=10) std::cout<<"           ide_vec size: "<<ide_vec.size()<<std::endl;

    //Fill True Energy Charge Deposited in ROI

    if( ide_vec.size() == 0 ) continue;

    //! This ROI has true signal
    roi.hasTrueSignal = true;

    //!pair<trkid, pair<numElectrons, energy>> , sum up n electrons and total energy within the ROI for each trkid

    //! fill vector of ide
    if (fLogLevel>=10) std::cout<<"           fill truth loop: begin"<<std::endl;
    for( auto ide : ide_vec ) //loop 5
    {
      //! check if track has been saved in the trkID_sum vector.
      auto it = trkID_sum.begin();
      for( ; it!=trkID_sum.end(); ++it ) {
        if( it->first == ide.trackID ) break; 
      }

      if (it == trkID_sum.end()) {
        //! if not, add the track id and initialize with the IDE data
        std::pair< int, std::pair<float,float> > data;
        data.first = ide.trackID;
        data.second.first =ide.numElectrons;
        data.second.second = ide.energy;
        trkID_sum.push_back(data);
      } else {
        //! if yes, sum nElectron and energy to the track id.
        it->second.first+=ide.numElectrons;
        it->second.second+=ide.energy;
      }
    }// end loop 5

    if (fLogLevel>=10) std::cout<<"           fill truth loop: end"<<std::endl;

    //! sort vector of (trkid,(#electron, energy)) by # electrons
    //! could also sort by energy but #e should be a more direct observable?
  }//end loop 1


  if (fLogLevel>=10) std::cout<<"   filled all trkID_sum"<<std::endl;
  // Once all ROI has been analyzed, sort trkID by #electrons
  std::sort( trkID_sum.begin(), trkID_sum.end(), [](auto &a, auto &b){ return a.second.first > b.second.first;} );
  if (fLogLevel>=10) std::cout<<"   trkID_sum sorted"<<std::endl;

  cluster.trkID_sum = trkID_sum;
  bool hasMCTrack =  trkID_sum.size() > 0;
  if( hasMCTrack ) cluster.label = this->trkid_to_label_map[ trkID_sum[0].first ];
  else cluster.label = "None";

  if( cluster.label == "None" ) cluster.label_code = kNone;
  else if( cluster.label == "NuEScatter" ) cluster.label_code = kElastic;
  else if( cluster.label == "marley" ) cluster.label_code = kMarley;
  else cluster.label_code = kOthers;

  if (fLogLevel>=10) std::cout<<"   label set"<<std::endl;
  //Grep the most deposited particle first

  if( !hasMCTrack ) return;
  
  int trkid = trkID_sum.front().first;
  const simb::MCParticle * part = PIS->TrackIdToMotherParticle_P(trkid);
  art::Ptr<simb::MCTruth> truth = PIS->TrackIdToMCTruth_P(trkid);

  if ( truth->NeutrinoSet() )
  {
    cluster.truthFromNeutrino = true;
    const simb::MCNeutrino neutrino = truth->GetNeutrino();
    cluster.momentum_neutrino = neutrino.Nu().Momentum();
    for( int ipart = 0; ipart < truth->NParticles(); ++ipart )
    {
      int pdg = truth->GetParticle( ipart ).PdgCode();
      if (abs(pdg) == 12 || abs(pdg) == 14 || abs(pdg) == 16) cluster.n_nu++;
      else if(abs(pdg) == 11 || abs(pdg) == 13 || abs(pdg) == 15) cluster.n_lepton++;
      else if (abs(pdg) == 2212) cluster.n_proton++;
      else if (abs(pdg) == 2112) cluster.n_neutron++;
      else if (abs(pdg) == 22 ) cluster.n_photon++;
      else if (abs(pdg)>=100 && abs(pdg)<600 ) cluster.n_meson++;
      else if (abs(pdg)>= 1e10 ) cluster.n_nucleus++;
    }
  }
  cluster.momentum_part = part->Momentum();

  for( auto &trkid : trkID_sum )
  {
    auto part1 = PIS->TrackIdToMotherParticle_P(trkid.first);
    int pdg = part1->PdgCode();
    cluster.pdg_energy_list.push_back( std::make_pair( pdg, trkid.second ));
  }
}



void 
wireana::WireAna::TagAllROITruth(const detinfo::DetectorClocksData &clock )
{
  for( auto &pvv : plane_view_roi_map )
  {
    for ( auto & vv : pvv.second )
    {
      for (auto & roi : vv.second )
      {
        TagROITruth( roi, clock );
      }
    }
  }
  if (fLogLevel >= 4 )
  {
    std::cout<<"========================================="<<std::endl;
    std::cout<<"Printing ROI Truth:"<<std::endl;
    for( auto &pvv : plane_view_roi_map )
    {
      for ( auto & vv : pvv.second )
      {
        std::cout<<Form("Plane %d, View %d: ", pvv.first,vv.first)<<std::endl;
        for (auto & roi: vv.second )
        {
          std::cout<<"-----------------"<<std::endl;
          std::cout<<
            Form("  roi channel: %d, tick range: (%d, %d)",
              roi.channel, roi.begin_index, roi.end_index)
            <<std::endl;
          std::cout<<
            Form("  has_truth: %d, true_electron_dep: %.2f",
              roi.hasTrueSignal, roi.true_leading_energy_deposit )
            <<std::endl;
        }
      }
    }
  }

}

void
wireana::WireAna::TagROITruth( wireana::roi &roi, const detinfo::DetectorClocksData &clock )
{
  if (fLogLevel>=10) std::cout<<"Entering TagROITruth"<<std::endl;

  std::vector< std::pair< int, std::pair<float,float> > > trkID_sum; 
  int roi_channel = roi.wire->Channel();
  if (fLogLevel>=10) 
  {
    std::cout<<"TagROITruth:look at channel "<<roi_channel<<std::endl;
    std::cout<<"            Map Size: "<<ch_w_sc.size()<<std::endl;
    std::cout<<"            Channel ID Exists? "<< (ch_w_sc.find(roi_channel) != ch_w_sc.end())<<std::endl;
    std::cout<<"            Wire Exist? "<<ch_w_sc[ roi_channel ].first->Channel()<<std::endl;
    std::cout<<"            SimChannel Exist? "<<ch_w_sc[ roi_channel ].second->Channel()<<std::endl;
  }
  art::Ptr<sim::SimChannel> &simchannel = ch_w_sc[ roi_channel ].second;

  //! Get the TPC time from the ticks in the ROI datarange_t
  double startTDC = clock.TPCTick2TDC( roi.begin_index );
  double endTDC   = clock.TPCTick2TDC( roi.end_index );
  if (fLogLevel>=10) 
  {
    std::cout<<"            SimChannel exist? "<<simchannel->Channel()<<std::endl;
    std::cout<<"            startTDC, endTDC: "<<startTDC<<", "<<endTDC<<std::endl;
  }
  //--------------------Tagging Energy/Electron Deposits------------
  sim::SimChannel::TDC_t totalTDC = sim::SimChannel::TDC_t(endTDC)-sim::SimChannel::TDC_t(startTDC)+1;
  sim::SimChannel::TDC_t nTDC = 0;

  for( sim::SimChannel::TDC_t i = (sim::SimChannel::TDC_t) startTDC; i <= (sim::SimChannel::TDC_t) endTDC; i++ )
  {
    //Fill TrueEnergyChargeDepositedInROI
    std::vector< sim::IDE >  ides = simchannel->TrackIDsAndEnergies(i,i);

    std::vector<float> energies(partTypes.size(),fECMin );
    std::vector<float> charges(partTypes.size(),fECMin );
    std::vector<float> energiesNeut(partTypes.size(),fECMin );
    std::vector<float> chargesNeut(partTypes.size(),fECMin );
    std::vector<float> energiesRad(partTypes.size(),fECMin );
    std::vector<float> chargesRad(partTypes.size(),fECMin );

    for( auto &ide: ides )
    {
      bool isSignal = trkid_to_label_map[ ide.trackID ] == "NuEScatter" || trkid_to_label_map[ ide.trackID ] == "marley";
      int pdg = PIS->TrackIdToParticle_P( ide.trackID )->PdgCode();
      float energy = ide.energy;
      float numElectrons = ide.numElectrons;
      energies[kAll]+=energy;
      charges[kAll]+=numElectrons;
      int partType = -1;
      if( abs(pdg) == 11 || abs(pdg) == 13 || abs(pdg) == 15 )
      {
        partType = kElectron;
      } else if( abs(pdg) == 2212)
      {
        partType = kProton;
      } else if(abs(pdg) == 2112)
      {
        partType = kNeutron;
      } else if(abs(pdg) == 22)
      {
        partType = kPhoton;
      } else
      {
        partType = kNuc;
      }
      //parsed particle, accumulate energy
      energies[partType]+=energy; charges[partType]+=numElectrons;
      if( isSignal )
      {
        energiesNeut[partType]+=energy; chargesNeut[partType]+=numElectrons;
      }
      else
      {
        energiesRad[partType]+=energy; chargesRad[partType]+=numElectrons;
      }
    }

    FillHistogram( energies, charges, kSAll, true );
    FillHistogram( energiesNeut, chargesNeut, kSNeutrino, true );
    FillHistogram( energiesRad, chargesRad, kSRad, true );

    TrueEnergyDepositedInROI->Fill(energies[0]);
    TrueChargeDepositedInROI->Fill(charges[0]);

    if (energies[0]>0 || charges[0]>0) 
    {
      ++nTDC;
      roi.true_energy_deposit+= energies[0];
      roi.true_electron_deposit+= charges[0];
    }
  }
  roi.tdc_hit_fraction = double(nTDC)/double(totalTDC);
  if(fLogLevel>=10)
  {
    std::cout<<Form("startTDC, endTDC: (%f,%f), nTDC: %d, totalTDC: %d",
        startTDC,endTDC, nTDC, totalTDC )<<std::endl;
  }


  //--------------------Tagging Truth Particles---------------------
  //! Get all the IDE (Ionization at a point of the TPC sensitive volume. ) that deposited energy in the datarange_t
  std::vector< sim::IDE >  ide_vec = simchannel->TrackIDsAndEnergies(startTDC, endTDC);
  if (fLogLevel>=10) std::cout<<"           ide_vec size: "<<ide_vec.size()<<std::endl;

  if( ide_vec.size() == 0 ) return;

  //! This ROI has true signal
  roi.hasTrueSignal = true;

  //!pair<trkid, pair<numElectrons, energy>> , sum up n electrons and total energy within the ROI for each trkid

  //! fill vector of ide
  if (fLogLevel>=10) std::cout<<"           fill truth loop: begin"<<std::endl;
  for( auto ide : ide_vec ) //loop 5
  {
    //! check if track has been saved in the trkID_sum vector.
    auto it = trkID_sum.begin();
    for( ; it!=trkID_sum.end(); ++it ) {
      if( it->first == ide.trackID ) break; 
    }

    if (it == trkID_sum.end()) {
      //! if not, add the track id and initialize with the IDE data
      std::pair< int, std::pair<float,float> > data;
      data.first = ide.trackID;
      data.second.first =ide.numElectrons;
      data.second.second = ide.energy;
      trkID_sum.push_back(data);
    } else {
      //! if yes, sum nElectron and energy to the track id.
      it->second.first+=ide.numElectrons;
      it->second.second+=ide.energy;
    }
  }// end loop 5

  if (fLogLevel>=10) std::cout<<"           fill truth loop: end"<<std::endl;

  //! sort vector of (trkid,(#electron, energy)) by # electrons
  //! could also sort by energy but #e should be a more direct observable?


  if (fLogLevel>=10) std::cout<<"   filled all trkID_sum"<<std::endl;
  // Once all ROI has been analyzed, sort trkID by #electrons
  std::sort( trkID_sum.begin(), trkID_sum.end(), [](auto &a, auto &b){ return a.second.first > b.second.first;} );
  if (fLogLevel>=10) std::cout<<"   trkID_sum sorted"<<std::endl;

  roi.hasTrueSignal =  trkID_sum.size() > 0;
  if( roi.hasTrueSignal )
  {
    roi.true_leading_pdg = PIS->TrackIdToParticle_P( trkID_sum.front().first )->PdgCode();
    roi.true_leading_electron_deposit = trkID_sum.front().second.first;
    roi.true_leading_energy_deposit = trkID_sum.front().second.second;
  }
}



std::vector<double>
wireana::WireAna::GetArrayFromWire( std::vector<art::Ptr<recob::Wire>> &wirelist, wireana::roicluster &cluster, int channel_width, int tick_width )
{
  if( fLogLevel >= 10 )
  {
    std::cout<<Form("GetArrayFromWire::Log: begin")<<std::endl;
  }
 
  double centroid_channel = cluster.abs_centroidChannel, centroid_tick = cluster.abs_centroidIndex;
  int c0,t0;
  c0= centroid_channel - channel_width/2. + 1;
  t0= centroid_tick - tick_width/2.;
  std::vector<double> ret( channel_width*tick_width, 0 );
  if( fLogLevel >= 3 )
  {
    std::cout<<Form("GetArrayFromWire::Log:  c(c,t):(%.2f, %.2f). (c0,t0):(%d,%d). cpos(%d,%d) ",centroid_channel, centroid_tick, c0, t0, int(centroid_channel-c0), int(centroid_tick-t0))<<std::endl;
  }
  auto ptr1 = wirelist.begin();
  while ( (int) (*ptr1)->Channel() < c0 ) ++ptr1; //increment ptr1 until it reaches the first channel within range
  //while (  ptr1 != wirelist.end() ) 
  //{
  //  if( (*ptr1)->Channel() < c0 ) ptr1++;
  //  else break;
  //}
  auto ptr_end = ptr1;
  while ((int)  (*ptr_end)->Channel() < c0+channel_width ) 
  {
    ++ptr_end; //increment ptr1 until it reaches the first channel within range
    if( ptr_end == wirelist.end() ) break;
  }
  //while (  ptr_end != wirelist.end() ) 
  //{
  //  if( (*ptr_end)->Channel() < c0+channel_width ) ++ptr_end; 
  //  else break;
  //}
  if( fLogLevel >= 3 )
  {
    std::cout<<Form("                        Print ptr:")<<std::endl;
    std::cout<<Form("                        (ch_begin,ch_end): (%d, %d)",(*ptr1)->Channel(),(*(ptr_end-1))->Channel() )<<std::endl;
  }


  //increment ptr_end until it reaches the first channel outside the range. 
  //i.e. cfirst, clast = 1,5 --> width = 5
  //1+5 = 6 is the first channel outside the range
  //so ptr_end will end when channel >= 6 
  //

  for ( auto ptr = ptr1; ptr != ptr_end; ++ptr )
  {
    int channel = (*ptr)->Channel();
    int dC = channel - c0;
    if( dC>=channel_width) break;
    for( int i = 0; i < tick_width; ++i )
    {
      int index = dC + i*channel_width ;
      ret[index] = (*ptr)->SignalROI()[t0+i];
    }
  }
   if( fLogLevel >= 3 )
  {
    std::cout<<Form("GetArrayFromWire::Log: processed")<<std::endl;
  }
 
  return ret;
  
}

std::vector<double> 
wireana::WireAna::CombineTicks( const std::vector<double> &input, int channel_width, int nticks)
{
  std::vector<double> output;
  std::string func_name = "CombineTicks: ";
  int n_totalticks = input.size()/channel_width;

  if( n_totalticks%nticks != 0 ) 
  {
    std::cout<<"Err: CombineTicks: nticks does not divide total number of ticks"<<std::endl;
    return output;
  }

  std::vector<unsigned int> index_list;
  for( int i = 0; i< n_totalticks/nticks; i++ )
  {
    unsigned int start_index = (i*nticks)*channel_width; //start index at each chunk
    unsigned int end_index = start_index+channel_width;  //end index of first line in each chunk
    if(fLogLevel>=10) std::cout<<func_name<<"start_index: "<<start_index<<" end_index: "<<end_index<<std::endl
                               <<func_name<<"summing: ";
    for( unsigned int index = start_index; index<end_index; ++index )
    {
      //for each pass, loop through channels once. At each channel sum up nticks together
      double v=0;
      for( int t = 0; t < nticks; ++t ) 
      {
        unsigned int sum_index=index+ t*channel_width;
        v+= input[sum_index];
        if(fLogLevel>=10) std::cout<<func_name<<"("<<sum_index<<", "<<input[sum_index]<<") ";
      }
      output.push_back( v );
      if( fLogLevel >= 10 ) std::cout<<func_name<<std::endl;;
    }
    if( fLogLevel >= 10 ) std::cout<<func_name<<std::endl;

  }
  return output;
}


std::vector<double> 
wireana::WireAna::ScaleArray( const std::vector<double> &input, double min, double max )
{
  std::vector<double> output;
  if( max-min <= 0 ) return output;
  double delta = max - min;
  double this_min = *std::min_element(input.begin(), input.end());
  double this_max = *std::max_element(input.begin(), input.end());
  double this_delta = this_max-this_min;
  for( auto v: input )
  {
    //linearly transform element to new values
    // map to (0,1): v1 = (v-this_min)/this_delta
    // map to (min, max) vnew = v1*delta + min 
    if(v == 0 )
    {
      //do not transform 0, 
      //TODO: initialize array with default value, i.e. -999
      output.push_back(0);
      continue;
    }
    double vnew = (v-this_min)/this_delta * delta + min;
    output.push_back(vnew);
  }
  return output;
}

int 
wireana::WireAna::CalculateIndex( int c, int t, int c_width, int t_width )
{
  return (c%c_width) + (t%t_width)*c_width;
}
std::pair<int,int> 
wireana::WireAna::CalculateCT( int index, int c_width )
{
  int c = index%c_width;
  int t = index/c_width;
  return std::make_pair( c, t);
}


void
wireana::WireAna::WriteNumPy( wireana::matchedroicluster& cluster,  std::vector<art::Ptr<recob::Wire>>& wirelist)
{
  if(cluster.clusters.size() != 3 ) return;

  c2numpy_uint16(&npywriter, (unsigned int) run );
  c2numpy_uint16(&npywriter, (unsigned int) subrun );
  c2numpy_uint16(&npywriter, (unsigned int) event );
  c2numpy_uint16(&npywriter, (unsigned int) MC );

  c2numpy_uint16(&npywriter, (unsigned int) image_channel_width );
  c2numpy_uint16(&npywriter, (unsigned int) image_tick_width );
  c2numpy_uint16(&npywriter, (unsigned int) image_rebin_tick );

  c2numpy_uint16(&npywriter, (unsigned int) cluster.gencode() );
  c2numpy_uint16(&npywriter, (unsigned int) cluster.labelmatch() );
  c2numpy_uint16(&npywriter, (unsigned int) cluster.trkmatch());

  c2numpy_uint16(&npywriter, (unsigned int) cluster.clusters[0].n_photon );
  c2numpy_uint16(&npywriter, (unsigned int) cluster.clusters[0].n_proton );
  c2numpy_uint16(&npywriter, (unsigned int) cluster.clusters[0].n_neutron);
  c2numpy_uint16(&npywriter, (unsigned int) cluster.clusters[0].n_meson  );
  c2numpy_uint16(&npywriter, (unsigned int) cluster.clusters[0].n_nucleus);

  c2numpy_float(&npywriter, cluster.clusters[0].momentum_neutrino.Px() );
  c2numpy_float(&npywriter, cluster.clusters[0].momentum_neutrino.Py() );
  c2numpy_float(&npywriter, cluster.clusters[0].momentum_neutrino.Pz() );
  c2numpy_float(&npywriter, cluster.clusters[0].momentum_part.Px() );
  c2numpy_float(&npywriter, cluster.clusters[0].momentum_part.Py() );
  c2numpy_float(&npywriter, cluster.clusters[0].momentum_part.Pz() );

  std::vector<double> u_vec = GetArrayFromWire( wirelist, cluster.clusters[0], image_channel_width, image_tick_width );
  std::vector<double> v_vec = GetArrayFromWire( wirelist, cluster.clusters[1], image_channel_width, image_tick_width );
  std::vector<double> z_vec = GetArrayFromWire( wirelist, cluster.clusters[2], image_channel_width, image_tick_width );

  std::shared_ptr<std::vector<double>> u_vecc = std::make_shared< std::vector<double> >( CombineTicks( u_vec, image_channel_width, image_rebin_tick )) ;
  std::shared_ptr<std::vector<double>> v_vecc = std::make_shared< std::vector<double> >( CombineTicks( v_vec, image_channel_width, image_rebin_tick )) ;
  std::shared_ptr<std::vector<double>> z_vecc = std::make_shared< std::vector<double> >( CombineTicks( z_vec, image_channel_width, image_rebin_tick )) ;

  std::vector< std::shared_ptr<std::vector<double>> > vecc({ u_vecc,
                                                             v_vecc,
                                                             z_vecc } ); 
  for( int i = 0; i < image_size; i++ )
  {
    for( int v = 0; v < 3; v++ )
    {
      std::string name=Form("%s_%08d", fViewMap[v].c_str(),i);
      c2numpy_float(&npywriter, vecc[v]->at(i) );
    }
  }
}

//============================ Histogram Functions ================================
void
wireana::WireAna::FillHistogram(std::vector<float> energies, std::vector<float>charges ,SType s, bool inROI)
{
  string hnamebase = ( (inROI)? "TrueEnergyChargeDepositedInROI" : "TrueEnergyChargeDeposited" )+selTypes[s];
  for( unsigned int i = 0; i< partTypes.size(); i++ )
  {
    string hname = (hnamebase+partTypes[i]);
    if( fLogLevel >= 10 ) std::cout<<hname<<std::endl;
    TH2FMap[hname]->Fill(energies[i],charges[i]);
  }
}



DEFINE_ART_MODULE(wireana::WireAna)

#include "ROIAna_module.h"
//Constructor for cnn struct


roiana::ROIAna::ROIAna(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset}  ,
  fWireProducerLabel(pset.get< art::InputTag >("InputWireProducerLabel", "caldata")),
  fRawProducerLabel(pset.get< art::InputTag >("InputRawProducerLabel", "tpcrawdecoder:daq")),
  fSimChannelLabel(pset.get< art::InputTag >("SimChannelLabel", "elecDrift")),
  fSimulationProducerLabel(pset.get< art::InputTag >("SimulationProducerLabel", "largeant"))
{
  fLogLevel           = pset.get<int>("LogLevel", 10);
  fNChanPerApa        = pset.get<int>("ChannelPerApa", 2560);
  fNTicksPerWire      = pset.get<int>("TicksPerWire", 6000);
  //auto const* geo = lar::providerFrom<geo::Geometry>();
  geo = lar::providerFrom<geo::Geometry>();
  //geo = art::ServiceHandle<geo::Geometry>::get();
  fNPlanes = geo->Nplanes();

  fROI_Peak  = pset.get<float>("ROIPeak", 20);
  fROI_Range = pset.get<int>("ROIRange", 120);
  fROI_CH    = pset.get<int>("ROICH", 10);

  fTreeName          = pset.get<std::string>("TREENAME", "wireana");
  fECMin             = pset.get<float>("fECMin", -1e-8); //setting energy and charge accumulation to start at this value --> ensures true background, i.e 0 energy/0 charge falls into the underflow bin
  fHistEnergyMax     = pset.get<float>("fHistEnergyMax", 100); 
  fHistChargeMax     = pset.get<float>("fHistChargeMax", 1000); 

}

void roiana::ROIAna::analyze(art::Event const & evt) {


  //get detector property
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

  //get event data
  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();
  std::cout<<"########## EvtNo."<<event<<std::endl;
  std::cout<<"ROI peak threshold: "<<fROI_Peak<<std::endl;
  std::cout<<"ROI timetick width: "<<fROI_Range<<std::endl;
  std::cout<<"ROI channel width: "<<fROI_CH<<std::endl;
  std::cout<<"Let's TRY!"<<std::endl;

  MC = !evt.isRealData();
  ////////////////////////////////////////////////////////////
  //FillTruthInfo to internal data objects
  /// i.e. trkid_to_label_map

  // This attempt used MCTruths_marley, but didn't give TrackIds that matched SimChannels
  
  /*
  auto MarleyTruthHandle = evt.getValidHandle<std::vector<simb::MCTruth>>("marley");
  //std::cout << "MarleyTruth.size()=" << MarleyTruth->size() << std::endl;
  for (size_t i=0; i< MarleyTruthHandle->size(); i++) {
    simb::MCTruth MarleyTruth = MarleyTruthHandle->at(i);
    //std::cout << "MarleyTruth.size()=" << MarleyTruth.NParticles() << std::endl;
    for (int j=0; j < MarleyTruth.NParticles(); j++) {
      const simb::MCParticle ThisPar = MarleyTruth.GetParticle(j);
      std::cout << "MarleyTruth.TrackId()=" << ThisPar.TrackId() << std::endl;
      trkid_to_label_map[ThisPar.TrackId()] = "marley";
    }
  }
  */

  // This isanother attempt at pulling out truth info for a sensitivity study
  /*
  auto LargeantAssnsHandle = evt.getValidHandle<art::Assns<simb::MCTruth,simb::MCParticle,sim::GeneratedParticleInfo>>("largeant");
  //std::cout << "LargeantAssnsHandle.size()=" << LargeantAssnsHandle->size() << std::endl;
  for (auto i=LargeantAssnsHandle->begin(); i!= LargeantAssnsHandle->end(); i++) {
    //auto LargeantAssn = LargeantAssnsHandle->at(i);
    //const simb::MCTruth* ThisParMCTruth = i->first.get();
    const simb::MCParticle* ThisParMCParticle = i->second.get();
    //const sim::GeneratedParticleInfo& ThisParData = data(i);
    //const simb::MCGeneratorInfo ThisParGenInfo = ThisParMCTruth->GeneratorInfo();
    
    //std::cout << "ThisParGenInfo.generatorConfig.size()=" << ThisParGenInfo.generatorConfig.size() << std::endl;
    //std::cout << "ThisParMCTruth.NParticles()=" << ThisParMCTruth->NParticles() << std::endl;
    if (ThisParMCParticle->PdgCode() == 12) {
      std::cout << "Found a neutrino - ThisParMCParticle.TrackId()=" << ThisParMCParticle->TrackId() << std::endl;
      trkid_to_label_map[ThisParMCParticle->TrackId()]="marley?";
    }
    //std::cout << "ThisParMCParticle.TrackId()=" << ThisParMCParticle->TrackId() << std::endl;
  }
  */

  ////////////////////////////////////////////////////////////
  //Build Wire SimChannel List

  art::Handle<std::vector<sim::SimChannel>> simChannelListHandle;
  std::vector<art::Ptr<sim::SimChannel>> channellist;
  if (evt.getByLabel(fSimChannelLabel, simChannelListHandle)) 
    art::fill_ptr_vector(channellist, simChannelListHandle);
  SortWirePtrByChannel( channellist, true );

  ////////////////////////////////////////////////////////////
  //Build RawDigit List
  art::Handle<std::vector<raw::RawDigit>> rawListHandle, noiselessRawListHandle;
  std::vector<art::Ptr<raw::RawDigit>> rawList, noiselessRawList;

  if (evt.getByLabel(fRawProducerLabel, rawListHandle)) art::fill_ptr_vector(rawList, rawListHandle);
  SortWirePtrByChannel( rawList, true );

  // //Get channel-> wire,simchannel map
  if( fLogLevel >= 3 ) std::cout<<"Fill ch_w_sc"<<std::endl;
  for( auto w: rawList ) 
    ch_w_sc[ w->Channel() ].first = w;
  for( auto w: channellist) 
    ch_w_sc[ w->Channel() ].second= w;

  if( fLogLevel >= 3 ) std::cout << "starting TruthFilter" << std::endl;
  TruthFilter();

  //Creating ROI
  if( fLogLevel >= 3 ) std::cout << "starting ProcessROI" << std::endl;
  std::map<int,bool> ret;
  ret = ProcessROIWide(rawList);

  if( fLogLevel >= 3 ) std::cout << "starting ROIFilter" << std::endl;
  ROIFilter(ret);

  //compute efficiency and data reduction
  if( fLogLevel >= 3 ) std::cout << "starting ROIEfficiencies" << std::endl;
  int n_channels = rawList.size();
  ROIEfficiencies(ret, n_channels);
}


std::map<int,bool> roiana::ROIAna::ProcessROIWide(std::vector<art::Ptr<raw::RawDigit>>& rawList)
{
  int num_channels = rawList.size();

  // ret maps channel number to bool, where True means inROI
  std::map<int,bool> ret;
  for(int i=0; i<num_channels; i++)
  {
    ret[i] = false;
  }

  const int adc_size = rawList.front()->NADC();

  for( auto rawdigit : rawList )
  {
    int channel =  rawdigit->Channel();
    std::vector<float> charges;
    for( auto adc : rawdigit->ADCs() )
    {
      charges.push_back( adc  );
    }
    std::vector<float> newcharge(adc_size,0.0);
    auto chargessort = charges;

    std::nth_element(chargessort.begin(), chargessort.begin()+chargessort.size()/2, chargessort.end());
    float median = chargessort[ chargessort.size()/2];

    auto max_it = std::max_element(charges.begin(), charges.end());
    auto min_it = std::min_element(charges.begin(), charges.end());

    float adc_max = *max_it - median;
    float adc_min = *min_it - median;

    // Add channels to ROI
    if(adc_min<-fROI_Peak or adc_max>fROI_Peak)
    {
      for (int update_channel = (channel-fROI_CH); update_channel <= (channel+fROI_CH); update_channel++)
      {
        ret[update_channel] = true;
      }
    }
  }

  return ret;
}

void roiana::ROIAna::beginJob() {

  gROOT->SetBatch(1);

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>(fTreeName.c_str() ,fTreeName.c_str() );

  fTree->Branch("run", &run);
  fTree->Branch("subrun", &subrun);
  fTree->Branch("event", &event);
  fTree->Branch("MC", &MC);

  std::string name1 = Form("TrueEnergyDeposited_%s",fTreeName.c_str() );
  std::string name2 = Form("TrueEnergyDepositedInROI_%s",fTreeName.c_str() );
  std::string name3 = Form("TrueEnergyDepositedRatio_%s",fTreeName.c_str() );
  TrueEnergyDeposited = tfs->make<TH1F>( name1.c_str(), "Energy; E (MeV)",100,0,fHistEnergyMax);
  TrueEnergyDepositedInROI =  tfs->make<TH1F>( name2.c_str(), "Energy; E (MeV)",100,0,fHistEnergyMax);
  TrueEnergyDepositedRatio =  tfs->make<TH1F>( name3.c_str(), "Energy; E (MeV)",100,0,fHistEnergyMax);

  std::string name4 = Form("TrueChargeDeposited_%s",fTreeName.c_str() );
  std::string name5 = Form("TrueChargeDepositedInROI_%s",fTreeName.c_str() );
  std::string name6 = Form("TrueChargeDepositedRatio_%s",fTreeName.c_str() );
  TrueChargeDeposited = tfs->make<TH1F>( name4.c_str(), "Charge; Number of ionization electrons",100,0,fHistChargeMax);
  TrueChargeDepositedInROI =  tfs->make<TH1F>( name5.c_str(), "Charge; Number of ionization electrons",100,0,fHistChargeMax);
  TrueChargeDepositedRatio =  tfs->make<TH1F>( name6.c_str(), "Charge; Number of ionization electrons",100,0,fHistChargeMax);

  std::string name7 = Form("DataReductionRate_%s",fTreeName.c_str() );
  DataReductionRate = tfs->make<TH1F>( name7.c_str(), "Data retention fraction; (channels in ROI)/(total channels)",100,0,1);
}

void roiana::ROIAna::endJob()
{
}

template<class T>
void roiana::ROIAna::SortWirePtrByChannel( std::vector<art::Ptr<T>> &vec, bool increasing )
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

void roiana::ROIAna::TruthFilter()
{
  for( auto m: ch_w_sc)
  {
    //auto channel = m.first;
    //auto rawdigit = m.second.first;
    auto sim = m.second.second;
    for( auto &tdcide: sim->TDCIDEMap() )
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
        //std::cout<<"TruthFilter: ide.trackID: "<<ide.trackID<<std::endl;
        //std::string label = trkid_to_label_map[ ide.trackID ];
        //if (label == "marley?") {
        //  std::cout << "label is marley(?)" << std::endl;
        //  std::cout<<"ide.trackID: "<<ide.trackID<<std::endl;
        //}
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
        } else if(abs(pdg) == 12)
        {
          partType = kNeutrino;
          //std::cout<<"TruthFilter: neutrino with trackID: "<<ide.trackID<<std::endl;
        }else
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
          //std::cout << "label: " << label << std::endl;
        }
      }
      
      if( fLogLevel >= 3 ) std::cout<<"FillHistogram: begin"<<std::endl;
      //Some function to fill histogram
      TrueEnergyDeposited->Fill(energies[0]);
      TrueChargeDeposited->Fill(charges[0]);
      if( fLogLevel >= 3 ) std::cout<<"FillHistogram: end"<<std::endl;
    }
  }

}


void roiana::ROIAna::ROIFilter( std::map<int,bool> ret )
{
  for( auto m: ch_w_sc)
  {
    auto channel = m.first;
    //auto rawdigit = m.second.first;
    auto sim = m.second.second;
    if (!ret[channel]) continue;

    //accumulate energy and charge for all channels
    for( auto &tdcide: sim->TDCIDEMap() )
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
        //std::string label = trkid_to_label_map[ ide.trackID ];
        //if (label == "marley?") {
        //  std::cout << "label is marley(?)" << std::endl;
        //  std::cout<<"ide.trackID: "<<ide.trackID<<std::endl;
        //}
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
        } else if(abs(pdg) == 12)
        {
          partType = kNeutrino;
          //std::cout<<"ROIFilter: neutrino with trackID: "<<ide.trackID<<std::endl;
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
      //if( fLogLevel >= 3 ) std::cout<<"adding charge to ROI: "<<charges[0]<<std::endl;
      TrueEnergyDepositedInROI->Fill(energies[0]);
      TrueChargeDepositedInROI->Fill(charges[0]);
      if( fLogLevel >= 3 ) std::cout<<"FillHistogram: end"<<std::endl;

    }
  }
}

void roiana::ROIAna::ROIEfficiencies( std::map<int,bool> ret, int n_channels )
{
  // True energy & charge ratio in ROI
  TrueEnergyDepositedRatio->Divide(TrueEnergyDepositedInROI, TrueEnergyDeposited);
  TrueChargeDepositedRatio->Divide(TrueChargeDepositedInROI, TrueChargeDeposited);

  // data reduction rate estimation
  int n_channels_in_ROI = 0;
  for (const auto& ROIChannel : ret)
  {
    bool ChInROI = ROIChannel.second;
    if (ChInROI) n_channels_in_ROI += 1;
  }
  double ROIDataFraction = static_cast<double>(n_channels_in_ROI)/n_channels;
  std::cout << "Fraction of data kept by ROI: " << ROIDataFraction << std::endl;
  DataReductionRate->Fill(ROIDataFraction);
}


DEFINE_ART_MODULE(roiana::ROIAna)

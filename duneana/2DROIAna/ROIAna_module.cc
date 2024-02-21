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

  fLabels      = pset.get<std::vector<std::string>>("ParticleLabelVector");
  fInteraction = pset.get<std::vector<std::string>>("InteractionLabelVector");

  fROI_Peak  = pset.get<float>("ROIPeak", 20);
  //fROI_Range = pset.get<int>("ROIRange", 120);
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
  //std::cout<<"ROI timetick width: "<<fROI_Range<<std::endl;
  std::cout<<"ROI channel width: "<<fROI_CH<<std::endl;
  std::cout<<"Let's TRY!"<<std::endl;

  MC = !evt.isRealData();
  std::map<int,simb::MCParticle> ThisGeneratorParts;

  //--------------------------------------------------------------------------------------------//
  //--------------------------------- Create maps for ID tracking ------------------------------//
  //--------------------------------------------------------------------------------------------//
  // -------- Fill MC Truth IDs to tracking vectors. Get a list of all of my particles ---------//
  mf::LogInfo lparticle("particles");
  std::string lparticlestr = "";
  const sim::ParticleList& PartList = pi_serv->ParticleList();
  lparticle << "\nTotal #particles: " << PartList.size();
  
  // Loop over all signal+bkg handles and collect track IDs
  for ( size_t i = 0; i < fLabels.size(); i++){
    Parts.push_back(ThisGeneratorParts); // For each label insert empty list
    
    art::Handle<std::vector<simb::MCTruth>> ThisHandle;
    evt.getByLabel(fLabels[i], ThisHandle);
    
    if(ThisHandle){
      lparticlestr = PrintInColor(lparticlestr,"\n"+fLabels[i]+" *is generated!",GetColor("green"));
      
      auto ThisValidHandle = evt.getValidHandle<std::vector<simb::MCTruth>>(fLabels[i]); 
      art::FindManyP<simb::MCParticle> Assn(ThisValidHandle,evt,fSimulationProducerLabel);            
      
      FillMyMaps( Parts[i], Assn, ThisValidHandle);                                                                           
      //FillMCInteractionTree(Parts[i], fInteraction, fLogLevel);
      TPart.push_back(Parts[i].size()); // Insert #signal+bkg particles generated
      lparticlestr = PrintInColor(lparticlestr,"\n-> #Particles: "+std::to_string(Parts[i].size()),GetColor("blue"));
    }
    else{
      TPart.push_back(0);
      lparticlestr = PrintInColor(lparticlestr,"\n"+fLabels[i]+": *not generated!",GetColor("yellow"));
    }
  }
  lparticle << lparticlestr;
  fMCTruthTree->Fill();  

  ////////////////////////////////////////////////////////////
  //Build Wire SimChannel List

  art::Handle<std::vector<sim::SimChannel>> simChannelListHandle;
  std::vector<art::Ptr<sim::SimChannel>> channellist;
  if (evt.getByLabel(fSimChannelLabel, simChannelListHandle)) 
    art::fill_ptr_vector(channellist, simChannelListHandle);
  SortWirePtrByChannel( channellist, true );

  // MarleyTrackIDs and MarleyChannels store TrackId and Channel+IDE of energy deposits from Marley produced signals
  std::vector<int> MarleyTrackIDs;
  std::map<int,sim::IDE> MarleyChannels;
  //std::vector<float> Marley
  for (const auto& [TrackID, MCPart] : Parts[0]) {
    //std::cout << "Marley particle TrackID: " << TrackID << std::endl;
    MarleyTrackIDs.push_back(TrackID);
  }
  for (auto mySimChannel : channellist) {
    std::vector<sim::IDE> simIDEList = mySimChannel->TrackIDsAndEnergies(0, fNTicksPerWire);
    for (auto simIDE : simIDEList) {
      if (std::find(MarleyTrackIDs.begin(), MarleyTrackIDs.end(), simIDE.trackID) != MarleyTrackIDs.end()) {
        //std::cout << "marley particle in channel: " << mySimChannel->Channel() << std::endl;
        MarleyChannels[mySimChannel->Channel()] = simIDE;
        //MarleyChannels.push_back(mySimChannel->Channel());
      }
    }
  }

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
  ROIEfficiencies(ret, n_channels, MarleyChannels);
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
  fMCTruthTree = tfs->make<TTree>("MCTruthTree","MC Truth Tree");
  //fInteractionTree = tfs->make<TTree>("MCInteraction","MC Event Tree");

  fTree->Branch("run", &run);
  fTree->Branch("subrun", &subrun);
  fTree->Branch("event", &event);
  fTree->Branch("MC", &MC);

  fMCTruthTree -> Branch("Event",                &event,          "Event/I");  // Event number.
  fMCTruthTree -> Branch("Flag",                 &flag,           "Flag/I");   // Flag used to match truth with reco tree entries.
  fMCTruthTree -> Branch("TruthPart",            &TPart);                      // Number particles per generator.

  /*
  fInteractionTree -> Branch("Event",            &event,          "Event/I");  // Event number.
  fInteractionTree -> Branch("Flag",             &flag,           "Flag/I");   // Flag used to match truth with reco tree entries.
  fInteractionTree -> Branch("PDG",              &PDG,            "PDG/I");    // Main interacting particle PDG.
  fInteractionTree -> Branch("Energy",           &Energy,         "Energy/F"); // Main interacting particle energy [GeV^2].
  fInteractionTree -> Branch("Interaction",      &Interaction);                // Type of interaction.
  fInteractionTree -> Branch("Momentum",         &Momentum);                   // Main interacting particle momentum [GeV^2].
  fInteractionTree -> Branch("StartVertex",      &StartVertex);                // Main interacting particle start vertex [cm].
  fInteractionTree -> Branch("EndVertex",        &EndVertex);                  // Main interacting particle end vertex [cm].
  fInteractionTree -> Branch("DaughterPDG",      &DaughterPDG);                // Main interacting particle daughter PDG.
  fInteractionTree -> Branch("DaughterE",        &DaughterE);                  // Main interacting particle daughter energy [GeV^2].
  fInteractionTree -> Branch("DaughterPx",       &DaughterPx);                 // Main interacting particle daughter momentum X [GeV^2].
  fInteractionTree -> Branch("DaughterPy",       &DaughterPy);                 // Main interacting particle daughter momentum Y [GeV^2].
  fInteractionTree -> Branch("DaughterPz",       &DaughterPz);                 // Main interacting particle daughter momentum Z [GeV^2].
  fInteractionTree -> Branch("DaughterStartVx",  &DaughterStartVx);            // Main interacting particle daughter start vertex X [cm].
  fInteractionTree -> Branch("DaughterStartVy",  &DaughterStartVy);            // Main interacting particle daughter start vertex Y [cm].
  fInteractionTree -> Branch("DaughterStartVz",  &DaughterStartVz);            // Main interacting particle daughter start vertex Z [cm].
  fInteractionTree -> Branch("DaughterEndVx",    &DaughterEndVx);              // Main interacting particle daughter end vertex X [cm].
  fInteractionTree -> Branch("DaughterEndVy",    &DaughterEndVy);              // Main interacting particle daughter end vertex Y [cm].
  fInteractionTree -> Branch("DaughterEndVz",    &DaughterEndVz);              // Main interacting particle daughter end vertex Z [cm].
  */

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

  std::string name8 = Form("MarleySignalSensitivity_%s",fTreeName.c_str() );
  MarleySignalSensitivity = tfs->make<TH1F>( name8.c_str(), "Marley event energy fraction in ROI; (energy in ROI)/(total energy)",100,0,1);
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

void roiana::ROIAna::ROIEfficiencies( std::map<int,bool> ret, int n_channels, std::map<int,sim::IDE> MarleyChannels )
/*
 * Fill ROI performance histograms:
 *   */

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

  // Marley sensitivity
  float MarleyEnergyTot = 0.0;
  float MarleyChargeTot = 0.0;
  float MarleyEnergyROI = 0.0;
  float MarleyChargeROI = 0.0;
  for (const auto& Ch_IDE : MarleyChannels) {
    int MarleyCh = Ch_IDE.first;
    sim::IDE MarleyIDE = Ch_IDE.second;

    MarleyEnergyTot += MarleyIDE.energy;
    MarleyChargeTot += MarleyIDE.numElectrons;
    if (ret[MarleyCh]) {
      //if signal channel is in ROI
      MarleyEnergyROI += MarleyIDE.energy;
      MarleyChargeROI += MarleyIDE.numElectrons;
    } 

  }//for MarleyChannels
  std::cout << "Marley signal energy total (MeV): " << MarleyEnergyTot << std::endl;
  std::cout << "Marley signal energy fraction in ROI: " << MarleyEnergyROI/MarleyEnergyTot << std::endl;
  std::cout << "Marley signal charge fraction in ROI: " << MarleyEnergyROI/MarleyEnergyTot << std::endl;
  MarleySignalSensitivity->Fill(MarleyEnergyROI/MarleyEnergyTot);
}

//......................................................
void roiana::ROIAna::FillMCInteractionTree( std::map< int, simb::MCParticle> &MCParticleList, std::vector<std::string> ProcessList, int fLogLevel )
/*
Fill MCInteraction Tree with information about the main interaction in the event:
- MCParticleList is the list of MCParticles with a given generator label in the event
- ProcessList is the list of processes to be considered as main interactions
- HeavDebug is a boolean to turn on/off debugging statements
*/
{ 
  mf::LogInfo lheader("header");
  std::string lheaderstr = "";
  // Make a copy of MCParticleList to be used for finding Daughter info
  std::map< int, simb::MCParticle> MCParticleListCopy = MCParticleList;
  bool FoundInteraction;

  if (ProcessList.empty()){
    lheaderstr = PrintInColor(lheaderstr,"\n-> No processes in the list!",GetColor("red"));
    lheader << lheaderstr;
    return;  
  }
  
  for (size_t j = 0; j < ProcessList.size(); j++){
    FoundInteraction = false;    
    for ( std::map<int,simb::MCParticle>::iterator mainiter = MCParticleList.begin(); mainiter != MCParticleList.end(); mainiter++ ){
      if ( mainiter->second.Process() != ProcessList[j] && mainiter->second.EndProcess() != ProcessList[j]){continue;}
      if( fLogLevel >= 3 ) lheaderstr = lheaderstr+"\nFound a main interaction "+mainiter->second.EndProcess();
      FoundInteraction = true;
      simb::MCParticle MCParticle = mainiter->second;
      Interaction =  MCParticle.EndProcess();
      PDG =          MCParticle.PdgCode();
      Energy =       MCParticle.E();
      Momentum =    {MCParticle.Px(),MCParticle.Py(),MCParticle.Pz()};
      StartVertex = {MCParticle.Vx(),MCParticle.Vy(),MCParticle.Vz()};
      EndVertex =   {MCParticle.EndX(),MCParticle.EndY(),MCParticle.EndZ()};
      
      std::vector<int> DaughterList = {};
      for (int i = 0; i < MCParticle.NumberDaughters(); i++){
        DaughterList.push_back(MCParticle.Daughter(i));
      }
      // Print nice output with all the main interaction info
      if( fLogLevel >= 3 ) {
        lheaderstr = PrintInColor(lheaderstr,"\nMain interacting particle for process "+mainiter->second.Process()+": ",GetColor("magenta"));
        lheaderstr = PrintInColor(lheaderstr,"\nPDG ->\t"         + std::to_string(PDG),GetColor("cyan"));
        lheaderstr = PrintInColor(lheaderstr,"\nEnergy ->\t"      + std::to_string(Energy),GetColor("cyan"));
        lheaderstr = PrintInColor(lheaderstr,"\nMomentum ->\t"    + std::to_string(Momentum[0]) + " " + std::to_string(Momentum[1]) + " " + std::to_string(Momentum[2]),GetColor("cyan"));
        lheaderstr = PrintInColor(lheaderstr,"\nStartVertex ->\t" + std::to_string(StartVertex[0]) + " " + std::to_string(StartVertex[1]) + " " + std::to_string(StartVertex[2]),GetColor("cyan"));
        lheaderstr = PrintInColor(lheaderstr,"\nEndVertex ->\t"   + std::to_string(EndVertex[0]) + " " + std::to_string(EndVertex[1]) + " " + std::to_string(EndVertex[2]),GetColor("cyan"));
      }

      for ( std::map<int,simb::MCParticle>::iterator daughteriter = MCParticleListCopy.begin(); daughteriter != MCParticleListCopy.end(); daughteriter++ ){
        for (size_t i = 0; i < DaughterList.size(); i++){
          if (daughteriter->first == MCParticle.Daughter(i)){
            DaughterPDG.push_back(daughteriter->second.PdgCode());
            DaughterE.push_back(daughteriter->second.E());
            DaughterPx.push_back(daughteriter->second.Px());
            DaughterPy.push_back(daughteriter->second.Py());
            DaughterPz.push_back(daughteriter->second.Pz());
            DaughterStartVx.push_back(daughteriter->second.Vx());
            DaughterStartVy.push_back(daughteriter->second.Vy());
            DaughterStartVz.push_back(daughteriter->second.Vz());
            DaughterEndVx.push_back(daughteriter->second.EndX());
            DaughterEndVy.push_back(daughteriter->second.EndY());
            DaughterEndVz.push_back(daughteriter->second.EndZ());
          } // If the particle is a daughter of the main interaction
        } // Loop over all daughters
      } // Loop over all particles in the map
      fInteractionTree -> Fill();
    } // Loop over all particles in the map
    if( fLogLevel >= 3 ) {
      if (!FoundInteraction) lheaderstr = PrintInColor(lheaderstr,"\n-> No main interaction found for process "+ProcessList[j]+"!",GetColor("yellow"));
      else lheaderstr = PrintInColor(lheaderstr,"\n-> Filled MCINteraction Tree for process "+ProcessList[j]+"!",GetColor("green"));
    }
  } // Loop over all processes in the list
  if( fLogLevel >= 3 ) { lheader << lheaderstr;}
  return;
} // FillMCInteractionTree

//......................................................
void roiana::ROIAna::FillMyMaps( std::map< int, simb::MCParticle> &MyMap, art::FindManyP<simb::MCParticle> Assn, art::ValidHandle< std::vector<simb::MCTruth> > Hand )
/*
 * This function fills a map with the MCParticles from a given MCTruth
 * */ 
{
  for ( size_t L1=0; L1 < Hand->size(); ++L1 ) {
    for ( size_t L2=0; L2 < Assn.at(L1).size(); ++L2 ) {
      const simb::MCParticle ThisPar = (*Assn.at(L1).at(L2));
      MyMap[ThisPar.TrackId()] = ThisPar;
      if (fLogLevel >= 3 ) std::cout << ThisPar.PdgCode() << " " << ThisPar.E() << std::endl;
    }
  }
  return;
}

//......................................................
// This function checks if a given TrackID is in a given map
bool roiana::ROIAna::InMyMap( int TrID, std::map< int, simb::MCParticle> ParMap ){
  std::map<int, simb::MCParticle>::iterator ParIt;
  ParIt = ParMap.find( TrID );
  if (ParIt != ParMap.end()) {return true;}
  else return false;
}

//......................................................
// This function creates a terminal color printout
std::string roiana::ROIAna::PrintInColor( std::string InputString, std::string MyString, int Color ){
  std::string OutputString = InputString + "\033[" + std::to_string(Color) + "m" + MyString + "\033[0m";
  return OutputString;
}

// ......................................................
// This function returns an integer that corresponds to a given color name
int roiana::ROIAna::GetColor( std::string ColorName ){
  if (ColorName == "black") return 30;
  else if (ColorName == "red") return 31;
  else if (ColorName == "green") return 32;
  else if (ColorName == "yellow") return 33;
  else if (ColorName == "blue") return 34;
  else if (ColorName == "magenta") return 35;
  else if (ColorName == "cyan") return 36;
  else if (ColorName == "white") return 37;
  else {std::cout << "Color " << ColorName << " not recognized. Returning white." << std::endl; return 37;}
  return 0;
}

DEFINE_ART_MODULE(roiana::ROIAna)

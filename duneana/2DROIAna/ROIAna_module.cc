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
  geo = lar::providerFrom<geo::Geometry>();
  fNPlanes = geo->Nplanes();

  fLabels      = pset.get<std::vector<std::string>>("ParticleLabelVector");
  fInteraction = pset.get<std::vector<std::string>>("InteractionLabelVector");

  fROI_Peak  = pset.get<float>("ROIPeak", 20);
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
      FillMCInteractionTree(Parts[i], fInteraction, fLogLevel);
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

  ////////////////////////////////////////////////////////////
  //Store TrackIDs for each generator
  // MarleyChannels stores Channel+IDE of energy deposits from Marley produced signals
  std::vector<int> MarleyTrackIDs;
  std::vector<std::vector<int>> AllPartsTrackIDs(fLabels.size()); // list of trackIDs for each generator label in Parts
  std::map<int,sim::IDE> MarleyChannels;
  
  for ( size_t i = 0; i < fLabels.size(); i++){
    for (const auto& [TrackID, MCPart] : Parts[i]) AllPartsTrackIDs[i].push_back(TrackID);
  }

  //TODO move this into a function after creating ROI
  for (auto mySimChannel : channellist) {
    std::vector<sim::IDE> simIDEList = mySimChannel->TrackIDsAndEnergies(0, fNTicksPerWire);
    for (auto simIDE : simIDEList) {
      if (std::find(AllPartsTrackIDs[0].begin(), AllPartsTrackIDs[0].end(), simIDE.trackID) != AllPartsTrackIDs[0].end()) {
        if( fLogLevel >= 3 ) std::cout << "marley particle in channel: " << mySimChannel->Channel() << std::endl;
        MarleyChannels[mySimChannel->Channel()] = simIDE;
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

  //Creating ROI
  if( fLogLevel >= 3 ) std::cout << "starting ProcessROI" << std::endl;
  std::map<int,bool> ret;
  ret = ProcessROIWide(rawList);
  
  // ROIMetrics
  if( fLogLevel >= 3 ) std::cout << "starting ROIMetrics" << std::endl;
  ROIMetrics(ret);

  // ROIFilter
  //if( fLogLevel >= 3 ) std::cout << "starting ROIFilter" << std::endl;
  //ROIFilter(rawList, ret);

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
  //std::vector<short> ADC_zeros(adc_size, 0);

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
  fROITree = tfs->make<TTree>("MCROITree" , "MC ROI Tree (for energy depo per timetick per TrackID)" );
  fMCTruthTree = tfs->make<TTree>("MCTruthTree","contains particle counts from each generator");
  fInteractionTree = tfs->make<TTree>("MCInteraction","MC Interaction Tree");

  fROITree -> Branch("Run",      &run);
  fROITree -> Branch("Subrun",   &subrun);
  fROITree -> Branch("Event",    &event);
  fROITree -> Branch("PDG",                    &PDGROITree,   "PDG/I");       // Energy depo per TrackID, PDG
  fROITree -> Branch("TrackID",                &TrackIDROITree);              // Energy depo per TrackID, TrackID
  fROITree -> Branch("Generator",              &Generator);                   // Energy depo per TrackID, generator
  fROITree -> Branch("MCPartEnergy",           &MCParticleEnergy,  "Energy/F");  // Energy depo per TrackID, MC particle energy
  fROITree -> Branch("TrueEnergyDepositedU",    &TEnergyDepositedU);            // Energy depo per TrackID, energy deposited, u plane
  fROITree -> Branch("TrueChargeDepositedU",    &TChargeDepositedU);            // Energy depo per TrackID, charge induced, u plane
  fROITree -> Branch("TrueEnergyDepositedV",    &TEnergyDepositedV);            // Energy depo per TrackID, energy deposited, v plane
  fROITree -> Branch("TrueChargeDepositedV",    &TChargeDepositedV);            // Energy depo per TrackID, charge induced, v plane
  fROITree -> Branch("TrueEnergyDepositedX",    &TEnergyDepositedX);            // Energy depo per TrackID, energy deposited, x plane (collection)
  fROITree -> Branch("TrueChargeDepositedX",    &TChargeDepositedX);            // Energy depo per TrackID, charge induced, x plane (collection)
  fROITree -> Branch("TrueEnergyDepositedROIU", &TEnergyDepositedROIU);         // Energy depo per TrackID, energy in ROI, u plane
  fROITree -> Branch("TrueChargeDepositedROIU", &TChargeDepositedROIU);         // Energy depo per TrackID, charge in ROI, u plane
  fROITree -> Branch("TrueEnergyDepositedROIV", &TEnergyDepositedROIV);         // Energy depo per TrackID, energy in ROI, v plane
  fROITree -> Branch("TrueChargeDepositedROIV", &TChargeDepositedROIV);         // Energy depo per TrackID, charge in ROI, v plane
  fROITree -> Branch("TrueEnergyDepositedROIX", &TEnergyDepositedROIX);         // Energy depo per TrackID, energy in ROI, x plane (collection)
  fROITree -> Branch("TrueChargeDepositedROIX", &TChargeDepositedROIX);         // Energy depo per TrackID, charge in ROI, x plane (collection)

  fMCTruthTree -> Branch("Event",                &event,          "Event/I");  // Event number.
  fMCTruthTree -> Branch("Flag",                 &flag,           "Flag/I");   // Flag used to match truth with reco tree entries.
  fMCTruthTree -> Branch("TruthPart",            &TPart);                      // Number particles per generator.
 
  fInteractionTree -> Branch("Event",            &event,          "Event/I");  // Event number.
  fInteractionTree -> Branch("Flag",             &flag,           "Flag/I");   // Flag used to match truth with reco tree entries.
  fInteractionTree -> Branch("PDG",              &PDG,            "PDG/I");    // Main interacting particle PDG.
  fInteractionTree -> Branch("Energy",           &Energy,         "Energy/F"); // Main interacting particle energy [GeV^2].
  fInteractionTree -> Branch("Interaction",      &Interaction);                // Type of interaction.
  fInteractionTree -> Branch("Momentum",         &Momentum);                   // Main interacting particle momentum [GeV^2].
  fInteractionTree -> Branch("StartVertex",      &StartVertex);                // Main interacting particle start vertex [cm].
  fInteractionTree -> Branch("EndVertex",        &EndVertex);                  // Main interacting particle end vertex [cm].
  //fInteractionTree -> Branch("DaughterPDG",      &DaughterPDG);                // Main interacting particle daughter PDG.
  //fInteractionTree -> Branch("DaughterE",        &DaughterE);                  // Main interacting particle daughter energy [GeV^2].
  //fInteractionTree -> Branch("DaughterPx",       &DaughterPx);                 // Main interacting particle daughter momentum X [GeV^2].
  //fInteractionTree -> Branch("DaughterPy",       &DaughterPy);                 // Main interacting particle daughter momentum Y [GeV^2].
  //fInteractionTree -> Branch("DaughterPz",       &DaughterPz);                 // Main interacting particle daughter momentum Z [GeV^2].
  //fInteractionTree -> Branch("DaughterStartVx",  &DaughterStartVx);            // Main interacting particle daughter start vertex X [cm].
  //fInteractionTree -> Branch("DaughterStartVy",  &DaughterStartVy);            // Main interacting particle daughter start vertex Y [cm].
  //fInteractionTree -> Branch("DaughterStartVz",  &DaughterStartVz);            // Main interacting particle daughter start vertex Z [cm].
  //fInteractionTree -> Branch("DaughterEndVx",    &DaughterEndVx);              // Main interacting particle daughter end vertex X [cm].
  //fInteractionTree -> Branch("DaughterEndVy",    &DaughterEndVy);              // Main interacting particle daughter end vertex Y [cm].
  //fInteractionTree -> Branch("DaughterEndVz",    &DaughterEndVz);              // Main interacting particle daughter end vertex Z [cm].
  
  TrueEnergyDeposited = tfs->make<TH1F>( "TrueEnergyDeposited", "Energy; E (MeV)",100,0,fHistEnergyMax);
  TrueEnergyDepositedInROI =  tfs->make<TH1F>( "TrueEnergyDepositedInROI", "Energy; E (MeV)",100,0,fHistEnergyMax);
  TrueEnergyDepositedRatio =  tfs->make<TH1F>( "TrueEnergyDepositedRatio", "Energy; E (MeV)",100,0,fHistEnergyMax);

  TrueChargeDeposited = tfs->make<TH1F>( "TrueChargeDeposited", "Charge; Number of ionization electrons",100,0,fHistChargeMax);
  TrueChargeDepositedInROI =  tfs->make<TH1F>( "TrueChargeDepositedInROI", "Charge; Number of ionization electrons",100,0,fHistChargeMax);
  TrueChargeDepositedRatio =  tfs->make<TH1F>( "TrueChargeDepositedRatio", "Charge; Number of ionization electrons",100,0,fHistChargeMax);

  DataReductionRate = tfs->make<TH1F>( "DataReductionRate", "Data retention fraction; (channels in ROI)/(total channels)",100,0,1.02);

  MarleySignalSensitivity = tfs->make<TH1F>( "SignalSensitivity", "Marley event energy fraction in ROI; (energy in ROI)/(total energy)",100,0,1.02);
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

void roiana::ROIAna::ROIMetrics( std::map<int,bool> ret )
/*
  Fill ROI charge and energy histograms:
*/
{
  for( auto m: ch_w_sc)
  {
    auto channel = m.first;
    auto channel_reduced = channel%fNChanPerApa;
    //auto rawdigit = m.second.first;
    auto sim = m.second.second;
    //for( auto &tdcide: sim->TDCIDEMap() )
    //These energies and charges are per channel
    float energies=fECMin;
    float charges=fECMin;
    float energies_roi=fECMin;
    float charges_roi=fECMin;
    
    //loop over trackIDs for each channel
    for( auto &trackide: sim->TrackIDsAndEnergies(0, fNTicksPerWire) )
    {
      //These energies and charges are per channel
      //float energies=fECMin;
      //float charges=fECMin;
      //float energies_roi=fECMin;
      //float charges_roi=fECMin;
      
      //for( auto &ide: tdcide.second )
      //{
      int TrackID = abs(trackide.trackID);
      float energy = trackide.energy;
      float numElectrons = trackide.numElectrons;
      if( fLogLevel >= 3 ) std::cout<<"trackID: "<<TrackID<<std::endl;
      //std::cout<<"trackID: "<<trackide.trackID<<std::endl;
      //if (trackide.trackID==41562) std::cout<<"found trackID 41562 "<<std::endl;

      //fill energy by TrackID and wire plane
      //TrackIDEnergyMap[TrackID] += energy;
      //TrackIDChargeMap[TrackID] += numElectrons;

      if (channel_reduced <800) {
        TrackIDEnergyMapU[TrackID] += energy;
        TrackIDChargeMapU[TrackID] += numElectrons;
      } 
      else if (channel_reduced <1600) {
        TrackIDEnergyMapV[TrackID] += energy;
        TrackIDChargeMapV[TrackID] += numElectrons;
      }
      else if (channel_reduced <2560) {
        TrackIDEnergyMapX[TrackID] += energy;
        TrackIDChargeMapX[TrackID] += numElectrons;
      }
      else{std::cout<<"WARNING: channel number above 2560"<<std::endl;}

      //fill ROI energies by TrackID and wire plane
      if (ret[channel] && channel_reduced <800) {
        TrackIDEnergyMapROIU[TrackID] += energy;
        TrackIDChargeMapROIU[TrackID] += numElectrons;
      }
      else if (ret[channel] && channel_reduced <1600) {
        TrackIDEnergyMapROIV[TrackID] += energy;
        TrackIDChargeMapROIV[TrackID] += numElectrons;
      }
      else if (ret[channel]){
        TrackIDEnergyMapROIX[TrackID] += energy;
        TrackIDChargeMapROIX[TrackID] += numElectrons;
      }

      energies+=energy;
      charges+=numElectrons;
      if (ret[channel]) {
        energies_roi+=energy;
        charges_roi+=numElectrons;
      }
      //}
    }  
    if( fLogLevel >= 3 ) std::cout<<"FillHistogram: begin"<<std::endl;
    //Some function to fill histogram
    TrueEnergyDeposited->Fill(energies);
    TrueChargeDeposited->Fill(charges);
    if (ret[channel]) {
      TrueEnergyDepositedInROI->Fill(energies_roi);
      TrueChargeDepositedInROI->Fill(charges_roi);
    }
    if( fLogLevel >= 3 ) std::cout<<"FillHistogram: end"<<std::endl;
    
  }
  // True energy & charge ratio in ROI
  TrueEnergyDepositedRatio->Divide(TrueEnergyDepositedInROI, TrueEnergyDeposited);
  TrueChargeDepositedRatio->Divide(TrueChargeDepositedInROI, TrueChargeDeposited);
}

void roiana::ROIAna::ROIFilter( std::vector<art::Ptr<raw::RawDigit>>& rawList, std::map<int,bool> ret )
/*
 * Apply ROI Filter on RawDigit output:
*/
{
  const int adc_size = rawList.front()->NADC();
  const std::vector<short> ADC_zeros(adc_size, 0);

  for( auto& rawdigit_ptr : rawList )
  {
    auto rawdigit = *rawdigit_ptr;
    int channel = rawdigit.Channel();
    
    if (!ret[channel])
    {
      // Construct RawDigits with ADC_zeros
      //raw::RawDigit empty_rawdigit = raw::RawDigit(channel, fNTicksPerWire, ADC_zeros);
      //art::Ptr<raw::RawDigit> empty_rawdigit_ptr = &empty_rawdigit
      //empty_rawdigit_ptr->RawDigit(channel, fNTicksPerWire, ADC_zeros, rawdigit->Compression());
      //rawdigit_ptr = &empty_rawdigit;
    }
  }
}


void roiana::ROIAna::ROIEfficiencies( std::map<int,bool> ret, int n_channels, std::map<int,sim::IDE> MarleyChannels )
/*
 * Fill ROI performance histograms:
*/
{
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
  float MarleyEnergyTot = fECMin;
  float MarleyChargeTot = fECMin;
  float MarleyEnergyROI = fECMin;
  float MarleyChargeROI = fECMin;
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
  //std::cout << "Marley signal charge fraction in ROI: " << MarleyEnergyROI/MarleyEnergyTot << std::endl;
  MarleySignalSensitivity->Fill(MarleyEnergyROI/MarleyEnergyTot);

  //populate MCROITree, with info per trackID
  for ( size_t i = 0; i < fLabels.size(); i++){
    for (const auto& [TrackID, MCPart] : Parts[i]){
      PDGROITree =          MCPart.PdgCode();
      TrackIDROITree =      TrackID;
      Generator =           fLabels[i];
      MCParticleEnergy =    MCPart.E();
      TEnergyDepositedU =    TrackIDEnergyMapU[TrackID];
      TChargeDepositedU =    TrackIDChargeMapU[TrackID];
      TEnergyDepositedV =    TrackIDEnergyMapV[TrackID];
      TChargeDepositedV =    TrackIDChargeMapV[TrackID];
      TEnergyDepositedX =    TrackIDEnergyMapX[TrackID];
      TChargeDepositedX =    TrackIDChargeMapX[TrackID];
      TEnergyDepositedROIU = TrackIDEnergyMapROIU[TrackID];
      TChargeDepositedROIU = TrackIDChargeMapROIU[TrackID];
      TEnergyDepositedROIV = TrackIDEnergyMapROIV[TrackID];
      TChargeDepositedROIV = TrackIDChargeMapROIV[TrackID];
      TEnergyDepositedROIX = TrackIDEnergyMapROIX[TrackID];
      TChargeDepositedROIX = TrackIDChargeMapROIX[TrackID];

      fROITree -> Fill();
    }
  }
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
  //mf::LogInfo lheader("header");
  //std::string lheaderstr = "";
  // Make a copy of MCParticleList to be used for finding Daughter info
  std::map< int, simb::MCParticle> MCParticleListCopy = MCParticleList;
  //bool FoundInteraction;

  if (ProcessList.empty()){
    //lheaderstr = PrintInColor(lheaderstr,"\n-> No processes in the list!",GetColor("red"));
    //lheader << lheaderstr;
    return;  
  }
  
  for (size_t j = 0; j < ProcessList.size(); j++){
    //FoundInteraction = false;    
    for ( std::map<int,simb::MCParticle>::iterator mainiter = MCParticleList.begin(); mainiter != MCParticleList.end(); mainiter++ ){
      if ( mainiter->second.Process() != ProcessList[j] && mainiter->second.EndProcess() != ProcessList[j]){continue;}
      //if( fLogLevel >= 3 ) lheaderstr = lheaderstr+"\nFound a main interaction "+mainiter->second.EndProcess();
      //FoundInteraction = true;
      simb::MCParticle MCParticle = mainiter->second;
      Interaction =  MCParticle.EndProcess();
      PDG =          MCParticle.PdgCode();
      Energy =       MCParticle.E();
      Momentum =    {MCParticle.Px(),MCParticle.Py(),MCParticle.Pz()};
      StartVertex = {MCParticle.Vx(),MCParticle.Vy(),MCParticle.Vz()};
      EndVertex =   {MCParticle.EndX(),MCParticle.EndY(),MCParticle.EndZ()};
      
      //std::vector<int> DaughterList = {};
      //for (int i = 0; i < MCParticle.NumberDaughters(); i++){
      //  DaughterList.push_back(MCParticle.Daughter(i));
      //}
      // Print nice output with all the main interaction info
      /*if( fLogLevel >= 3 ) {
        lheaderstr = PrintInColor(lheaderstr,"\nMain interacting particle for process "+mainiter->second.Process()+": ",GetColor("magenta"));
        lheaderstr = PrintInColor(lheaderstr,"\nPDG ->\t"         + std::to_string(PDG),GetColor("cyan"));
        lheaderstr = PrintInColor(lheaderstr,"\nEnergy ->\t"      + std::to_string(Energy),GetColor("cyan"));
        lheaderstr = PrintInColor(lheaderstr,"\nMomentum ->\t"    + std::to_string(Momentum[0]) + " " + std::to_string(Momentum[1]) + " " + std::to_string(Momentum[2]),GetColor("cyan"));
        lheaderstr = PrintInColor(lheaderstr,"\nStartVertex ->\t" + std::to_string(StartVertex[0]) + " " + std::to_string(StartVertex[1]) + " " + std::to_string(StartVertex[2]),GetColor("cyan"));
        lheaderstr = PrintInColor(lheaderstr,"\nEndVertex ->\t"   + std::to_string(EndVertex[0]) + " " + std::to_string(EndVertex[1]) + " " + std::to_string(EndVertex[2]),GetColor("cyan"));
      }*/

      /* for ( std::map<int,simb::MCParticle>::iterator daughteriter = MCParticleListCopy.begin(); daughteriter != MCParticleListCopy.end(); daughteriter++ ){
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
      */
      fInteractionTree -> Fill();
    } // Loop over all particles in the map
    /*if( fLogLevel >= 3 ) {
      if (!FoundInteraction) lheaderstr = PrintInColor(lheaderstr,"\n-> No main interaction found for process "+ProcessList[j]+"!",GetColor("yellow"));
      else lheaderstr = PrintInColor(lheaderstr,"\n-> Filled MCINteraction Tree for process "+ProcessList[j]+"!",GetColor("green"));
    }*/
  } // Loop over all processes in the list
  //if( fLogLevel >= 3 ) { lheader << lheaderstr;}
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

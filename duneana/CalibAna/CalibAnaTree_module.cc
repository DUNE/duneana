////////////////////////////////////////////////////////////////////////
// Class:       CalibAnaTree
// Plugin Type: analyzer
// File:        CalibAnaTree_module.cc
//
// Based on TrackCaloSkimmer.h by Gray Putnam,
// https://github.com/SBNSoftware/sbncode/blob/develop/sbncode/Calibration/TrackCaloSkimmer.h.
// To be used as a builder for simpler ntuples for calibration analyses.
//
////////////////////////////////////////////////////////////////////////

#include "art/Utilities/make_tool.h"

#include "CalibAnaTree.h"
#include "G4processUtil.h"

#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"

// Global functions / data for fitting
const size_t MAX_N_FIT_DATA = 30;

static int N_FIT_DATA = -1;
static double FIT_RR[MAX_N_FIT_DATA];
static double FIT_DQDX[MAX_N_FIT_DATA];

void ConstResiduals(int &npar, double *g, double &result, double *par, int flag) {
  double ret = 0;

  double C = *par;

  for (int i = 0; i < N_FIT_DATA; i++) {
    double diff = FIT_DQDX[i] - C;
    ret += diff*diff;
  }

  result = sqrt(ret);
}

void ExpResiduals(int &npar, double *g, double &result, double *par, int flag) {
  double ret = 0;

  double A = par[0];
  double R = par[1];

  for (int i = 0; i < N_FIT_DATA; i++) {
    double diff = FIT_DQDX[i] - A*exp(-FIT_RR[i]/R);
    ret += diff*diff;
  }

  result = sqrt(ret);
}

dune::CalibAnaTree::~CalibAnaTree() {
  delete fTrack;
  delete fCluster;

}

dune::CalibAnaTree::CalibAnaTree(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fFitExp(2),
    fFitConst(1)
{
  // Grab config
  fLowEnergyClusterAnalysis     = p.get<bool>("LowEnergyClusterAnalysis", true);
  fRDTLabel                     = p.get<std::string>("RDTLabel","tpcrawdecoder:daq");

  fRadiusInt                    = p.get<float>("RadiusInt");
  fRadiusExt                    = p.get<float>("RadiusExt");

  fCoincidenceWd1_left          = p.get<float>("CoincidenceWindow1_left"  ,3); //in mus,
  fCoincidenceWd1_right         = p.get<float>("CoincidenceWindow1_right",3); //in mus,
  fCoincidenceWd2_left          = p.get<float>("CoincidenceWindow2_left"  ,3); //in mus,
  fCoincidenceWd2_right         = p.get<float>("CoincidenceWindow2_right",3); //in mus,

  fPitch                        = p.get<float>("Pitch",0.5); //cm
  fPitchMultiplier              = p.get<float>("PitchMultiplier",1.2); // 20% error    

  bIs3ViewsCoincidence          = p.get<bool>("Is3ViewsCoincidence");
  bIsPDVD                       = p.get<bool>("IsPDVD",false);
  bIsPDHD                       = p.get<bool>("IsPDHD",false);

  fNumberInitClusters           = p.get<int>("NumberInitClusters",25);
  fMaxSizeCluster               = p.get<float>("MaxSizeCluster",0);
  fMinSizeCluster               = p.get<float>("MinSizeCluster",0);

  if (fMaxSizeCluster == 0) fMaxSizeCluster = fRadiusInt*1.75;
  if (fMinSizeCluster == 0) fMinSizeCluster = fRadiusInt;
  
  fClusterSizeMulti             = p.get<float>("ClusterSizeMulti",1.2);
  fNumberConvStep               = p.get<int>("NumberConvStep",300);
  fCovering                     = p.get<float>("Covering",0.99);

  fTrackAnalysis                = p.get<bool>("TrackAnalysis" , true);
  fPFPproducer                  = p.get< art::InputTag > ("PFPproducer","pandoraGausCryo0");
  fT0producers                  = p.get< std::vector<art::InputTag> > ("T0producers", {"pandoraGausCryo0"} );
  fCALOproducer                 = p.get< art::InputTag > ("CALOproducer");
  fTRKproducer                  = p.get< art::InputTag > ("TRKproducer" );
  fTRKHMproducer                = p.get< art::InputTag > ("TRKHMproducer", "");
  fHITproducer                  = p.get< art::InputTag > ("HITproducer" );
  fG4producer                   = p.get< std::string > ("G4producer" );
  fSimChannelproducer           = p.get< std::string > ("SimChannelproducer" );
  fRequireT0                    = p.get<bool>("RequireT0", false);
  fDoTailFit                    = p.get<bool>("DoTailFit", true);
  fVerbose                      = p.get<bool>("Verbose", false);
  fSilenceMissingDataProducts   = p.get<bool>("SilenceMissingDataProducts", false);
  fHitRawDigitsTickCollectWidth = p.get<double>("HitRawDigitsTickCollectWidth", 50.);
  fHitRawDigitsWireCollectWidth = p.get<int>("HitRawDigitsWireCollectWidth", 5);
  fTailFitResidualRange         = p.get<double>("TailFitResidualRange", 5.);
  fFillTrackEndHits             = p.get<bool>("FillTrackEndHits", true);
  fTrackEndHitWireBox           = p.get<float>("TrackEndHitWireBox", 60); // 30 cm in the plane projection
  fTrackEndHitTimeBox           = p.get<float>("TrackEndHitTimeBox", 300); // 150 us, about 25 cm

  fRawDigitproducers            = p.get<std::vector<art::InputTag>>("RawDigitproducers", {}); 
 
  if (fTailFitResidualRange > 10.) {
    std::cout << "dune::CalibAnaTree: Bad tail fit residual range config :(" << fTailFitResidualRange << "). Fits will not be meaningful.\n";
  }
  fFillTrackEndHits = p.get<bool>("FillTrackEndHits", true);
  fTrackEndHitWireBox = p.get<float>("TrackEndHitWireBox", 60); // 30 cm in the plane projection
  fTrackEndHitTimeBox = p.get<float>("TrackEndHitTimeBox", 300); // 150 us, about 25 cm

  fRawDigitproducers = p.get<std::vector<art::InputTag>>("RawDigitproducers", {});

  std::vector<fhicl::ParameterSet> selection_tool_configs(p.get<std::vector<fhicl::ParameterSet>>("SelectionTools"), {});
  for (const fhicl::ParameterSet &p: selection_tool_configs) {
    fSelectionTools.push_back(art::make_tool<dune::ICATSelectionTool>(p));
  }

  // Setup meta info
  fMeta.iproc = -1;
  fMeta.ifile = -1;
  const char *process_str = std::getenv("PROCESS");
  if (process_str) {
    try {
      fMeta.iproc = std::stoi(process_str);
    }
    catch (...) {}
  }

  //geometry stuff
  fGeom    = &*art::ServiceHandle<geo::Geometry>();
  //Get detector Boundaries
  unsigned fNtpcs = fGeom->NTPC();

  for(unsigned t_tpc_id=0;t_tpc_id<fNtpcs;t_tpc_id++)
  {
    geo::TPCID tpcid{0, t_tpc_id};
    if(fgeoXmin > fGeom->TPC(tpcid).BoundingBox().MinX()) fgeoXmin = fGeom->TPC(tpcid).BoundingBox().MinX();
    if(fgeoXmax < fGeom->TPC(tpcid).BoundingBox().MaxX()) fgeoXmax = fGeom->TPC(tpcid).BoundingBox().MaxX();
    if(fgeoYmin > fGeom->TPC(tpcid).BoundingBox().MinY()) fgeoYmin = fGeom->TPC(tpcid).BoundingBox().MinY();
    if(fgeoYmax < fGeom->TPC(tpcid).BoundingBox().MaxY()) fgeoYmax = fGeom->TPC(tpcid).BoundingBox().MaxY();
    if(fgeoZmin > fGeom->TPC(tpcid).BoundingBox().MinZ()) fgeoZmin = fGeom->TPC(tpcid).BoundingBox().MinZ();
    if(fgeoZmax < fGeom->TPC(tpcid).BoundingBox().MaxZ()) fgeoZmax = fGeom->TPC(tpcid).BoundingBox().MaxZ();
  }
  if(fVerbose)
  {
    std::cout << " -- detector boundaries -- " << std::endl;
    std::cout << "  " << fgeoXmin << " < X < " << fgeoXmax << std::endl;
    std::cout << "  " << fgeoYmin << " < Y < " << fgeoYmax << std::endl;
    std::cout << "  " << fgeoZmin << " < Z < " << fgeoZmax << std::endl;
  }

  // Output stuff
  fTrack = new dune::TrackInfo();
  fCluster = new dune::ClusterInfo();

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("CalibAnaTree", "CalibAna Tree");
  //if (fTrackAnalysis)            
  btrk              = fTree->Branch("trk", &fTrack);
  //if (fLowEnergyClusterAnalysis)
  bLowEnergyCluster = fTree->Branch("LowEnergyClusters", &fCluster);
}

void dune::CalibAnaTree::analyze(art::Event const& e)
{
  unsigned evt = e.event();
  unsigned sub = e.subRun();
  unsigned run = e.run();
  if (fVerbose) {
    std::cout << "[CalibAnaTree::analyzeEvent] Run: " << run << ", SubRun: " << sub << ", Event: "<< evt << ", Is Data: " << e.isRealData() << std::endl;
  }

  fMeta.evt = evt;
  fMeta.subrun = sub;
  fMeta.run = run;
  fMeta.time = e.time().value();

  // Services
  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();
  const geo::WireReadoutGeom *wireReadout = &art::ServiceHandle<geo::WireReadout>()->Get();
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const dprop =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);

  float fTickTimeInMus  = detinfo::sampling_rate(clock_data)/1000; //caution sampling_rate(clock_data) is in ns
  float fCoincidenceWd1_right_inmus = fCoincidenceWd1_right/fTickTimeInMus; // in tt
  float fCoincidenceWd1_left_inmus  = fCoincidenceWd1_left/fTickTimeInMus; // in tt
  float fCoincidenceWd2_left_inmus  = fCoincidenceWd2_left/fTickTimeInMus; // in tt
  float fCoincidenceWd2_right_inmus = fCoincidenceWd2_right/fTickTimeInMus; // in tt

  float fElectronVelocity   = dprop.DriftVelocity();

  if(fVerbose){
    std::cout << " -- timing parameters -- "                             << std::endl;
    std::cout << "   sampling_rate         " << detinfo::sampling_rate(clock_data) << std::endl;
    std::cout << "   fTickTimeInMus:       " << fTickTimeInMus           << std::endl;
    std::cout << "   CoincidenceWd1:       " << fCoincidenceWd1_right_inmus + fCoincidenceWd1_left_inmus << std::endl;
    std::cout << "   CoincidenceWd2:       " << fCoincidenceWd2_right_inmus + fCoincidenceWd2_left_inmus << std::endl;
  }

  fMeta.ttinmus          = fTickTimeInMus;
  fMeta.electronvelocity = fElectronVelocity;

  // Identify which detector: can only detect either sbnd or icarus
  // Viktor Pec copied this for DUNE..  gdml seems ignored
  
  //std::string gdml = geometry->GDMLFile();
  //gdml = basename(gdml.c_str());
  //for(unsigned int i = 0; i <gdml.size(); ++i) gdml[i] = std::tolower(gdml[i]);

  EDet det = kNOTDEFINED;

  // Setup the volumes
  std::vector<std::vector<geo::BoxBoundedGeo>> TPCVols;
  std::vector<geo::BoxBoundedGeo> AVs;

  // First the TPC
  for (auto const &cryo: geometry->Iterate<geo::CryostatGeo>()) {
    std::vector<geo::BoxBoundedGeo> this_tpc_volumes;
    for (auto const& TPC : geometry->Iterate<geo::TPCGeo>(cryo.ID())) {
      this_tpc_volumes.push_back(TPC.ActiveBoundingBox());
    }
    TPCVols.push_back(std::move(this_tpc_volumes));
  }

  for (const std::vector<geo::BoxBoundedGeo> &tpcs: TPCVols) {
    double XMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinX() < rhs.MinX(); })->MinX();
    double YMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinY() < rhs.MinY(); })->MinY();
    double ZMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinZ() < rhs.MinZ(); })->MinZ();

    double XMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxX() < rhs.MaxX(); })->MaxX();
    double YMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxY() < rhs.MaxY(); })->MaxY();
    double ZMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxZ() < rhs.MaxZ(); })->MaxZ();

    AVs.emplace_back(XMin, XMax, YMin, YMax, ZMin, ZMax);
  }

  // Truth information
  std::vector<art::Ptr<simb::MCParticle>> mcparticles;
  std::vector<art::Ptr<sim::SimChannel>> simchannels;
  if (!e.isRealData())
  {
    if (fG4producer.size()) {
      art::ValidHandle<std::vector<simb::MCParticle>> mcparticle_handle = e.getValidHandle<std::vector<simb::MCParticle>>(fG4producer);
      art::fill_ptr_vector(mcparticles, mcparticle_handle);
    }

    if (fSimChannelproducer.size()) {
      art::ValidHandle<std::vector<sim::SimChannel>> simchannel_handle = e.getValidHandle<std::vector<sim::SimChannel>>(fSimChannelproducer);
      art::fill_ptr_vector(simchannels, simchannel_handle);
    }
  }

  if (fVerbose) {
      std::cout<<"Received "<<simchannels.size()<<" SimChannels for this event."<<std::endl;
  }

  if (fLowEnergyClusterAnalysis)
  {
    if (fVerbose)
    {
      std::cout<< "Starting Low Energy clusters analysis for event " << e.event() << std::endl;
    }

    std::vector<dune::ClusterInfo*> vCluster = SingleHitAnalysis(e, fRDTLabel, fG4producer, fHITproducer,
					       	          bIsPDVD, bIsPDHD, fCoincidenceWd1_left_inmus,
		      					  fCoincidenceWd1_right_inmus, fCoincidenceWd2_left_inmus,
       							  fCoincidenceWd2_right_inmus, bIs3ViewsCoincidence,
					       		  fPitch, fPitchMultiplier, fVerbose, fMinSizeCluster,
						      	  fMaxSizeCluster, fNumberInitClusters, fRadiusInt,
						      	  fRadiusExt, fgeoYmin, fgeoYmax, fgeoZmin, fgeoZmax,
							  fElectronVelocity , fTickTimeInMus);
    if (!vCluster.empty()) 
    {
      if (fVerbose) std::cout<< "There are " << vCluster.size() << " clusters in event " << fMeta.evt << std::endl;

      for( int i = 0 ; i < (int) vCluster.size() ; i++ )
      {
	fCluster = vCluster[i];
	fCluster->meta = fMeta;
	if(fVerbose) std::cout << "y = " << fCluster->y << "  z = " << fCluster->z << " charge = " << fCluster->ChargeCollection << std::endl;
        bLowEnergyCluster->Fill();
      }
    }

  }// end if Low Energy Clusters

  if (!fTrackAnalysis) return;

  // Reconstructed Information
  std::vector<art::Ptr<recob::PFParticle>> PFParticleList;
  try {
    art::ValidHandle<std::vector<recob::PFParticle>> pfparticles = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
    art::fill_ptr_vector(PFParticleList, pfparticles);
  }
  catch(...) {
    std::cout << "PFP's with tag: " << fPFPproducer << " not present.\n";
    // Real data may have missing products -- just ignore the event
    if (fSilenceMissingDataProducts) return;
    else throw;
  }

  // PFP-associated data
  std::vector<art::FindManyP<anab::T0>> fmT0;
  for (unsigned i_label = 0; i_label < fT0producers.size(); i_label++) {
    fmT0.emplace_back(PFParticleList, e, fT0producers[i_label]);
  }
  art::FindManyP<recob::SpacePoint> PFParticleSPs(PFParticleList, e, fPFPproducer);

  // Now we don't need to guard access to further data. If this is an empty event it should be caught by PFP's or Hit's
  art::ValidHandle<std::vector<recob::Track>> tracks = e.getValidHandle<std::vector<recob::Track>>(fTRKproducer);

  // Track - associated data
  art::FindManyP<recob::Track> fmTracks(PFParticleList, e, fTRKproducer);

  art::InputTag thm_label = fTRKHMproducer.empty() ? fTRKproducer : fTRKHMproducer;
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmtrkHits(tracks, e, thm_label);
  art::FindManyP<anab::Calorimetry> fmCalo(tracks, e, fCALOproducer);

  // Collect raw digits for saving hits
  std::vector<art::Ptr<raw::RawDigit>> rawdigitlist;
  for (const art::InputTag &t: fRawDigitproducers) {
    try {
      art::ValidHandle<std::vector<raw::RawDigit>> thisdigits = e.getValidHandle<std::vector<raw::RawDigit>>(t);
      art::fill_ptr_vector(rawdigitlist, thisdigits);
    }
    catch(...) {
      if (!fSilenceMissingDataProducts) throw;
      else {} // Allow Raw Digits to not be present
    }
  }

  // The raw digit list is not sorted, so make it into a map on the WireID
  std::map<geo::WireID, art::Ptr<raw::RawDigit>> rawdigits;
  for (const art::Ptr<raw::RawDigit> &d: rawdigitlist) {

    std::vector<geo::WireID> wids;
    // Handle bad channel ID
    try {
      wids = wireReadout->ChannelToWire(d->Channel());
    }
    catch(...) {
      //Fail if something...fails...
      throw cet::exception("CalibAnaTree_module") << "Raw digit " <<
            d->Channel() << " has no wires";
    }

    //We shouldn't ever have a raw digit NOT mapped to a wire
    if (wids.size() == 0) {
      throw cet::exception("CalibAnaTree_module") << "Raw digit " <<
            d->Channel() << " has no wires";
    }

    //A raw digit can come from multiple wire segments because some wires
    //are wrapped
    for (const auto & w : wids) {
      // Ignore wires that are already mapped
      if (rawdigits.count(w)) continue;
      rawdigits[w] = d;
    }
  }

  // Collect all hits
  art::ValidHandle<std::vector<recob::Hit>> allhit_handle = e.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
  std::vector<art::Ptr<recob::Hit>> allHits;
  art::fill_ptr_vector(allHits, allhit_handle);

  // And lookup the SP's
  art::FindManyP<recob::SpacePoint> allHitSPs(allHits, e, fPFPproducer);

  // Prep truth-to-reco-matching info
  //
  std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE*>>> id_to_ide_map;
  std::map<int, std::vector<art::Ptr<recob::Hit>>> id_to_truehit_map;
  const cheat::BackTrackerService *bt = NULL;

  if (simchannels.size()) {
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    id_to_ide_map = PrepSimChannels(simchannels, *wireReadout);
    id_to_truehit_map = PrepTrueHits(allHits, clock_data, *bt_serv.get());
    bt = bt_serv.get();
  }

  // service data

  // Build global track info
  std::vector<GlobalTrackInfo> track_infos;
  for (const recob::Track &t: *tracks) {
    track_infos.push_back({
	t.Start(), t.End(), t.StartDirection(), t.EndDirection(), t.ID()
      });
  }

  for (art::Ptr<recob::PFParticle> p_pfp: PFParticleList) {
    const recob::PFParticle &pfp = *p_pfp;

    const std::vector<art::Ptr<recob::Track>> thisTrack = fmTracks.at(p_pfp.key());
    if (thisTrack.size() != 1)
      continue;

    art::Ptr<recob::Track> trkPtr = thisTrack.at(0);

    std::vector<art::Ptr<anab::Calorimetry>> emptyCaloVector;
    const std::vector<art::Ptr<anab::Calorimetry>> &calo = fmCalo.isValid() ? fmCalo.at(trkPtr.key()) : emptyCaloVector;

    std::vector<art::Ptr<recob::Hit>> emptyHitVector;
    const std::vector<art::Ptr<recob::Hit>> &trkHits  = fmtrkHits.isValid() ? fmtrkHits.at(trkPtr.key()) : emptyHitVector;

    art::FindManyP<recob::SpacePoint> fmtrkHitSPs(trkHits, e, fPFPproducer);

    std::vector<const recob::TrackHitMeta*> emptyTHMVector;
    const std::vector<const recob::TrackHitMeta*> &trkHitMetas = fmtrkHits.isValid() ? fmtrkHits.data(trkPtr.key()) : emptyTHMVector;

    art::Ptr<recob::SpacePoint> nullSP;
    std::vector<art::Ptr<recob::SpacePoint>> trkHitSPs;
    if (fmtrkHitSPs.isValid()) {
      for (unsigned i_hit = 0; i_hit < trkHits.size(); i_hit++) {
        const std::vector<art::Ptr<recob::SpacePoint>> &h_sp = fmtrkHitSPs.at(i_hit);
        if (h_sp.size()) {
          trkHitSPs.push_back(h_sp.at(0));
        }
        else {
          trkHitSPs.push_back(nullSP);
        }
      }
    }

    int whicht0 = -1;
    float t0 = std::numeric_limits<float>::signaling_NaN();
    for (unsigned i_t0 = 0; i_t0 < fmT0.size(); i_t0++) {
      if (fmT0[i_t0].isValid() && fmT0[i_t0].at(p_pfp.key()).size()) {
        t0 = fmT0[i_t0].at(p_pfp.key()).at(0)->Time();
        whicht0 = i_t0;

        if (fVerbose) std::cout << "Track: " << trkPtr->ID() << " Has T0 (" << fT0producers[i_t0] << ")\n";

        break;
      }
    }

    if (fRequireT0 && whicht0 < 0) {
      continue;
    }

    if (fVerbose) std::cout << "Processing new track! ID: " << trkPtr->ID() << " time: " << t0 << std::endl;

    // Reset the track object
    *fTrack = dune::TrackInfo();

    // Reset other persistent info
    fSnippetCount.clear();
    fWiresToSave.clear();

    // Fill the track!
    FillTrack(*trkPtr, pfp, t0, trkHits, trkHitMetas, trkHitSPs, calo, rawdigits, track_infos, wireReadout, clock_data, bt, det);
    fTrack->whicht0 = whicht0;

    FillTrackDaughterRays(*trkPtr, pfp, PFParticleList, PFParticleSPs);

    if (fFillTrackEndHits) FillTrackEndHits(geometry, wireReadout, dprop, *trkPtr, allHits, allHitSPs);

    // Fill the truth information if configured
    if (simchannels.size()) FillTrackTruth(clock_data, trkHits, mcparticles, AVs, TPCVols, id_to_ide_map, id_to_truehit_map, dprop, geometry, wireReadout);

    // Save?
    bool select = false;
    if (!fSelectionTools.size()) select = true;

    // Take the OR of each selection tool
    int i_select = 0;
    for (const std::unique_ptr<dune::ICATSelectionTool> &t: fSelectionTools) {
      if (t->DoSelect(*fTrack)) {
        select = true;
        fTrack->selected = i_select;
        fTrack->nprescale = t->GetPrescale();
        break;
      }
      i_select ++;
    }

    // Save!
    if (select) {
      if (fVerbose) std::cout << "Track Selected! By tool: " << i_select << std::endl;
      btrk->Fill();
    }
  }
}

// helpers

dune::Vector3D ConvertTVector(const TVector3 &tv) {
  dune::Vector3D v;
  v.x = tv.X();
  v.y = tv.Y();
  v.z = tv.Z();

  return v;
}

// Turn a particle position to a space-charge induced position
geo::Point_t TrajectoryToWirePosition(const geo::Point_t &loc, const geo::Vector_t& driftdir) {
  auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  geo::Point_t ret = loc;

  // Returned X is the drift -- multiply by the drift direction to undo this
  int corr = driftdir.X();

  if (sce && sce->EnableSimSpatialSCE()) {
    geo::Vector_t offset = sce->GetPosOffsets(ret);

    ret.SetX(ret.X() + corr * offset.X());
    ret.SetY(ret.Y() + offset.Y());
    ret.SetZ(ret.Z() + offset.Z());
  }

  return ret;
}

// Turn a space-charge induced position to a trajectory Position
geo::Point_t WireToTrajectoryPosition(const geo::Point_t &loc, const geo::TPCID &tpc) {
  auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  geo::Point_t ret = loc;

  if (sce && sce->EnableSimSpatialSCE()) {
    geo::Vector_t offset = sce->GetCalPosOffsets(ret, tpc.TPC);

    ret.SetX(ret.X() + offset.X());
    ret.SetY(ret.Y() + offset.Y());
    ret.SetZ(ret.Z() + offset.Z());
  }

  return ret;

}

// Collect MCParticle information
dune::TrueParticle TrueParticleInfo(const simb::MCParticle &particle,
				    const std::vector<geo::BoxBoundedGeo> &active_volumes,
				    const std::vector<std::vector<geo::BoxBoundedGeo>> &tpc_volumes,
				    const std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE *>>> &id_to_ide_map,
				    const std::map<int, std::vector<art::Ptr<recob::Hit>>> &id_to_truehit_map,
				    const detinfo::DetectorPropertiesData &dprop,
				    const geo::GeometryCore *geo,
				    const geo::WireReadoutGeom *wireReadout) {

  std::vector<std::pair<geo::WireID, const sim::IDE *>> empty;
  const std::vector<std::pair<geo::WireID, const sim::IDE *>> &particle_ides = id_to_ide_map.count(particle.TrackId()) ? id_to_ide_map.at(particle.TrackId()) : empty;

  std::vector<art::Ptr<recob::Hit>> emptyHits;
  const std::vector<art::Ptr<recob::Hit>> &particle_hits = id_to_truehit_map.count(particle.TrackId()) ? id_to_truehit_map.at(particle.TrackId()) : emptyHits;

  dune::TrueParticle trueparticle;

  trueparticle.length = 0.;
  trueparticle.crosses_tpc = false;
  trueparticle.plane0VisE = 0.;
  trueparticle.plane1VisE = 0.;
  trueparticle.plane2VisE = 0.;
  trueparticle.plane0nhit = 0;
  trueparticle.plane1nhit = 0;
  trueparticle.plane2nhit = 0;
  for (auto const &ide_pair: particle_ides) {
    const geo::WireID &w = ide_pair.first;
    const sim::IDE *ide = ide_pair.second;

    if (w.Plane == 0) {
      trueparticle.plane0VisE += ide->energy / 1000. /* MeV -> GeV*/;
    }
    else if (w.Plane == 1) {
      trueparticle.plane1VisE += ide->energy / 1000. /* MeV -> GeV*/;
    }
    else if (w.Plane == 2) {
      trueparticle.plane2VisE += ide->energy / 1000. /* MeV -> GeV*/;
    }
  }

  for (const art::Ptr<recob::Hit> h: particle_hits) {
    const geo::WireID &w = h->WireID();

    if (w.Plane == 0) {
      trueparticle.plane0nhit ++;
    }
    else if (w.Plane == 1) {
      trueparticle.plane1nhit ++;
    }
    else if (w.Plane == 2) {
      trueparticle.plane2nhit ++;
    }

  }

  // if no trajectory points, then assume outside AV
  trueparticle.cont_tpc = particle.NumberTrajectoryPoints() > 0;
  trueparticle.contained = particle.NumberTrajectoryPoints() > 0;

  // Get the entry and exit points
  int entry_point = -1;

  int cryostat_index = -1;
  int tpc_index = -1;

  for (unsigned j = 0; j < particle.NumberTrajectoryPoints(); j++) {
    for (unsigned i = 0; i < active_volumes.size(); i++) {
      if (active_volumes.at(i).ContainsPosition(particle.Position(j).Vect())) {
	entry_point = j;
	cryostat_index = i;
	break;
      }
    }
    if (entry_point != -1) break;
  }

  int exit_point = -1;

  // now setup the cryostat the particle is in
  std::vector<geo::BoxBoundedGeo> volumes;
  if (entry_point >= 0) {
    volumes = tpc_volumes.at(cryostat_index);
    for (unsigned i = 0; i < volumes.size(); i++) {
      if (volumes[i].ContainsPosition(particle.Position(entry_point).Vect())) {
	tpc_index = i;
	trueparticle.cont_tpc = entry_point == 0;
	break;
      }
    }
    trueparticle.contained = entry_point == 0;
  }
  // if we couldn't find the initial point, set not contained
  else {
    trueparticle.contained = false;
  }

  if (tpc_index < 0) {
    trueparticle.cont_tpc = false;
  }

  // Get the length and determine if any point leaves the active volume
  // Use every trajectory point if possible
  if (entry_point >= 0) {
    // particle trajectory
    const simb::MCTrajectory &trajectory = particle.Trajectory();
    TVector3 pos = trajectory.Position(entry_point).Vect();
    for (unsigned i = entry_point+1; i < particle.NumberTrajectoryPoints(); i++) {
      TVector3 this_point = trajectory.Position(i).Vect();
      // get the exit point
      // update if particle is contained
      // check if particle has crossed TPC
      if (!trueparticle.crosses_tpc) {
	for (unsigned j = 0; j < volumes.size(); j++) {
	  if (volumes[j].ContainsPosition(this_point) && tpc_index >= 0 && j != ((unsigned)tpc_index)) {
	    trueparticle.crosses_tpc = true;
	    break;
	  }
	}
      }
      // check if particle has left tpc
      if (trueparticle.cont_tpc) {
	trueparticle.cont_tpc = volumes[tpc_index].ContainsPosition(this_point);
      }

      if (trueparticle.contained) {
	trueparticle.contained = active_volumes.at(cryostat_index).ContainsPosition(this_point);
      }

      trueparticle.length += (this_point - pos).Mag();

      if (!active_volumes.at(cryostat_index).ContainsPosition(this_point) && active_volumes.at(cryostat_index).ContainsPosition(pos)) {
	exit_point = i-1;
      }

      pos = trajectory.Position(i).Vect();
    }
  }
  if (exit_point < 0 && entry_point >= 0) {
    exit_point = particle.NumberTrajectoryPoints() - 1;
  }

  // other truth information
  trueparticle.pdg = particle.PdgCode();

  trueparticle.gen = ConvertTVector(particle.NumberTrajectoryPoints() ? particle.Position().Vect() : TVector3(-9999, -9999, -9999));
  trueparticle.genT = particle.NumberTrajectoryPoints() ? particle.Position().T() / 1000. /* ns -> us*/: -9999;
  trueparticle.genp = ConvertTVector(particle.NumberTrajectoryPoints() ? particle.Momentum().Vect(): TVector3(-9999, -9999, -9999));
  trueparticle.genE = particle.NumberTrajectoryPoints() ? particle.Momentum().E(): -9999;

  trueparticle.start = ConvertTVector((entry_point >= 0) ? particle.Position(entry_point).Vect(): TVector3(-9999, -9999, -9999));
  trueparticle.startT = (entry_point >= 0) ? particle.Position(entry_point).T() / 1000. /* ns-> us*/: -9999;
  trueparticle.end = ConvertTVector((exit_point >= 0) ? particle.Position(exit_point).Vect(): TVector3(-9999, -9999, -9999));
  trueparticle.endT = (exit_point >= 0) ? particle.Position(exit_point).T() / 1000. /* ns -> us */ : -9999;

  trueparticle.startp = ConvertTVector((entry_point >= 0) ? particle.Momentum(entry_point).Vect() : TVector3(-9999, -9999, -9999));
  trueparticle.startE = (entry_point >= 0) ? particle.Momentum(entry_point).E() : -9999.;
  trueparticle.endp = ConvertTVector((exit_point >= 0) ? particle.Momentum(exit_point).Vect() : TVector3(-9999, -9999, -9999));
  trueparticle.endE = (exit_point >= 0) ? particle.Momentum(exit_point).E() : -9999.;

  trueparticle.start_process = (int)GetG4ProcessID(particle.Process());
  trueparticle.end_process = (int)GetG4ProcessID(particle.EndProcess());

  trueparticle.G4ID = particle.TrackId();
  trueparticle.parent = particle.Mother();

  // Organize deposition info into per-wire true "Hits" -- key is the Channel Number
  std::map<unsigned, dune::TrueHit> truehits;

  for (auto const &ide_pair: particle_ides) {
    const geo::WireID &w = ide_pair.first;
    unsigned c = wireReadout->PlaneWireToChannel(w);
    const sim::IDE *ide = ide_pair.second;

    // Set stuff
    truehits[c].cryo = w.Cryostat;
    truehits[c].tpc = w.TPC;
    truehits[c].plane = w.Plane;
    truehits[c].wire = w.Wire;
    truehits[c].channel = c;

    // Average stuff using charge-weighting
    float old_elec = truehits[c].nelec;
    float new_elec = old_elec + ide->numElectrons;
    truehits[c].p.x = (truehits[c].p.x*old_elec + ide->x*ide->numElectrons) / new_elec;
    truehits[c].p.y = (truehits[c].p.y*old_elec + ide->y*ide->numElectrons) / new_elec;
    truehits[c].p.z = (truehits[c].p.z*old_elec + ide->z*ide->numElectrons) / new_elec;

    // Also get the position with space charge un-done
    geo::Point_t ide_p(ide->x, ide->y, ide->z);
    geo::Point_t ide_p_scecorr = WireToTrajectoryPosition(ide_p, w);

    truehits[c].p_scecorr.x = (truehits[c].p_scecorr.x*old_elec + ide_p_scecorr.x()*ide->numElectrons) / new_elec;
    truehits[c].p_scecorr.y = (truehits[c].p_scecorr.y*old_elec + ide_p_scecorr.y()*ide->numElectrons) / new_elec;
    truehits[c].p_scecorr.z = (truehits[c].p_scecorr.z*old_elec + ide_p_scecorr.z()*ide->numElectrons) / new_elec;

    // Sum stuff
    truehits[c].nelec += ide->numElectrons;
    truehits[c].e += ide->energy;
    truehits[c].ndep += 1;
  }

  // Compute widths
  for (auto const &ide_pair: particle_ides) {
    const geo::WireID &w = ide_pair.first;
    unsigned c = wireReadout->PlaneWireToChannel(w);
    const sim::IDE *ide = ide_pair.second;

    geo::Point_t ide_p(ide->x, ide->y, ide->z);
    geo::Point_t ide_p_scecorr = WireToTrajectoryPosition(ide_p, w);

    // Average stuff using charge-weighting
    float this_elec = ide->numElectrons;

    truehits[c].p_width.x += (ide_p.x() - truehits[c].p.x) * (ide_p.x() - truehits[c].p.x) * this_elec / truehits[c].nelec;
    truehits[c].p_width.y += (ide_p.y() - truehits[c].p.y) * (ide_p.y() - truehits[c].p.y) * this_elec / truehits[c].nelec;
    truehits[c].p_width.z += (ide_p.z() - truehits[c].p.z) * (ide_p.z() - truehits[c].p.z) * this_elec / truehits[c].nelec;

    truehits[c].p_scecorr_width.x += (ide_p_scecorr.x() - truehits[c].p_scecorr.x) * (ide_p_scecorr.x() - truehits[c].p_scecorr.x) * this_elec / truehits[c].nelec;
    truehits[c].p_scecorr_width.y += (ide_p_scecorr.y() - truehits[c].p_scecorr.y) * (ide_p_scecorr.y() - truehits[c].p_scecorr.y) * this_elec / truehits[c].nelec;
    truehits[c].p_scecorr_width.z += (ide_p_scecorr.z() - truehits[c].p_scecorr.z) * (ide_p_scecorr.z() - truehits[c].p_scecorr.z) * this_elec / truehits[c].nelec;
  }

  // Convert to vector
  std::vector<dune::TrueHit> truehits_v;
  for (auto const &p: truehits) {
    truehits_v.push_back(p.second);
  }

  // Compute the time of each hit
  for (dune::TrueHit &h: truehits_v) {
    h.time = dprop.ConvertXToTicks(h.p.x, h.plane, h.tpc, h.cryo);

    double xdrift = abs(h.p.x - wireReadout->Plane(geo::PlaneID(h.cryo, h.tpc, 0)).GetCenter().X());
    h.tdrift = xdrift / dprop.DriftVelocity();
  }

  // Compute the pitch of each hit and order it in the trajectory
  for (dune::TrueHit &h: truehits_v) {
    // Use the SCE-undone hit since this matches to the Trajectory
    TVector3 h_p(h.p_scecorr.x, h.p_scecorr.y, h.p_scecorr.z);

    TVector3 direction;
    float closest_dist = -1.;
    int traj_index = -1;
    for (unsigned i_traj = 0; i_traj < particle.NumberTrajectoryPoints(); i_traj++) {
      if (closest_dist < 0. || (particle.Position(i_traj).Vect() - h_p).Mag() < closest_dist) {
	direction = particle.Momentum(i_traj).Vect().Unit();
	closest_dist = (particle.Position(i_traj).Vect() - h_p).Mag();
	traj_index = i_traj;
      }
    }

    // If we got a direction, get the pitch
    if (closest_dist >= 0. && direction.Mag() > 1e-4) {
      geo::PlaneID plane(h.cryo, h.tpc, h.plane);
      geo::PlaneGeo const& planeGeo = wireReadout->Plane(plane);
      float angletovert = wireReadout->WireAngleToVertical(planeGeo.View(), plane) - 0.5*::util::pi<>();
      float cosgamma = abs(cos(angletovert) * direction.Z() + sin(angletovert) * direction.Y());
      float pitch = planeGeo.WirePitch() / cosgamma;
      h.pitch = pitch;
    }
    else {
      h.pitch = -1.;
    }
    // And the pitch induced by SCE
    if (closest_dist >= 0. && direction.Mag() > 1e-4) {
      geo::PlaneID plane(h.cryo, h.tpc, h.plane);
      geo::PlaneGeo const& planeGeo = wireReadout->Plane(plane);
      float angletovert = wireReadout->WireAngleToVertical(planeGeo.View(), plane) - 0.5*::util::pi<>();
      
      TVector3 loc_mdx_v = h_p - direction * (planeGeo.WirePitch() / 2.);
      TVector3 loc_pdx_v = h_p + direction * (planeGeo.WirePitch() / 2.);

      // Convert types for helper functions
      geo::Point_t loc_mdx(loc_mdx_v.X(), loc_mdx_v.Y(), loc_mdx_v.Z());
      geo::Point_t loc_pdx(loc_pdx_v.X(), loc_pdx_v.Y(), loc_pdx_v.Z());
      geo::Point_t h_p_point(h_p.X(), h_p.Y(), h_p.Z());

      auto const driftdir = geo->TPC(plane).DriftDir();
      loc_mdx = TrajectoryToWirePosition(loc_mdx, driftdir);
      loc_pdx = TrajectoryToWirePosition(loc_pdx, driftdir);

      // Direction at wires
      geo::Vector_t dir = (loc_pdx - loc_mdx) /  (loc_mdx - loc_pdx).r();

      // Pitch at wires
      double cosgamma = std::abs(std::sin(angletovert)*dir.Y() + std::cos(angletovert)*dir.Z());
      double pitch;
      if (cosgamma) {
	pitch = planeGeo.WirePitch()/cosgamma;
      }
      else {
	pitch = 0.;
      }

      // Now bring that back to the particle trajectory
      geo::Point_t loc_w = TrajectoryToWirePosition(h_p_point, driftdir);

      geo::Point_t locw_pdx_traj = WireToTrajectoryPosition(loc_w + pitch*dir, plane);
      geo::Point_t loc = WireToTrajectoryPosition(loc_w, plane);

      h.pitch_sce = (locw_pdx_traj - loc).R();
    }
    else {
      h.pitch_sce = -1.;
    }

    // And the trajectory location
    h.itraj = traj_index;

    // And the residual range of the hit
    h.rr = 0.;
    if (traj_index >= 0) {
      for (int i_traj = traj_index+1; i_traj < (int)particle.NumberTrajectoryPoints(); i_traj++) {
	h.rr += (particle.Position(i_traj).Vect() - particle.Position(i_traj-1).Vect()).Mag();
      }

      // Also account for the distance from the Hit point to the matched trajectory point
      double hit_distance_along_particle = (h_p - particle.Position(traj_index).Vect()).Dot(particle.Momentum(traj_index).Vect().Unit());
      h.rr += -hit_distance_along_particle;

    }
  }

  // Order the hits by their location along the trajectory, start to end
  std::sort(truehits_v.begin(), truehits_v.end(),
	    [](auto const &lhs, auto const &rhs) {
	      return lhs.itraj < rhs.itraj;
	    });

  // Save depositions into the True Particle
  for (dune::TrueHit &h: truehits_v) {
    if (h.plane == 0) {
      trueparticle.truehits0.push_back(h);
    }
    else if (h.plane == 1) {
      trueparticle.truehits1.push_back(h);
    }
    else if (h.plane == 2) {
      trueparticle.truehits2.push_back(h);
    }
  }

  // Save the true trajectory
  for (unsigned i_traj = 0; i_traj < particle.NumberTrajectoryPoints(); i_traj++) {
    // Get trajectory point
    TVector3 traj = particle.Position(i_traj).Vect();
    geo::Point_t traj_p(traj.X(), traj.Y(), traj.Z());

    // lookup TPC
    geo::TPCGeo const* tpc{nullptr}; // invalid by default
    for (auto const &cryo: geo->Iterate<geo::CryostatGeo>()) {
      for (auto const& TPC : geo->Iterate<geo::TPCGeo>(cryo.ID())) {
	if (TPC.ActiveBoundingBox().ContainsPosition(traj_p)) {
	  tpc = &TPC;
	  break;
	}
      }
      if (tpc && tpc->ID().isValid) break;
    }

    // add in space-charge-deflected position if applicable
    geo::Point_t traj_p_sce = tpc ? TrajectoryToWirePosition(traj_p, tpc->DriftDir()) : traj_p;

    dune::Vector3D traj_v;
    traj_v.x = traj_p.x();
    traj_v.y = traj_p.y();
    traj_v.z = traj_p.z();

    dune::Vector3D traj_v_sce;
    traj_v_sce.x = traj_p_sce.x();
    traj_v_sce.y = traj_p_sce.y();
    traj_v_sce.z = traj_p_sce.z();

    trueparticle.traj.push_back(traj_v);
    trueparticle.traj_sce.push_back(traj_v_sce);
  }

  return trueparticle;
}

void dune::CalibAnaTree::FillTrackEndHits(const geo::GeometryCore *geometry,
					  const geo::WireReadoutGeom *wireReadout,
					  const detinfo::DetectorPropertiesData &dprop,
					  const recob::Track &track,
					  const std::vector<art::Ptr<recob::Hit>> &allHits,
					  const art::FindManyP<recob::SpacePoint> &allHitSPs) {

  (void) dprop; // TODO: use??

  geo::TPCID tpc_end = geometry->FindTPCAtPosition(track.End());
  if (!tpc_end) return;

  geo::PlaneID plane_end(tpc_end, 2 /* collection */);

  float end_w = wireReadout->Plane(plane_end).WireCoordinate(track.End());

  float end_t = -1000.;
  float closest_wire_dist = -1.;
  // Get the hit closest to the end to get the end time
  for (const TrackHitInfo &h: fTrack->hits2) {
    if (h.oncalo && (closest_wire_dist < 0. || abs(h.h.wire - end_w) < closest_wire_dist)) {
      closest_wire_dist = abs(h.h.wire - end_w);
      end_t = h.h.time;
    }
  }

  for (const art::Ptr<recob::Hit> &hit: allHits) {
    geo::PlaneID h_p = hit->WireID();
    if (h_p != plane_end) continue;

    // Inside the box?
    float h_w = (float)hit->WireID().Wire;
    float h_t = hit->PeakTime();

    if (abs(h_w - end_w) < fTrackEndHitWireBox &&
	abs(h_t - end_t) < fTrackEndHitTimeBox) {

      // HitInfo to save
      dune::HitInfo hinfo;

      // information from the hit object
      hinfo.integral = hit->Integral();
      hinfo.sumadc = hit->ROISummedADC();
      hinfo.width = hit->RMS();
      hinfo.time = hit->PeakTime();
      hinfo.mult = hit->Multiplicity();
      hinfo.wire = hit->WireID().Wire;
      hinfo.plane = hit->WireID().Plane;
      hinfo.tpc = hit->WireID().TPC;
      hinfo.end = hit->EndTick();
      hinfo.start = hit->StartTick();
      hinfo.id = hit.key();

      const std::vector<art::Ptr<recob::SpacePoint>> &h_sp = allHitSPs.at(hit.key());
      if (h_sp.size()) {
	const recob::SpacePoint &sp = *h_sp[0];
	hinfo.sp.x = sp.position().x();
	hinfo.sp.y = sp.position().y();
	hinfo.sp.z = sp.position().z();

	hinfo.hasSP = true;
      }
      else {
	hinfo.hasSP = false;
      }

      fTrack->endhits.push_back(hinfo);

    }
  }

}

void dune::CalibAnaTree::FillTrackTruth(const detinfo::DetectorClocksData &clock_data,
					const std::vector<art::Ptr<recob::Hit>> &trkHits,
					const std::vector<art::Ptr<simb::MCParticle>> &mcparticles,
					const std::vector<geo::BoxBoundedGeo> &active_volumes,
					const std::vector<std::vector<geo::BoxBoundedGeo>> &tpc_volumes,
					const std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE*>>> id_to_ide_map,
					const std::map<int, std::vector<art::Ptr<recob::Hit>>> id_to_truehit_map,
					const detinfo::DetectorPropertiesData &dprop,
					const geo::GeometryCore *geo,
					const geo::WireReadoutGeom *wireReadout) {

  // Lookup the true-particle match -- use utils ported from SBNCode CAF
  std::vector<std::pair<int, float>> matches = AllTrueParticleIDEnergyMatches(clock_data, trkHits, true);
  float total_energy = TotalHitEnergy(clock_data, trkHits);

  if (fVerbose) {
      std::cout<<"Matched " << matches.size() << " true particles to the track. "<<std::endl;
      std::cout<<"Total track's hit true energy: "<<total_energy<<std::endl;
  }

  fTrack->truth.depE = total_energy / 1000. /* MeV -> GeV */;

  // sort highest energy match to lowest
  std::sort(matches.begin(), matches.end(),
	    [](const auto &a, const auto &b) {
	      return a.second > b.second;
	    }
	    );

  // Save the best match
  if (matches.size()) {
    std::pair<int, float> bestmatch = matches[0];

    fTrack->truth.pur = bestmatch.second / total_energy;

    for (const art::Ptr<simb::MCParticle> &p_mcp: mcparticles) {
      if (p_mcp->TrackId() == bestmatch.first) {
	if (fVerbose)
	    std::cout << "Matched! G4 Track ID: " << p_mcp->TrackId()
		      << " pdg: " << p_mcp->PdgCode()
		      << " process: " << p_mcp->EndProcess() << std::endl;
	fTrack->truth.p = TrueParticleInfo(*p_mcp, active_volumes, tpc_volumes, id_to_ide_map, id_to_truehit_map, dprop, geo, wireReadout);
	fTrack->truth.eff = fTrack->truth.depE / (fTrack->truth.p.plane0VisE + fTrack->truth.p.plane1VisE + fTrack->truth.p.plane2VisE);

	// Lookup any Michel
	for (const art::Ptr<simb::MCParticle> &d_mcp: mcparticles) {
	  if (d_mcp->Mother() == p_mcp->TrackId() && // correct parent
	      (d_mcp->Process() == "Decay" || d_mcp->Process() == "muMinusCaptureAtRest") && // correct process
	      abs(d_mcp->PdgCode()) == 11) { // correct PDG code

	    fTrack->truth.michel = TrueParticleInfo(*d_mcp, active_volumes, tpc_volumes, id_to_ide_map, id_to_truehit_map, dprop, geo, wireReadout);
	    break;
	  }
	}

	break;
      }
    }
  }
}


void dune::CalibAnaTree::FillTrackDaughterRays(const recob::Track &trk,
					       const recob::PFParticle &pfp,
					       const std::vector<art::Ptr<recob::PFParticle>> &PFParticleList,
					       const art::FindManyP<recob::SpacePoint> &PFParticleSPs) {

  for (unsigned d: pfp.Daughters()) {
    const recob::PFParticle &d_pfp = *PFParticleList[d];

    fTrack->daughter_pdg.push_back(d_pfp.PdgCode());

    unsigned nsp = 0;
    float min_distance = -1.;
    for (const art::Ptr<recob::SpacePoint> &sp: PFParticleSPs.at(d)) {
      if (min_distance < 0. || (sp->position() - trk.End()).r() < min_distance) {
	min_distance = (sp->position() - trk.End()).r();
      }
      nsp++;
    }

    fTrack->daughter_sp_toend_dist.push_back(min_distance);
    fTrack->daughter_nsp.push_back(nsp);
  }

}

void dune::CalibAnaTree::FillTrack(const recob::Track &track,
				   const recob::PFParticle &pfp, float t0,
				   const std::vector<art::Ptr<recob::Hit>> &hits,
				   const std::vector<const recob::TrackHitMeta*> &thms,
				   const std::vector<art::Ptr<recob::SpacePoint>> &sps,
				   const std::vector<art::Ptr<anab::Calorimetry>> &calo,
				   const std::map<geo::WireID, art::Ptr<raw::RawDigit>> &rawdigits,
				   const std::vector<GlobalTrackInfo> &tracks,
				   const geo::WireReadoutGeom *wireReadout,
				   const detinfo::DetectorClocksData &clock_data,
				   const cheat::BackTrackerService *bt_serv,
				   const dune::EDet det) {

  // Fill top level stuff
  fTrack->meta = fMeta;
  fTrack->t0 = t0;
  fTrack->id = track.ID();
  fTrack->clear_cosmic_muon = pfp.Parent() == recob::PFParticle::kPFParticlePrimary;

  fTrack->length = track.Length();
  fTrack->start.x = track.Start().X();
  fTrack->start.y = track.Start().Y();
  fTrack->start.z = track.Start().Z();
  fTrack->end.x = track.End().X();
  fTrack->end.y = track.End().Y();
  fTrack->end.z = track.End().Z();
  fTrack->dir.x = track.StartDirection().X();
  fTrack->dir.y = track.StartDirection().Y();
  fTrack->dir.z = track.StartDirection().Z();

  if (hits.size() > 0) {
    fTrack->cryostat = hits[0]->WireID().Cryostat;
  }

  // Fill each hit
  for (unsigned i_hit = 0; i_hit < hits.size(); i_hit++) {
    dune::TrackHitInfo hinfo = MakeHit(*hits[i_hit], hits[i_hit].key(), *thms[i_hit], track, sps[i_hit], calo, wireReadout, clock_data, bt_serv);
    if (hinfo.h.plane == 0) {
      fTrack->hits0.push_back(hinfo);
    }
    else if (hinfo.h.plane == 1) {
      fTrack->hits1.push_back(hinfo);
    }
    else if (hinfo.h.plane == 2) {
      fTrack->hits2.push_back(hinfo);
    }
  }

  // Save information on a fit to the end of the track
  if (fDoTailFit) DoTailFit();

  // Save the Wire ADC values we need to
  for (auto const &w_pair: fWiresToSave) {
    geo::WireID wire = w_pair.first;

    if (rawdigits.count(wire)) {
      const raw::RawDigit &thisdigit = *rawdigits.at(wire);
      int min_tick = std::max(0, w_pair.second.first);
      int max_tick = std::min((int)thisdigit.NADC(), w_pair.second.second);

      // collect the adcs
      std::vector<short> adcs;
      for (int t = min_tick; t < max_tick; t++) {
	adcs.push_back(thisdigit.ADC(t));
      }

      WireInfo winfo;
      winfo.wire = wire.Wire;
      winfo.plane = wire.Plane;
      winfo.tpc = wire.TPC;
      winfo.channel = wireReadout->PlaneWireToChannel(wire);
      winfo.tdc0 = min_tick;
      winfo.adcs = adcs;

      if (winfo.plane == 0) {
	fTrack->wires0.push_back(winfo);
      }
      else if (winfo.plane == 1) {
	fTrack->wires1.push_back(winfo);
      }
      else if (winfo.plane == 2) {
	fTrack->wires2.push_back(winfo);
      }
    }
  }
  // Sort the ADC values by wire
  std::sort(fTrack->wires0.begin(), fTrack->wires0.end(), [](auto const &lhs, auto const &rhs) {return lhs.wire < rhs.wire;});
  std::sort(fTrack->wires1.begin(), fTrack->wires1.end(), [](auto const &lhs, auto const &rhs) {return lhs.wire < rhs.wire;});
  std::sort(fTrack->wires2.begin(), fTrack->wires2.end(), [](auto const &lhs, auto const &rhs) {return lhs.wire < rhs.wire;});

  // get information on nearby tracks
  for (const GlobalTrackInfo &othr: tracks) {
    if (othr.ID == track.ID()) continue;

    if ((track.End() - othr.start).r() < 50. || (track.End() - othr.end).r() < 50.) {
      fTrack->tracks_near_end_dist.push_back(std::min((track.End() - othr.start).r(), (track.End() - othr.end).r()));
      fTrack->tracks_near_end_costh.push_back(
					      (track.End() - othr.start).r() < (track.End() - othr.end).r() ?
					      track.EndDirection().Dot(othr.dir) : track.EndDirection().Dot(othr.enddir));
    }
  }

  for (const GlobalTrackInfo &othr: tracks) {
    if (othr.ID == track.ID()) continue;

    if ((track.Start() - othr.start).r() < 50. || (track.Start() - othr.end).r() < 50.) {
      fTrack->tracks_near_start_dist.push_back(std::min((track.Start() - othr.start).r(), (track.Start() - othr.end).r()));
      fTrack->tracks_near_start_costh.push_back(
						(track.Start() - othr.start).r() < (track.Start() - othr.end).r() ?
						track.StartDirection().Dot(othr.dir) : track.StartDirection().Dot(othr.enddir));
    }
  }

}

void dune::CalibAnaTree::DoTailFit() {
  // Try fitting the constant and exponentials to the tail of dQ/dx v. RR on the collection plane
  std::vector<double> fit_rr;
  std::vector<double> fit_dqdx;

  for (const TrackHitInfo &h: fTrack->hits2) {
    if (h.oncalo && h.rr > 0. && h.rr < fTailFitResidualRange) {
      fit_rr.push_back(h.rr);
      fit_dqdx.push_back(h.dqdx);
    }
  }

  // TODO: should we throw an exception here??
  //
  // If there is too much data to fit in the array, throw exception
  if (fit_rr.size() > MAX_N_FIT_DATA) {
    throw cet::exception("dune::CalibAnaTree::DoTailFit: More fitting points required ("
			 + std::to_string(fit_rr.size()) + ") than available in fit array (" + std::to_string(MAX_N_FIT_DATA) + ").\n");
  }

  // Copy the fit data to the global array
  for (unsigned i = 0; i < fit_rr.size() && i < MAX_N_FIT_DATA; i++) {
    FIT_RR[i] = fit_rr[i];
    FIT_DQDX[i] = fit_dqdx[i];
  }
  N_FIT_DATA = std::min(fit_rr.size(), MAX_N_FIT_DATA);

  fTrack->n_fit_point = N_FIT_DATA;
  if (fTrack->n_fit_point > 2) { // need more points than params
    // Fit the Exponential
    fFitExp.SetFCN(ExpResiduals);
    fFitExp.SetParameter(0, "A", *std::max_element(fit_dqdx.begin(), fit_dqdx.end()), 200, 0, 5000);
    fFitExp.SetParameter(1, "R", 10., 0.5, 0, 1000);
    fFitExp.ExecuteCommand("MIGRAD", 0, 0);

    double A = fFitExp.GetParameter(0);
    double R = fFitExp.GetParameter(1);

    int nparam;
    double param[2] {A, R};
    double residuals = -1;
    ExpResiduals(nparam, NULL, residuals, param, 0);

    fTrack->exp_fit_A = A;
    fTrack->exp_fit_R = R;
    fTrack->exp_fit_residuals = residuals;

    // Fit the Constant
    fFitConst.SetFCN(ConstResiduals);
    fFitConst.SetParameter(0, "C", std::accumulate(fit_dqdx.begin(), fit_dqdx.end(), 0.), 200, 0, 5000);
    fFitConst.ExecuteCommand("MIGRAD", 0, 0);

    double C = fFitConst.GetParameter(0);

    double cresiduals = -1;
    ConstResiduals(nparam, NULL, cresiduals, &C, 0);

    fTrack->const_fit_C = C;
    fTrack->const_fit_residuals = cresiduals;
  }
}


dune::TrackHitInfo dune::CalibAnaTree::MakeHit(const recob::Hit &hit,
					       unsigned hkey,
					       const recob::TrackHitMeta &thm,
					       const recob::Track &trk,
					       const art::Ptr<recob::SpacePoint> &sp,
					       const std::vector<art::Ptr<anab::Calorimetry>> &calo,
					       const geo::WireReadoutGeom *wireReadout,
					       const detinfo::DetectorClocksData &dclock,
					       const cheat::BackTrackerService *bt_serv) {

  // TrackHitInfo to save
  dune::TrackHitInfo hinfo;

  // information from the hit object
  hinfo.h.integral = hit.Integral();
  hinfo.h.sumadc = hit.ROISummedADC();
  hinfo.h.width = hit.RMS();
  hinfo.h.time = hit.PeakTime();
  hinfo.h.mult = hit.Multiplicity();
  hinfo.h.wire = hit.WireID().Wire;
  hinfo.h.plane = hit.WireID().Plane;
  hinfo.h.channel = wireReadout->PlaneWireToChannel(hit.WireID());
  hinfo.h.tpc = hit.WireID().TPC;
  hinfo.h.end = hit.EndTick();
  hinfo.h.start = hit.StartTick();
  hinfo.h.id = (int)hkey;

  // Do back-tracking on each hit
  if (bt_serv) {
    // The default BackTracking function goes from (peak - width, peak + width).
    //
    // This time range does not match well hits with a non-Gaussian shape where
    // the Gaussian-fit-width does not replicate the width of the pulse.
    //
    // Instead, we use the Hit (start, end) time range. This is also more relevant
    // for (e.g.) the SummedADC charge extraction method.
    //
    // Don't use this:
    // std::vector<sim::TrackIDE> ides = bt_serv->HitToTrackIDEs(dclock, hit);
    //
    // Use this:
    std::vector<sim::TrackIDE> ides = bt_serv->ChannelToTrackIDEs(dclock, hit.Channel(), hit.StartTick(), hit.EndTick());

    hinfo.h.truth.e = 0.;
    hinfo.h.truth.nelec = 0.;

    for (const sim::TrackIDE &ide: ides) {
      hinfo.h.truth.e += ide.energy;
      hinfo.h.truth.nelec += ide.numElectrons;
    }
  }
  else {
    hinfo.h.truth.e = -1.;
    hinfo.h.truth.nelec = -1.;
  }

  // look up the snippet
  dune::CalibAnaTree::Snippet snippet {hit.WireID(), hit.StartTick(), hit.EndTick()};
  if (!fSnippetCount.count(snippet)) {
    fSnippetCount[snippet] = 0;
    hinfo.i_snippet = 0;
  }
  else {
    fSnippetCount[snippet] ++;
    hinfo.i_snippet = fSnippetCount[snippet];
  }

  // Which wires to save
  int min_tick = (int)std::floor(hit.PeakTime() - fHitRawDigitsTickCollectWidth);
  int max_tick = (int)std::ceil(hit.PeakTime() + fHitRawDigitsTickCollectWidth);
  for (int wire = hinfo.h.wire - fHitRawDigitsWireCollectWidth; wire <= hinfo.h.wire + fHitRawDigitsWireCollectWidth; wire++) {
    geo::WireID w(hit.WireID(), wire);

    if (fWiresToSave.count(w)) {
      fWiresToSave.at(w).first = std::min(fWiresToSave.at(w).first, min_tick);
      fWiresToSave.at(w).second = std::max(fWiresToSave.at(w).second, max_tick);
    }
    else {
      fWiresToSave[w] = {min_tick, max_tick};
    }
  }

  // Information from the TrackHitMeta
  bool badhit = (thm.Index() == std::numeric_limits<unsigned int>::max()) ||
    (!trk.HasValidPoint(thm.Index()));

  hinfo.ontraj = !badhit;

  // Save trajectory information if we can
  if (!badhit) {
    geo::Point_t loc = trk.LocationAtPoint(thm.Index());
    hinfo.tp.x = loc.X();
    hinfo.tp.y = loc.Y();
    hinfo.tp.z = loc.Z();

    geo::Vector_t dir = trk.DirectionAtPoint(thm.Index());
    hinfo.dir.x = dir.X();
    hinfo.dir.y = dir.Y();
    hinfo.dir.z = dir.Z();

    // And determine if the Hit is on a Calorimetry object
    for (const art::Ptr<anab::Calorimetry> &c: calo) {
      if (c->PlaneID().Plane != hinfo.h.plane) continue;

      // Found the plane! Now find the hit:
      for (unsigned i_calo = 0; i_calo < c->dQdx().size(); i_calo++) {
	if (c->TpIndices()[i_calo] == hkey) { // "TpIndices" match to the hit key
	  // Fill the calo information associated with the hit
	  hinfo.oncalo = true;
	  hinfo.pitch = c->TrkPitchVec()[i_calo];
	  hinfo.dqdx = c->dQdx()[i_calo];
	  hinfo.rr = c->ResidualRange()[i_calo];
	  break;
	}
      }
      break;
    }
  }

  // Save SpacePoint information
  if (sp) {
    hinfo.h.sp.x = sp->position().x();
    hinfo.h.sp.y = sp->position().y();
    hinfo.h.sp.z = sp->position().z();

    hinfo.h.hasSP = true;
  }
  else {
    hinfo.h.hasSP = false;
  }

  return hinfo;
}



/// Code taken from sbncode/CAFMaker/FillTrue.cxx
std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE*>>> dune::CalibAnaTree::PrepSimChannels(const std::vector<art::Ptr<sim::SimChannel>> &simchannels,
												   const geo::WireReadoutGeom &geo)
/// Creates map of WireID and sim::IDE by backtracked track ID
{
  std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE*>>> ret;

  for (const art::Ptr<sim::SimChannel> sc : simchannels) {
    // Lookup the wire of this channel
    raw::ChannelID_t channel = sc->Channel();
    std::vector<geo::WireID> maybewire = geo.ChannelToWire(channel);
    geo::WireID thisWire; // Default constructor makes invalid wire
    if (maybewire.size()) thisWire = maybewire[0];

    for (const auto &item : sc->TDCIDEMap()) {
      for (const sim::IDE &ide: item.second) {
	// indexing initializes empty vector
	ret[abs(ide.trackID)].push_back({thisWire, &ide});
      }
    }
  }

  if (fVerbose) {
      std::cout<<"Prepared IDEs for "<<ret.size()<<" G4 tracks."<<std::endl;
  }
  return ret;
}

std::map<int, std::vector<art::Ptr<recob::Hit>>> dune::CalibAnaTree::PrepTrueHits(const std::vector<art::Ptr<recob::Hit>> &allHits,
							      const detinfo::DetectorClocksData &clock_data,
							      const cheat::BackTrackerService &backtracker)
{
  std::map<int, std::vector<art::Ptr<recob::Hit>>> ret;
  for (const art::Ptr<recob::Hit> h: allHits) {
    for (int ID: backtracker.HitToTrackIds(clock_data, *h)) {
      ret[abs(ID)].push_back(h);
    }
  }
  return ret;
}


std::vector<std::pair<int, float>> dune::CalibAnaTree::AllTrueParticleIDEnergyMatches(const detinfo::DetectorClocksData &clock_data, const std::vector<art::Ptr<recob::Hit> >& hits, bool rollup_unsaved_ids) {
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  std::map<int, float> trackIDToEDepMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;
    std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(clock_data, hit);
    for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
      int id = trackIDs[idIt].trackID;
      if (rollup_unsaved_ids) id = std::abs(id);
      id = GetShowerPrimary(id);
      trackIDToEDepMap[id] += trackIDs[idIt].energy;
    }
  }

  std::vector<std::pair<int, float>> ret;
  for (auto const &pair: trackIDToEDepMap) {
    ret.push_back(pair);
  }
  return ret;
}

int dune::CalibAnaTree::GetShowerPrimary(const int g4ID)
{
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
  const sim::ParticleList& particles = particleInventory->ParticleList();
  const sim::ParticleList::const_iterator part_iter = particles.find(g4ID);
  if(part_iter == particles.end()) return g4ID;

  auto temp_iter = part_iter;
  int primary_id = part_iter->second->TrackId();

  while (std::abs(temp_iter->second->PdgCode()) == 11 || temp_iter->second->PdgCode() == 22)
    {
      primary_id = temp_iter->second->TrackId();
      temp_iter = particles.find(temp_iter->second->Mother());
      if(temp_iter == particles.end()) break;
    }

  return primary_id;
}

float dune::CalibAnaTree::TotalHitEnergy(const detinfo::DetectorClocksData &clock_data, const std::vector<art::Ptr<recob::Hit> >& hits) {
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  float ret = 0.;

  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;
    std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(clock_data, hit);
    for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
      ret += trackIDs[idIt].energy;
    }
  }
  return ret;
}


                   //////////////////////////////
		   //                          //
		   //     CLUSTER FUNCTION     //
		   //                          //
		   //////////////////////////////



bool dune::CalibAnaTree::Inside( int k , std::list<int> list){
  return (std::find(list.begin(), list.end(), k) != list.end());
}

bool dune::CalibAnaTree::AllSame( std::vector<int> v)
{
  if (v.size() == 0) return false;
  bool allsame = true;
  for( int i = 0; i < (int) v.size() ; i++)
  {
    if( v[i] == v[0] ) continue;
    allsame = false;
    break;
  }
  return allsame;
}
int dune::CalibAnaTree::NearOrFar( bool IsPDVD , bool IsPDHD , const recob::Hit & hit)
{
  if (bIsPDHD)
  {
    if ( (hit.WireID().TPC == 2)|| (hit.WireID().TPC == 6) || (hit.WireID().TPC == 3) || (hit.WireID().TPC == 7) ) return  -1; // 3 and 7 are dumie TPC
    if ( (hit.WireID().TPC == 1)|| (hit.WireID().TPC == 5) || (hit.WireID().TPC == 0) || (hit.WireID().TPC == 4) ) return  1;  // 0 and 4 are dumie TPC
  }
  if (bIsPDVD)
  {
    if (hit.WireID().TPC <= 7 ) return -1;
    if (hit.WireID().TPC  > 7 ) return 1;
  }
  return -999;
}


void dune::CalibAnaTree::GetListOfTimeCoincidenceHit( bool IsPDVD , bool IsPDHD , art::Event const & ev, art::InputTag HitLabel, float const CoincidenceWd1_l , float const CoincidenceWd1_r, float const CoincidenceWd2_l , float const CoincidenceWd2_r, const recob::Hit & HitCol, 
                                                                                  std::list<geo::WireID> & WireInd1,
                                                                                  std::list<geo::WireID> & WireInd2,
                                                                                  std::list<int>   & ChannelInd1,
                                                                                  std::list<int>   & ChannelInd2,
                                                                                  std::list<float> & EInd1,
                                                                                  std::list<float> & EInd2,
                                                                                  std::list<float> & PTInd1,
                                                                                  std::list<float> & PTInd2,
										  std::list<float> & PAInd1,
                                                                                  std::list<float> & PAInd2)
{
  auto const hitlist = ev.getValidHandle<std::vector<recob::Hit>>(HitLabel);

  recob::Hit hit = hitlist->at(0);

  float PeakTimeCol    = HitCol.PeakTime();
  //float RMSPeakTimeCol = HitCol.RMS();

  float EndTime1   = PeakTimeCol + CoincidenceWd1_r;
  float EndTime2   = PeakTimeCol + CoincidenceWd2_r;
  float StartTime1 = PeakTimeCol - CoincidenceWd1_l;
  float StartTime2 = PeakTimeCol - CoincidenceWd2_l;

  float PeakTime = -999;
  int   Plane    = -999;
  int NoFCol = NearOrFar(IsPDVD,IsPDHD,HitCol);
  int NoF = -4;

  for (int i=0, sz=hitlist->size(); i!=sz; ++i)
  { 
    hit   = hitlist->at(i);
    Plane = hit.WireID().Plane;
    if (Plane == 2) continue;

    NoF = NearOrFar(IsPDVD,IsPDHD,hit);
    if (NoF != NoFCol) continue;

    PeakTime = hit.PeakTime();
    if (Plane == 0)
    {
      if ((PeakTime < StartTime1)||(PeakTime > EndTime1)) continue;

      WireInd1.push_back(hit.WireID());
      ChannelInd1.push_back(hit.Channel());
      //EInd1.push_back(hit.ROISummedADC());
      EInd1.push_back(hit.Integral());
      PTInd1.push_back(PeakTime);
      PAInd1.push_back(hit.PeakAmplitude());
      continue;
    }
    if (Plane == 1)
    {
      if ((PeakTime < StartTime2)||(PeakTime > EndTime2)) continue;

      WireInd2.push_back(hit.WireID());
      ChannelInd2.push_back(hit.Channel());
      //EInd2.push_back(hit.ROISummedADC());
      EInd2.push_back(hit.Integral());
      PTInd2.push_back(PeakTime);
      PAInd2.push_back(hit.PeakAmplitude());
    }
  }
}

bool dune::CalibAnaTree::IntersectOutsideOfTPC( float Ymin , float Ymax , float Zmin , float Zmax ,
		                                double ChInd_start_y , double ChInd_start_z ,
					       	double ChInd_end_y ,double ChInd_end_z ,
						double ChCol_start_y , double ChCol_start_z ,
					    	double ChCol_end_y , double ChCol_end_z ,
		     				double& y , double& z )
{
  // Equation from http://en.wikipedia.org/wiki/Line%E2%80%93line_intersection

  double const denom = (ChInd_start_y - ChInd_end_y) * (ChCol_start_z - ChCol_end_z) - (ChInd_start_z - ChInd_end_z) * (ChCol_start_y - ChCol_end_y);

  if (denom == 0) return false;

  double const A = (ChInd_start_y * ChInd_end_z - ChInd_start_z * ChInd_end_y) / denom;
  double const B = (ChCol_start_y * ChCol_end_z - ChCol_start_z * ChCol_end_y) / denom;

  y = (ChCol_start_y - ChCol_end_y) * A - (ChInd_start_y - ChInd_end_y) * B;
  z = (ChCol_start_z - ChCol_end_z) * A - (ChInd_start_z - ChInd_end_z) * B;

  bool drap = ( y > Ymin ) && ( y < Ymax ) && ( z > Zmin ) && ( z < Zmax ) ;
  
  return drap;
}

void dune::CalibAnaTree::GetListOfCrossingChannel(   bool IsPDVD , bool IsPDHD , float Ymin , float Ymax , float Zmin , float Zmax ,
		                                    geo::WireID & WireCol , std::list<geo::WireID> & WireInd1 , std::list<geo::WireID> & WireInd2 , 
						    std::list<int>  & ChInd1 , std::list<float> & EInd1 , std::list<float> & YInd1 , std::list<float> & ZInd1 , 
						    std::list<int>  & ChIntersectInd1 , std::list<float> & EIntersectInd1 ,
                                                    std::list<int>  & ChInd2 , std::list<float> & EInd2 , std::list<float> & YInd2 , std::list<float> & ZInd2 , 
						    std::list<int>  & ChIntersectInd2 , std::list<float> & EIntersectInd2)
{
  geo::Point_t point = geo::Point_t(-999,-999,-999);
  bool drap;

  double y = -999. , z = -999.;
  //auto const wcol = fGeom->WireEndPoints(WireCol);
  auto const [wcolstart, wcolend] = fWireReadout.WireEndPoints(WireCol);

  ChIntersectInd1.clear();
  EIntersectInd1.clear();
  ChIntersectInd2.clear();
  EIntersectInd2.clear();
  
  std::list<int>::iterator ch1  = ChInd1.begin();
  std::list<float>::iterator e1   = EInd1.begin();

  for (auto const elementInd1 : WireInd1)
  {
    if (WireCol.TPC != elementInd1.TPC)
    {
    
      if (bIsPDVD)
      {  
        //auto const wind1 = fGeom->WireEndPoints(elementInd1);
        auto const [wind1start, wind1end] = fWireReadout.WireEndPoints(elementInd1);
        bool flag = IntersectOutsideOfTPC( Ymin , Ymax , Zmin ,Zmax , wind1start.Y() , wind1start.Z() , wind1end.Y() , wind1end.Z() , wcolstart.Y() , wcolstart.Z() , wcolend.Y() , wcolend.Z() , y , z);
        if (flag)
        {
	  YInd1.push_back(y);
          ZInd1.push_back(z);
          EIntersectInd1.push_back(*e1);
	  ChIntersectInd1.push_back(*ch1);
          ++e1; 
	  ++ch1;
	  continue;
        }
      }
      ++e1;
      ++ch1;
      continue ;
    }

    //drap = fGeom->WireIDsIntersect( WireCol , elementInd1 , point);
    drap = fWireReadout.WireIDsIntersect( WireCol , elementInd1 , point);
    if ( drap )
    {
      YInd1.push_back(point.Y());
      ZInd1.push_back(point.Z());
      ChIntersectInd1.push_back(*ch1);
      EIntersectInd1.push_back(*e1);
    }
    ++e1;
    ++ch1;
  }
  std::list<int>::iterator ch2  = ChInd2.begin();
  std::list<float>::iterator e2   = EInd2.begin();

  for (auto const elementInd2 : WireInd2)
  { 
    if (WireCol.TPC != elementInd2.TPC ) 
    { 

      if (bIsPDVD)
      {
        //auto const wind2 = fGeom->WireEndPoints(elementInd2);
        auto const [wind2start, wind2end] = fWireReadout.WireEndPoints(elementInd2);
        bool flag = IntersectOutsideOfTPC( Ymin , Ymax , Zmin ,Zmax , wind2start.Y() , wind2start.Z() , wind2end.Y() , wind2end.Z() , wcolstart.Y() , wcolstart.Z() , wcolend.Y() , wcolend.Z() , y , z);
        if (flag)
        {
          YInd2.push_back(y);
          ZInd2.push_back(z);
          ChIntersectInd2.push_back(*ch2);
	  EIntersectInd2.push_back(*e2);
	  ++e2;
          ++ch2;
          continue;
        }
      }
      ++e2;
      ++ch2; 
      continue ;
    }
    //drap = fGeom->WireIDsIntersect( WireCol , elementInd2 , point);
    drap = fWireReadout.WireIDsIntersect( WireCol , elementInd2 , point);
    if ( drap ) 
    { 
      YInd2.push_back(point.Y());
      ZInd2.push_back(point.Z());
      ChIntersectInd2.push_back(*ch2);
      EIntersectInd2.push_back(*e2);
    }
    ++e2;
    ++ch2;
  }
}

void dune::CalibAnaTree::GetListOf3ViewsPoint( float pitch , float alpha , 
                                        std::list<int> & ChIntersectInd1 , std::list<float> YInd1 , std::list<float> ZInd1 , std::list<float> EIntersectInd1, 
                                        std::list<int> & ChIntersectInd2 , std::list<float> YInd2 , std::list<float> ZInd2 , std::list<float> EIntersectInd2 , 
                                        std::list<float> & listYSP       , std::list<float> & listZSP , 
                                        std::list<float> & listEind1SP   , std::list<float> & listEind2SP , 
                                        std::list<int> & listCh1SP       , std::list<int> & listCh2SP)
{

  std::list<int>::iterator  ch1t = ChIntersectInd1.begin();
  std::list<float>::iterator z1t = ZInd1.begin();
  std::list<float>::iterator e1t = EIntersectInd1.begin();

  float dy, dz, dr;

  for( auto const yind1 : YInd1)
  {
    std::list<int>::iterator  ch2t = ChIntersectInd2.begin();
    std::list<float>::iterator z2t = ZInd2.begin();
    std::list<float>::iterator e2t = EIntersectInd2.begin();

    for ( auto const yind2 : YInd2)
    {
      dy = yind1 - yind2;
      dz = *z1t - *z2t  ;
      dr = TMath::Sqrt( dy*dy + dz*dz );

      if ( dr <= pitch*alpha )
      {
        float y = ( (*e1t)*(yind1) + (*e2t)*(yind2) )/( *e1t + *e2t );
        float z = ( (*e1t)*(*z1t) + (*e2t)*(*z2t) )/( *e1t + *e2t );

        listYSP.push_back( y );     
        listZSP.push_back( z );
        listEind1SP.push_back( *e1t );
        listEind2SP.push_back( *e2t );
        listCh1SP.push_back( *ch1t );
        listCh2SP.push_back( *ch2t );
      }
      ++e2t;
      ++z2t;
      ++ch2t;
    }
    ++e1t;
    ++z1t;
    ++ch1t;
  }
}


std::vector<int> dune::CalibAnaTree::GetXYZIsolatedPoint( std::vector<float> vYPoint , std::vector<float> vZPoint , std::vector<float> vPeakTimeCol , std::vector<int> vNOF ,
					                  float fElectronVelocity , float fTickToMus , float radiusInt , float radiusExt )
{

  if (vYPoint.size() != vZPoint.size())  throw std::invalid_argument( "BIG PROBLEM" );

  std::vector<int> vIso;
  vIso.clear();

  if (radiusInt >= radiusExt) 
  {
    if (fVerbose) std::cout << " Size of DONUT NON PHYSIQUE " << std::endl;   
    return vIso;
  }
  
  int npoint = vYPoint.size();
  std::vector<int> vIsIsolated( npoint , -1 );
  vIso.resize( npoint, -1);

  float zIs = 0;
  float yIs = 0;
  float xIs = 0;
  float xDiff = 0;
  float yDiff = 0;
  float zDiff = 0;

  int cellX = 0;
  int cellY = 0;
  int cellZ = 0;

  bool IsIsolated = true;

  float distSq = 0;
  int iso_count = 0;

  float electronDriftScale = fElectronVelocity * fTickToMus;
  float radiusIntSq = radiusInt*radiusInt;
  float radiusExtSq = radiusExt*radiusExt;

  int nof = 0;

  float gridCellSize = 2*radiusExt; // The size of each grid cell is based on the external radius
  std::unordered_map<std::tuple<int, int, int>, std::vector<int>, GridHasher> grid;

  // Populate the grid with points
  for (int k = 0; k < npoint; k++) 
  {	  
    xIs = vPeakTimeCol[k] * electronDriftScale;
    yIs = vYPoint[k];
    zIs = vZPoint[k];

    // Skip invalid points
    if ((yIs == -999) || (zIs == -999)) 
    {
      vIsIsolated[k] = 0;
      continue;
    }

    // Determine which cell the point belongs to
    cellX = (int)(xIs / gridCellSize);
    cellY = (int)(yIs / gridCellSize);
    cellZ = (int)(zIs / gridCellSize);

    // Store the point in the appropriate grid cell
    grid[std::make_tuple(cellX, cellY, cellZ)].push_back(k);
  }

  for( int k = 0 ; k<npoint ; k++)
  {
    if (vIsIsolated[k] == 0)
    {
      continue;
    }

    zIs = vZPoint[k];
    yIs = vYPoint[k];
    xIs = vPeakTimeCol[k]*electronDriftScale;

    IsIsolated = true;

    nof = vNOF[k];
    if (( yIs == -999) || (zIs == -999))
    {
      vIsIsolated[k] = 0;
      continue;
    }

    // Determine the grid cell of the point
    cellX = (int)(xIs / gridCellSize);
    cellY = (int)(yIs / gridCellSize);
    cellZ = (int)(zIs / gridCellSize);

    for (int dx = -1; dx <= 1; dx++) 
    {
      for (int dy = -1; dy <= 1; dy++) 
      {
        for (int dz = -1; dz <= 1; dz++) 
	{
          auto neighborCell = std::make_tuple(cellX + dx, cellY + dy, cellZ + dz);

          // If the neighboring cell contains points
          if (grid.find(neighborCell) != grid.end()) 
	  {
            for (int i : grid[neighborCell]) // i indice in the vector of indices of point in neighborCell
	    {
              if (i == k) continue; // Skip comparing the point with itself
              if (nof != vNOF[i]) continue; 

              xDiff = xIs - vPeakTimeCol[i] * electronDriftScale;
              yDiff = yIs - vYPoint[i];
              zDiff = zIs - vZPoint[i];

              distSq = xDiff * xDiff + yDiff * yDiff + zDiff * zDiff;

              if (distSq > radiusIntSq && distSq < radiusExtSq) 
	      {
                vIsIsolated[k] = 0; // Mark both points as non-isolated
                vIsIsolated[i] = 0;
                IsIsolated = false;
                break;
              }
            }
          }
          if (!IsIsolated) break;
        }
        if (!IsIsolated) break;
      }
      if (!IsIsolated) break;
    }

    if (IsIsolated)
    {
      vIsIsolated[k] = 1;
      vIso[iso_count] = k;
      iso_count++;
    }
  }//end point loop

  vIso.resize(iso_count);
  return vIso;
}

///////////////////////////////////////////////////////////////////////////////////
// CLUSTER FUNCTIONS

float dune::CalibAnaTree::dist2(dune::CalibAnaTree::point a, dune::CalibAnaTree::point b)
{
    float z = a->z - b->z, y = a->y - b->y;
    return z*z + y*y;
}

float dune::CalibAnaTree::randf(float m)
{
    return m * rand() / (RAND_MAX - 1.);
}

dune::CalibAnaTree::point dune::CalibAnaTree::gen_yz(int size , std::vector<int> vIndex , std::vector<float> vY , std::vector<float> vZ , std::vector<int> vNOF)
{
  int i = 0;
  dune::CalibAnaTree::point p, pt = (dune::CalibAnaTree::point) malloc(sizeof(dune::CalibAnaTree::point_t) * size);

  for (p = pt + size; p-- > pt;)
  {
    p->y = vY[vIndex[i]];
    p->index = vIndex[i];

    if (vNOF[vIndex[i]] == -1)
    {
      p->z = -1*vZ[vIndex[i]] - fgeoZmax;
    }
    else p->z = vZ[vIndex[i]];
    i++;
  }

  return pt;
}


int dune::CalibAnaTree::nearest(dune::CalibAnaTree::point pt, dune::CalibAnaTree::point cent, int n_cluster, float *d2)
{
    int i = 0;
    int  min_i = 0;
    dune::CalibAnaTree::point c;
    float d = 0.0;
    float  min_d = 0.0;

#       define for_n for (c = cent, i = 0; i < n_cluster; i++, c++)
    for_n {
        min_d = HUGE_VAL;
        min_i = pt->group;
        for_n {
            if (min_d > (d = dist2(c, pt))) {
                min_d = d; min_i = i;
            }
        }
    }
    if (d2) *d2 = min_d;
    return min_i;
}

int dune::CalibAnaTree::reallocate(dune::CalibAnaTree::point pt, std::vector<std::vector<float>> ClusterPosition , float threshold)
{
    int  min_i = pt->group;
    float  min_d = HUGE_VAL;

    for( int k = 0 ; k < (int) ClusterPosition[0].size() ; k++) 
    {

      float dist = sqrt(GetDist2D(pt->z,pt->y,ClusterPosition[0][k],ClusterPosition[1][k]));
      if (min_d > dist ) 
      {
         min_d = dist; 
	 min_i = k;
      }
     
    }
    if ( min_d > threshold ) return -1;
    return min_i;
}


float dune::CalibAnaTree::GetDist2D(float y0,float z0,float y1,float z1){
    float z = z0-z1;
    float y = y0-y1;
    return z*z+y*y;
}

float dune::CalibAnaTree::mean(float y,float z){
    return (z+y)/2.;
}

void dune::CalibAnaTree::kpp(dune::CalibAnaTree::point pts, int len, dune::CalibAnaTree::point cent, int n_cent)
{
#       define for_len for (j = 0, p = pts; j < len; j++, p++)
    int j;
    int n_cluster;
    float sum, *d = (float*)malloc(sizeof(float) * len);

    dune::CalibAnaTree::point p;
    cent[0] = pts[ rand() % len ];
    for (n_cluster = 1; n_cluster < n_cent; n_cluster++) {
        sum = 0;
        for_len {
            nearest(p, cent, n_cluster, d + j);
            sum += d[j];
        }
        sum = randf(sum);
        for_len {
            if ((sum -= d[j]) > 0) continue;
            cent[n_cluster] = pts[j];
            break;
        }
    }
    for_len p->group = nearest(p, cent, n_cluster, 0);
    free(d);
}

std::vector<std::vector<float>> dune::CalibAnaTree::lloyd(dune::CalibAnaTree::point pts, int len, int n_cluster)
{
    int i, j, min_i;
    int changed;

    dune::CalibAnaTree::point cent = (dune::CalibAnaTree::point)malloc(sizeof(dune::CalibAnaTree::point_t) * n_cluster), p, c;

    /* assign init grouping randomly */
    //for_len p->group = j % n_cluster;

    /* or call k++ init */
    kpp(pts, len, cent, n_cluster);

    do {
        /* group element for centroids are used as counters */
        for_n { c->group = 0; c->z = c->y = 0; }
        for_len {
            c = cent + p->group;
            c->group++;
            c->z += p->z; c->y += p->y;
        }
        for_n { c->z /= c->group; c->y /= c->group; }

        changed = 0;
        /* fInd closest centroid of each point */
        for_len {
            min_i = nearest(p, cent, n_cluster, 0);
            if (min_i != p->group) {
                changed++;
                p->group = min_i;
            }
        }
    } while (changed > (len >> 10)); /* stop when 99.9% of points are good */

    for_n { c->group = i; }

    std::vector<std::vector<float> > clusterPos;
    std::vector<float> clusterPosY,clusterPosZ;

    dune::CalibAnaTree::point result;

    for(i = 0, result = cent; i < n_cluster; i++, result++) {
        clusterPosZ.push_back(result->z);
        clusterPosY.push_back(result->y);
    }
    clusterPos.push_back(clusterPosZ);
    clusterPos.push_back(clusterPosY);

    return clusterPos;
}


std::vector<std::vector<float>> dune::CalibAnaTree::GetData(int len , dune::CalibAnaTree::point data){

    std::vector<std::vector<float> > dataPos;
    std::vector<float> dataPosZ,dataPosY;


  for(int i = 0; i < len; i++, data++) {
        dataPosZ.push_back(data->z);
        dataPosY.push_back(data->y);
    }
    dataPos.push_back(dataPosZ);
    dataPos.push_back(dataPosY);

    return dataPos;
}

std::vector<int> dune::CalibAnaTree::CheckCompletude(std::vector<std::vector<float> > &data,std::vector<std::vector<float> > &cluster, float RMS , float mult )
{
    int Npts = data[0].size();
    int Ncls = cluster[0].size();

    int Nin = 0, Nin2 = 0, Nout = 0;
    float dist;
    std::vector<int> IDin(Npts,0);

    for(int i = 0 ; i < Ncls ; i++ )
    {
        Nout = 0;
        for(int j = 0;j<Npts;j++)
        {
            if(IDin[j] == 0)
            {
                dist = sqrt(GetDist2D(data[0][j],data[1][j],cluster[0][i],cluster[1][i]));
                // printf("Distance to cluster : %f %f \n",i,Nin,Nout);
                if(dist <= RMS){ Nin++; IDin[j] = 1;}
                else if(dist <= mult*RMS )
                {
                  Nin2++;
                  IDin[j] = 1;
                }
                else{ Nout++; }
            }
        }
    }

    //std::cout << "!!!!!!!!!!!!!!! COMPLETUDE : " << Nin << " " << Nin2 << " " << Nout << std::endl;
    std::vector<int> N(3,0);
    N[0] = Nin;
    N[1] = Nin2;
    N[2] = Nout;

    return N;
}

std::vector<int> dune::CalibAnaTree::CheckClusters(std::vector<std::vector<float> > &data,std::vector<std::vector<float> > &cluster, float RMS , float mult , float tmp)
{

    std::vector<int> v(2,0);
    int Npts = data[0].size();
    int Ncls = cluster[0].size();
    //std::cout << " nombre of point for check cluster : " << Npts << std::endl;

    int Nin = 0, Nin2 = 0, Nout = 0;

    std::vector<int> NComp(3,0);
    NComp  =  CheckCompletude( data, cluster , RMS , mult);
    Nin  = NComp[0];
    Nin2 = NComp[1];
    Nout = NComp[2];

    //std::cout<< " Nin " << Nin << " Nin2 " << Nin2 << " Nout " << Nout << " Npst " << Npts << " Ncls " << Ncls << std::endl;

    if ( fVerbose) printf("Counting : %.03f %.03f %.03f sum : %.02f \n",float(Nin)/float(Npts),float(Nin2)/float(Npts),float(Nout)/float(Npts),float(Nin+Nin2+Nout)/float(Npts));


    if((float(Nin+Nin2)/float(Npts) > tmp -0.04) && float(Nin)/float(Npts) < tmp)
    {
      v[0] = 2;
      v[1] = cluster[0].size();
      return v;
    }
    else if(float(Nin)/float(Npts) < tmp)
    {
      v[0] = 0;
      v[1] = cluster[0].size();
      return v;
    }
    //if(float(Nin)/float(Npts) < 0.99) return 0;


    float dist;
    std::vector<std::vector<int>> IDoverlap( Ncls , std::vector<int> (Ncls ,0));
    int overlap = 0;

    for(int i = 0 ; i < Ncls ; i++)
    {
        for(int j = i ; j < Ncls ; j++)
        {
            if(j > i)
            {
                dist = sqrt(GetDist2D(cluster[0][j],cluster[1][j],cluster[0][i],cluster[1][i]));
                if(dist < 2.*RMS)
                {
                    IDoverlap[i][j] = 1;
                    overlap++;
                }
            }
        }
    }
    int overlap_counter = 0;

    std::vector<float> newclusterZ, newclusterY;
    float meanZ, meanY;
    while( (overlap_counter<10)&&(overlap>0) )
    {
      newclusterZ.clear(); 
      newclusterY.clear();
      meanZ = 0;
      meanY = 0;

      for(int i = 0;i<Ncls;i++)
      {
        meanZ = 0.;
        meanY = 0.;
        overlap = 0;

        for(int j = i;j<Ncls;j++)
        {
          if(IDoverlap[i][j] == 1 && cluster[0][j] != -999)
          {
            meanZ += mean(cluster[0][i],cluster[0][j]);
            meanY += mean(cluster[1][i],cluster[1][j]);
            cluster[0][j] = -999;
            cluster[1][j] = -999;
            overlap++;
          }
        }
        if(overlap == 0 && cluster[0][i] != -999)
        {
          newclusterZ.push_back(cluster[0][i]);
          newclusterY.push_back(cluster[1][i]);
        }
        else if(cluster[0][i] != -999)
        {
          newclusterZ.push_back(meanZ/float(overlap));
          newclusterY.push_back(meanY/float(overlap));
        }
      }

      cluster.clear();
      cluster.push_back(newclusterZ);
      cluster.push_back(newclusterY);

      if ( fVerbose) printf("%lu clusters has been removed at iteration %d \n",Ncls-cluster[0].size(),overlap_counter);

      dist = 0;
      Ncls = cluster[0].size();
      IDoverlap.clear();
      overlap = 0;

      for(int i = 0 ; i < Ncls ; i++)
      {
	std::vector<int> v(Ncls , 0); 
        for(int j = i ; j < Ncls ; j++)
        {
          if(j > i)
          {
            dist = sqrt(GetDist2D(cluster[0][j],cluster[1][j],cluster[0][i],cluster[1][i]));
            if(dist < 2.*RMS)
            {
              v[j] = 1;
              overlap++;
            }
          }
        }
	IDoverlap.push_back(v);
      }

      overlap_counter++;
    }

    NComp = CheckCompletude( data, cluster , RMS , mult );
    Nin  = NComp[0];
    Nin2 = NComp[1];
    Nout = NComp[2];

    if ( fVerbose) printf("Counting : %.03f %.03f %.03f sum : %.02f \n",float(Nin)/float(Npts),float(Nin2)/float(Npts),float(Nout)/float(Npts),float(Nin+Nin2+Nout)/float(Npts));

    if(float(Nin+Nin2)/float(Npts) > tmp && float(Nin)/float(Npts) < tmp)
    {
      v[0] = 2;
      v[1] = cluster[0].size();
      return v;
    }
    else if(float(Nin)/float(Npts) < tmp)
    {
      v[0] = 0;
      v[1] = cluster[0].size();
      return v;
    }

    v[0] = 1;
    v[1] = cluster[0].size();
    return v;

}

std::vector<dune::ClusterInfo*> dune::CalibAnaTree::GetCluster( bool fAsConverged , float CRP_T0, int n_point , int n_cluster , dune::CalibAnaTree::point p , 
                                                     std::vector<float> vEInd1PointByEvent , std::vector<float> vEInd2PointByEvent , 
                                                     std::vector<int> vChInd1PointByEvent , std::vector<int> vChInd2PointByEvent ,  
                                                     std::vector<float> vEnergyColByEvent , std::vector<float> vPeakTimeColByEvent , std::vector<int> vChannelColByEvent ,  
						     bool truth,
                                                     std::vector<int> vMCPDGByEvent , std::vector<int> vMCMOMpdgByEvent ,std::vector<float> vMCWeightByEvent ,  
                                                     std::vector<std::string> vGeneratorTagByEvent ,  
                                                     std::vector<float> vMCXByEvent ,  std::vector<float> vMCYByEvent , std::vector<float> vMCZByEvent , 
                                                     std::vector<int> vNoFByEvent, std::vector<float> vMCEByEvent , std::vector<int> vMCNeByEvent)
{

  int k;
  dune::CalibAnaTree::point vp;
  TempCluster NullCluster;
  NullCluster.Sumz = 0;
  NullCluster.Sumy = 0;
  NullCluster.Npoint = 0;
  NullCluster.ECol = 0;
  NullCluster.PeakTime = 0;
  NullCluster.NCol = 0;
  NullCluster.vNOF.clear();
  NullCluster.vMCPDG.clear();
  NullCluster.vMCMOMpdg.clear();
  NullCluster.vMCWEI.clear();
  NullCluster.vMCGenTag.clear();
  NullCluster.vMCX.clear();
  NullCluster.vMCY.clear();
  NullCluster.vMCZ.clear();
  NullCluster.vMCE.clear();
  NullCluster.vMCNe.clear();
  NullCluster.lChannelCol.clear();
  NullCluster.lChannelInd1.clear();
  NullCluster.EInd1 = 0;
  NullCluster.NInd1 = 0;
  NullCluster.lChannelInd2.clear();
  NullCluster.EInd2 = 0;
  NullCluster.NInd2= 0;
    
  std::vector<TempCluster> vTempCluster(n_cluster, NullCluster);
 
  int out = 0;

  if( fVerbose) std::cout << "there are " << n_cluster << " clusters to check" << std::endl;
  for( k = 0 , vp=p ; k<n_point ; k++ , vp++)
  {

    int index = vp->index;
    int ClusterID   = vp->group;


    int NoF         = vNoFByEvent[index];
    int ChannelCol  = vChannelColByEvent[index];
    int ChannelInd1 = vChInd1PointByEvent[index];
    int ChannelInd2 = vChInd2PointByEvent[index];
    float ECol      = vEnergyColByEvent[index];
    float PeakTime  = vPeakTimeColByEvent[index];

    int MCPart_pdg = vMCPDGByEvent[index];
    int MCPart_mompdg = vMCMOMpdgByEvent[index];
    float MCPart_weight = vMCWeightByEvent[index];
    float MCPart_x = vMCXByEvent[index];
    float MCPart_y = vMCYByEvent[index];
    float MCPart_z = vMCZByEvent[index];
    int MCPart_Ne  = vMCNeByEvent[index];
    int MCPart_E   = vMCEByEvent[index];

    std::string Generator_tag = vGeneratorTagByEvent[index];

    if (ClusterID >=  n_cluster)
    {
      out++;
      continue;
    }

    if (ClusterID ==  -1)
    {
      out++;
      continue;
    }
    vTempCluster[ClusterID].Npoint += 1;

    if ( !Inside( ChannelCol , vTempCluster[ClusterID].lChannelCol ) )
    {
      if (NoF == -1) vTempCluster[ClusterID].Sumz += ECol * ( -1*(vp->z)-fgeoZmax );
      else vTempCluster[ClusterID].Sumz += ECol * ( vp->z );

      vTempCluster[ClusterID].Sumy += ECol * ( vp->y );
      vTempCluster[ClusterID].ECol += ECol;
      vTempCluster[ClusterID].PeakTime += PeakTime;
      vTempCluster[ClusterID].NCol += 1;
      if (truth)
      {
        (vTempCluster[ClusterID].vMCPDG).push_back( MCPart_pdg );
        (vTempCluster[ClusterID].vMCMOMpdg).push_back( MCPart_mompdg );
        (vTempCluster[ClusterID].vMCWEI).push_back( MCPart_weight );
        (vTempCluster[ClusterID].vMCGenTag).push_back( Generator_tag );
        (vTempCluster[ClusterID].vMCX).push_back( MCPart_x );
        (vTempCluster[ClusterID].vMCY).push_back( MCPart_y );
        (vTempCluster[ClusterID].vMCZ).push_back( MCPart_z );
	(vTempCluster[ClusterID].vMCE).push_back( MCPart_E );
	(vTempCluster[ClusterID].vMCNe).push_back( MCPart_Ne );
      }
      (vTempCluster[ClusterID].lChannelCol).push_back( ChannelCol );
      (vTempCluster[ClusterID].vNOF).push_back( NoF );
    }
    if ( (ChannelInd1 != -1) && (!Inside( ChannelInd1 , vTempCluster[ClusterID].lChannelInd1 ) ) )
    {
      vTempCluster[ClusterID].EInd1 += vEInd1PointByEvent[index];
      vTempCluster[ClusterID].NInd1 += 1;
      (vTempCluster[ClusterID].lChannelInd1).push_back( ChannelInd1 );
    }
    if ( (ChannelInd2 != -1) && (!Inside( ChannelInd2 , vTempCluster[ClusterID].lChannelInd2 ) ) )
    {
      vTempCluster[ClusterID].EInd2 += vEInd2PointByEvent[index];
      vTempCluster[ClusterID].NInd2 += 1;
      (vTempCluster[ClusterID].lChannelInd2).push_back( ChannelInd2 );
    }

  }
  if ( fVerbose) std::cout << "temporary clusters done" << std::endl;

  std::vector<dune::ClusterInfo*> vCluster(n_cluster);

  for( int j = 0 ; j < n_cluster ; j++ )
  {

    vCluster[j] = new dune::ClusterInfo();
    vCluster[j]->AsConverged           = fAsConverged;
    vCluster[j]->truth                 = truth;
    vCluster[j]->CRP_T0                = CRP_T0;

    if (vTempCluster[j].ECol) 
    {
      vCluster[j]->z = ( vTempCluster[j].Sumz )/(vTempCluster[j].ECol);
      vCluster[j]->y = ( vTempCluster[j].Sumy )/(vTempCluster[j].ECol);
    }
    else 
    {
      vCluster[j]->z = -999.;
      vCluster[j]->y = -999.;
      if ( fVerbose) std::cout << "cluster with null energy on colection" << std::endl;
    }

    vCluster[j]->NumberOfPoint       = vTempCluster[j].Npoint;
    vCluster[j]->NumberOfCollection  = vTempCluster[j].NCol;
    vCluster[j]->NumberOfInduction1  = vTempCluster[j].NInd1;
    vCluster[j]->NumberOfInduction2  = vTempCluster[j].NInd2;

    std::vector<int> vtemp_NOF = vTempCluster[j].vNOF;
    if ( AllSame( vtemp_NOF ) )  vCluster[j]->nof = vTempCluster[j].vNOF[0];
    else 
    {
      vCluster[j]->nof = -999;
      if (fVerbose) std::cout << "CLUSTER WITH HITS FROM MULTIPLES TPC" << std::endl;
    }
	  
    vCluster[j]->ChargeCollection    = vTempCluster[j].ECol;
    vCluster[j]->ChargeInduction1    = vTempCluster[j].EInd1;
    vCluster[j]->ChargeInduction2    = vTempCluster[j].EInd2;
    if (vTempCluster[j].NCol) vCluster[j]->peaktime = vTempCluster[j].PeakTime/vTempCluster[j].NCol;
    else
    {
      vCluster[j]->peaktime = -999.;
      if ( fVerbose) std::cout << "cluster with no hits on colection" << std::endl;
    }

    if (truth)
    {

      vCluster[j]->vMC_pdg    = vTempCluster[j].vMCPDG;
      vCluster[j]->vMC_mompdg = vTempCluster[j].vMCMOMpdg;
      vCluster[j]->vMC_weight = vTempCluster[j].vMCWEI;
      vCluster[j]->vMC_gentag = vTempCluster[j].vMCGenTag;

      vCluster[j]->vMC_x = vTempCluster[j].vMCX;
      vCluster[j]->vMC_y = vTempCluster[j].vMCY;
      vCluster[j]->vMC_z = vTempCluster[j].vMCZ;

      vCluster[j]->vMC_energy = vTempCluster[j].vMCE;
      vCluster[j]->vMC_nelec  = vTempCluster[j].vMCNe;
    }
  }

  if ( fVerbose) std::cout << "WE HAVE " << (float) out/n_point << " POINTS OUTSIDE OF CLUSTERS" << std::endl;
  return vCluster;
}



std::vector<std::string> dune::CalibAnaTree::GetGeneratorTag( art::Event const &e , art::InputTag fG4producer , art::ServiceHandle<cheat::BackTrackerService> bt_serv )
{
    //int MCPartcounter = 0;
    std::vector<std::pair<int, std::string>> vTrackIdToLabelPair;

    // get all MC truth object
    std::vector<art::Handle<std::vector<simb::MCTruth>>> mcTruths;
    mcTruths = e.getMany<std::vector<simb::MCTruth>>();

    for (auto const &mcTruth : mcTruths)
    {
      // get generator tag (ie name)
      const std::string &sModuleLabel = mcTruth.provenance()->moduleLabel();

      // MC truth (gen) Mc particle association (g4)
      art::FindManyP<simb::MCParticle> MCTruthsToMCParticles( mcTruth , e , fG4producer );
      std::vector<art::Ptr<simb::MCParticle>> mcParts = MCTruthsToMCParticles.at(0);

      //MCPartcounter += (int) mcParts.size();

      for (const art::Ptr<simb::MCParticle> ptr : mcParts)
      {
        int track_id = ptr->TrackId();
        //creation trackID -> Label association
        vTrackIdToLabelPair.push_back(std::make_pair(track_id, sModuleLabel));
      }

      if( fVerbose) std::cout << "THERE ARE " << (int) mcParts.size() << " MCPARTICLES FROM GENERATOR " << sModuleLabel << std::endl;

    }// end for MCtruth

    //sort to but greates trackID in first position
    std::sort(vTrackIdToLabelPair.begin(), vTrackIdToLabelPair.end(), [](std::pair<int,std::string> a, std::pair<int,std::string> b){ return a.first > b.first;});

    // reassociation for quick access
    std::string noTrackID = "no association";
    std::vector<std::string> vGeneratorLabels( vTrackIdToLabelPair[0].first +1 , noTrackID);

    for(int j = 0 ; j < (int) vTrackIdToLabelPair.size() ; j++)
    {
      if (vGeneratorLabels[vTrackIdToLabelPair[j].first] == noTrackID )
      {
        vGeneratorLabels[vTrackIdToLabelPair[j].first] = vTrackIdToLabelPair[j].second;
      }
      else
      {
        if ( fVerbose) std::cout << "EXCEPTION ERROR : ISSUE WITH ASSOCIATION " << vTrackIdToLabelPair[j].first << std::endl;
        vGeneratorLabels[vTrackIdToLabelPair[j].first] = vTrackIdToLabelPair[j].second;
      }

    }// end for pair(trackID,tag)

    return vGeneratorLabels;
}


std::vector<dune::ClusterInfo*> dune::CalibAnaTree::SingleHitAnalysis(
    art::Event const& e,
    art::InputTag fRDTLabel,
    art::InputTag fG4producer,
    art::InputTag fHITproducer,
    bool  bIsPDVD,
    bool  bIsPDHD,
    float fCoincidenceWd1_left,
    float fCoincidenceWd1_right,
    float fCoincidenceWd2_left,
    float fCoincidenceWd2_right,
    bool  bIs3ViewsCoincidence,
    float fPitch,
    float fPitchMultiplier,
    bool  fVerbose,
    float fMinSizeCluster,
    float fMaxSizeCluster,
    float fNumberInitClusters,
    float fRadiusInt,
    float fRadiusExt,
    float fgeoYmin,
    float fgeoYmax,
    float fgeoZmin,
    float fgeoZmax,
    float fElectronVelocity, 
    float fTickTimeInMus)
{
  //definition vector
  std::vector<ClusterInfo*> vec;

  //Set event ID
  int fEventID = e.id().event();

  //clock 
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clock_data);
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  
  //Get DAQ Time Stamp
  bool truth  = true;
  float CRP_T0 = 0;

  if ( e.isRealData())
  {
    auto const& rdts = *e.getValidHandle<std::vector<raw::RDTimeStamp>>(fRDTLabel);
    CRP_T0 = rdts[0].GetTimeStamp();
    truth = false;
  }
 
  // definition of by event
  std::vector<float> vYPointByEvent;
  std::vector<float> vZPointByEvent;

  std::vector<float> vEnergyColByEvent;
  std::vector<float> vEInd1PointByEvent;
  std::vector<float> vEInd2PointByEvent;

  std::vector<float> vPeakTimeColByEvent;

  std::vector<int>   vChannelColByEvent;
  std::vector<int>   vChInd1PointByEvent;
  std::vector<int>   vChInd2PointByEvent;
  std::vector<int>   vNoFByEvent;

  std::vector<int>   vMCMOMpdgByEvent;
  std::vector<int>   vMCPDGByEvent;
  std::vector<float> vMCWeightByEvent;
  std::vector<float> vMCXByEvent;
  std::vector<float> vMCYByEvent;
  std::vector<float> vMCZByEvent;
  std::vector<float> vMCEByEvent;
  std::vector<int>   vMCNeByEvent;
  std::vector<std::string> vGeneratorTagByEvent;

  //retrieve map trackID MC particle to genrator tag
  std::vector<std::string> vTrackIDToGeneratorTag;
  if (!e.isRealData()) vTrackIDToGeneratorTag = GetGeneratorTag( e , fG4producer , bt_serv );
 
  //retrieve hit list
  auto const HitList = e.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
  int fNHits = HitList->size();
 
  if (fNHits == 0)
  {
    return vec;
  }

  for(int index =0 ; index<fNHits; index++)
  {

    const recob::Hit& hit = HitList->at(index);

    geo::WireID fWire   = hit.WireID();
    int fChannel        = hit.Channel();
    int fPlane          = hit.WireID().Plane;

    int fNearOrFarToTheBeam = NearOrFar(bIsPDVD , bIsPDHD , hit);

    //float fEnergy         = hit.ROISummedADC();///fADCtoEl;
    float fEnergy         = hit.Integral();///fADCtoEl;
    float fPeakTime       = hit.PeakTime();//*ftick_in_mus;

    if (fPlane != 2) continue;

    std::vector<int> vMCPart_pdgCode, vMCPart_motherPdg;
    std::vector<float> vMCPart_weight, vMCPart_Endx, vMCPart_Endy, vMCPart_Endz;
    std::vector<std::string> vGenerator_tag;

    float sumMCEnergy = 0;
    int   totElec     = 0;

    if ( truth ) 
    {

      using weightedMCPair = std::pair<const simb::MCParticle*, float>;
      std::vector<weightedMCPair> tempMCPair;

      for (const auto & ide : bt_serv->HitToEveTrackIDEs(clock_data, hit))
      {
        const simb::MCParticle* curr_part = pi_serv->TrackIdToParticle_P(ide.trackID);
        tempMCPair.push_back(std::make_pair(curr_part,ide.energy));
        sumMCEnergy += ide.energy;
	totElec     += ide.numElectrons;
      }
      
      std::sort(tempMCPair.begin(), tempMCPair.end(), [](weightedMCPair a, weightedMCPair b){ return a.second > b.second;});
      
      for (weightedMCPair& p : tempMCPair ) 
      {   
        vMCPart_pdgCode.push_back( (p.first)->PdgCode()   );
        vMCPart_weight.push_back(  (p.second)/sumMCEnergy );
	vGenerator_tag.push_back( vTrackIDToGeneratorTag[(p.first)->TrackId()] );

	vMCPart_Endx.push_back( (float) (p.first)->EndX() );
        vMCPart_Endy.push_back( (float) (p.first)->EndY() );
        vMCPart_Endz.push_back( (float) (p.first)->EndZ() );

        if ( (p.first)->Mother() == 0 )
        {
          vMCPart_motherPdg.push_back( 0 ); //primary particle
          continue;
        }
        const simb::MCParticle* curr_part_mom = pi_serv->TrackIdToParticle_P((p.first)->Mother());
        vMCPart_motherPdg.push_back( curr_part_mom->PdgCode() );
      }
      
      
    }// end if event != real data

    //definition of list 
    std::list<geo::WireID>   lWireInd1; //working variable
    std::list<int>           lChannelInd1;
    std::list<float>         lEnergyInd1;
    std::list<float>         lPeakTimeInd1;
    std::list<float>         lPeakAmpInd1;
    std::list<float>         lYInd1;
    std::list<float>         lZInd1;
    std::list<int>           lChIntersectInd1;
    std::list<float>         lEIntersectInd1;

    std::list<geo::WireID>   lWireInd2; //working variable
    std::list<int>           lChannelInd2;
    std::list<float>         lEnergyInd2;
    std::list<float>         lPeakTimeInd2;
    std::list<float>         lPeakAmpInd2;
    std::list<float>         lYInd2;
    std::list<float>         lZInd2;
    std::list<int>           lChIntersectInd2;
    std::list<float>         lEIntersectInd2;

    std::list<float>         lYPoint;
    std::list<float>         lZPoint;
    std::list<int>           lChInd1Point;
    std::list<int>           lChInd2Point;
    std::list<float>         lEInd1Point;
    std::list<float>         lEInd2Point;

    // Coincidence research

    GetListOfTimeCoincidenceHit( bIsPDVD, bIsPDHD , e, fHITproducer, fCoincidenceWd1_left, fCoincidenceWd1_right ,fCoincidenceWd2_left, fCoincidenceWd2_right , hit, lWireInd1, lWireInd2, lChannelInd1, lChannelInd2, lEnergyInd1, lEnergyInd2, lPeakTimeInd1, lPeakTimeInd2, lPeakAmpInd1, lPeakAmpInd2);

    int fCoincidence = 0;
    if ( !lWireInd1.empty() || !lWireInd2.empty() ) fCoincidence += 1;
    if ( !lWireInd1.empty() && !lWireInd2.empty() ) fCoincidence += 1;

    if ( fCoincidence > 0 )
    {
      GetListOfCrossingChannel(  bIsPDVD, bIsPDHD , fgeoYmin , fgeoYmax , fgeoZmin , fgeoZmax , fWire , lWireInd1 , lWireInd2 , 
			          lChannelInd1 , lEnergyInd1 , lYInd1 , lZInd1 , lChIntersectInd1 , lEIntersectInd1 , 
				  lChannelInd2 , lEnergyInd2 , lYInd2 , lZInd2 , lChIntersectInd2 , lEIntersectInd2); 
      if ( bIs3ViewsCoincidence ) GetListOf3ViewsPoint( fPitch , fPitchMultiplier , lChIntersectInd1 , lYInd1 , lZInd1 , lEIntersectInd1 , lChIntersectInd2 , lYInd2 , lZInd2 , lEIntersectInd2 , lYPoint , lZPoint , lEInd1Point , lEInd2Point , lChInd1Point , lChInd2Point);
      else
      {
        //induction 1
    	lYPoint.insert( lYPoint.end() , lYInd1.begin() , lYInd1.end() );
	lZPoint.insert( lZPoint.end() , lZInd1.begin() , lZInd1.end() );
	lEInd1Point.insert( lEInd1Point.end() , lEIntersectInd1.begin() , lEIntersectInd1.end() );
	lChInd1Point.insert( lChInd1Point.end() , lChIntersectInd1.begin() , lChIntersectInd1.end() );
	int N1 = lYInd1.size();
	lEInd2Point.insert( lEInd2Point.end() , N1 , 0 );
	lChInd2Point.insert( lChInd2Point.end() , N1 , -1 );
 
	//induction 2
	lYPoint.insert( lYPoint.end() , lYInd2.begin() , lYInd2.end() );
	lZPoint.insert( lZPoint.end() , lZInd2.begin() , lZInd2.end() );
	lEInd2Point.insert( lEInd2Point.end() , lEIntersectInd2.begin() , lEIntersectInd2.end() );
	lChInd2Point.insert( lChInd2Point.end() , lChIntersectInd2.begin() , lChIntersectInd2.end() );
	int N2 = lYInd2.size();
	lEInd1Point.insert( lEInd1Point.end() , N2 , 0 );
	lChInd1Point.insert( lChInd1Point.end() , N2 , -1 );
      }
	  
    }
    

    // Retreiving hit info by event

    std::list<int>::iterator ch1  = lChInd1Point.begin();
    std::list<int>::iterator ch2  = lChInd2Point.begin();
    std::list<float>::iterator e1 = lEInd1Point.begin();
    std::list<float>::iterator e2 = lEInd2Point.begin();
    std::list<float>::iterator z  = lZPoint.begin();

    
    for ( auto const y : lYPoint)
    {
      if(( y == -999) || (*z == -999) || (*e1 == -999) || (*e2 == -999) || (*ch1 == -999) || (*ch2 == -999))
      {
        z++;
        ch1++;
        ch2++;
        e1++;
        e2++;
	continue;
      }

      vYPointByEvent.push_back( y );
      vZPointByEvent.push_back( *z );
      vEInd1PointByEvent.push_back( *e1 );
      vEInd2PointByEvent.push_back( *e2 );
      vChInd1PointByEvent.push_back( *ch1 );
      vChInd2PointByEvent.push_back( *ch2 );
      vNoFByEvent.push_back( fNearOrFarToTheBeam );
      vEnergyColByEvent.push_back( fEnergy   );
      vPeakTimeColByEvent.push_back(  fPeakTime );
      vChannelColByEvent.push_back( fChannel );

      // take first origin MC truth for now
      if (( truth )&&( !vMCPart_pdgCode.empty() ))
      {
        vMCPDGByEvent.push_back( vMCPart_pdgCode[0] );
        vMCMOMpdgByEvent.push_back( vMCPart_motherPdg[0] );
        vMCWeightByEvent.push_back( vMCPart_weight[0] );
        vMCXByEvent.push_back( vMCPart_Endx[0] );
        vMCYByEvent.push_back( vMCPart_Endy[0] );
        vMCZByEvent.push_back( vMCPart_Endz[0] );
        vMCNeByEvent.push_back( totElec );
	vMCEByEvent.push_back( sumMCEnergy );

        vGeneratorTagByEvent.push_back( vGenerator_tag[0] );
      }
      else
      {
	vMCPDGByEvent.push_back( -999);
        vMCMOMpdgByEvent.push_back( -999);
        vMCWeightByEvent.push_back(-999 );
        vMCXByEvent.push_back( -999 );
        vMCYByEvent.push_back( -999 );
        vMCZByEvent.push_back( -999 );
        vMCNeByEvent.push_back(-999 );
        vMCEByEvent.push_back(-999);

        vGeneratorTagByEvent.push_back( "data" );
      }

      z++;
      ch1++;
      ch2++;
      e1++;
      e2++;

    }

  }// end hit loop

  std::vector<int> vIso = GetXYZIsolatedPoint( vYPointByEvent , vZPointByEvent , vPeakTimeColByEvent , vNoFByEvent , fElectronVelocity , fTickTimeInMus , fRadiusInt , fRadiusExt );
  int PTSIsolated = (int) vIso.size();

  if (PTSIsolated == 0)
  {
    if ( fVerbose) std::cout << "EXEPTION ERROR : THERE IS NO ISOLATED POINT IN EVENT " << fEventID << std::endl;
    return vec;
  }

  if( fVerbose)
  {
  std::cout << " THERE ARE " << vYPointByEvent.size() << " POINTS IN EVENT " << fEventID << std::endl;
  std::cout << " THERE ARE " << PTSIsolated << " ISOLATED POINTS IN EVENT " << fEventID << std::endl;
  }

  dune::CalibAnaTree::point v;
  v = gen_yz( PTSIsolated , vIso , vYPointByEvent , vZPointByEvent , vNoFByEvent );


  std::vector<std::vector<float> > dataPos = GetData(PTSIsolated,v);
  std::vector<std::vector<float> > clustersPos;

  std::vector<int> vchecks(2,0);

  int K = fNumberInitClusters;

  int check = 0;
  float threshold = fMinSizeCluster;

  int fAsConverged = false;
  for( int i = 0 ; i < fNumberConvStep ; i++ )
  {
    if ( fVerbose) printf("%d %d %d ",i,K,check);

    clustersPos = lloyd(v, PTSIsolated, K);
    vchecks = CheckClusters( dataPos , clustersPos , threshold , fClusterSizeMulti, fCovering);
    check = vchecks[0];

    if(check == 1)  
    {
      fAsConverged = true;	     
      break;
    }
    if(check == 2)
    {
      if (threshold < fMaxSizeCluster)
      {
        threshold *= fClusterSizeMulti;
        if ( fVerbose) printf("Threshold increased \n");
        K = vchecks[1];
      }
      else
      {
        K = vchecks[1] + 1;
        if ( fVerbose) printf("Threshold Max reached \n");
      }
    }
    else K = vchecks[1] + 5;
  }

  if ( fVerbose) printf("Data size : %lu x %lu \n",dataPos.size(),dataPos[0].size());
  if ( fVerbose) printf("Cluster size : %lu x %lu \n",clustersPos.size(),clustersPos[0].size());

  if ( fVerbose) printf("Data clustering ended successfully \n");

  int j = 0;
  dune::CalibAnaTree::point p;
  for (j = 0, p = v; j < PTSIsolated ; j++, p++)
  {
    p->group = reallocate( p , clustersPos , threshold);
  }
  if (fVerbose) std::cout<< " reallocation done " << std::endl;

  vec = GetCluster( fAsConverged, CRP_T0, PTSIsolated , clustersPos[0].size() , v , vEInd1PointByEvent , vEInd2PointByEvent , vChInd1PointByEvent , vChInd2PointByEvent , vEnergyColByEvent , vPeakTimeColByEvent , vChannelColByEvent , truth , vMCPDGByEvent , vMCMOMpdgByEvent , vMCWeightByEvent , vGeneratorTagByEvent , vMCXByEvent , vMCYByEvent , vMCZByEvent , vNoFByEvent , vMCEByEvent , vMCNeByEvent);
  //NCluster = vCluster.size();

  return vec;
}

DEFINE_ART_MODULE(dune::CalibAnaTree)

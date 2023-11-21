#ifndef WIREANA_UTILS_H
#define WIREANA_UTILS_H

#include <vector>
#include <cmath>
#include <cmath>

#define UNCLASSIFIED -1
#define CORE_POINT 1
#define BORDER_POINT 2
#define NOISE -2
#define SUCCESS 0
#define FAILURE -3

namespace wireana{
  using std::string;
  using std::vector;
  using std::pair;
  using std::map;
  using std::set;

  using art::Ptr;
  using art::Event;
  using art::ServiceHandle;

  using cheat::BackTrackerService;
  using cheat::ParticleInventoryService;

  using simb::MCParticle;
  using sim::SimChannel;

  using detinfo::DetectorClocksData;

  using recob::Wire;
  using recob::Hit;

  enum PType{ kUnknown=0, kMarl, kAPA, kCPA, kAr39, kAr42, kNeutron, kKryp, kPlon, kRdon, kNPTypes };

  enum EType{
    kES,
    kCC,
    kRad
  };

  enum GType{
    kNone,kElastic,kMarley,kOthers,kMismatch
  };

  

  struct wirecluster{
    double totalCharge=0.;
    double totalAbsCharge=0.;
    int nWires=0;
    std::vector<art::Ptr<recob::Wire> > wire_pointers;
    double avgCharge(){ return totalCharge/nWires; }
    double avgAbsCharge() {return totalAbsCharge/nWires; }
    void push_back( art::Ptr<recob::Wire> &p )
    {
      this->wire_pointers.push_back( p );
      this->nWires++;
    }
    art::Ptr<recob::Wire> back()
    {
      return this->wire_pointers.back();
    }
  };

  struct roi{
    roi( const art::Ptr<recob::Wire> &wire, recob::Wire::RegionsOfInterest_t::datarange_t r)
      :wire(wire)
    {
      updateroi(wire, r);
    }
    //roi( const roi &r )
    //{
    //  roi( r.wire, r.range );
    //}
    roi(){};

    const art::Ptr<recob::Wire> wire;
    recob::Wire::RegionsOfInterest_t::datarange_t range;
    int channel=-1;
    int begin_index=-1;
    int end_index=-1;
    int width=-1;
    int clusterID=-1;
    double sum=0;
    double abs_sum=0;
    double centroid=0;
    double abs_centroid=0;
    bool hasTrueSignal = false;
    int true_leading_pdg=-999999;
    double true_leading_energy_deposit = 0;
    double true_leading_electron_deposit = 0;
    double true_energy_deposit = 0;
    double true_electron_deposit = 0;
    double tdc_hit_fraction = 0;

    void updateroi( const art::Ptr<recob::Wire> &wire, recob::Wire::RegionsOfInterest_t::datarange_t r)
    {
      this->range = r;
      this->channel = wire->Channel();
      this->begin_index = r.begin_index();
      this->end_index = r.end_index();
      this->width = this->end_index - this->begin_index;
      this->sum = this->Sum();
      this->abs_sum = this->Sum(true);
      this->centroid = this->Centroid();
      this->abs_centroid = this->Centroid(true);
    }

    double Sum(bool useabs=false)
    {
      double ret=0;
      for( auto it = this->range.begin(); it!=this->range.end(); ++it ) ret+= (useabs)? abs(*it):*it;
      return ret;
    }
    double Centroid(bool useabs=false)
    {
      double denum = 0;
      double num = 0;
      for( size_t i = this->range.begin_index(); i != this->range.end_index(); i++ )
      {
        double v = this->range[i];
        if (useabs) v = abs(v);
        denum+=v;
        num+=(i*v);
      }
      if(denum!=0) return num/denum;
      else return -999;
    }

    friend bool operator == ( const roi &lhs, const roi &rhs );
  };
  bool operator == ( const roi &lhs, const roi &rhs )
  {
    return ( lhs.wire == rhs.wire &&  lhs.range == rhs.range && lhs.hasTrueSignal == rhs.hasTrueSignal );
  }



  struct roicluster{
    int clusterID = -1;
    int nWires = -1;
    int channel_min = -1;
    int channel_max = -1;
    int begin_index = -1;
    int end_index = -1;
    int planeid = -1;
    int view = -1;
    bool truthFromNeutrino=false;
    double sum = 0;
    double abs_sum = 0;
    double centroidChannel=0;
    double centroidIndex=0;
    double abs_centroidChannel=0;
    double abs_centroidIndex=0;

    int label_code = -1 ;

    int n_nu = 0;
    int n_lepton = 0;
    int n_proton = 0;
    int n_neutron = 0;
    int n_photon = 0;
    int n_nucleus = 0;
    int n_meson = 0;

    std::vector< roi > ROIs;
    std::vector< std::pair< int, std::pair<float,float> > > pdg_energy_list;
    std::vector< std::pair< int, std::pair<float,float> > > trkID_sum; 
    std::string label;
    TLorentzVector momentum_part;
    TLorentzVector momentum_neutrino;
    std::vector<double> array;
    int width_channels;
    int width_ticks;


    void AddROI( roi &roi, int planeid )
    {
      if(ROIs.size() == 0)
      {
        this->nWires=1;
        this->clusterID=roi.clusterID;
        this->channel_min = roi.channel;
        this->channel_max = roi.channel;
        this->begin_index = roi.begin_index;
        this->end_index = roi.end_index;
        this->view = roi.wire->View();
        this->planeid = planeid;
      }
      else
      {
        if ( this->clusterID != roi.clusterID || this->view != roi.wire->View() || this->planeid != planeid ) 
        {
          std::cout<<"ROI did not match ROICluster!"<<std::endl;
          return;
        }
        ++nWires;
        channel_min = (roi.channel < channel_min)? roi.channel: channel_min;
        channel_max = (roi.channel > channel_max)? roi.channel: channel_max;
        begin_index = (roi.begin_index < begin_index)? roi.begin_index : begin_index;
        end_index = (roi.end_index > end_index)? roi.end_index : end_index;
      }
      ROIs.push_back(roi);
      // centroid = Sum( index*energy )/sum(energy)
      this->abs_centroidChannel*= this->abs_sum;
      this->abs_centroidChannel+= roi.channel * roi.abs_sum;
      this->abs_centroidIndex*= this->abs_sum;
      this->abs_centroidIndex+= roi.abs_centroid * roi.abs_sum;
      this->abs_sum+= roi.abs_sum;
      this->abs_centroidChannel/= this->abs_sum;
      this->abs_centroidIndex/= this->abs_sum;

      this->centroidChannel*= this->sum;
      this->centroidChannel+= roi.channel * roi.sum;
      this->centroidIndex*= this->sum;
      this->centroidIndex+= roi.centroid * roi.sum;
      this->sum+= roi.sum;
      this->centroidChannel/= this->sum;
      this->centroidIndex/= this->sum;

      return;
    }

    double TotalSignal(bool useabs=false)
    {
      double ret = 0;
      for( auto &roi : this->ROIs ) ret+=roi.Sum(useabs);
      return ret;
    }
    double GetNROIS(){ return this->ROIs.size(); }
    int GetWidthTick(){ return (this->end_index - this->begin_index + 1);}
    int GetWidthChannel(){ return (this->channel_max- this->channel_min + 1);}

    std::vector< art::Ptr<recob::Wire> > GetWires()
    {
      std::vector< art::Ptr<recob::Wire> > ret;
      for( auto &roi: this->ROIs ) 
      {
        if( std::find( ret.begin(), ret.end(), roi.wire ) == ret.end() )
        {
          ret.push_back( roi.wire );
        }
      }
      return ret;
    }

    void CalculateArray( double channel_width, double tick_width )
    {
      this->width_channels = channel_width;
      this->width_ticks = tick_width;
      std::vector<double> ret( channel_width*tick_width, 0 );

      int c0,t0;//c1,t0,t1;
      c0=this->abs_centroidChannel - channel_width/2.;
      t0=this->abs_centroidIndex - tick_width/2.;

      for( auto& roi : this->ROIs )
      {
        int c = roi.channel - c0;
        for ( int index = roi.begin_index; index<=roi.end_index; ++index )
        {
          int t = index - t0;
          int arr_index = c + t*channel_width;
          ret[arr_index] = roi.wire->SignalROI()[index];
        }
      }
      this->array = ret;
    }
    
    unsigned int GetIndex(unsigned int c, unsigned int t)
    {
      return c%this->width_channels + t*this->width_channels;
    }

    std::pair<unsigned int, unsigned int> GetCT( unsigned index )
    {
      unsigned int c = index%this->width_channels;
      unsigned int t = index/this->width_channels;
      return std::make_pair(c,t);
    }

    double GetArray(unsigned int c, unsigned int t)
    {
      unsigned int index = GetIndex(c,t);
      if (this->array.size() == 0 || index >= (this->array.size()) ) return -1;
      return this->array[index];
    }

    friend bool operator==(const roicluster&l, const roicluster& r );
  };
  bool operator==(const roicluster& l ,const roicluster& r )
  {
    return (
        l.clusterID == r.clusterID &&
        l.nWires == r.nWires &&
        l.channel_min        ==      r.channel_min           &&
        l.channel_max        ==      r.channel_max           &&
        l.begin_index        ==      r.begin_index           &&
        l.end_index          ==      r.end_index             &&
        l.planeid            ==      r.planeid               &&
        l.view               ==      r.view                  &&
        l.truthFromNeutrino  ==      r.truthFromNeutrino     &&
        l.sum                ==      r.sum                   &&
        l.abs_sum            ==      r.abs_sum               &&
        l.centroidChannel    ==      r.centroidChannel       &&
        l.centroidIndex      ==      r.centroidIndex         &&
        l.abs_centroidChannel==      r.abs_centroidChannel   &&
        l.abs_centroidIndex  ==      r.abs_centroidIndex     &&
        l.label_code         ==      r.label_code            &&
        l.n_nu               ==      r.n_nu                  &&
        l.n_lepton           ==      r.n_lepton              &&
        l.n_proton           ==      r.n_proton              &&
        l.n_neutron          ==      r.n_neutron             &&
        l.n_photon           ==      r.n_photon              &&
        l.n_nucleus          ==      r.n_nucleus             &&
        l.n_meson            ==      r.n_meson          &&
        l.ROIs               ==      r.ROIs      );
  }

  struct matchedroicluster{
    matchedroicluster( int planeid, double metric, roicluster &u, roicluster &v, roicluster &z )
    {
      this->planeid = planeid;
      this->metric = metric;
      clusters.push_back(u);
      clusters.push_back(v);
      clusters.push_back(z);
      totalROIs = u.ROIs.size() + v.ROIs.size() + z.ROIs.size();
    }
    int planeid;
    double metric;
    int totalROIs;
    //u=0, v=1, z=2
    std::vector<roicluster> clusters;

    unsigned int gencode()
    {
      std::set<unsigned int> labelset;
      for( auto &cluster:clusters )
      {
        labelset.insert( cluster.label_code );
      }
      if ( labelset.size() == 1 ) return clusters.front().label_code;
      return kMismatch;
    }

    bool labelmatch()
    {
      bool matched = false;
      if(clusters.size()==3)
      {
        std::string l0 = clusters[0].label;
        std::string l1 = clusters[1].label;
        std::string l2 = clusters[2].label;
        matched=((l0==l1)&&(l1==l2));
      }
      return matched;
    }
    bool trkmatch()
    {
      bool matched = false;
      if(clusters.size()==3)
      {
        for( auto &cluster : clusters )
        {
          if ( cluster.trkID_sum.size() == 0 ) return false;
        }
        int l0 = clusters[0].trkID_sum.front().first;
        int l1 = clusters[1].trkID_sum.front().first;
        int l2 = clusters[2].trkID_sum.front().first;
        matched=((l0==l1)&&(l1==l2));
      }
      return matched;
    }

  };

  typedef std::map<int, std::map< geo::View_t, std::vector<roi> > > PlaneViewROIMap;
  typedef std::map<int, std::map< geo::View_t, std::vector<roicluster> > > PlaneViewROIClusterMap;

  struct ApaROIConstainer
  {
    int PlaneID;
    std::map<int, std::vector<roi>> ViewROIMap;
  };


  class WireAnaDBSCAN;
  class ROIMatcher;
  class numpywriter;
}


class wireana::WireAnaDBSCAN {
public:    
    WireAnaDBSCAN(){
    }
    ~WireAnaDBSCAN(){}

    void SetParameters(std::vector<wireana::roi> &ROIs, unsigned int minPts=2, float eps=4.5, float drift=1.6, float pitch=3.0 ){
        m_minPoints = minPts;
        m_epsilon = eps;
        m_ROIs = &ROIs;
        m_roiSize = ROIs.size();
        m_drift = drift;
        m_pitch = pitch;
    }

    int run();
    void promoteSingleROI( int minWidth=150, double minAbsSum = 80 );

    vector<int> calculateCluster(wireana::roi &r);
    int expandCluster(wireana::roi &r, int clusterID);

    float calculateTickDistance( const wireana::roi& r1, const wireana::roi&r2 );
    float calculateDistance( const wireana::roi &r1, const wireana::roi &r2, float driftspeed=1.6, float pitch=3.0 );

    void setDriftAndPitch( float drift, float pitch ){ m_drift=drift; m_pitch = pitch ; }
    int getTotalROISize() {return m_roiSize;}
    int getMinimumClusterSize() {return m_minPoints;}
    int getEpsilonSize() {return m_epsilon;}

    std::vector<wireana::roi> &GetROIs(){ return *m_ROIs; }

    std::vector<wireana::roi> *m_ROIs;
    
private:    
    unsigned int m_roiSize;
    unsigned int m_minPoints;
    float m_epsilon;
    float m_drift;
    float m_pitch;
};




int 
wireana::WireAnaDBSCAN::run()
{
  int clusterID = 1;
  std::vector<roi>::iterator iter;
  if (!m_ROIs) return -1;
  for(iter = m_ROIs->begin(); iter != m_ROIs->end(); ++iter)
  {
    if ( iter->clusterID == UNCLASSIFIED )
    {
      if ( expandCluster(*iter, clusterID) != FAILURE )
      {
        clusterID += 1;
      }
    }
  }
  return 0;
}


std::vector<int> 
wireana::WireAnaDBSCAN::calculateCluster(wireana::roi &r)
{
  int index = 0;
  std::vector<wireana::roi>::iterator iter;
  std::vector<int> clusterIndex;
  for( iter = m_ROIs->begin(); iter != m_ROIs->end(); ++iter)
  {
    if ( calculateDistance(r, *iter, m_drift, m_pitch) <= m_epsilon )
    {
      clusterIndex.push_back(index);
    }
    index++;
  }
  return clusterIndex;
}

int 
wireana::WireAnaDBSCAN::expandCluster(wireana::roi &r, int clusterID)
{  
  std::vector<int> clusterSeeds = calculateCluster(r);

  if ( clusterSeeds.size() < m_minPoints )
  {
    r.clusterID = NOISE;
    return FAILURE;
  }
  else
  {
    int index = 0, indexCorePoint = 0;
    std::vector<int>::iterator iterSeeds;
    for( iterSeeds = clusterSeeds.begin(); iterSeeds != clusterSeeds.end(); ++iterSeeds)
    {
      m_ROIs->at(*iterSeeds).clusterID = clusterID;
      if (m_ROIs->at(*iterSeeds).channel == r.channel && m_ROIs->at(*iterSeeds).begin_index == r.begin_index && m_ROIs->at(*iterSeeds).end_index == r.end_index )
      {
        indexCorePoint = index;
      }
      ++index;
    }
    clusterSeeds.erase(clusterSeeds.begin()+indexCorePoint);

    for( std::vector<int>::size_type i = 0, n = clusterSeeds.size(); i < n; ++i )
    {
      std::vector<int> clusterNeighors = calculateCluster(m_ROIs->at(clusterSeeds[i]));

      if ( clusterNeighors.size() >= m_minPoints )
      {
        std::vector<int>::iterator iterNeighors;
        for ( iterNeighors = clusterNeighors.begin(); iterNeighors != clusterNeighors.end(); ++iterNeighors )
        {
          if ( m_ROIs->at(*iterNeighors).clusterID == UNCLASSIFIED || m_ROIs->at(*iterNeighors).clusterID == NOISE )
          {
            if ( m_ROIs->at(*iterNeighors).clusterID == UNCLASSIFIED )
            {
              clusterSeeds.push_back(*iterNeighors);
              n = clusterSeeds.size();
            }
            m_ROIs->at(*iterNeighors).clusterID = clusterID;
          }
        }
      }
    }
    return SUCCESS;
  }
}

float 
wireana::WireAnaDBSCAN::calculateTickDistance(const wireana::roi &r1, const wireana::roi &r2)
  {
    roi rr1 = (r1.begin_index < r2.begin_index)? r1 : r2;
    roi rr2 = (r1.begin_index < r2.begin_index)? r2 : r1;
    float ret = rr2.begin_index - rr1.end_index;
    return (ret>0)? ret: 0;
  }
float 
wireana::WireAnaDBSCAN::calculateDistance( const roi &r1, const roi &r2, float driftspeed, float pitch )
{
  float dx = this->calculateTickDistance(r1,r2)*driftspeed/2.;
  float dz = (r1.channel - r2.channel)*pitch;
  return pow(dx*dx+dz*dz,0.5);
}


class wireana::ROIMatcher
{
  public:
  ROIMatcher(double dist = 20 /*mm*/){
    fGeometry = &*(art::ServiceHandle<geo::Geometry const>());
    SetMatchDistance(dist);
  }
  ~ROIMatcher(){}

  void Reset();
  void SetData( PlaneViewROIClusterMap pvrm ) { m_planeViewRoillusterMap = pvrm; }
  void SetMatchDistance( float dist = 20 /*mm*/ ) { m_minMatchDistance = dist ;}
  void MatchROICluster();
  void CleanDuplicates( float delta = 0. );
  void SortROIClustersBySize();
  const std::vector<matchedroicluster> &GetMatchedClusters() { return m_matchedclusters; }

  private:
  bool OverlapInTime(const roicluster& r1, const roicluster& r2, double frac=.8);
  bool OverlapInSpace(const roicluster& r1, const roicluster& r2);

  double DeltaT(const roicluster& r1, const roicluster& r2);
  double DeltaC(const roicluster& r1, const roicluster& r2);

  std::vector<geo::Point_t> GetIntersectionPoints(const roicluster& r1, const roicluster& r2);
  double DeltaDist(const roicluster& r1, const roicluster& r2, const roicluster& r3);

  double DeltaN(const roicluster& r1, const roicluster& r2);
  double DeltaCT( const wireana::roicluster& r1, const wireana::roicluster &r2 );

  geo::GeometryCore const* fGeometry;

  PlaneViewROIClusterMap m_planeViewRoillusterMap;
  std::vector<matchedroicluster> m_matchedclusters; 

  double m_minMatchDistance;
};

void 
wireana::ROIMatcher::Reset()
{
  m_planeViewRoillusterMap.clear();
  m_matchedclusters.clear();
}

bool 
wireana::ROIMatcher::OverlapInTime(const roicluster& r1, const roicluster& r2, double frac)
{
  roicluster rr1 = (r1.begin_index < r2.begin_index)? r1 : r2;
  roicluster rr2 = (r1.begin_index < r2.begin_index)? r2 : r1;
  bool overlap = (rr2.begin_index - rr1.end_index)<=0;
  //must overlap and similar sized
  int dt1 = rr1.end_index - rr1.begin_index;
  int dt2 = rr2.end_index - rr2.begin_index;
  bool similar = false;
  if (dt1 > dt2) similar = ((dt1-dt2)/dt1 > frac);
  else similar = ((dt2-dt1)/dt2 > frac);
  return (overlap && similar);
}

double
wireana::ROIMatcher::DeltaT(const wireana::roicluster& r1, const wireana::roicluster& r2)
{
  return abs( r1.begin_index - r2.begin_index ) + abs( r1.end_index - r2.end_index );
}

double
wireana::ROIMatcher::DeltaN(const wireana::roicluster& r1, const wireana::roicluster& r2)
{
  return abs( r1.nWires - r2.nWires );
}

double
wireana::ROIMatcher::DeltaC(const roicluster& r1, const roicluster& r2)
{
  double ret = 0;
  std::map<raw::ChannelID_t, std::vector<geo::WireID>> cwidmap;
  for( auto &rr1 : r1.ROIs ){
    raw::ChannelID_t c1 = rr1.wire->Channel();
    auto c1towid = fGeometry->ChannelToWire(c1);
    for( auto &rr2: r2.ROIs )
    {
      //double y,z;
      raw::ChannelID_t c2 = rr2.wire->Channel();
      if( cwidmap.find(c2) == cwidmap.end() ) cwidmap[c2] = fGeometry->ChannelToWire(c2);
      //loop through each wid, set intersect to true if any has intersection

      bool intersect = false;
      for( auto &w1: c1towid )
      {
        for( auto &w2: cwidmap[c2] )
        {
          if( w1.asTPCID() != w2 ) continue;
          geo::Point_t intersection_point;
          intersect = fGeometry->WireIDsIntersect(w1,w2,intersection_point);
          if( intersect ) break;
        }
      }
      ret+= (!intersect);
    };
  }
  return ret;
}

std::vector<geo::Point_t>
wireana::ROIMatcher::GetIntersectionPoints( const roicluster& r1, const roicluster& r2)
{

  std::vector<geo::Point_t> points;
  int c1 = std::round(r1.centroidChannel), c2 = std::round(r2.centroidChannel);
  std::vector< geo::WireID > c1wid = fGeometry->ChannelToWire( c1 );
  std::vector< geo::WireID > c2wid = fGeometry->ChannelToWire( c2 );
  for( auto &w1: c1wid )
  {
    for( auto &w2:c2wid )
    {
      geo::Point_t point(-99999,-99999,-99999);
      if ( w1.asTPCID() != w2.asTPCID() ) continue;
      if ( w1.asPlaneID() == w2.asPlaneID() ) continue;
      bool intersect = fGeometry->WireIDsIntersect(w1,w2,point);
      if (intersect) points.push_back( point );
    }
  }
  return points;
}

double 
wireana::ROIMatcher::DeltaDist( const roicluster &u, const roicluster &v, const roicluster&z)
{
  auto uv = wireana::ROIMatcher::GetIntersectionPoints(u,v);
  auto uz = wireana::ROIMatcher::GetIntersectionPoints(u,z);
  auto vz = wireana::ROIMatcher::GetIntersectionPoints(v,z);
  if (uv.size() == 0 || uz.size() == 0 || vz.size() == 0 ) return false;
  //then find the maximum delta distance between each triplet
  std::vector<double> distances;

  std::vector<double> p12d, p13d, p23d;
  for( auto &p1 : uv )
  {
    for (auto &p2: uz )
    {
      double d12 = (p1-p2).R();
      for (auto &p3: vz )
      {
        double d13 = (p1-p3).R();
        double d23 = (p2-p3).R();
        double max = d12;
        if( d13>max ) max = d13;
        if( d23>max ) max = d23;
        distances.push_back(max);
      }
    }
  }
  return *std::min_element( distances.begin(), distances.end() );
}

double 
wireana::ROIMatcher::DeltaCT( const wireana::roicluster& r1, const wireana::roicluster &r2 )
{
  double dT=0,dN=0,dC=0;
  dT = DeltaT(r1,r2);
  dN = DeltaN(r1,r2);
  //dC = DeltaC(r1,r2);
  return dT+dN+dC;
}

void 
wireana::ROIMatcher::SortROIClustersBySize()
{
  std::sort(m_matchedclusters.begin(), m_matchedclusters.end(),
      [](const auto &a, const auto &b){ return a.totalROIs > b.totalROIs; });
}
void 
wireana::ROIMatcher::MatchROICluster()
{
  m_matchedclusters.clear();

  for( auto& pvv: m_planeViewRoillusterMap )//plane: (view, vector<roicluster>)
  {
    int planeid = pvv.first;
    //u=0, v=1, z=2
    for( auto & rcu: pvv.second[geo::kU] )
    {
      if( rcu.clusterID < 0 ) continue;
      int n_rcu = rcu.ROIs.size();
      for( auto & rcv: pvv.second[geo::kV] )
      {
        if( rcv.clusterID < 0 ) continue;
        int n_rcv = rcv.ROIs.size();
        double uv = DeltaCT(rcu,rcv);
        for( auto & rcz: pvv.second[geo::kZ] )
        {
          if( rcz.clusterID < 0 ) continue;
          double dist = DeltaDist( rcu,rcv,rcz );
          if( dist > m_minMatchDistance ) continue;
          int n_rcz = rcz.ROIs.size();
          double uz = DeltaCT(rcu,rcz);
          double vz = DeltaCT(rcv,rcz);
          double metric = pow( uv*uv+uz*uz+vz*vz, 0.5 )/pow(n_rcu+n_rcv+n_rcz,3);
          m_matchedclusters.push_back( matchedroicluster(planeid, metric, rcu, rcv, rcz ));
        }//loop z
      }//loop v
    }//loop u 
    std::sort(m_matchedclusters.begin(), m_matchedclusters.end(),
      [](const auto &a, const auto &b){ return a.metric < b.metric; });
  }
}
void 
wireana::ROIMatcher::CleanDuplicates(float delta)
{
  if ( m_matchedclusters.size() == 0 ) return;
  std::cout<<"Initial Matched Clusters Size: "<<m_matchedclusters.size()<<std::endl;
  auto it = m_matchedclusters.begin(); 
  while ( it != m_matchedclusters.end() )
  {
    //it->clusters;
    auto it2 = it+1;
    while ( it2 != m_matchedclusters.end() )
    {
      //loop over rest of the clusters and remove any that contains duplicate clusters. 
      bool found_duplicate_cluster = false;
      if( std::abs(it->metric - it2->metric) > delta )
      {
        for( auto clus = it2->clusters.begin(); clus != it2->clusters.end(); ++clus )
        {
          if ( std::find( it->clusters.begin(), it->clusters.end(), *clus ) != it->clusters.end() )
          {
            //found roicluster
            found_duplicate_cluster = true;
            break;
          }
        }
      }
      if ( found_duplicate_cluster ) 
        m_matchedclusters.erase( it2 ); 
      // since current it2 erased, it2 automatically becomes the "next" element
      else ++it2;
    }
    ++it;
  }
  std::cout<<"Final Matched Clusters Size: "<<m_matchedclusters.size()<<std::endl;
}

#endif

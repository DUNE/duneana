#include "CalibAnaTree.h"
#include "larcore/Geometry/WireReadout.h"
#include "art/Framework/Principal/Event.h"

typedef struct //!< 2D point for clustering : - group (number of the associated cluster)
                 //!<                           - index (index to retreive the info like energy to the associated hit)
{
    float y, z, t; //!<  [cm]
    int group, index;
} point_t, *point;

typedef struct //!< Temporary Clusters used to construct the Clusters (ClusterInfo)
  {
    float Sumy , Sumz , ECol , EInd1 , EInd2 , PeakTime;
    int Npoint, NCol , NInd1 , NInd2;
    std::list<int> lChannelCol , lChannelInd1 , lChannelInd2;
    std::vector<int> vMCPDG , vMCMOMpdg , vNOF;
    std::vector<float> vMCWEI;
    std::vector<float> vMCX , vMCY , vMCZ;
    std::vector<std::string> vMCGenTag;
    std::vector<float> vMCE;
    std::vector<int> vMCNe;
} TempCluster ;

//!< stucture used for association between a grid coordinate and a int
//!< this structure is used for voxelisation and isolation search
struct GridHasher
{
    std::size_t operator()(const std::tuple<int, int, int>& key) const {
        // Use a combination of prime numbers to create a unique hash from 3D coordinates
        return std::get<0>(key) * 73856093 ^ std::get<1>(key) * 19349663 ^ std::get<2>(key) * 83492791;
    }
};

bool Inside( int k , std::list<int> list){
  return (std::find(list.begin(), list.end(), k) != list.end());
}

bool AllSame( std::vector<int> v)
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
int NearOrFar( bool IsPDVD , bool IsPDHD , bool IsFDVD , bool IsFDHD , const recob::Hit & hit)
{
  if (IsPDHD)
  {
    if ( (hit.WireID().TPC == 2)|| (hit.WireID().TPC == 6) || (hit.WireID().TPC == 3) || (hit.WireID().TPC == 7) ) return  -1; // 3 and 7 are dumie TPC
    if ( (hit.WireID().TPC == 1)|| (hit.WireID().TPC == 5) || (hit.WireID().TPC == 0) || (hit.WireID().TPC == 4) ) return  1;  // 0 and 4 are dumie TPC
  }
  if (IsPDVD)
  {
    if (hit.WireID().TPC <= 7 ) return -1;
    if (hit.WireID().TPC  > 7 ) return 1;
  }
  if (IsFDHD)
  {
    if( hit.WireID().TPC %2 ==0 ) return +1;
    else return -1;
  }
  if (IsFDVD)
  {
    if(hit.WireID().TPC % 8 < 4) return +1;
    else return -1;
  }
  return -999;
}


void GetListOfTimeCoincidenceHit( bool IsPDVD , bool IsPDHD , bool IsFDVD , bool IsFDHD , art::Event const & ev, art::InputTag HitLabel, float const CoincidenceWd1_l , float const CoincidenceWd1_r, float const CoincidenceWd2_l , float const CoincidenceWd2_r, const recob::Hit & HitCol, 
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
  int NoFCol = NearOrFar(IsPDVD,IsPDHD,IsFDVD,IsFDHD,HitCol);
  int NoF = -4;

  for (int i=0, sz=hitlist->size(); i!=sz; ++i)
  { 
    hit   = hitlist->at(i);
    Plane = hit.WireID().Plane;
    if (Plane == 2) continue;

    NoF = NearOrFar(IsPDVD,IsPDHD,IsFDVD,IsFDHD,hit);
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

bool IntersectOutsideOfTPC( float Ymin , float Ymax , float Zmin , float Zmax ,
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

void GetListOfCrossingChannel(   bool IsPDVD , bool IsPDHD , float Ymin , float Ymax , float Zmin , float Zmax ,
		                                    geo::WireID & WireCol , std::list<geo::WireID> & WireInd1 , std::list<geo::WireID> & WireInd2 , 
						    std::list<int>  & ChInd1 , std::list<float> & EInd1 , std::list<float> & YInd1 , std::list<float> & ZInd1 , 
						    std::list<int>  & ChIntersectInd1 , std::list<float> & EIntersectInd1 ,
                                                    std::list<int>  & ChInd2 , std::list<float> & EInd2 , std::list<float> & YInd2 , std::list<float> & ZInd2 , 
						    std::list<int>  & ChIntersectInd2 , std::list<float> & EIntersectInd2,
						    const geo::WireReadoutGeom& fWireReadout)
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
    
      if (IsPDVD)
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

      if (IsPDVD)
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

void GetListOf3ViewsPoint( float pitch , float alpha , 
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


std::vector<int> GetXYZIsolatedPoint( std::vector<float> vYPoint , 
		                      std::vector<float> vZPoint , 
				      std::vector<float> vPeakTimeCol , 
				      std::vector<int> vNOF ,
				      float fElectronVelocity , 
				      float fTickToMus , 
				      float radiusInt , 
				      float radiusExt ,
				      bool fVerbose )
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

float dist2( point a, point b)
{
    float z = a->z - b->z, y = a->y - b->y , t = a->t - b->t;
    return z*z + y*y + t*t;
}

float randf(float m)
{
    return m * rand() / (RAND_MAX - 1.);
}

point gen_yzt(int size , std::vector<int> vIndex , std::vector<float> vY , std::vector<float> vZ , std::vector<float> vT , std::vector<int> vNOF, float fElectronVelocity , float fTickToMus ,float fgeoZmax)
{
  int i = 0;
  point p, pt = ( point) malloc(sizeof( point_t) * size);

  for (p = pt + size; p-- > pt;)
  {
    p->y = vY[vIndex[i]];
    p->index = vIndex[i];

    float time_in_cm = vT[vIndex[i]]*electronDriftScale;
    p->t = time_in_cm;
    if (vNOF[vIndex[i]] == -1)
    {
      p->z = -1*vZ[vIndex[i]] - fgeoZmax;
    }
    else p->z = vZ[vIndex[i]];
    i++;
  }

  return pt;
}


int nearest(point pt, point cent, int n_cluster, float *d2)
{
    int i = 0;
    int  min_i = 0;
    point c;
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
/*
float GetDist2D(float y0,float z0,float y1,float z1){
    float z = z0-z1;
    float y = y0-y1;
    return z*z+y*y;
}
*/
int reallocate(point pt, std::vector<std::vector<float>> ClusterPosition , float threshold)
{
    int  min_i = pt->group;
    float  min_d = HUGE_VAL;

    for( int k = 0 ; k < (int) ClusterPosition[0].size() ; k++) 
    {

      float dist = GetDist(pt->z,pt->y,pt->t,ClusterPosition[0][k],ClusterPosition[1][k],ClusterPosition[2][k]);
      if (min_d > dist ) 
      {
         min_d = dist; 
	 min_i = k;
      }
     
    }
    if ( min_d > threshold ) return -1;
    return min_i;
}


float mean(float y,float z){
    return (z+y)/2.;
}

void kpp(point pts, int len, point cent, int n_cent)
{
#       define for_len for (j = 0, p = pts; j < len; j++, p++)
    int j;
    int n_cluster;
    float sum, *d = (float*)malloc(sizeof(float) * len);

    point p;
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

std::vector<std::vector<float>> lloyd(point pts, int len, int n_cluster)
{
    int i, j, min_i;
    int changed;

     point cent = ( point)malloc(sizeof( point_t) * n_cluster), p, c;

    /* assign init grouping randomly */
    //for_len p->group = j % n_cluster;

    /* or call k++ init */
    kpp(pts, len, cent, n_cluster);

    do {
        /* group element for centroids are used as counters */
        for_n { c->group = 0; c->z = c->y = c->t = 0; }
        for_len {
            c = cent + p->group;
            c->group++;
            c->z += p->z;
	    c->y += p->y;
	    c->t += p->t;
        }
        for_n 
	{ 
	    c->z /= c->group; 
            c->y /= c->group;
            c->t /= c->group;
        }

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
    std::vector<float> clusterPosY,clusterPosZ,clusterPosT;

     point result;

    for(i = 0, result = cent; i < n_cluster; i++, result++) {
        clusterPosZ.push_back(result->z);
        clusterPosY.push_back(result->y);
        clusterPosT.push_back(result->t);
    }
    clusterPos.push_back(clusterPosZ);
    clusterPos.push_back(clusterPosY);
    clusterPos.push_back(clusterPosT);

    return clusterPos;
}


std::vector<std::vector<float>> GetData(int len ,  point data){

    std::vector<std::vector<float> > dataPos;
    std::vector<float> dataPosZ,dataPosY;


  for(int i = 0; i < len; i++, data++) {
        dataPosZ.push_back(data->z);
        dataPosY.push_back(data->y);
	dataPosT.push_back(data->t);
    }
    dataPos.push_back(dataPosZ);
    dataPos.push_back(dataPosY);
    dataPos.push_back(dataPosT);
    return dataPos;
}

std::vector<int> CheckCompletude(std::vector<std::vector<float> > &data,std::vector<std::vector<float> > &cluster, float RMS , float mult )
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
                dist = GetDist(data[0][j],data[1][j],data[2][j],cluster[0][i],cluster[1][i],cluster[2][i]);
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

std::vector<int> CheckClusters(std::vector<std::vector<float> > &data,std::vector<std::vector<float> > &cluster, float RMS , float mult , float tmp , bool fVerbose)
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


    if((float(Nin+Nin2)/float(Npts) > tmp ) && float(Nin)/float(Npts) < tmp)
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
                dist = GetDist(cluster[0][j],cluster[1][j],cluster[2][j],cluster[0][i],cluster[1][i],cluster[2][i]);
                if(dist < 2.*RMS)
                {
                    IDoverlap[i][j] = 1;
                    overlap++;
                }
            }
        }
    }
    int overlap_counter = 0;

    std::vector<float> newclusterZ, newclusterY, newclusterT;
    float meanZ, meanY, meanT;
    while( (overlap_counter<10)&&(overlap>0) )
    {
      newclusterZ.clear(); 
      newclusterY.clear();
      newclusterT.clear();	    
      meanZ = 0;
      meanY = 0;
      meanT = 0;
	    
      for(int i = 0;i<Ncls;i++)
      {
        meanZ = 0.;
        meanY = 0.;
	meanT = 0.;      
        overlap = 0;

        for(int j = i;j<Ncls;j++)
        {
          if(IDoverlap[i][j] == 1 && cluster[0][j] != -999)
          {
            meanZ += mean(cluster[0][i],cluster[0][j]);
            meanY += mean(cluster[1][i],cluster[1][j]);
	    meanT += mean(cluster[2][i],cluster[2][j]);
            cluster[0][j] = -999;
            cluster[1][j] = -999;
	    cluster[2][j] =-999;  
            overlap++;
          }
        }
        if(overlap == 0 && cluster[0][i] != -999)
        {
          newclusterZ.push_back(cluster[0][i]);
          newclusterY.push_back(cluster[1][i]);
	  newclusterT.push_back(cluster[2][i]);
        }
        else if(cluster[0][i] != -999)
        {
          newclusterZ.push_back(meanZ/float(overlap));
          newclusterY.push_back(meanY/float(overlap));
	  newclusterT.push_back(meanT/float(overlap));
        }
      }

      cluster.clear();
      cluster.push_back(newclusterZ);
      cluster.push_back(newclusterY);
      cluster.push_back(newclusterT);

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
            dist = GetDist(cluster[0][j],cluster[1][j],cluster[2][j],cluster[0][i],cluster[1][i],cluster[2][i]);
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

std::vector<dune::ClusterInfo*> GetCluster( bool fAsConverged , float CRP_T0, int n_point , 
		                            int n_cluster , point p , 
                                            std::vector<float> vEInd1PointByEvent , std::vector<float> vEInd2PointByEvent , 
                                            std::vector<int> vChInd1PointByEvent , std::vector<int> vChInd2PointByEvent ,  
                                            std::vector<float> vEnergyColByEvent , std::vector<float> vPeakTimeColByEvent ,
					    std::vector<int> vChannelColByEvent ,  
					    bool truth,
                                            std::vector<int> vMCPDGByEvent , std::vector<int> vMCMOMpdgByEvent ,
					    std::vector<float> vMCWeightByEvent ,  
                                            std::vector<std::string> vGeneratorTagByEvent ,  
                                            std::vector<float> vMCXByEvent ,  std::vector<float> vMCYByEvent ,
					    std::vector<float> vMCZByEvent , 
                                            std::vector<int> vNoFByEvent, std::vector<float> vMCEByEvent , 
					    std::vector<int> vMCNeByEvent , bool fVerbose,
					    float fgeoZmax)
{

  int k;
  point vp;
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



std::vector<std::string> GetGeneratorTag( art::Event const &e , art::InputTag fG4producer , art::ServiceHandle<cheat::BackTrackerService> bt_serv , bool fVerbose)
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


std::vector<dune::ClusterInfo*> SingleHitAnalysis(
    art::Event const& e,
    art::InputTag fRDTLabel,
    art::InputTag fG4producer,
    art::InputTag fHITproducer,
    bool  IsPDVD,
    bool  IsPDHD,
    bool  IsFDVD,
    bool  IsFDHD,
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
    int   fNumberConvStep,
    float fClusterSizeMulti,
    float fCovering,
    float fRadiusInt,
    float fRadiusExt,
    float fgeoYmin,
    float fgeoYmax,
    float fgeoZmin,
    float fgeoZmax,
    float fElectronVelocity, 
    float fTickTimeInMus,
    const geo::WireReadoutGeom& fWireReadout)
{
  //definition vector
  std::vector<dune::ClusterInfo*> vec;

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
  if (!e.isRealData()) vTrackIDToGeneratorTag = GetGeneratorTag( e , fG4producer , bt_serv , fVerbose);
 
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

    int fNearOrFarToTheBeam = NearOrFar(IsPDVD , IsPDHD ,IsFDVD , IsFDHD, hit);

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

    GetListOfTimeCoincidenceHit( IsPDVD, IsPDHD ,IsFDVD , IsFDHD, e, fHITproducer, fCoincidenceWd1_left, fCoincidenceWd1_right ,fCoincidenceWd2_left, fCoincidenceWd2_right , hit, lWireInd1, lWireInd2, lChannelInd1, lChannelInd2, lEnergyInd1, lEnergyInd2, lPeakTimeInd1, lPeakTimeInd2, lPeakAmpInd1, lPeakAmpInd2);

    int fCoincidence = 0;
    if ( !lWireInd1.empty() || !lWireInd2.empty() ) fCoincidence += 1;
    if ( !lWireInd1.empty() && !lWireInd2.empty() ) fCoincidence += 1;

    if ( fCoincidence > 0 )
    {
      GetListOfCrossingChannel(  IsPDVD, IsPDHD , fgeoYmin , fgeoYmax , fgeoZmin , fgeoZmax , fWire , lWireInd1 , lWireInd2 , 
			          lChannelInd1 , lEnergyInd1 , lYInd1 , lZInd1 , lChIntersectInd1 , lEIntersectInd1 , 
				  lChannelInd2 , lEnergyInd2 , lYInd2 , lZInd2 , lChIntersectInd2 , lEIntersectInd2,
				  fWireReadout); 
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

  std::vector<int> vIso = GetXYZIsolatedPoint( vYPointByEvent , vZPointByEvent , vPeakTimeColByEvent , vNoFByEvent , fElectronVelocity , fTickTimeInMus , fRadiusInt , fRadiusExt , fVerbose);
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

  point v;
  v = gen_yzt( PTSIsolated , vIso , vYPointByEvent , vZPointByEvent , vPeakTimeColByEvent ,vNoFByEvent , fElectronVelocity , fTickTimeInMus, fgeoZmax);


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
    vchecks = CheckClusters( dataPos , clustersPos , threshold , fClusterSizeMulti, fCovering , fVerbose);
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
        threshold = fMaxSizeCluster;
        K = vchecks[1] + 5;
        if ( fVerbose) printf("Threshold Max reached \n");
      }
    }
    else K = vchecks[1] + 5;
  }

  if ( fVerbose) printf("Data size : %lu x %lu \n",dataPos.size(),dataPos[0].size());
  if ( fVerbose) printf("Cluster size : %lu x %lu \n",clustersPos.size(),clustersPos[0].size());

  if ( fVerbose) printf("Data clustering ended successfully \n");

  int j = 0;
  point p;
  for (j = 0, p = v; j < PTSIsolated ; j++, p++)
  {
    p->group = reallocate( p , clustersPos , threshold);
  }
  if (fVerbose) std::cout<< " reallocation done " << std::endl;

  vec = GetCluster( fAsConverged, CRP_T0, PTSIsolated , clustersPos[0].size() , v , vEInd1PointByEvent , vEInd2PointByEvent , vChInd1PointByEvent , vChInd2PointByEvent , vEnergyColByEvent , vPeakTimeColByEvent , vChannelColByEvent , truth , vMCPDGByEvent , vMCMOMpdgByEvent , vMCWeightByEvent , vGeneratorTagByEvent , vMCXByEvent , vMCYByEvent , vMCZByEvent , vNoFByEvent , vMCEByEvent , vMCNeByEvent, fVerbose, fgeoZmax);
  //NCluster = vCluster.size();

  return vec;
}


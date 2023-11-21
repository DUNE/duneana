#ifndef CLASS_CHANNELCLUSTER_H
#define CLASS_CHANNELCLUSTER_H

#include "class_RecoHit.h"
#include "class_Cluster.h"

class channelCluster
{
  public:
    channelCluster(int cEvent, std::vector<recoHit> cHitVector, unsigned int cConfig);
    channelCluster();

    int getNClusters();
    std::vector<cluster> getClusterVector();

  private:
    int fNClusters       = 0;
    unsigned int fConfig = 0;
    std::vector<recoHit> fHitVector;
    std::vector<cluster> fVecClusters;
};

#endif


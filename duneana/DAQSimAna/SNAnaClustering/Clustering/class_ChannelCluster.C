#include "Module_SNClustering_Config.h"
#include "class_ChannelCluster.h"

channelCluster::channelCluster(int cEvent, std::vector<recoHit> cHitVector, unsigned int cConfig)
{
  fConfig = cConfig;
  fHitVector.clear();
  fHitVector = cHitVector;

  //ORDER IN THE HIT VECTOR BY TIME.
  std::sort(fHitVector.begin(), fHitVector.end(),
            [](const recoHit& lhs, const recoHit& rhs){return lhs.fHitTime < rhs.fHitTime;});

  for(unsigned int i = 0; i < fHitVector.size()-1; i++)
  {
    std::vector<recoHit> vec_TempHits;
    if(std::abs(fHitVector.at(i).getHitTime()-fHitVector.at(i+1).getHitTime())<=cut_TimeWindowSize.at(cConfig))
    {
      int  timeCount = 1;

      vec_TempHits.push_back(fHitVector.at(i));
      vec_TempHits.push_back(fHitVector.at(i+1));

      while((i+timeCount+1)<fHitVector.size() 
          && std::abs(fHitVector.at(i+timeCount).getHitTime()-fHitVector.at(i+timeCount+1).getHitTime())
          <=cut_TimeWindowSize.at(cConfig))
      {
        vec_TempHits.push_back(fHitVector.at(i + timeCount + 1));
        timeCount++;
      }

      i = i + timeCount;
      cluster temp(fHitVector.at(0).getEvent(), vec_TempHits);
      fVecClusters.push_back(temp);
    }
  }

  fNClusters = fVecClusters.size();
}

channelCluster::channelCluster(){}

int channelCluster::getNClusters(){return fNClusters;}
std::vector<cluster> channelCluster::getClusterVector(){return fVecClusters;}

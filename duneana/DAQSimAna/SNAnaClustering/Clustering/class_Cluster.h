#ifndef CLASS_CLUSTER_H
#define CLASS_CLUSTER_H

#include "class_RecoHit.h"

class cluster
{
  public:
    cluster(int cEvent, std::vector<recoHit> cHitVector);
    cluster();

    int getEvent();
    int getStartChan();
    int getEndChan();
    int getNChan();
    int getChanWidth();
    int getNHits();
    int getType();
    int getTriggerFlag();
    float getHitSADC();
    float getFirstTimeHit();
    float getLastTimeHit();
    float getTimeWidth();
    double getMC_EnergyNu();
    double getMC_EnergyLep();
    double getMC_MarlTime();
    void  printCluster();
    void  setHitSADC(float cHitSADC);
    void  setTriggerFlag(int cTriggerFlag);
    void  setMC_EnergyNu(double cMC_EnergyNu);
    void  setMC_EnergyLep(double cMC_EnergyLep);
    void  setMC_MarlTime(double cMC_Marltime);
    std::vector<recoHit> getHits();

  private:
    int fEvent       = 0;
    int fStartChan   = 0;
    int fEndChan     = 0;
    int fNChan       = 0;
    int fChanWidth   = 0;
    int fNHits       = 0;
    int fType        = 0;
    int fTriggerFlag = 0;
    float fHitSADC   = 0;
    float fFirstHitTime = 0;
    float fLastHitTime  = 0;
    float fTimeWidth    = 0;
    double fMC_EnergyNu  = 0;
    double fMC_EnergyLep = 0;
    double fMC_MarlTime  = 0;
    std::vector<recoHit> fHitVector;
};

#endif

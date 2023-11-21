#include "class_Cluster.h"

cluster::cluster(int cEvent, std::vector<recoHit> cHitVector)
{
  int type(0);
  std::set<int> channels;
  fHitVector.clear();
  fStartChan    = 50000;
  fEndChan      = 0;
  fFirstHitTime = 50000;
  fLastHitTime  = 0;
  fEvent        = cEvent;
  fHitVector    = cHitVector;
  fNHits        = fHitVector.size();
  for(int i = 0; i < fNHits; i++)
  {
    fHitSADC  += fHitVector.at(i).getHitSADC(); 
    int   chan = fHitVector.at(i).getHitChan();
    float time = fHitVector.at(i).getHitTime(); 

    channels.insert(chan);

    if(chan < fStartChan)
    {
      fStartChan = chan;
    }
    if(chan > fEndChan)
    {
      fEndChan = chan;
    }
    if(time < fFirstHitTime)
    {
      fFirstHitTime = time;
    }
    if(time > fLastHitTime)
    {
      fLastHitTime = time;
    }
    if(fHitVector.at(i).getGenType()==1)
    {
      type++;
    }
  }

  fNChan = channels.size();

  //CALL THE CLUSTER MARLEY IF THERE ARE MORE THAN TWO MARLEY HITS IN IT.
  if(type>=2)
  {
    fType = 1;
  }
  else
  {
    fType = 0;
  }

  fChanWidth = fEndChan - fStartChan;
  fTimeWidth = fLastHitTime - fFirstHitTime;

  //ORDER IN THE HIT VECTOR BY TIME.
  std::sort(fHitVector.begin(), fHitVector.end(),
            [](const recoHit& lhs, const recoHit& rhs){return lhs.fHitTime < rhs.fHitTime;});
}

cluster::cluster(){}

int cluster::getEvent(){return fEvent;}
int cluster::getStartChan(){return fStartChan;}
int cluster::getEndChan(){return fEndChan;}
int cluster::getNChan(){return fNChan;}
int cluster::getChanWidth(){return fChanWidth;}
int cluster::getNHits(){return fNHits;}
int cluster::getType(){return fType;}
int cluster::getTriggerFlag(){return fTriggerFlag;}
float cluster::getHitSADC(){return fHitSADC;}
float cluster::getFirstTimeHit(){return fFirstHitTime;}
float cluster::getLastTimeHit(){return fLastHitTime;}
float cluster::getTimeWidth(){return fTimeWidth;}
double cluster::getMC_EnergyLep(){return fMC_EnergyLep;}
double cluster::getMC_EnergyNu(){return fMC_EnergyNu;}
double cluster::getMC_MarlTime(){return fMC_MarlTime;}
void  cluster::setMC_EnergyNu(double cMC_EnergyNu){ fMC_EnergyNu = 1000*cMC_EnergyNu; return;}
void  cluster::setMC_EnergyLep(double cMC_EnergyLep){ fMC_EnergyLep = 1000*cMC_EnergyLep; return;}
void  cluster::setMC_MarlTime(double cMC_Marltime){ fMC_MarlTime = cMC_Marltime; return;}
void  cluster::setHitSADC(float cHitSADC){fHitSADC = cHitSADC; return;}
void  cluster::setTriggerFlag(int cTriggerFlag){fTriggerFlag = cTriggerFlag; return;}
std::vector<recoHit> cluster::getHits(){return fHitVector;}

void cluster::printCluster()
{
  std::cout << "********************************************************************************************************" << std::endl;
  std::cout << "THE CLUSTER BELONGS TO EVENT: "    << fEvent  << ", CONTAINS " <<  fNChan << " CHANNELS AND "
            << fNHits << " HITS WITH " << fHitSADC << "ADC. TYPE IS " << fType << ", THE HITS ARE:" << std::endl;

  for(int i = 0; i < fNHits; i++)
  {
    fHitVector.at(i).printHit();
  }

  std::cout << "********************************************************************************************************" << std::endl;

  return;
}



#ifndef CLASS_RECOHIT_H
#define CLASS_RECOHIT_H

#include <iostream>
#include <vector>
#include <set>

class recoHit
{
  public:
    recoHit(int cEvent, int cHitView, int cGenType, int cHitChan, float cHitTime, float cHitSADC, float cHitRMS);

    int   getEvent();
    int   getHitView();
    int   getGenType();
    int   getHitChan();
    float getHitTime();
    float getHitSADC();
    float getHitRMS();
    float fHitTime;
    void  printHit();

  private:
    int   fEvent = 0;
    int   fHitView = 0;
    int   fGenType = 0;
    int   fHitChan = 0;
    float fHitSADC = 0;
    float fHitRMS  = 0;
};

#endif

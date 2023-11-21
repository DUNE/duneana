//
// - Sets the navigation state to a value that the event display service
//   interprets as an action to perform when TApplication returns after 
//   receiving a ROOT signal.  Based on Mark Messier's code for Nova's
//   event display service.
//

#include "NavState.h"
#include "TROOT.h"
#include "TApplication.h"

namespace dune{

  static int gsNavState    = 0;
  static int gsTargetRun   = 0;
  static int gsTargetEvent = 0;

  //......................................................................

  int NavState::Which() { return gsNavState; }

  //......................................................................

  void NavState::Set(int which)
  {
    gsNavState = which;
    gROOT->GetApplication()->Terminate();
  }

  //......................................................................

  void NavState::SetTarget(int run, int event)
  {
    gsTargetRun = run;
    gsTargetEvent = event;
  }

  //......................................................................

  int NavState::TargetRun() { return gsTargetRun; }

  //......................................................................

  int NavState::TargetEvent() { return gsTargetEvent; }

}// namespace

// DuneEventFilter_module.cc
//
// Module to select or reject events based on event number.
//
// David Adams
// October 2018
// June 2019: Update using event range "events" from the job IndexRangeTool.
//
// Configuration:
//   LogLevel: 0=quiet, 1=ctor message, 2=status of each event
//   SelectEvents: If this has entries, the event must be among them
//   RejectEvents: Event is rejected if it is in this list.
//   EventBegin, EventEnd: if EventEnd > EventBegin, event must be in (EventBegin, EventEnd].
//   EventModFreq, EventModVal: keep events for which ievt%EventModFreq == EventModVal
//   JobIndexRangeTool: Name of the job IndexRangeTool, e.g. jobRanges.
//   SkipEventTool: Name of a tool (IndexVectorTool) that lists events to be skipped
//                  indexed by run number. Blank means no tool.
//
// The range of event to process is (EventBegin, EventEnd] if EventEnd > EventBegin.
// Otherwise the range is taken from the index range "events"  in the job IndexRangeTool.
// If that range does not exist or is not valid, all events are processed.

#include <iostream>
#include <sstream>
#include <vector>
#include <set>

#include "dunecore/DuneInterface/Tool/IndexRangeTool.h"
#include "dunecore/DuneInterface/Tool/IndexVectorMapTool.h"
#include "dunecore/ArtSupport/DuneToolManager.h"
#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Run.h" 
#include "art/Framework/Principal/SubRun.h" 

class DuneEventFilter : public art::EDFilter {

public:

  using Index = unsigned int;
  using IndexVector = std::vector<Index>;
  using IndexSet = std::set<Index>;
  using Name = std::string;

  explicit DuneEventFilter(fhicl::ParameterSet const & pset);
  virtual ~DuneEventFilter();
  virtual bool filter(art::Event& e);

private:

  // Configuration parameters.
  Index m_LogLevel;
  IndexVector m_SelectEventVector;
  IndexVector m_RejectEventVector;
  Index m_EventBegin;
  Index m_EventEnd;
  Index m_EventModFreq;
  Index m_EventModVal;
  Name m_JobIndexRangeTool;
  Name m_SkipEventTool;

  // Derived from configuration.
  Index m_beginEvent;
  Index m_endEvent;
  IndexSet m_SelectEvents;
  IndexSet m_RejectEvents;
  const IndexVectorMapTool* m_pSkipEventTool;

  // Counters.
  Index m_nproc;
  Index m_nsel;

};

//*****************************************************************************

DuneEventFilter::DuneEventFilter(fhicl::ParameterSet const & pset)
: EDFilter(pset),
  m_LogLevel(pset.get<Index>("LogLevel")),
  m_SelectEventVector(pset.get<IndexVector>("SelectEvents")),
  m_RejectEventVector(pset.get<IndexVector>("RejectEvents")),
  m_EventBegin(pset.get<Index>("EventBegin")),
  m_EventEnd(pset.get<Index>("EventEnd")),
  m_EventModFreq(pset.get<Index>("EventModFreq")),
  m_EventModVal(pset.get<Index>("EventModVal")),
  m_JobIndexRangeTool(pset.get<Name>("JobIndexRangeTool")),
  m_SkipEventTool(pset.get<Name>("SkipEventTool")),
  m_beginEvent(0), m_endEvent(0),
  m_pSkipEventTool(nullptr),
  m_nproc(0), m_nsel(0) {
  using std::cout;
  using std::endl;
  using std::string;
  const string myname = "DuneEventFilter::ctor: ";
  for ( Index ievt : m_SelectEventVector ) m_SelectEvents.insert(ievt);
  for ( Index ievt : m_RejectEventVector ) m_RejectEvents.insert(ievt);
  if ( m_LogLevel >= 1 ) {
    cout << myname << "      LogLevel: " << m_LogLevel << endl;
    cout << myname << "  SelectEvents: [";
    bool first = true;
    for ( unsigned int ievt : m_SelectEvents ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << ievt;
    }
    cout << "]" << endl;
    cout << myname << "  RejectEvents: [";
    first = true;
    for ( unsigned int ievt : m_RejectEvents ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << ievt;
    }
    cout << "]" << endl;
    cout << myname << "         EventBegin: " << m_EventBegin << endl;
    cout << myname << "           EventEnd: " << m_EventEnd << endl;
    cout << myname << "        EventModVal: " << m_EventModVal << endl;
    cout << myname << "       EventModFreq: " << m_EventModFreq << endl;
    cout << myname << "  JobIndexRangeTool: " << m_JobIndexRangeTool << endl;
    cout << myname << "      SkipEventTool: " << m_SkipEventTool << endl;
  }
  if ( m_EventEnd > m_EventBegin ) {
    m_beginEvent = m_EventBegin;
    m_endEvent = m_EventEnd;
  } else if ( m_JobIndexRangeTool.size() ) {
    DuneToolManager* ptm = DuneToolManager::instance();
    const IndexRangeTool* pjrt = ptm->getShared<IndexRangeTool>(m_JobIndexRangeTool);
    if ( pjrt == nullptr ) {
      cout << myname << "ERROR: Job index range tool not found: " << m_JobIndexRangeTool << endl;
    } else {
      IndexRange ran = pjrt->get("events");
      if ( ! ran.isValid() ) {
        cout << myname << "Event range is unrestricted." << endl;
      } else {
        m_beginEvent = ran.begin;
        m_endEvent = ran.end;
      }
    }
  }
  if ( m_endEvent > m_beginEvent ) {
    cout << myname << "Event selection range is [" << m_beginEvent << ", "
         << m_endEvent << ")." << endl;
  } else {
    cout << myname << "No event selection range." << endl;
  }
  if ( m_SkipEventTool.size() ) {
    DuneToolManager* ptm = DuneToolManager::instance();
    m_pSkipEventTool = ptm->getShared<IndexVectorMapTool>(m_SkipEventTool);
    if ( m_pSkipEventTool == nullptr ) {
      cout << "WARNING: Unable to find SkipEventTool " << m_SkipEventTool << endl;
    } else {
      cout << myname << "Using SkipEventTool @" << m_pSkipEventTool << endl;
    }
  }
}

//*****************************************************************************

bool DuneEventFilter::filter(art::Event & evt) {
  using std::cout;
  using std::endl;
  using std::string;
  const string myname = "DuneEventFilter::filter: ";
  std::ostringstream ssrev;
  Index irun = evt.getRun().id().run();
  Index isru = evt.getSubRun().id().subRun();
  ssrev << "run " << irun;
  if ( isru > 0 ) {
    ssrev << "-" << isru;
  }
  Index ievt = evt.event();
  ssrev << " event " << ievt;
  string srev = ssrev.str();
  string srevc = srev;
  srevc[0] = 'R';
  ++m_nproc;
  bool keep = true;
  if ( keep && m_SelectEvents.size() ) keep = m_SelectEvents.count(ievt);
  if ( keep ) keep = m_RejectEvents.count(ievt) == 0;
  if ( keep && m_endEvent > m_beginEvent ) keep = ievt >= m_beginEvent && ievt < m_endEvent;
  if ( keep && m_EventModFreq ) keep = (ievt % m_EventModFreq) == m_EventModVal;
  if ( keep && m_pSkipEventTool != nullptr ) {
    const IndexVector& skipEvents = m_pSkipEventTool->get(evt.run());
    if ( find(skipEvents.begin(), skipEvents.end(), ievt) != skipEvents.end() ) {
      keep = false;
      if ( m_LogLevel >= 3 ) {
        cout << myname << "  " << srevc << " rejected by SkipEventTool" << endl;
      }
    }
  }
  if ( m_LogLevel >= 2 ) {
    cout << myname << (keep ? "Sel" : "Rej") << "ecting " << srev << endl;
  }
  if ( keep ) ++m_nsel;
  return keep;
}

//*****************************************************************************

DuneEventFilter::~DuneEventFilter() {
  using std::cout;
  using std::endl;
  using std::string;
  const string myname = "DuneEventFilter::dtor: ";
  if ( m_LogLevel >= 1 ) {
    float fsel = double(m_nsel)/double(m_nproc);
    cout << myname << " Events processed: " << m_nproc << endl;
    cout << myname << "  Events selected: " << m_nsel << " (" << fsel << ")" << endl;
  }
}

//*****************************************************************************

DEFINE_ART_MODULE(DuneEventFilter)

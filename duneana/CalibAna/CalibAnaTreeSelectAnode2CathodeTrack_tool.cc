#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "ICATSelectionTool.h"

namespace dune {

class CalibAnaTreeSelectAnode2CathodeTrack: public ICATSelectionTool {
public:

  CalibAnaTreeSelectAnode2CathodeTrack(const fhicl::ParameterSet &p);
  ~CalibAnaTreeSelectAnode2CathodeTrack() {}

  bool Select(const TrackInfo &t) override;

private:
     // config
     double fTickCut;
     };
  
     CalibAnaTreeSelectAnode2CathodeTrack::CalibAnaTreeSelectAnode2CathodeTrack(const fhicl::ParameterSet &p):
       ICATSelectionTool(p),
         fTickCut(p.get<double>("TickCut"))
         {}
  
        bool CalibAnaTreeSelectAnode2CathodeTrack::Select(const TrackInfo &t) {

		// use the collection plane
		double minBeam = -1;
		double maxBeam = -1;
		double minNoBeam = -1;
		double maxNoBeam = -1;
		double hit_is_BeamSide = -1;
		double hit_is_NoBeamSide = -1;

		for (const dune::TrackHitInfo &h: t.hits2) {
			// In PDHD, TPC 1 and 5 Beam-side
		        hit_is_BeamSide = (h.h.tpc == 1 or h.h.tpc == 5);		
			hit_is_NoBeamSide = (h.h.tpc == 2 or h.h.tpc == 6);		                

			if (h.oncalo && hit_is_BeamSide == 1) {
                        	if (maxBeam < 0. || h.h.time > maxBeam) maxBeam = h.h.time;
				if (minBeam < 0. || h.h.time < minBeam) minBeam = h.h.time;
			}
			if (h.oncalo && hit_is_NoBeamSide == 1) {
			        if (maxNoBeam < 0. || h.h.time > maxNoBeam) maxNoBeam = h.h.time;
                                if (minNoBeam < 0. || h.h.time < minNoBeam) minNoBeam = h.h.time;
			}
		}		

        	return std::max(abs(maxBeam - minBeam), abs(maxNoBeam - minNoBeam)) > fTickCut;
  
        }
  
	DEFINE_ART_CLASS_TOOL(CalibAnaTreeSelectAnode2CathodeTrack)
  
} // end namespace dune

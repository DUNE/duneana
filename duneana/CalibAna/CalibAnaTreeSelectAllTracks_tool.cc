// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "ICATSelectionTool.h"

namespace dune {

class CalibAnaTreeSelectAllTracks: public ICATSelectionTool {
public:

  CalibAnaTreeSelectAllTracks(const fhicl::ParameterSet &p);
  ~CalibAnaTreeSelectAllTracks() {}

  bool Select(const TrackInfo &t) override;

private:
};

CalibAnaTreeSelectAllTracks::CalibAnaTreeSelectAllTracks(const fhicl::ParameterSet &p):
  ICATSelectionTool(p)
{}

bool CalibAnaTreeSelectAllTracks::Select(const TrackInfo &t) {
  return true;
}

DEFINE_ART_CLASS_TOOL(CalibAnaTreeSelectAllTracks)

} // end namespace dune

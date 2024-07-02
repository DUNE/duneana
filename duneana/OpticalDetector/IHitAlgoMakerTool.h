/**
 * @file   larana/OpticalDetector/IHitAlgoMakerTool.h
 * @brief  Tool interface for creating a hit finder algorithm.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   February 12, 2023
 *
 * This library is header-only.
 */

#ifndef LARANA_OPTICALDETECTOR_IHITALGOMAKERTOOL_H
#define LARANA_OPTICALDETECTOR_IHITALGOMAKERTOOL_H

// LArSoft libraries
#include "larana/OpticalDetector/OpHitFinder/PMTPulseRecoBase.h"

// C/C++ standard libraries
#include <memory>

// -----------------------------------------------------------------------------
namespace opdet {
  struct IHitAlgoMakerTool;
}
/**
 * @brief Tool interface for creating a hit finder algorithm.
 *
 * The hit finder algorithms in `larana` implement a common abstraction and
 * interface (`pmtana::PMTPulseRecoBase`).
 * In principle, wrappers may be written as _art_ tools that implement that same
 * interface and wrap the actual algorithms. In this case, users will actually
 * use the tool classes as hit finder algorithms.
 * In alternative, _art_ tools may be written to _create_ the current hit finder
 * algorithms. In this case, users will use the tool only at setup stage to
 * create the hit finder algorithms, and then the tools will have no further
 * role.
 *
 * This class provides the interface for a tool following this second design:
 * a tool following this interface will be able to create hit finder algorithm
 * objects by executing `makeAlgo()`, and that will be the only function of the
 * tool.
 */
struct opdet::IHitAlgoMakerTool {

  virtual ~IHitAlgoMakerTool() = default;

  /**
   * @brief Creates and returns a new instance of hit finder algorithm.
   * @return the newly created algorithm instance
   *
   * The returned object is completely independent of this tool: after calling
   * this function, the tool can in principle be discarded.
   *
   * Note that all the information necessary to the creation of the algorithm
   * must have already been passed to the tool (and stored) in the FHiCL
   * configuration on construction.
   */
  virtual std::unique_ptr<pmtana::PMTPulseRecoBase> makeAlgo() = 0;

}; // opdet::IHitAlgoMakerTool

// -----------------------------------------------------------------------------

#endif // LARANA_OPTICALDETECTOR_IHITALGOMAKERTOOL_H

/**
 * @file   larana/OpticalDetector/PedAlgoUBMaker_tool.cc
 * @brief  _art_ tool to create a `pmtana::PedAlgoUB` algorithm.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   February 12, 2023
 */

// LArSoft libraries
#include "larana/OpticalDetector/OpHitFinder/PedAlgoUB.h"
#include "larana/OpticalDetector/PedAlgoMakerToolBase.h"

// framework libraries
#include "art/Utilities/ToolMacros.h"

// -----------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(opdet::PedAlgoMakerToolBase<pmtana::PedAlgoUB>)

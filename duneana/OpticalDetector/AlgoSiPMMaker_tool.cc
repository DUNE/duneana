/**
 * @file   larana/OpticalDetector/AlgoSiPMMaker_tool.cc
 * @brief  _art_ tool to create a `pmtana::AlgoSiPM` algorithm.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   February 12, 2023
 */

// LArSoft libraries
#include "larana/OpticalDetector/HitAlgoMakerToolBase.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSiPM.h"

// framework libraries
#include "art/Utilities/ToolMacros.h"

// -----------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(opdet::HitAlgoMakerToolBase<pmtana::AlgoSiPM>)

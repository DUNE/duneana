// ========================================================================================
// AdjOpHits.h
// This library is based on Michael Baird's DAQSimAna_module.
// It is used to find adjacent hits in time and space to create clusters in the context of
// the SolarNuAna module and DUNE's solar neutrino analysis.
// 
// @authors     : Sergio Manthey Corchado
// @created     : Apr, 2024 
//=========================================================================================

#ifndef AdjOpHitsTool_h
#define AdjOpHitsTool_h

#include <cmath>
#include <iostream>
#include <vector>
#include <fcntl.h>

#include "art/Persistency/Common/PtrMaker.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/WireReadout.h"
#include "fhiclcpp/ParameterSet.h"

#include "TH1I.h"
#include "TH1F.h"

namespace solar
{
    class AdjOpHitsUtils
    {
        public:
            struct FlashInfo
            {
                int NHit;
                double Time;
                double TimeWidth;
                double PE;
                double MaxPE;
                std::vector<double> PEperOpDet;
                double FastToTotal;
                double X;
                double Y;
                double Z;
                double YWidth;
                double ZWidth;
            };
            explicit AdjOpHitsUtils( fhicl::ParameterSet const& p);
            void MakeFlashVector(std::vector<FlashInfo> &FlashVec, std::vector<std::vector<art::Ptr<recob::OpHit>>> &Clusters, art::Event const &evt);
            void CalcAdjOpHits(std::vector<art::Ptr<recob::OpHit>> Vec, std::vector<std::vector<art::Ptr<recob::OpHit>>> &Clusters, bool HeavDebug);
            void CalcAdjOpHitsFast(std::vector<art::Ptr<recob::OpHit>> Vec, std::vector<std::vector<art::Ptr<recob::OpHit>>> &Clusters, bool HeavDebug);
            void CalcCentroid(std::vector<art::Ptr<recob::OpHit>> Hits, double &x, double &y, double &z);
            double GaussianPDF(double x, double mean, double sigma);
            // Write a struct to store the flash information
        
        private:
            geo::WireReadoutGeom const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
            // From fhicl configuration
            const float fOpFlashAlgoTime;
            const float fOpFlashAlgoRad;
            const float fOpFlashAlgoPE;
            const float fOpFlashAlgoTriggerPE;
            const bool fOpFlashAlgoCentroid;
            const bool fOpFlashAlgoDebug;
    };
}
#endif

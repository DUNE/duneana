// -*- mode: c++; c-basic-offset: 2; -*-
/*!
 * Title:   OpFlash Algorithims
 * Authors, editors:  Ben Jones, MIT
 *                    Wes Ketchum wketchum@lanl.gov
 *                    Gleb Sinev  gleb.sinev@duke.edu
 *                    Alex Himmel ahimmel@fnal.gov
 *
 * Description:
 * These are the algorithms used by OpFlashFinder to produce flashes.
 */

#include "OpFlashAlg.h"

#include "TFile.h"
#include "TH1.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/DetectorInfo/ElecClock.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric> // std::iota()

namespace opdet {

  //----------------------------------------------------------------------------
  void writeHistogram(std::vector<double> const& binned)
  {
    TH1D* binned_histogram = new TH1D("binned_histogram",
                                      "Collection of All OpHits;Time (ms);PEs",
                                      binned.size(),
                                      0,
                                      binned.size());
    for (size_t i = 0; i < binned.size(); ++i)
      binned_histogram->SetBinContent(i, binned.at(i));

    TFile f_out("output_hist.root", "RECREATE");
    binned_histogram->Write();
    f_out.Close();

    delete binned_histogram;
  }

  //----------------------------------------------------------------------------
  void checkOnBeamFlash(std::vector<recob::OpFlash> const& FlashVector)
  {
    for (auto const& flash : FlashVector)
      if (flash.OnBeamTime() == 1)
        std::cout << "OnBeamFlash with time " << flash.Time() << std::endl;
  }

  //----------------------------------------------------------------------------
  void RunFlashFinder(std::vector<recob::OpHit> const& HitVector,
                      std::vector<recob::OpFlash>& FlashVector,
                      std::vector<std::vector<int>>& AssocList,
                      double const BinWidth,
                      geo::GeometryCore const& geom,
                      float const FlashThreshold,
                      float const WidthTolerance,
                      detinfo::DetectorClocksData const& ClocksData,
                      float const TrigCoinc)
  {
    // Initial size for accumulators - will be automatically extended if needed
    int initialsize = 6400;

    // These are the accumulators which will hold broad-binned light yields
    std::vector<double> Binned1(initialsize);
    std::vector<double> Binned2(initialsize);

    // These will keep track of which pulses put activity in each bin
    std::vector<std::vector<int>> Contributors1(initialsize);
    std::vector<std::vector<int>> Contributors2(initialsize);

    // These will keep track of where we have met the flash condition
    // (in order to prevent second pointless loop)
    std::vector<int> FlashesInAccumulator1;
    std::vector<int> FlashesInAccumulator2;

    double minTime = std::numeric_limits<float>::max();
    for (auto const& hit : HitVector)
      if (hit.PeakTime() < minTime) minTime = hit.PeakTime();

    for (auto const& hit : HitVector) {

      double peakTime = hit.PeakTime();

      unsigned int AccumIndex1 = GetAccumIndex(peakTime, minTime, BinWidth, 0.0);

      unsigned int AccumIndex2 = GetAccumIndex(peakTime, minTime, BinWidth, BinWidth / 2.0);

      // Extend accumulators if needed (2 always larger than 1)
      if (AccumIndex2 >= Binned1.size()) {
        std::cout << "Extending vectors to " << AccumIndex2 * 1.2 << std::endl;
        Binned1.resize(AccumIndex2 * 1.2);
        Binned2.resize(AccumIndex2 * 1.2);
        Contributors1.resize(AccumIndex2 * 1.2);
        Contributors2.resize(AccumIndex2 * 1.2);
      }

      size_t const hitIndex = &hit - &HitVector[0];

      FillAccumulator(AccumIndex1,
                      hitIndex,
                      hit.PE(),
                      FlashThreshold,
                      Binned1,
                      Contributors1,
                      FlashesInAccumulator1);

      FillAccumulator(AccumIndex2,
                      hitIndex,
                      hit.PE(),
                      FlashThreshold,
                      Binned2,
                      Contributors2,
                      FlashesInAccumulator2);

    } // End loop over hits

    // Now start to create flashes.
    // First, need vector to keep track of which hits belong to which flashes
    std::vector<std::vector<int>> HitsPerFlash;

    //if (Frame == 1) writeHistogram(Binned1);

    AssignHitsToFlash(FlashesInAccumulator1,
                      FlashesInAccumulator2,
                      Binned1,
                      Binned2,
                      Contributors1,
                      Contributors2,
                      HitVector,
                      HitsPerFlash,
                      FlashThreshold);

    // Now we do the fine grained part.
    // Subdivide each flash into sub-flashes with overlaps within hit widths
    // (assumed wider than photon travel time)
    std::vector<std::vector<int>> RefinedHitsPerFlash;
    for (auto const& HitsThisFlash : HitsPerFlash)
      RefineHitsInFlash(
        HitsThisFlash, HitVector, RefinedHitsPerFlash, WidthTolerance, FlashThreshold);

    // Now we have all our hits assigned to a flash.
    // Make the recob::OpFlash objects
    for (auto const& HitsPerFlashVec : RefinedHitsPerFlash)
      ConstructFlash(HitsPerFlashVec, HitVector, FlashVector, geom, ClocksData, TrigCoinc);

    RemoveLateLight(FlashVector, RefinedHitsPerFlash);

    //checkOnBeamFlash(FlashVector);

    // Finally, write the association list.
    // back_inserter tacks the result onto the end of AssocList
    for (auto& HitIndicesThisFlash : RefinedHitsPerFlash)
      AssocList.push_back(HitIndicesThisFlash);

  } // End RunFlashFinder

  //----------------------------------------------------------------------------
  unsigned int GetAccumIndex(double const PeakTime,
                             double const MinTime,
                             double const BinWidth,
                             double const BinOffset)
  {
    return static_cast<unsigned int>((PeakTime - MinTime + BinOffset) / BinWidth);
  }

  //----------------------------------------------------------------------------
  void FillAccumulator(unsigned int const& AccumIndex,
                       unsigned int const& HitIndex,
                       double const PE,
                       float const FlashThreshold,
                       std::vector<double>& Binned,
                       std::vector<std::vector<int>>& Contributors,
                       std::vector<int>& FlashesInAccumulator)
  {

    Contributors.at(AccumIndex).push_back(HitIndex);

    Binned.at(AccumIndex) += PE;

    // If this wasn't a flash already, add it to the list
    if (Binned.at(AccumIndex) >= FlashThreshold && (Binned.at(AccumIndex) - PE) < FlashThreshold)
      FlashesInAccumulator.push_back(AccumIndex);
  }

  //----------------------------------------------------------------------------
  void FillFlashesBySizeMap(
    std::vector<int> const& FlashesInAccumulator,
    std::vector<double> const& BinnedPE,
    int const& Accumulator,
    std::map<double, std::map<int, std::vector<int>>, std::greater<double>>& FlashesBySize)
  {
    for (auto const& flash : FlashesInAccumulator)
      FlashesBySize[BinnedPE.at(flash)][Accumulator].push_back(flash);
  }

  //----------------------------------------------------------------------------
  void FillHitsThisFlash(std::vector<std::vector<int>> const& Contributors,
                         int const& Bin,
                         std::vector<int> const& HitClaimedByFlash,
                         std::vector<int>& HitsThisFlash)
  {
    // For each hit in the flash
    for (auto const& HitIndex : Contributors.at(Bin))
      // If unclaimed, claim it
      if (HitClaimedByFlash.at(HitIndex) == -1) HitsThisFlash.push_back(HitIndex);
  }

  //----------------------------------------------------------------------------
  void ClaimHits(std::vector<recob::OpHit> const& HitVector,
                 std::vector<int> const& HitsThisFlash,
                 float const FlashThreshold,
                 std::vector<std::vector<int>>& HitsPerFlash,
                 std::vector<int>& HitClaimedByFlash)
  {
    // Check for newly claimed hits
    double PE = 0;
    for (auto const& Hit : HitsThisFlash)
      PE += HitVector.at(Hit).PE();

    if (PE < FlashThreshold) return;

    // Add the flash to the list
    HitsPerFlash.push_back(HitsThisFlash);

    // And claim all the hits
    for (auto const& Hit : HitsThisFlash)
      if (HitClaimedByFlash.at(Hit) == -1) HitClaimedByFlash.at(Hit) = HitsPerFlash.size() - 1;
  }

  //----------------------------------------------------------------------------
  void AssignHitsToFlash(std::vector<int> const& FlashesInAccumulator1,
                         std::vector<int> const& FlashesInAccumulator2,
                         std::vector<double> const& Binned1,
                         std::vector<double> const& Binned2,
                         std::vector<std::vector<int>> const& Contributors1,
                         std::vector<std::vector<int>> const& Contributors2,
                         std::vector<recob::OpHit> const& HitVector,
                         std::vector<std::vector<int>>& HitsPerFlash,
                         float const FlashThreshold)
  {
    // Sort all the flashes found by size. The structure is:
    // FlashesBySize[flash size][accumulator_num] = [flash_index1, flash_index2...]
    std::map<double, std::map<int, std::vector<int>>, std::greater<double>> FlashesBySize;

    // Sort the flashes by size using map
    FillFlashesBySizeMap(FlashesInAccumulator1, Binned1, 1, FlashesBySize);
    FillFlashesBySizeMap(FlashesInAccumulator2, Binned2, 2, FlashesBySize);

    // This keeps track of which hits are claimed by which flash
    std::vector<int> HitClaimedByFlash(HitVector.size(), -1);

    // Walk from largest to smallest, claiming hits.
    // The biggest flash always gets dibbs,
    // but we keep track of overlaps for re-merging later (do we? ---WK)
    for (auto const& itFlash : FlashesBySize)
      // If several with same size, walk through accumulators
      for (auto const& itAcc : itFlash.second) {

        int Accumulator = itAcc.first;

        // Walk through flash-tagged bins in this accumulator
        for (auto const& Bin : itAcc.second) {

          std::vector<int> HitsThisFlash;

          if (Accumulator == 1)
            FillHitsThisFlash(Contributors1, Bin, HitClaimedByFlash, HitsThisFlash);
          else if (Accumulator == 2)
            FillHitsThisFlash(Contributors2, Bin, HitClaimedByFlash, HitsThisFlash);

          ClaimHits(HitVector, HitsThisFlash, FlashThreshold, HitsPerFlash, HitClaimedByFlash);

        } // End loop over this accumulator

      } // End loops over accumulators
    // End of loops over sorted flashes

  } // End AssignHitsToFlash

  //----------------------------------------------------------------------------
  void FindSeedHit(std::map<double, std::vector<int>, std::greater<double>> const& HitsBySize,
                   std::vector<bool>& HitsUsed,
                   std::vector<recob::OpHit> const& HitVector,
                   std::vector<int>& HitsThisRefinedFlash,
                   double& PEAccumulated,
                   double& FlashMaxTime,
                   double& FlashMinTime)
  {
    for (auto const& itHit : HitsBySize)
      for (auto const& HitID : itHit.second) {

        if (HitsUsed.at(HitID)) continue;

        PEAccumulated = HitVector.at(HitID).PE();
        FlashMaxTime = HitVector.at(HitID).PeakTime() + 0.5 * HitVector.at(HitID).Width();
        FlashMinTime = HitVector.at(HitID).PeakTime() - 0.5 * HitVector.at(HitID).Width();

        HitsThisRefinedFlash.clear();
        HitsThisRefinedFlash.push_back(HitID);

        HitsUsed.at(HitID) = true;
        return;

      } // End loop over inner vector
    // End loop over HitsBySize map

  } // End FindSeedHit

  //----------------------------------------------------------------------------
  void AddHitToFlash(int const& HitID,
                     std::vector<bool>& HitsUsed,
                     recob::OpHit const& currentHit,
                     double const WidthTolerance,
                     std::vector<int>& HitsThisRefinedFlash,
                     double& PEAccumulated,
                     double& FlashMaxTime,
                     double& FlashMinTime)
  {
    if (HitsUsed.at(HitID)) return;

    double HitTime = currentHit.PeakTime();
    double HitWidth = 0.5 * currentHit.Width();
    double FlashTime = 0.5 * (FlashMaxTime + FlashMinTime);
    double FlashWidth = 0.5 * (FlashMaxTime - FlashMinTime);

    if (std::abs(HitTime - FlashTime) > WidthTolerance * (HitWidth + FlashWidth)) return;

    HitsThisRefinedFlash.push_back(HitID);
    FlashMaxTime = std::max(FlashMaxTime, HitTime + HitWidth);
    FlashMinTime = std::min(FlashMinTime, HitTime - HitWidth);
    PEAccumulated += currentHit.PE();
    HitsUsed[HitID] = true;

  } // End AddHitToFlash

  //----------------------------------------------------------------------------
  void CheckAndStoreFlash(std::vector<std::vector<int>>& RefinedHitsPerFlash,
                          std::vector<int> const& HitsThisRefinedFlash,
                          double const PEAccumulated,
                          float const FlashThreshold,
                          std::vector<bool>& HitsUsed)
  {
    // If above threshold, we just add hits to the flash vector, and move on
    if (PEAccumulated >= FlashThreshold) {
      RefinedHitsPerFlash.push_back(HitsThisRefinedFlash);
      return;
    }

    // If there is only one hit in collection, we can immediately move on
    if (HitsThisRefinedFlash.size() == 1) return;

    // We need to release all other hits (allow possible reuse)
    for (std::vector<int>::const_iterator hitIterator = std::next(HitsThisRefinedFlash.begin());
         hitIterator != HitsThisRefinedFlash.end();
         ++hitIterator)
      HitsUsed.at(*hitIterator) = false;

  } // End CheckAndStoreFlash

  //----------------------------------------------------------------------------
  void RefineHitsInFlash(std::vector<int> const& HitsThisFlash,
                         std::vector<recob::OpHit> const& HitVector,
                         std::vector<std::vector<int>>& RefinedHitsPerFlash,
                         float const WidthTolerance,
                         float const FlashThreshold)
  {
    // Sort the hits by their size using map
    // HitsBySize[HitSize] = [hit1, hit2 ...]
    std::map<double, std::vector<int>, std::greater<double>> HitsBySize;
    for (auto const& HitID : HitsThisFlash)
      HitsBySize[HitVector.at(HitID).PE()].push_back(HitID);

    // Heres what we do:
    //  1.Start with the biggest remaining hit
    //  2.Look for any within one width of this hit
    //  3.Find the new upper and lower bounds of the flash
    //  4.Collect again
    //  5.Repeat until no new hits collected
    //  6.Remove these hits from consideration and repeat

    std::vector<bool> HitsUsed(HitVector.size(), false);
    double PEAccumulated, FlashMaxTime, FlashMinTime;
    std::vector<int> HitsThisRefinedFlash;

    while (true) {

      HitsThisRefinedFlash.clear();
      PEAccumulated = 0;
      FlashMaxTime = 0;
      FlashMinTime = 0;

      FindSeedHit(HitsBySize,
                  HitsUsed,
                  HitVector,
                  HitsThisRefinedFlash,
                  PEAccumulated,
                  FlashMaxTime,
                  FlashMinTime);

      if (HitsThisRefinedFlash.size() == 0) return;

      // Start this at zero to do the while at least once
      size_t NHitsThisRefinedFlash = 0;

      // If size of HitsThisRefinedFlash is same size,
      // that means we're not adding anymore
      while (NHitsThisRefinedFlash < HitsThisRefinedFlash.size()) {
        NHitsThisRefinedFlash = HitsThisRefinedFlash.size();

        for (auto const& itHit : HitsBySize)
          for (auto const& HitID : itHit.second)
            AddHitToFlash(HitID,
                          HitsUsed,
                          HitVector.at(HitID),
                          WidthTolerance,
                          HitsThisRefinedFlash,
                          PEAccumulated,
                          FlashMaxTime,
                          FlashMinTime);
      }

      // We did our collecting, now check if the flash is
      // still good and push back
      CheckAndStoreFlash(
        RefinedHitsPerFlash, HitsThisRefinedFlash, PEAccumulated, FlashThreshold, HitsUsed);

    } // End while there are hits left

  } // End RefineHitsInFlash

  //----------------------------------------------------------------------------
  void AddHitContribution(recob::OpHit const& currentHit,
                          double& MaxTime,
                          double& MinTime,
                          double& AveTime,
                          double& FastToTotal,
                          double& AveAbsTime,
                          double& TotalPE,
                          std::vector<double>& PEs)
  {
    double PEThisHit = currentHit.PE();
    double TimeThisHit = currentHit.PeakTime();
    if (TimeThisHit > MaxTime) MaxTime = TimeThisHit;
    if (TimeThisHit < MinTime) MinTime = TimeThisHit;

    // These quantities for the flash are defined
    // as the weighted averages over the hits
    AveTime += PEThisHit * TimeThisHit;
    FastToTotal += PEThisHit * currentHit.FastToTotal();
    AveAbsTime += PEThisHit * currentHit.PeakTimeAbs();

    // These are totals
    TotalPE += PEThisHit;
    PEs.at(static_cast<unsigned int>(currentHit.OpChannel())) += PEThisHit;
  }

  //----------------------------------------------------------------------------
  void GetHitGeometryInfo(recob::OpHit const& currentHit,
                          geo::GeometryCore const& geom,
                          std::vector<double>& sumw,
                          std::vector<double>& sumw2,
                          double& sumy,
                          double& sumy2,
                          double& sumz,
                          double& sumz2)
  {
    auto const xyz = geom.OpDetGeoFromOpChannel(currentHit.OpChannel()).GetCenter();
    double PEThisHit = currentHit.PE();

    geo::TPCID tpc = geom.FindTPCAtPosition(xyz);
    // if the point does not fall into any TPC,
    // it does not contribute to the average wire position
    if (tpc.isValid) {
      for (size_t p = 0; p != geom.Nplanes(); ++p) {
        geo::PlaneID const planeID(tpc, p);
        unsigned int w = geom.NearestWireID(xyz, planeID).Wire;
        sumw.at(p) += PEThisHit * w;
        sumw2.at(p) += PEThisHit * w * w;
      }
    } // if we found the TPC
    sumy += PEThisHit * xyz.Y();
    sumy2 += PEThisHit * xyz.Y() * xyz.Y();
    sumz += PEThisHit * xyz.Z();
    sumz2 += PEThisHit * xyz.Z() * xyz.Z();
  }

  //----------------------------------------------------------------------------
  double CalculateWidth(double const sum, double const sum_squared, double const weights_sum)
  {
    if (sum_squared * weights_sum - sum * sum < 0)
      return 0;
    else
      return std::sqrt(sum_squared * weights_sum - sum * sum) / weights_sum;
  }

  //----------------------------------------------------------------------------
  void ConstructFlash(std::vector<int> const& HitsPerFlashVec,
                      std::vector<recob::OpHit> const& HitVector,
                      std::vector<recob::OpFlash>& FlashVector,
                      geo::GeometryCore const& geom,
                      detinfo::DetectorClocksData const& ClocksData,
                      float const TrigCoinc)
  {
    double MaxTime = -std::numeric_limits<double>::max();
    double MinTime = std::numeric_limits<double>::max();

    std::vector<double> PEs(geom.MaxOpChannel() + 1, 0.0);
    unsigned int Nplanes = geom.Nplanes();
    std::vector<double> sumw(Nplanes, 0.0);
    std::vector<double> sumw2(Nplanes, 0.0);

    double TotalPE = 0;
    double AveTime = 0;
    double AveAbsTime = 0;
    double FastToTotal = 0;
    double sumy = 0;
    double sumz = 0;
    double sumy2 = 0;
    double sumz2 = 0;

    for (auto const& HitID : HitsPerFlashVec) {
      AddHitContribution(
        HitVector.at(HitID), MaxTime, MinTime, AveTime, FastToTotal, AveAbsTime, TotalPE, PEs);
      GetHitGeometryInfo(HitVector.at(HitID), geom, sumw, sumw2, sumy, sumy2, sumz, sumz2);
    }

    AveTime /= TotalPE;
    AveAbsTime /= TotalPE;
    FastToTotal /= TotalPE;

    double meany = sumy / TotalPE;
    double meanz = sumz / TotalPE;

    double widthy = CalculateWidth(sumy, sumy2, TotalPE);
    double widthz = CalculateWidth(sumz, sumz2, TotalPE);

    std::vector<double> WireCenters(Nplanes, 0.0);
    std::vector<double> WireWidths(Nplanes, 0.0);

    for (size_t p = 0; p != Nplanes; ++p) {
      WireCenters.at(p) = sumw.at(p) / TotalPE;
      WireWidths.at(p) = CalculateWidth(sumw.at(p), sumw2.at(p), TotalPE);
    }

    // Emprical corrections to get the Frame right.
    // Eventual solution - remove frames
    int Frame = ClocksData.OpticalClock().Frame(AveAbsTime - 18.1);
    if (Frame == 0) Frame = 1;

    int BeamFrame = ClocksData.OpticalClock().Frame(ClocksData.TriggerTime());
    bool InBeamFrame = false;
    if (!(ClocksData.TriggerTime() < 0)) InBeamFrame = (Frame == BeamFrame);

    double TimeWidth = (MaxTime - MinTime) / 2.0;

    int OnBeamTime = 0;
    if (InBeamFrame && (std::abs(AveTime) < TrigCoinc)) OnBeamTime = 1;

    FlashVector.emplace_back(AveTime,
                             TimeWidth,
                             AveAbsTime,
                             Frame,
                             PEs,
                             InBeamFrame,
                             OnBeamTime,
                             FastToTotal,
                             meany,
                             widthy,
                             meanz,
                             widthz,
                             WireCenters,
                             WireWidths);
  }

  //----------------------------------------------------------------------------
  double GetLikelihoodLateLight(double const iPE,
                                double const iTime,
                                double const iWidth,
                                double const jPE,
                                double const jTime,
                                double const jWidth)
  {
    if (iTime > jTime) return 1e6;

    // Calculate hypothetical PE if this were actually a late flash from i.
    // Argon time const is 1600 ns, so 1.6.
    double HypPE = iPE * jWidth / iWidth * std::exp(-(jTime - iTime) / 1.6);
    double nsigma = (jPE - HypPE) / std::sqrt(HypPE);
    return nsigma;
  }

  //----------------------------------------------------------------------------
  void MarkFlashesForRemoval(std::vector<recob::OpFlash> const& FlashVector,
                             size_t const BeginFlash,
                             std::vector<bool>& MarkedForRemoval)
  {
    for (size_t iFlash = BeginFlash; iFlash != FlashVector.size(); ++iFlash) {

      double iTime = FlashVector.at(iFlash).Time();
      double iPE = FlashVector.at(iFlash).TotalPE();
      double iWidth = FlashVector.at(iFlash).TimeWidth();

      for (size_t jFlash = iFlash + 1; jFlash != FlashVector.size(); ++jFlash) {

        if (MarkedForRemoval.at(jFlash - BeginFlash)) continue;

        double jTime = FlashVector.at(jFlash).Time();
        double jPE = FlashVector.at(jFlash).TotalPE();
        double jWidth = FlashVector.at(jFlash).TimeWidth();

        // If smaller than, or within 2sigma of expectation,
        // attribute to late light and toss out
        if (GetLikelihoodLateLight(iPE, iTime, iWidth, jPE, jTime, jWidth) < 3.0)
          MarkedForRemoval.at(jFlash - BeginFlash) = true;
      }
    }
  }

  //----------------------------------------------------------------------------
  void RemoveFlashesFromVectors(std::vector<bool> const& MarkedForRemoval,
                                std::vector<recob::OpFlash>& FlashVector,
                                size_t const BeginFlash,
                                std::vector<std::vector<int>>& RefinedHitsPerFlash)
  {
    for (int iFlash = MarkedForRemoval.size() - 1; iFlash != -1; --iFlash)
      if (MarkedForRemoval.at(iFlash)) {
        RefinedHitsPerFlash.erase(RefinedHitsPerFlash.begin() + iFlash);
        FlashVector.erase(FlashVector.begin() + BeginFlash + iFlash);
      }
  }

  //----------------------------------------------------------------------------
  void RemoveLateLight(std::vector<recob::OpFlash>& FlashVector,
                       std::vector<std::vector<int>>& RefinedHitsPerFlash)
  {
    std::vector<bool> MarkedForRemoval(RefinedHitsPerFlash.size(), false);

    size_t const BeginFlash = FlashVector.size() - RefinedHitsPerFlash.size();

    recob::OpFlashSortByTime sort_flash_by_time;

    // Determine the sort of FlashVector starting at BeginFlash
    auto sort_order = sort_permutation(FlashVector, BeginFlash, sort_flash_by_time);

    // Sort the RefinedHitsPerFlash in the same way as tail end of FlashVector
    apply_permutation(RefinedHitsPerFlash, sort_order);

    std::sort(FlashVector.begin() + BeginFlash, FlashVector.end(), sort_flash_by_time);

    MarkFlashesForRemoval(FlashVector, BeginFlash, MarkedForRemoval);

    RemoveFlashesFromVectors(MarkedForRemoval, FlashVector, BeginFlash, RefinedHitsPerFlash);

  } // End RemoveLateLight

  //----------------------------------------------------------------------------
  template <typename T, typename Compare>
  std::vector<int> sort_permutation(std::vector<T> const& vec, int offset, Compare compare)
  {

    std::vector<int> p(vec.size() - offset);
    std::iota(p.begin(), p.end(), 0);
    std::sort(
      p.begin(), p.end(), [&](int i, int j) { return compare(vec[i + offset], vec[j + offset]); });
    return p;
  }

  //----------------------------------------------------------------------------
  template <typename T>
  void apply_permutation(std::vector<T>& vec, std::vector<int> const& p)
  {

    std::vector<T> sorted_vec(p.size());
    std::transform(p.begin(), p.end(), sorted_vec.begin(), [&](int i) { return vec[i]; });
    vec = sorted_vec;
  }

} // End namespace opdet

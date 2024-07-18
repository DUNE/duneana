#include "AdjOpHitsUtils.h"

#include "larcorealg/Geometry/OpDetGeo.h"

namespace solar
{
  AdjOpHitsUtils::AdjOpHitsUtils(fhicl::ParameterSet const &p)
      : fOpFlashAlgoTime(p.get<float>("OpFlashAlgoTime")),
        fOpFlashAlgoRad(p.get<float>("OpFlashAlgoRad")),
        fOpFlashAlgoPE(p.get<float>("OpFlashAlgoPE")),
        fOpFlashAlgoTriggerPE(p.get<float>("OpFlashAlgoTriggerPE")),
        fOpFlashAlgoCentroid(p.get<bool>("OpFlashAlgoCentroid")),
        fOpFlashAlgoDebug(p.get<bool>("OpFlashAlgoDebug"))
  {
  }
  void AdjOpHitsUtils::MakeFlashVector(std::vector<FlashInfo> &FlashVec, std::vector<std::vector<art::Ptr<recob::OpHit>>> &Clusters, art::Event const &evt){
    // This is the constructor for the OpFlash object
    // OpFlash (double time, double timewidth, double abstime, unsigned int frame, std::vector< double > PEperOpDet, bool InBeamFrame=0, int OnBeamTime=0, double FastToTotal=1, double yCenter=0, double yWidth=0, double zCenter=0, double zWidth=0, std::vector< double > WireCenters=std::vector< double >(0), std::vector< double > WireWidths=std::vector< double >(0))
    auto assns = std::make_unique<art::Assns<recob::OpFlash, recob::OpHit>>(); 
    
    for (std::vector<art::Ptr<recob::OpHit>> Cluster : Clusters)
    {
      std::vector<art::Ptr<recob::OpHit>> ClusterCopy = Cluster;
      std::sort(ClusterCopy.begin(), ClusterCopy.end(), [](art::Ptr<recob::OpHit> a, art::Ptr<recob::OpHit> b) { return a->PeakTime() < b->PeakTime(); });
      int NHit = 0;
      double Time = ClusterCopy[0]->PeakTime();
      double TimeWidth = 0;
      double TimeSum = 0;
      double PE = 0;
      double MaxPE = 0;
      std::vector<double> PEperOpDet = {};
      double FastToTotal = 1;
      double X = 0;
      double Y = 0;
      double Z = 0;
      double YWidth = 0;
      double ZWidth = 0;
      double XSum = 0;
      double YSum = 0;
      double ZSum = 0;
      
      for (art::Ptr<recob::OpHit> PDSHit : ClusterCopy)
      {
        NHit++;
        PE += PDSHit->PE();
        if (PDSHit->PE() > MaxPE)
          MaxPE = PDSHit->PE();
        PEperOpDet.push_back(PDSHit->PE());
        TimeSum += PDSHit->PeakTime() * PDSHit->PE();
        auto OpHitXYZ = wireReadout.OpDetGeoFromOpChannel(PDSHit->OpChannel()).GetCenter();
        XSum += OpHitXYZ.X() * PDSHit->PE();
        YSum += OpHitXYZ.Y() * PDSHit->PE();
        ZSum += OpHitXYZ.Z() * PDSHit->PE();
      }
      Time = TimeSum / PE;
      X = XSum / PE;
      Y = YSum / PE;
      Z = ZSum / PE;
      // Alternatively compute the centroid of the flash in 3D space
      if (fOpFlashAlgoCentroid) 
      { 
        CalcCentroid(ClusterCopy, X, Y, Z);
      }

      for (art::Ptr<recob::OpHit> PDSHit : ClusterCopy)
      {
        auto OpHitXYZ = wireReadout.OpDetGeoFromOpChannel(PDSHit->OpChannel()).GetCenter();
        TimeWidth += (PDSHit->PeakTime() - Time) * (PDSHit->PeakTime() - Time);
        YWidth += (OpHitXYZ.Y() - Y) * (OpHitXYZ.Y() - Y);
        ZWidth += (OpHitXYZ.Z() - Z) * (OpHitXYZ.Z() - Z);
      }
      TimeWidth = sqrt(TimeWidth / ClusterCopy.size());
      YWidth = sqrt(YWidth / ClusterCopy.size());
      ZWidth = sqrt(ZWidth / ClusterCopy.size());
      
      // Compute FastToTotal according to the #PEs arriving within the first 10% of the time window wrt the total #PEs
      for (art::Ptr<recob::OpHit> PDSHit : ClusterCopy)
      {
        if (PDSHit->PeakTime() < Time + TimeWidth / 10)
          FastToTotal += PDSHit->PE();
      }
      FastToTotal /= PE;
      FlashVec.push_back(FlashInfo{NHit, Time, TimeWidth, PE, MaxPE, PEperOpDet, FastToTotal, X, Y, Z, YWidth, ZWidth});
    }
    return;
  }

  void AdjOpHitsUtils::CalcAdjOpHitsFast(std::vector<art::Ptr<recob::OpHit>> Vec, std::vector<std::vector<art::Ptr<recob::OpHit>>> &Clusters, bool HeavDebug)
  {
    const float TimeRange = fOpFlashAlgoTime;       // Time in ns
    const float RadRange  = fOpFlashAlgoRad;        // Range in cm
    const float MinPE     = fOpFlashAlgoPE;         // Minimum PE for a hit to be clustered
    const float TriggerPE = fOpFlashAlgoTriggerPE;  // Minimum PE for a hit to trigger a flash

    // Define MyVec as a copy of the input vector but only with hits with PE > MinPE
    std::vector<art::Ptr<recob::OpHit>> MyVec;
    for (const auto &hit : Vec)
    {
      if (hit->PE() >= MinPE)
        MyVec.push_back(hit);
    }

    // Sort hits according to time
    std::sort(MyVec.begin(), MyVec.end(), [](art::Ptr<recob::OpHit> a, art::Ptr<recob::OpHit> b) { return a->PeakTime() < b->PeakTime(); });
    if (HeavDebug) std::cout << "Selected ophits " << MyVec.size() << " from " << Vec.size() << std::endl;

    // Pre-calculate OpDet center coordinates
    std::unordered_map<int, TVector3> opDetCenters;
    for (const auto &hit : MyVec)
    {
      int opChannel = hit->OpChannel();
      if (opDetCenters.find(opChannel) == opDetCenters.end())
      {
        auto opDetXYZ = wireReadout.OpDetGeoFromOpChannel(opChannel).GetCenter();
        opDetCenters[opChannel] = TVector3(opDetXYZ.X(), opDetXYZ.Y(), opDetXYZ.Z());
      }
    }
    // Create a vector of bools to track if a hit has been clustered or not
    std::vector<bool> ClusteredHits(MyVec.size(), false);
    // Don't need cluster all the hits, only those with PE > MinPE that are close to a big hit
    for (auto it = MyVec.begin(); it != MyVec.end(); ++it)
    {
      const auto &hit = *it;
      if (hit->PE() < TriggerPE) {continue;}
      bool main_hit = true;
      // If a trigger hit is found, start a new cluster with the hits around it that are within the time and radius range
      std::vector<art::Ptr<recob::OpHit>> AdjHitVec = {};
      AdjHitVec.push_back(hit);
      if (HeavDebug) std::cout << "Trigger hit found: CH " << hit->OpChannel() << " Time " << hit->PeakTime() << std::endl;

      // Make use of the fact that the hits are sorted in time to only consider the hits that are adjacent in the vector up to a certain time range
      for (auto it2 = it + 1; it2 != MyVec.end(); ++it2)
      {
        const auto &adjHit = *it2; // Update adjHit here
        if (std::abs(adjHit->PeakTime() - hit->PeakTime()) > TimeRange)
          break;
        // If sign of x is the same, then the two hits are in the same drift volume and can be clustered, else skip
        if (opDetCenters[hit->OpChannel()].X() * opDetCenters[adjHit->OpChannel()].X() < 0)
          continue;
        // If hit has already been clustered, skip
        if (ClusteredHits[std::distance(MyVec.begin(), it2)])
          continue;
        if ((opDetCenters[hit->OpChannel()] - opDetCenters[adjHit->OpChannel()]).Mag() < RadRange)
        {
          if (adjHit->PE() > hit->PE())
          {
            if (HeavDebug) std::cout << "Hit with PE > TriggerPE found: CH " << adjHit->OpChannel() << " Time " << adjHit->PeakTime() << std::endl;
            main_hit = false;
            break;
          }
          if (HeavDebug) std::cout << "Adding hit: CH " << adjHit->OpChannel() << " Time " << adjHit->PeakTime() << std::endl;
          AdjHitVec.push_back(adjHit);
          ClusteredHits[std::distance(MyVec.begin(), it2)] = true;
        }
      }
      for (auto it3 = it - 1; it3 != MyVec.begin(); --it3)
      {
        const auto &adjHit = *it3;
        if (std::abs(adjHit->PeakTime() - hit->PeakTime()) > TimeRange)
          break;
        if (opDetCenters[hit->OpChannel()].X() * opDetCenters[adjHit->OpChannel()].X() < 0)
          continue;
        // if hit has already been clustered, skip
        if (ClusteredHits[std::distance(MyVec.begin(), it3)])
          continue;
        if ((opDetCenters[hit->OpChannel()] - opDetCenters[adjHit->OpChannel()]).Mag() < RadRange)
        {
          if (adjHit->PE() > hit->PE())
          {
            if (HeavDebug) std::cout << "Hit with PE > TriggerPE found: CH " << adjHit->OpChannel() << " Time " << adjHit->PeakTime() << std::endl;
            main_hit = false;
            break;
          }
          AdjHitVec.push_back(adjHit);
          ClusteredHits[std::distance(MyVec.begin(), it3)] = true;
          if (HeavDebug) std::cout << "Adding hit: CH " << adjHit->OpChannel() << " Time " << adjHit->PeakTime() << std::endl;
        }
      }
      if (main_hit)
      {
        Clusters.push_back(std::move(AdjHitVec));
        if (HeavDebug) std::cout << "Cluster size: " << Clusters.back().size() << "\n" << std::endl;
      }
    }
    return;
  }

  void AdjOpHitsUtils::CalcAdjOpHits(std::vector<art::Ptr<recob::OpHit>> Vec, std::vector<std::vector<art::Ptr<recob::OpHit>>> &Clusters, bool HeavDebug)
  {
    const float TimeRange = fOpFlashAlgoTime; // Time in ns
    const float RadRange = fOpFlashAlgoRad;   // Range in cm
    const float MinPE = fOpFlashAlgoPE;       // Minimum PE for a hit to be considered
    unsigned int FilledHits = 0;
    
    // Define MyVec as a copy of the input vector but only with hits with PE > MinPE
    std::vector<art::Ptr<recob::OpHit>> MyVec;
    for (const auto &hit : Vec)
    {
      if (hit->PE() >= MinPE)
        MyVec.push_back(hit);
    }
    unsigned int NumOriHits = MyVec.size();

    // Pre-calculate OpDet center coordinates
    std::unordered_map<int, TVector3> opDetCenters;
    for (const auto &hit : MyVec)
    {
      int opChannel = hit->OpChannel();
      if (opDetCenters.find(opChannel) == opDetCenters.end())
      {
        auto opDetXYZ = wireReadout.OpDetGeoFromOpChannel(opChannel).GetCenter();
        opDetCenters[opChannel] = TVector3(opDetXYZ.X(), opDetXYZ.Y(), opDetXYZ.Z());
      }
    }

    while (NumOriHits != FilledHits)
    {
      if (HeavDebug)
        std::cerr << "\nStart of my while loop" << std::endl;

      std::vector<art::Ptr<recob::OpHit>> AdjHitVec;
      AdjHitVec.reserve(MyVec.size());
      AdjHitVec.push_back(MyVec[0]);
      MyVec.erase(MyVec.begin() + 0);
      int LastSize = 0;
      int NewSize = AdjHitVec.size();

      while (LastSize != NewSize)
      {
        std::vector<int> AddNow;
        for (const auto &adjHit : AdjHitVec)
        {
          for (auto it = MyVec.begin(); it != MyVec.end();)
          {
            const auto &myHit = *it;
            if (HeavDebug)
            {
              std::cerr << "Looping though AdjVec and MyVec: AdjHitVec - " << adjHit->OpChannel() << " & " << adjHit->PeakTime() << std::endl
              << "MVec - " << myHit->OpChannel() << " & " << myHit->PeakTime() << std::endl
              << "Time " << std::abs(adjHit->PeakTime() - myHit->PeakTime()) << " bool " << (std::abs(adjHit->PeakTime() - myHit->PeakTime()) <= TimeRange)
              << std::endl;
            }

            int adjOpChannel = adjHit->OpChannel();
            int myOpChannel = myHit->OpChannel();
            double adjPeakTime = adjHit->PeakTime();
            double myPeakTime = myHit->PeakTime();

            auto adjHitXYZ = opDetCenters[adjOpChannel];
            auto myHitXYZ = opDetCenters[myOpChannel];

            if ((adjHitXYZ - myHitXYZ).Mag() <= RadRange &&
                std::abs(adjPeakTime - myPeakTime) <= TimeRange)
            {
              // --- Check that this element isn't already in AddNow.
              if (std::find(AddNow.begin(), AddNow.end(), std::distance(MyVec.begin(), it)) == AddNow.end())
                AddNow.push_back(std::distance(MyVec.begin(), it));
            }
            else
            {
              ++it;
            }
          }
        }

        // --- Now loop through AddNow and remove from MyVec whilst adding to AdjHitVec
        std::sort(AddNow.begin(), AddNow.end());
        for (auto it = AddNow.rbegin(); it != AddNow.rend(); ++it)
        {
          if (HeavDebug)
          {
            std::cerr << "\tRemoving element from MyVec ===> "
                      << MyVec[*it]->OpChannel() << " & " << MyVec[*it]->PeakTime()
                      << std::endl;
          }

          AdjHitVec.push_back(MyVec[*it]);
          MyVec.erase(MyVec.begin() + *it);
        }

        LastSize = NewSize;
        NewSize = AdjHitVec.size();
        if (HeavDebug)
        {
          std::cerr << "\t---After that pass, AddNow was size " << AddNow.size() << " ==> LastSize is " << LastSize << ", and NewSize is " << NewSize
                    << "\nLets see what is in AdjHitVec...." << std::endl;
          for (size_t aL = 0; aL < AdjHitVec.size(); ++aL)
          {
            std::cout << "\tElement " << aL << " is ===> " << AdjHitVec[aL]->OpChannel() << " & " << AdjHitVec[aL]->PeakTime() << std::endl;
          }
        }
      } // while ( LastSize != NewSize )

      int NumAdjColHits = AdjHitVec.size();

      if (HeavDebug)
        std::cerr << "After that loop, I had " << NumAdjColHits << " adjacent collection plane hits." << std::endl;

      FilledHits += NumAdjColHits;

      if (AdjHitVec.size() > 0)
        Clusters.push_back(std::move(AdjHitVec));
    }

    if (HeavDebug)
    {
      std::vector<double> avgChannel;
      std::vector<double> avgTick;
      std::vector<double> summedADCInt;

      for (const auto &hits : Clusters)
      {
        double adcInt = 0;
        double channel = 0;
        double tick = 0;

        for (const auto &PDSHit : hits)
        {
          tick += PDSHit->PE() * PDSHit->PeakTime();
          channel += PDSHit->PE() * PDSHit->OpChannel();
          adcInt += PDSHit->PE();
        }
        tick /= adcInt;
        channel /= adcInt;
        summedADCInt.push_back(adcInt);
        avgTick.push_back(tick);
        avgChannel.push_back(channel);
      }

      for (int i = 0; i < int(avgTick.size() - 1); i++)
      {
        for (int j = i + 1; j < int(avgTick.size()); j++)
        {
          std::cout << avgChannel[i] << " " << avgChannel[j] << "  " << std::abs(avgChannel[i] - avgChannel[j]) << std::endl;
          std::cout << avgTick[i] << " " << avgTick[j] << "  " << std::abs(avgTick[i] - avgTick[j]) << std::endl;
          std::cout << summedADCInt[i] << " " << summedADCInt[j] << std::endl;
        }
      }
    }
    return;
  }

  // Function to calculate the Gaussian probability density function
  double AdjOpHitsUtils::GaussianPDF(double x, double mean, double sigma)
  {
    return exp(-0.5 * pow((x - mean) / sigma, 2)) / (sqrt(2 * M_PI) * sigma);
  }

  void AdjOpHitsUtils::CalcCentroid(std::vector<art::Ptr<recob::OpHit>> Hits, double &x, double &y, double &z)
  {
    const double sigma = fOpFlashAlgoRad; // Gaussian sigma (range in cm)

    // Initialize variables
    double maxLikelihood = 0.0;
    double bestY = 0.0;
    double bestZ = 0.0;

    // Loop over possible x positions
    double firstHitX = wireReadout.OpDetGeoFromOpChannel(Hits[0]->OpChannel()).GetCenter().X();
    double firstHitY = wireReadout.OpDetGeoFromOpChannel(Hits[0]->OpChannel()).GetCenter().Y();
    double firstHitZ = wireReadout.OpDetGeoFromOpChannel(Hits[0]->OpChannel()).GetCenter().Z();

    for (double yPos = firstHitY - fOpFlashAlgoRad; yPos <= firstHitY + fOpFlashAlgoRad; yPos += 5)
    {
      for (double zPos = firstHitZ - fOpFlashAlgoRad; zPos <= firstHitZ + fOpFlashAlgoRad; zPos += 5)
      {
        // Skipt the yPos and zPos that are outside the circle of radius sigma around the first hit
        if (pow(yPos - firstHitY, 2) + pow(zPos - firstHitZ, 2) > pow(sigma, 2))
          continue;
        double likelihood = 0.0;
        double sumY = 0.0;
        double sumZ = 0.0;
        int count = 0;

        // Loop over hits
        for (const auto &hit : Hits)
        {
          double hitYPos = wireReadout.OpDetGeoFromOpChannel(hit->OpChannel()).GetCenter().Y();
          double hitZPos = wireReadout.OpDetGeoFromOpChannel(hit->OpChannel()).GetCenter().Z();

          // Calculate the likelihood for the hit
          double hitLikelihood = GaussianPDF(hitYPos, yPos, sigma) * GaussianPDF(hitZPos, zPos, sigma);

          // Accumulate the likelihood and calculate the weighted sum of Y and Z coordinates
          likelihood += log(hitLikelihood);
          sumY += hitLikelihood * hitYPos;
          sumZ += hitLikelihood * hitZPos;
          count++;
        }

        // Check if the current likelihood is the maximum
        if (likelihood > maxLikelihood)
        {
          maxLikelihood = likelihood;
          bestY = sumY / count;
          bestZ = sumZ / count;
        }
      }
    }

    // Set the best 3D spacepoint
    x = firstHitX;
    y = bestY;
    z = bestZ;
  }
} // namespace solar

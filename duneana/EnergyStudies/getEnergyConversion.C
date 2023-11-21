/////////////////////////////////////////////////////////////////////////////////////////////////////
// Mike Wallbank (m.wallbank@sheffield.ac.uk), August 2015
//
// Macro to produce energy-charge conversion for the DUNE FD/35t.
//
// Usage:
//   root getEnergyConversion.C [ran in the same directory as EMEnergyCalib.root]
//
// Description of intended use:
//
//   Used to provide conversion factors between collected charge and deposited energy in MC showers.
//   Relationship is linear and must be determined separately for each plane.
//
//   -- runs over the output file (EMEnergyCalib.root) produced from the EMEnergyCalib_module art
//      analyser
//   -- produces an output root file which contains all of the plots used to find the conversion
//      from deposited charge to true energy
//   -- function is of the form (straight line)
//          charge = A + (B * energy) and A and B are determined from this script
//   -- linear plots showing relation are found in the output root file
//   -- the values of A and B for each plane are also printed upon completion of the script and can
//      be used in reconstruction; in DUNE the values are placed in the fhicl file:
//          larreco/larreco/RecoAlg/showeralgorithms.fcl
//
//   It is intended to be used on a sample of PG showering particles of 10 different energies.
//   e.g. 1000 electrons at each of the energies 0.5 GeV, 1.0 GeV, .. , 5.0 GeV.
//   EMEnergyCalib must be run on the samples first before passing (as a single root file) into this
//   macro.
//   List the energies in the vector at the start (kParticleEnergies).
//   The macro will plot collected charge vs true deposited energy and provide conversion factors.
/////////////////////////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
#include <utility>
#include <map>
#include <algorithm>
#include <vector>
#include "TTree.h"
#include "TMath.h"

const int kMaxHits = 10000;
const double kEnergy[10] = { 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 };
const double kdE = 0.5;

class EMEnergyConversion {
public:

  EMEnergyConversion(TTree* tree);
  ~EMEnergyConversion();
  void MakeFits();
  std::pair<double,double> MakeFit(int num, TString plane);
  void ProcessEvent();
  void Run();
  void SaveHists();
  void SetBranchAddresses();

private:

  TTree* fTree;

  double TrueEnergy;
  double DepositU;
  double DepositV;
  double DepositZ;
  double CorrectedChargeU;
  double CorrectedChargeV;
  double CorrectedChargeZ;
  double VertexDetectorDist;
  int    NHits;
  int    Hit_TPC      [kMaxHits];
  int    Hit_Plane    [kMaxHits];
  int    Hit_Wire     [kMaxHits];
  int    Hit_Channel  [kMaxHits];
  double Hit_PeakT    [kMaxHits];
  double Hit_Charge   [kMaxHits];
  int    Hit_ClusterID[kMaxHits];

  TH2D* ChargeDepositEnergyU;
  TH2D* ChargeDepositEnergyV;
  TH2D* ChargeDepositEnergyZ;
  TH2D* EnergyDepositUDistance;
  TH2D* EnergyDepositVDistance;
  TH2D* EnergyDepositZDistance;

  TFile *outFile;

};

EMEnergyConversion::EMEnergyConversion(TTree* tree) {
  fTree = tree;
  ChargeDepositEnergyU = new TH2D("ChargeDepositEnergyU","Charge v Depositied Energy on U;Charge (ADC);Deposited Energy (GeV);",100,0,15e9,50,0,5);
  ChargeDepositEnergyV = new TH2D("ChargeDepositEnergyV","Charge v Depositied Energy on V;Charge (ADC);Deposited Energy (GeV);",100,0,15e9,50,0,5);
  ChargeDepositEnergyZ = new TH2D("ChargeDepositEnergyZ","Charge v Depositied Energy on Z;Charge (ADC);Deposited Energy (GeV);",100,0,15e9,50,0,5);
  EnergyDepositUDistance = new TH2D("EnergyDepositUDistance","Deposited Energy on U vs Distance from Detector Edge;Distance (cm);Fraction of Energy Deposited;",100,0,220,50,0,1.1);
  EnergyDepositVDistance = new TH2D("EnergyDepositVDistance","Deposited Energy on V vs Distance from Detector Edge;Distance (cm);Fraction of Energy Deposited;",100,0,220,50,0,1.1);
  EnergyDepositZDistance = new TH2D("EnergyDepositZDistance","Deposited Energy on Z vs Distance from Detector Edge;Distance (cm);Fraction of Energy Deposited;",100,0,220,50,0,1.1);
  EnergyCompleteness = new TH1D("EnergyCompleteness","Fraction of Energy Reconstructed in Largest Cluster;Energy Completeness;",101,0,1.01);
  outFile = TFile::Open("EnergyFits.root","RECREATE");
}

EMEnergyConversion::~EMEnergyConversion() {
  outFile->Close();
  delete outFile;
}

void EMEnergyConversion::SetBranchAddresses() {
  fTree->SetBranchAddress("TrueEnergy",&TrueEnergy);
  fTree->SetBranchAddress("DepositU",&DepositU);
  fTree->SetBranchAddress("DepositV",&DepositV);
  fTree->SetBranchAddress("DepositZ",&DepositZ);
  fTree->SetBranchAddress("CorrectedChargeU",&CorrectedChargeU);
  fTree->SetBranchAddress("CorrectedChargeV",&CorrectedChargeV);
  fTree->SetBranchAddress("CorrectedChargeZ",&CorrectedChargeZ);
  fTree->SetBranchAddress("VertexDetectorDist",&VertexDetectorDist);
  fTree->SetBranchAddress("NHits",&NHits);
  fTree->SetBranchAddress("Hit_TPC",&Hit_TPC);
  fTree->SetBranchAddress("Hit_Plane",&Hit_Plane);
  fTree->SetBranchAddress("Hit_Wire",&Hit_Wire);
  fTree->SetBranchAddress("Hit_Channel",&Hit_Channel);
  fTree->SetBranchAddress("Hit_PeakT",&Hit_PeakT);
  fTree->SetBranchAddress("Hit_Charge",&Hit_Charge);
  fTree->SetBranchAddress("Hit_ClusterID",&Hit_ClusterID);
}

void EMEnergyConversion::MakeFits() {

  // Fit for U,V,Z plane
  std::pair<double,double> fitParamsU = this->MakeFit(0, TString("U"));
  std::pair<double,double> fitParamsV = this->MakeFit(1, TString("V"));
  std::pair<double,double> fitParamsZ = this->MakeFit(2, TString("Z"));

  std::cout << "U plane fit has intercept " << fitParamsU.first << " and gradient " << fitParamsU.second << std::endl;
  std::cout << "V plane fit has intercept " << fitParamsV.first << " and gradient " << fitParamsV.second << std::endl;
  std::cout << "Z plane fit has intercept " << fitParamsZ.first << " and gradient " << fitParamsZ.second << std::endl;

}

std::pair<double,double> EMEnergyConversion::MakeFit(int num, TString plane) {

  // Find the highest and lowest charge
  std::map<int,double> hCharge;
  std::map<int,double> lCharge;
  for (int i = 1; i <= 10; ++i) {
    hCharge[i] = 0;
    lCharge[i] = 1e10;
  }

  for (unsigned int event = 0; event < fTree->GetEntriesFast(); ++event) {
    fTree->GetEntry(event);
    double deposit;
    double planeCharge;
    switch (num) {
    case 0:
      deposit = DepositU;
      planeCharge = CorrectedChargeU;
      break;
    case 1:
      deposit = DepositV;
      planeCharge = CorrectedChargeV;
      break;
    case 2:
      deposit = DepositZ;
      planeCharge = CorrectedChargeZ;
      break;
    }
    if (deposit > kEnergy[0]-kdE && deposit <= kEnergy[0]) { if (planeCharge > hCharge[1]) hCharge[1] = planeCharge; if (planeCharge < lCharge[1]) lCharge[1] = planeCharge; }
    if (deposit > kEnergy[1]-kdE && deposit <= kEnergy[1]) { if (planeCharge > hCharge[2]) hCharge[2] = planeCharge; if (planeCharge < lCharge[2]) lCharge[2] = planeCharge; }
    if (deposit > kEnergy[2]-kdE && deposit <= kEnergy[2]) { if (planeCharge > hCharge[3]) hCharge[3] = planeCharge; if (planeCharge < lCharge[3]) lCharge[3] = planeCharge; }
    if (deposit > kEnergy[3]-kdE && deposit <= kEnergy[3]) { if (planeCharge > hCharge[4]) hCharge[4] = planeCharge; if (planeCharge < lCharge[4]) lCharge[4] = planeCharge; }
    if (deposit > kEnergy[4]-kdE && deposit <= kEnergy[4]) { if (planeCharge > hCharge[5]) hCharge[5] = planeCharge; if (planeCharge < lCharge[5]) lCharge[5] = planeCharge; }
    if (deposit > kEnergy[5]-kdE && deposit <= kEnergy[5]) { if (planeCharge > hCharge[6]) hCharge[6] = planeCharge; if (planeCharge < lCharge[6]) lCharge[6] = planeCharge; }
    if (deposit > kEnergy[6]-kdE && deposit <= kEnergy[6]) { if (planeCharge > hCharge[7]) hCharge[7] = planeCharge; if (planeCharge < lCharge[7]) lCharge[7] = planeCharge; }
    if (deposit > kEnergy[7]-kdE && deposit <= kEnergy[7]) { if (planeCharge > hCharge[8]) hCharge[8] = planeCharge; if (planeCharge < lCharge[8]) lCharge[8] = planeCharge; }
    if (deposit > kEnergy[8]-kdE && deposit <= kEnergy[8]) { if (planeCharge > hCharge[9]) hCharge[9] = planeCharge; if (planeCharge < lCharge[9]) lCharge[9] = planeCharge; }
    if (deposit > kEnergy[9]-kdE && deposit <= kEnergy[9]) { if (planeCharge > hCharge[10]) hCharge[10] = planeCharge; if (planeCharge < lCharge[10]) lCharge[10] = planeCharge; }
  }

  // Make a load of histograms!
  long dCharge = 100;
  TH1D* ChargeDist1 = new TH1D(TString("ChargeDist1Plane")+plane,";Total ADC;",200,lCharge[1]-dCharge,hCharge[1]+dCharge);
  TH1D* ChargeDist2 = new TH1D(TString("ChargeDist2Plane")+plane,";Total ADC;",200,lCharge[2]-dCharge,hCharge[2]+dCharge);
  TH1D* ChargeDist3 = new TH1D(TString("ChargeDist3Plane")+plane,";Total ADC;",200,lCharge[3]-dCharge,hCharge[3]+dCharge);
  TH1D* ChargeDist4 = new TH1D(TString("ChargeDist4Plane")+plane,";Total ADC;",200,lCharge[4]-dCharge,hCharge[4]+dCharge);
  TH1D* ChargeDist5 = new TH1D(TString("ChargeDist5Plane")+plane,";Total ADC;",200,lCharge[5]-dCharge,hCharge[5]+dCharge);
  TH1D* ChargeDist6 = new TH1D(TString("ChargeDist6Plane")+plane,";Total ADC;",200,lCharge[6]-dCharge,hCharge[6]+dCharge);
  TH1D* ChargeDist7 = new TH1D(TString("ChargeDist7Plane")+plane,";Total ADC;",200,lCharge[7]-dCharge,hCharge[7]+dCharge);
  TH1D* ChargeDist8 = new TH1D(TString("ChargeDist8Plane")+plane,";Total ADC;",200,lCharge[8]-dCharge,hCharge[8]+dCharge);
  TH1D* ChargeDist9 = new TH1D(TString("ChargeDist9Plane")+plane,";Total ADC;",200,lCharge[9]-dCharge,hCharge[9]+dCharge);
  TH1D* ChargeDist10 = new TH1D(TString("ChargeDist10Plane")+plane,";Total ADC;",200,lCharge[10]-dCharge,hCharge[10]+dCharge);
  TH1D* EnergyDist1 = new TH1D(TString("EnergyDist1Plane")+plane,";Energy (GeV);",100,kEnergy[0]-kdE,kEnergy[0]);
  TH1D* EnergyDist2 = new TH1D(TString("EnergyDist2Plane")+plane,";Energy (GeV);",100,kEnergy[1]-kdE,kEnergy[1]);
  TH1D* EnergyDist3 = new TH1D(TString("EnergyDist3Plane")+plane,";Energy (GeV);",100,kEnergy[2]-kdE,kEnergy[2]);
  TH1D* EnergyDist4 = new TH1D(TString("EnergyDist4Plane")+plane,";Energy (GeV);",100,kEnergy[3]-kdE,kEnergy[3]);
  TH1D* EnergyDist5 = new TH1D(TString("EnergyDist5Plane")+plane,";Energy (GeV);",100,kEnergy[4]-kdE,kEnergy[4]);
  TH1D* EnergyDist6 = new TH1D(TString("EnergyDist6Plane")+plane,";Energy (GeV);",100,kEnergy[5]-kdE,kEnergy[5]);
  TH1D* EnergyDist7 = new TH1D(TString("EnergyDist7Plane")+plane,";Energy (GeV);",100,kEnergy[6]-kdE,kEnergy[6]);
  TH1D* EnergyDist8 = new TH1D(TString("EnergyDist8Plane")+plane,";Energy (GeV);",100,kEnergy[7]-kdE,kEnergy[7]);
  TH1D* EnergyDist9 = new TH1D(TString("EnergyDist9Plane")+plane,";Energy (GeV);",100,kEnergy[8]-kdE,kEnergy[8]);
  TH1D* EnergyDist10 = new TH1D(TString("EnergyDist10Plane")+plane,";Energy (GeV);",100,kEnergy[9]-kdE,kEnergy[9]);

  for (unsigned int event = 0; event < fTree->GetEntriesFast(); ++event) {
    fTree->GetEntry(event);
    double deposit;
    double planeCharge;
    switch (num) {
    case 0:
      deposit = DepositU;
      planeCharge = CorrectedChargeU;
      break;
    case 1:
      deposit = DepositV;
      planeCharge = CorrectedChargeV;
      break;
    case 2:
      deposit = DepositZ;
      planeCharge = CorrectedChargeZ;
      break;
    }
    if (deposit > kEnergy[0]-kdE && deposit <= kEnergy[0]) { ChargeDist1->Fill(planeCharge); EnergyDist1->Fill(deposit); }
    if (deposit > kEnergy[1]-kdE && deposit <= kEnergy[1]) { ChargeDist2->Fill(planeCharge); EnergyDist2->Fill(deposit); }
    if (deposit > kEnergy[2]-kdE && deposit <= kEnergy[2]) { ChargeDist3->Fill(planeCharge); EnergyDist3->Fill(deposit); }
    if (deposit > kEnergy[3]-kdE && deposit <= kEnergy[3]) { ChargeDist4->Fill(planeCharge); EnergyDist4->Fill(deposit); }
    if (deposit > kEnergy[4]-kdE && deposit <= kEnergy[4]) { ChargeDist5->Fill(planeCharge); EnergyDist5->Fill(deposit); }
    if (deposit > kEnergy[5]-kdE && deposit <= kEnergy[5]) { ChargeDist6->Fill(planeCharge); EnergyDist6->Fill(deposit); }
    if (deposit > kEnergy[6]-kdE && deposit <= kEnergy[6]) { ChargeDist7->Fill(planeCharge); EnergyDist7->Fill(deposit); }
    if (deposit > kEnergy[7]-kdE && deposit <= kEnergy[7]) { ChargeDist8->Fill(planeCharge); EnergyDist8->Fill(deposit); }
    if (deposit > kEnergy[8]-kdE && deposit <= kEnergy[8]) { ChargeDist9->Fill(planeCharge); EnergyDist9->Fill(deposit); }
    if (deposit > kEnergy[9]-kdE && deposit <= kEnergy[9]) { ChargeDist10->Fill(planeCharge); EnergyDist10->Fill(deposit); }
  }

  TF1* fit;
  double charge[10], energy[10];
  outFile->cd();

  ChargeDist1->Fit("gaus"); fit = ChargeDist1->GetFunction("gaus"); charge[0] = fit->GetParameter(1); ChargeDist1->Write();
  ChargeDist2->Fit("gaus"); fit = ChargeDist2->GetFunction("gaus"); charge[1] = fit->GetParameter(1); ChargeDist2->Write();
  ChargeDist3->Fit("gaus"); fit = ChargeDist3->GetFunction("gaus"); charge[2] = fit->GetParameter(1); ChargeDist3->Write();
  ChargeDist4->Fit("gaus"); fit = ChargeDist4->GetFunction("gaus"); charge[3] = fit->GetParameter(1); ChargeDist4->Write();
  ChargeDist5->Fit("gaus"); fit = ChargeDist5->GetFunction("gaus"); charge[4] = fit->GetParameter(1); ChargeDist5->Write();
  ChargeDist6->Fit("gaus"); fit = ChargeDist6->GetFunction("gaus"); charge[5] = fit->GetParameter(1); ChargeDist6->Write();
  ChargeDist7->Fit("gaus"); fit = ChargeDist7->GetFunction("gaus"); charge[6] = fit->GetParameter(1); ChargeDist7->Write();
  ChargeDist8->Fit("gaus"); fit = ChargeDist8->GetFunction("gaus"); charge[7] = fit->GetParameter(1); ChargeDist8->Write();
  ChargeDist9->Fit("gaus"); fit = ChargeDist9->GetFunction("gaus"); charge[8] = fit->GetParameter(1); ChargeDist9->Write();
  ChargeDist10->Fit("gaus"); fit = ChargeDist10->GetFunction("gaus"); charge[9] = fit->GetParameter(1); ChargeDist10->Write();
  energy[0] = EnergyDist1->GetMean(); EnergyDist1->Write();
  energy[1] = EnergyDist2->GetMean(); EnergyDist2->Write();
  energy[2] = EnergyDist3->GetMean(); EnergyDist3->Write();
  energy[3] = EnergyDist4->GetMean(); EnergyDist4->Write();
  energy[4] = EnergyDist5->GetMean(); EnergyDist5->Write();
  energy[5] = EnergyDist6->GetMean(); EnergyDist6->Write();
  energy[6] = EnergyDist7->GetMean(); EnergyDist7->Write();
  energy[7] = EnergyDist8->GetMean(); EnergyDist8->Write();
  energy[8] = EnergyDist9->GetMean(); EnergyDist9->Write();
  energy[9] = EnergyDist10->GetMean(); EnergyDist10->Write();

  delete fit;

  TGraph *graph = new TGraph(10, charge, energy);
  graph->SetName(TString("FitPlane")+plane);
  graph->GetXaxis()->SetTitle("Charge");
  graph->GetYaxis()->SetTitle("Deposited Energy (GeV)");
  graph->SetTitle(plane);
  graph->SetMarkerStyle(8);
  graph->SetMarkerSize(1);
  graph->Fit("pol1");
  graph->Write();
  fit = graph->GetFunction("pol1");
  std::pair<double,double> fitParameters = std::make_pair(fit->GetParameter(0), fit->GetParameter(1));
  delete fit;
  delete graph;

  return fitParameters;

}

void EMEnergyConversion::Run() {

  this->SetBranchAddresses();

  for (unsigned int event = 0; event < fTree->GetEntriesFast(); ++event) {
    if (event % 1000 == 0) std::cout << "Processing event " << event << std::endl;
    fTree->GetEntry(event);
    this->ProcessEvent();
  }

  this->MakeFits();
  this->SaveHists();

}

void EMEnergyConversion::ProcessEvent() {

  std::map<int,double> clusterChargeU, clusterChargeV, clusterChargeZ;

  for (int hit = 0; hit < NHits; ++hit) {
    switch (Hit_Plane[hit]) {
    case 0:
      clusterChargeU[Hit_ClusterID[hit]] += (Hit_Charge[hit] * TMath::Exp((500 * Hit_PeakT[hit])/3e6));
      break;
    case 1:
      clusterChargeV[Hit_ClusterID[hit]] += (Hit_Charge[hit] * TMath::Exp((500 * Hit_PeakT[hit])/3e6));
      break;
    case 2:
      clusterChargeZ[Hit_ClusterID[hit]] += (Hit_Charge[hit] * TMath::Exp((500 * Hit_PeakT[hit])/3e6));
      break;
    }

  } // hit loop

  // Find the highest charge cluster
  double highChargeU = 0;
  for (std::map<int,double>::iterator chargeIt = clusterChargeU.begin(); chargeIt != clusterChargeU.end(); ++chargeIt)
    if (chargeIt->second > highChargeU) highChargeU = chargeIt->second;
  double highChargeV = 0;
  for (std::map<int,double>::iterator chargeIt = clusterChargeV.begin(); chargeIt != clusterChargeV.end(); ++chargeIt)
    if (chargeIt->second > highChargeV) highChargeV = chargeIt->second;
  double highChargeZ = 0;
  for (std::map<int,double>::iterator chargeIt = clusterChargeZ.begin(); chargeIt != clusterChargeZ.end(); ++chargeIt)
    if (chargeIt->second > highChargeZ) highChargeZ = chargeIt->second;

  ChargeDepositEnergyU->Fill(DepositU, CorrectedChargeU);
  ChargeDepositEnergyV->Fill(DepositV, CorrectedChargeV);
  ChargeDepositEnergyZ->Fill(DepositZ, CorrectedChargeZ);
  EnergyDepositUDistance->Fill(VertexDetectorDist, (double)DepositU/(double)TrueEnergy);
  EnergyDepositVDistance->Fill(VertexDetectorDist, (double)DepositV/(double)TrueEnergy);
  EnergyDepositZDistance->Fill(VertexDetectorDist, (double)DepositZ/(double)TrueEnergy);

}

void EMEnergyConversion::SaveHists() {

  outFile->cd();

  TCanvas* cChargeDepositEnergyU = new TCanvas("cChargeDepositEnergyU","",800,600);
  ChargeDepositEnergyU->Draw("colz");
  cChargeDepositEnergyU->Write("ChargeDepositEnergyU");
  TCanvas* cChargeDepositEnergyV = new TCanvas("cChargeDepositEnergyV","",800,600);
  ChargeDepositEnergyV->Draw("colz");
  cChargeDepositEnergyV->Write("ChargeDepositEnergyV");
  TCanvas* cChargeDepositEnergyZ = new TCanvas("cChargeDepositEnergyZ","",800,600);
  ChargeDepositEnergyZ->Draw("colz");
  cChargeDepositEnergyZ->Write("ChargeDepositEnergyZ");

  EnergyDepositUDistance->Write("EnergyDepositUDistance");
  EnergyDepositVDistance->Write("EnergyDepositVDistance");
  EnergyDepositZDistance->Write("EnergyDepositZDistance");

  delete cChargeDepositEnergyU;
  delete cChargeDepositEnergyV;
  delete cChargeDepositEnergyZ;

}

void getEnergyConversion() {

  //gStyle->SetOptStat(0);
  //TFile* inFile = TFile::Open("/pnfs/lbne/persistent/users/tjyang/v04_20_00/mergeana/prod_gamma_0.1-1.0GeV_dune35t/2884577_0/EMEnergyCalib.root");
  TFile* inFile = TFile::Open("EMEnergyCalib.root");
  TTree* tree = (TTree*)inFile->Get("energyCalib/EMEnergyCalib");

  EMEnergyConversion emenergycalib = EMEnergyConversion(tree);
  emenergycalib.Run();

  inFile->Close();
  inFile->Delete();

}

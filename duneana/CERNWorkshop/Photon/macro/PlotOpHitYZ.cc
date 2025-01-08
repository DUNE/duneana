void PlotOpHitYZ(TString file_name, int event)
{
  TFile * file = new TFile(file_name,"READ");
  TTree * tree = (TTree *)file->Get("hitdumper/hitdumpertree");

  TCut c1 = Form("event == %d", event);

  //2D Histogram with DUNE PDS dimensions
  TH2F * hPhotocatode = new TH2F("hPhotocatode", ";Z [cm];Y [cm]; #PE", 80, 0, 1390, 70, -640, 640);
  hPhotocatode->SetStats(0);
  tree->Draw("ophit_opdet_y:ophit_opdet_z >> hPhotocatode", (c1)*"ophit_pe", "COLZ");

  gPad->SetRightMargin(0.125);
}

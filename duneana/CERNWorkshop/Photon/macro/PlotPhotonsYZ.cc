void PlotPhotonsYZ(TString file_name, int event)
{
  TFile * file = new TFile(file_name,"READ");
  TTree * tree = (TTree *)file->Get("opanalyzer/PhotonsPerOpDet");

  // Cuts to select only particular event
  TCut c2 = Form("EventID == %d", event);
  
  //2D Histogram with DUNE PDS dimensions
  TH2F * hPhotocatode = new TH2F("hPhotocathode", ";Z [cm];Y [cm]; #PE",80, 0, 1390, 70, -640, 640);
  hPhotocatode->SetStats(0);
  tree->Draw("OpDetY:OpDetZ >> hPhotocathode", (c2)*"CountAll", "COLZ");

  gPad->SetRightMargin(0.125);
}

using namespace RooFit;
void PlotFit(int sigmu = 0, double lumi = 1) { 
  double xlow(0.07), xhigh(0.95);
  double ylow(0.05), yhigh(0.90);
  double hratio(0.25);
  TString wspacePath = TString::Format("workspaces/weightedMC_mu%d_lumi%.0f.root",sigmu,lumi);
  TFile *infile = TFile::Open(wspacePath);
  RooWorkspace *w = (RooWorkspace*)infile->Get("w");
  RooRealVar *m = w->var("llphoton_m");
  RooAbsPdf *pdf = w->pdf("pdf");
  RooAbsData *data = w->data("data");
  TCanvas *can = new TCanvas("can","canvas",800,800);
  can->Divide(1,2);
  TVirtualPad *top = can->cd(1);
  TVirtualPad *bot = can->cd(2);
  // Margin: L, R, B, T
  // Pad xl, yl, xu, yu
  top->SetMargin(xlow,1-xhigh,0.0 ,1-yhigh);
  bot->SetMargin(xlow,1-xhigh,0.30,0.0);
  top->SetPad(0.0,hratio,1.0,1.0   );
  bot->SetPad(0.0,   0.0,1.0,hratio);
  bot->SetGridy(1);
  top->cd();
  RooPlot *plot = m->frame();
  int nbins(80);
  data->plotOn(plot,Binning(nbins));
//   pdf->SetLineWidth(1);
  pdf->plotOn(plot,RooFit::Components("background"),RooFit::LineColor(kGreen),RooFit::LineWidth(2));
  pdf->plotOn(plot,RooFit::LineColor(kRed),RooFit::LineWidth(2));
  pdf->paramOn(plot,data,"",2,"NELU",0.6,xhigh,yhigh);
  plot->GetYaxis()->SetTitleOffset(0.9);
  plot->GetXaxis()->SetLimits(100,180);
  plot->Draw("same");
  bot->cd();
  RooHist *pull = plot->pullHist();
  pull->SetMarkerStyle(21);
  pull->GetXaxis()->SetLimits(100,180);
  pull->SetTitle(";m_{ll#gamma} [GeV]; Pull");
  pull->GetXaxis()->SetTitleSize(0.1);
  pull->GetYaxis()->SetTitleSize(0.1);
  pull->GetXaxis()->SetLabelSize(0.08);
  pull->GetYaxis()->SetRangeUser(-3.2,3.2);
  pull->GetYaxis()->SetLabelSize(0.08);
  pull->GetYaxis()->CenterTitle();
  pull->GetYaxis()->SetTickLength(0);
  pull->GetYaxis()->SetTitleOffset(0.2);
  pull->Draw("AP same");
  TString outfile = TString::Format("plots/sb-weight-fit_mu%d_%0.fifb.pdf",sigmu,lumi);
  can->SaveAs(outfile);
  TCanvas *can2 = new TCanvas("can2","pull canvas",800,800);
  can2->cd();
  TH1D* hpull = new TH1D("hpull","m_{ll#gamma} pull distribution",8,-4,4);
  Double_t ax[nbins], ay[nbins];
  for(int i = 0; i < nbins; i++) {
    pull->GetPoint(i,ax[i],ay[i]);
    hpull->Fill(ay[i]);
  }
  hpull->SetLineColor(kBlack);
  hpull->SetFillColor(kAzure);
  hpull->Draw("B");
  outfile = TString::Format("plots/sb-weight-pulls_mu%d_%0.fifb.pdf",sigmu,lumi);
  can2->SaveAs(outfile);
  return;
}

void PlotFits() {
//   for(int i = 0; i < 11; i++)
  PlotFit(0, 42);
//   PlotFit(0 ,36);
//   PlotFit(1 ,36);
//   PlotFit(2 ,36);
//   PlotFit(3 ,36);
//   PlotFit(4 ,36);
//   PlotFit(5 ,36);
//   PlotFit(6 ,36);
//   PlotFit(7 ,36);
//   PlotFit(8 ,36);
//   PlotFit(9 ,36);
//   PlotFit(10,36);
//   PlotFit(50);
//   PlotFit(100);
//   PlotFit(150);
//   PlotFit(200);
//   PlotFit(250);
//   PlotFit(300);
  return;
}

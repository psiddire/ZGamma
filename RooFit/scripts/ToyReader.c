double GetFitYield(int mu = 0, double lumi = 1, int ntoys = 100, bool useMC = false) {
  TString inpath = TString::Format("toys/PDF_mu%d_%.0fifb_ntoys-%d.root",mu,lumi,ntoys);
  if(useMC) inpath.ReplaceAll("PDF","MC");
  TFile *infile = TFile::Open(inpath);
  TH1D *hsigmu = (TH1D*)infile->Get("hsigmu");
  return hsigmu->GetMean();
}
double GetFitError(int mu = 0, double lumi = 1, int ntoys = 100, bool useMC = false) {
  TString inpath = TString::Format("toys/PDF_mu%d_%.0fifb_ntoys-%d.root",mu,lumi,ntoys);
  if(useMC) inpath.ReplaceAll("PDF","MC");
  TFile *infile = TFile::Open(inpath);
  TH1D *hsigmu = (TH1D*)infile->Get("hsigmu");
//   return hsigmu->GetMeanError();
  return hsigmu->GetRMS();
}
void ToyReader() {
  double lumi(42);
  int ntoys(100);
  bool useMC(true);
  vector<int> mus = {0,1,2,3,4,5,6,7,8,9,10};
  double ax[mus.size()], ay[mus.size()],aey[mus.size()];
  // Fill arrays with yields and errors
  for(size_t imu = 0; imu < mus.size(); imu++) {
    ax[imu]  = mus.at(imu);
    ay[imu]  = GetFitYield(mus.at(imu), lumi, ntoys, useMC);
    aey[imu] = GetFitError(mus.at(imu), lumi, ntoys, useMC);
  }
  TGraphErrors *ivf = new TGraphErrors((int)mus.size(),ax,ay,NULL,aey);
  TString title = TString::Format("Injected Signal vs Fitted Signal (%.0f fb^{-1}, %d toys each); #mu_{Injected}; #hat{#mu}_{Fitted}",lumi,ntoys);
  ivf->SetTitle(title);
  TCanvas *can = new TCanvas("can","canvas",800,800);
  can->cd();
  // Fit points with line
  TF1 *line = new TF1("line","[0]*x+[1]",-50,300);
  line->SetParameters(1,0);
  ivf->Fit(line);
  ivf->SetMarkerStyle(20);
//   double lower = mus.at(mus.size()-1)*-0.1;
  double lower = ay[0] - 2;
//   double upper = mus.at(mus.size()-1)*1.3;
  double upper = ay[10]*1.2;
//   ivf->GetXaxis()->SetLimits(lower,upper);
  ivf->GetXaxis()->SetLimits(-0.5,10.5);
  ivf->GetYaxis()->SetRangeUser(lower,upper);
//   if(useMC) ivf->GetYaxis()->SetRangeUser(lower-4,upper-2);
  ivf->GetYaxis()->SetTitleOffset(1.3);
  ivf->GetYaxis()->SetLabelSize(0.03);
  ivf->GetXaxis()->SetLabelSize(0.03);
  ivf->Draw("AP");
  // Add legend with linear fit parameters 
  TLegend *leg = new TLegend(0.1,0.9,0.4,0.76);
  leg->SetTextSize(.032);
  TString fitresult1 = TString::Format("m = %.2f #pm %.2f",line->GetParameter(0),line->GetParError(0));
  TString fitresult2 = TString::Format("b = %.2f #pm %.2f",line->GetParameter(1),line->GetParError(1));
  leg->AddEntry(line,"y = m*x + b");
  leg->AddEntry((TObject*)0,fitresult1,"");
  leg->AddEntry((TObject*)0,fitresult2,"");
  leg->Draw();
  TString outpath = TString::Format("plots/PDF-ToyIvF_%.0fifb_%dtoys.pdf",lumi,ntoys);
  if(useMC) outpath.ReplaceAll("/PDF","/MC");
  can->SaveAs(outpath);
  return; 
}

double GetFitYield(int mu = 0, double lumi = 1) {
  TString wspacePath = TString::Format("workspaces/weightedMC_mu%d_lumi%.0f.root",mu,lumi);
  TFile *infile = TFile::Open(wspacePath);
  RooWorkspace *w = (RooWorkspace*)infile->Get("w");
  RooRealVar *n = w->var("norm_s");
  double yield = n->getVal();
  infile->Close();
  return yield;
}
double GetFitError(int mu = 0, double lumi = 1) {
  TString wspacePath = TString::Format("workspaces/weightedMC_mu%d_lumi%.0f.root",mu,lumi);
  TFile *infile = TFile::Open(wspacePath);
  RooWorkspace *w = (RooWorkspace*)infile->Get("w");
  RooRealVar *n = w->var("norm_s");
  double error = n->getError();
  infile->Close();
  return error;
}

// Function with retrieves fit signal yield (& uncertainties) from
// workspaces, then plots them vs injected signal. Assumes a workspace
// files exists in workspace/ for each signal strength considered.
void InjectedvsFit() {
  double lumi(42);
  vector<int> mus = {0,1,2,3,4,5,6,7,8,9,10};
  double ax[mus.size()], ay[mus.size()],aey[mus.size()];
  // Fill arrays with yields and errors
  for(size_t imu = 0; imu < mus.size(); imu++) {
//     ax[imu]  = mus.at(imu)*1.7624*lumi; // 1.7624 events expected/ifb in 2016
    ax[imu]  = mus.at(imu)*2.195*lumi; // 1.7624 events expected/ifb in 2017
    ay[imu]  = GetFitYield(mus.at(imu), lumi);
    aey[imu] = GetFitError(mus.at(imu), lumi);
  }
  TGraphErrors *ivf = new TGraphErrors((int)mus.size(),ax,ay,NULL,aey);
  ivf->SetTitle("Injected Signal vs Fitted Signal; N_{Injected}; N_{Fitted}");
  TCanvas *can = new TCanvas("can","canvas",800,800);
  can->cd();
  // Fit points with line
  TF1 *line = new TF1("line","[0]*x+[1]",-50,300);
  line->SetParameters(1,0);
  ivf->Fit(line);
  ivf->SetMarkerStyle(20);
  double lower = mus.at(mus.size()-1)*-0.3*1.762*lumi;
  double upper = mus.at(mus.size()-1)*1.3*1.762*lumi;
  ivf->GetXaxis()->SetLimits(lower,upper);
  ivf->GetYaxis()->SetRangeUser(lower,upper);
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
  can->SaveAs("plots/IVF.pdf");
  return; 
}


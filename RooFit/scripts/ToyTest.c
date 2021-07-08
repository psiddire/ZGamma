#include <iostream>
#include <iomanip>
#include <getopt.h>
#include "ModGaus.h"
#include "GausCB.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"

void GenerateToy(TH1D *parent, TH1D* toy);
void GenerateToy(TF1 *parent, TH1D* mllg, int nevents);
void PlotToy(RooDataHist data, RooAddPdf pdf, RooRealVar m, int sigmu, double lumi, int i, bool useMC = false, bool failed = false);
double GetSignalYield(int year = 2016);
double GetBackgroundYield(int year = 2016);

using namespace RooFit;

namespace {
  int year(2017);
  // 2016
//   double parent_vl(2.2), parent_vr(0.7), parent_s0(11),parent_sl(8.2), parent_sh(50), parent_m0(112.6);
  // 2017
  double parent_vl(2.6), parent_vr(1.7), parent_s0(9), parent_sl(9.0), parent_sh(56), parent_m0(113.1);
  double backYield = GetBackgroundYield(year);
  double sigYield = GetSignalYield(year);
}

double sb_pdf(double *x, double *par) {
  double norm_b(backYield);
  Double_t width, power;
  power = parent_vl + parent_vr*(x[0] - 100)/(180-100);
  if(x[0] <= parent_m0) 
    width = parent_sl + (parent_s0 - parent_sl)*(x[0] - 100)/(parent_m0-100);
  else
    width = parent_s0 + (parent_sh - parent_s0)*(x[0] - parent_m0)/(180 - parent_m0);
  double bcomp = exp(-1*pow(fabs((x[0]-parent_m0)/width), power));
  double mean(124.69), gaus_std(11.8), cb_std(1.66), cb_alpha(0.74), cb_n(10.3), cb_frac(0.979);
  double norm_s = sigYield*par[0]*par[1];
  double g = ROOT::Math::gaussian_pdf(x[0],gaus_std,mean);
  double cb = ROOT::Math::crystalball_pdf(x[0],cb_alpha, cb_n, cb_std, mean);
  double scomp = (1-cb_frac)*g+cb_frac*cb;
  return norm_b/33.3362*bcomp+norm_s*scomp;
}


void ToyTest(int injectedSignal, double lumi, int ntoys, bool useMC = false) {
// void ToyTest() {
  // Setup parent distribution for sampling
  TChain chain("tree");
  TString infile = TString::Format("MC/combined-MC-mu%d-%.0fifb.root",injectedSignal,lumi);
  chain.Add(infile);
  TH1D *hmllg = new TH1D("hmllg","ll#gamma mass",80,100,180);
  chain.Draw("llphoton_m[0]>>hmllg","weight");
  // Declare background fit pdf and variables
  RooRealVar m("llphoton_m","ll#gamma mass [GeV]"    ,125,100,180);                   
  RooRealVar m0("m_{0}","mass peak value"            ,parent_m0 ,105,125); m0.setError(0.1); 
  RooRealVar vl("#nu_{L}","low-end exponential"      ,parent_vl ,0.1,6 );  vl.setError(0.1); 
  RooRealVar vr("#nu_{r}","high-end exponential"     ,parent_vr ,  0,5 );  vr.setError(0.1); 
  RooRealVar s0("#sigma_{0}","peak width value"      ,parent_s0 ,  0,30);  s0.setError(0.1); 
  RooRealVar sl("#sigma_{L}","low-end width value"   ,parent_sl ,  0,15);  sl.setError(0.1); 
  RooRealVar sh("#sigma_{H}","high-end width value"  ,parent_sh , 20,80);  sh.setError(0.1); 
  ModGaus background("background","Modified Gaussian",m,m0,vl,vr,s0,sl,sh);
  // Declare signal fit pdf and variables (then set as constant)
  // parameters[1] = 124.69, 3.3, 1.40, 0.85, 2.9, 0.871
  RooRealVar mu      ("#mu"          ,"mean"          ,124.66 );// 124.66
  RooRealVar gaus_std("#sigma_{gaus}","Gaussian stdev", 11.8  );// 11.8
  RooRealVar cb_std  ("#sigma_{cb}"  ,"CB stdev"      ,  1.66 );// 1.66
  RooRealVar cb_alpha("#alpha_{cb}"  ,"CB alpha"      ,  0.74 );// 0.74
  RooRealVar cb_n    ("n_{cb}"       ,"CB n"          , 10.3  );// 10.3
  RooRealVar cb_frac ("f_{cb}"       ,"CB fraction"   ,  0.979);// 0.979
  mu.setConstant(); gaus_std.setConstant(); cb_std.setConstant();
  cb_alpha.setConstant(); cb_n.setConstant(); cb_frac.setConstant();
  GausCB signal("signal","Gaussian+Crystal Ball",m,mu,gaus_std,cb_std,cb_alpha,cb_n,cb_frac);
  // Combine signal and background pdfs
//   int siglimit(lumi*(200+injectedSignal*5));
  int range(400);
  if(lumi == 36) range = 3600;
  if(lumi == 42) range = 4200;
  RooRealVar norm_s("norm_s","N_{s}",0,-100*lumi,100*lumi); norm_s.setError(1);
  RooRealVar norm_b("norm_b","N_{b}",(int)backYield,(int)lumi*1000,(int)6000*lumi); norm_b.setError(1);
  const RooArgList components(signal,background);
  const RooArgList coeffs(norm_s,norm_b);
  RooAddPdf pdf("pdf","f_{s+b}",components,coeffs);
  // Minimize printing
  // Declare histograms for storing per-toy fit values
  int murange(100); if(lumi == 36) murange = 20;
  TH1D *hsigmu = new TH1D("hsigmu","Fitted signal strength",50,injectedSignal - murange,injectedSignal + murange);
  TH1D *hpull = new TH1D("hpull","Fit pull",20,-5,5);
  TH1D *fstatus = new TH1D("fstatus","Fit status",10,-0.5,9.5);
  TH1D *fcovQual = new TH1D("fcovQual","Fit covariance quality",5,-0.5,4.5);
  bool isSumW2(0);
  // Declare function for generating toys
  TF1 *fsb = new TF1("fsb",sb_pdf,100,180,2);
  fsb->SetParameters(injectedSignal, lumi);
  gRandom->SetSeed(0);
  TH1D *mllg = new TH1D("mllg","mllg",80,100,180);
  int nfails(0);
  for(int i = 0; i < ntoys; i++) {
    // Reset parameter values and errors each time
    // (Found for large ntoys variables would drift into
    //  bad values and subsequent fits would fail)
    m0.setVal(parent_m0); vl.setVal(parent_vl); vr.setVal(parent_vr);
    s0.setVal(parent_s0);  sl.setVal(parent_sl); sh.setVal(parent_sh);
    m0.setError(0.1); vl.setError(0.1); vr.setError(0.1);
    s0.setError(0.1); sl.setError(0.1); sh.setError(0.1);
    norm_b.setVal(backYield*lumi); norm_s.setVal(0);
    // Fill mllg with toy events
    int nev = (backYield+sigYield*injectedSignal)*lumi;
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    if(useMC) GenerateToy(hmllg,mllg);  // Fill from MC
    else      GenerateToy(fsb,mllg,nev);// Fill from best-fit function
    RooDataHist data("data","data",RooArgSet(m),mllg);
    RooFitResult *fit; // Use to check status + covariance quality in saved RooFit
    fit = pdf.fitTo(data,Extended(),Save(), Verbose(0),PrintLevel(-1), Offset(1), Strategy(2));
//     fit = pdf.fitTo(data,Extended(),Save(), Offset(1), Strategy(2));
    if(i%10 == 0) {
//       cout << "Toy #"<< i << " (mu = " << mu << ", status = " << fit->status() << ", covQual = " << fit->covQual() << ")" << endl;
      PlotToy(data, pdf, m, injectedSignal, lumi, i, useMC);
    }
    double Nexpected = sigYield*lumi*injectedSignal;
    double mu = norm_s.getVal()/(sigYield*lumi);
    if(fit->status() == 0) {
      fstatus->Fill(fit->status());
      fcovQual->Fill(fit->covQual());
      hsigmu->Fill(mu);
      hpull->Fill((norm_s.getVal() - Nexpected)/norm_s.getError());
      cout << "FIT SUCCEEDED. RETRYING." << endl;
    }
    else {
      i = i-1;
      cout << "FIT FAILED. RETRYING." << endl;
      PlotToy(data, pdf, m, injectedSignal, lumi, i, useMC, true);
      nfails++;
    }
    mllg->Reset();
  }
  cout << "Fit failure \%" << 100*double(nfails)/ntoys << endl;
  TString outname;
  if(useMC) outname = TString::Format("toys/MC_mu%d_%.0fifb_ntoys-%d.root",injectedSignal,lumi,ntoys);
  else      outname = TString::Format("toys/PDF_mu%d_%.0fifb_ntoys-%d.root",injectedSignal,lumi,ntoys);
  TFile outfile(outname,"RECREATE");
  fstatus->Write();
  fcovQual->Write();
  hsigmu->Write();
  hpull->Write();
  outfile.Close();
  return;
}

void GenerateToy(TH1D *parent, TH1D* toy) {
  for(int i = 0; i < parent->Integral(); i++) 
    toy->Fill(parent->GetRandom());
  return;
}

void GenerateToy(TF1 *parent, TH1D *mllg, int nevents) {
  for(int i = 0; i < nevents; i++)
    mllg->Fill(parent->GetRandom());
  return;
}

void PlotToy(RooDataHist data, RooAddPdf pdf, RooRealVar m, int sigmu, double lumi, int i, bool useMC = false, bool failed = false) {
  double xlow(0.07), xhigh(0.95);
  double ylow(0.05), yhigh(0.90);
  double hratio(0.25);
  gROOT->SetBatch(kTRUE);
  TCanvas *can = new TCanvas("can","canvas",800,800);
  can->Divide(1,2);
  TVirtualPad *top = can->cd(1);
  TVirtualPad *bot = can->cd(2);
  top->SetMargin(xlow,1-xhigh,0.0 ,1-yhigh);
  bot->SetMargin(xlow,1-xhigh,0.30,0.0);
  top->SetPad(0.0,hratio,1.0,1.0   );
  bot->SetPad(0.0,   0.0,1.0,hratio);
  bot->SetGridy(1);
  top->cd();
  RooPlot *plot = m.frame();
  int nbins(80);
  data.plotOn(plot,Binning(nbins));
  pdf.plotOn(plot,RooFit::Components("background"),RooFit::LineColor(kGreen),RooFit::LineWidth(2));
  pdf.plotOn(plot,RooFit::LineColor(kRed),RooFit::LineWidth(2));
  pdf.paramOn(plot,&data,"",2,"NELU",0.6,xhigh,yhigh);
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
  TString outfile = TString::Format("plots/PDF-Toy-SB-Fit_mu%d_%0.fifb/",sigmu,lumi);
  if(useMC) outfile.ReplaceAll("/PDF","/MC");
  if(failed) outfile += "failed_";
  if(i == 0) gSystem->MakeDirectory(outfile);
  outfile += TString::Format("%d.pdf",i);
  can->SaveAs(outfile);
  delete can;
  return;
}

double GetSignalYield(int year) {
  TChain chain("tree");
  TString inpath = "/Users/adorsett/Desktop/CERNbox/pico/NanoAODv5/zgamma_channelIslandsv3/"
                   + to_string(year) + "/HToZG/merged_zgmc_zgsr/*.root";
  chain.Add(inpath);
  TH1D *mllg = new TH1D("mllg","signal mllg",80,100,180);
  chain.Draw("llphoton_m[0]>>mllg","weight");
  double yield = mllg->Integral();
  mllg->Delete();
  return yield;
}
double GetBackgroundYield(int year) {
  TChain chain("tree");
  TString inpath = "/Users/adorsett/Desktop/CERNbox/pico/NanoAODv5/zgamma_channelIslandsv3/"
                   + to_string(year) + "/mc/merged_zgmc_zgsr/*.root";
  chain.Add(inpath);
  TH1D *mllg = new TH1D("mllg","signal mllg",80,100,180);
  chain.Draw("llphoton_m[0]>>mllg","weight");
  double yield = mllg->Integral();
  mllg->Delete();
  return yield;
}
 

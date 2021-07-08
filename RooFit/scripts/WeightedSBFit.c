#include "ModGaus.h"
#include "GausCB.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "CombinePicos.c"
#include "FitPlotter.c"
using namespace RooFit;
void WeightedFit(int year = 2016, int sigmu = 0, double lumi = 1) {
  double defaults[6] = {112.6, 2.2, 0.7, 11, 8.2, 50};
  if(year == 2017) {
    defaults[0] = 113.1; defaults[1] = 2.2; defaults[2] = 1.7;
    defaults[3] =     9; defaults[4] = 9.0; defaults[5] =  56;
  }
  // Declare background fit pdf and variables
  RooRealVar m("llphoton_m","ll#gamma mass [GeV]"     ,125  ,100,180);
  RooRealVar m0("m_{0}","mass peak value"             ,defaults[0],100,180); m0.setError(0.1);
  RooRealVar vl("#nu_{L}","low-end exponential"       ,defaults[1],0,20);    vl.setError(0.1);
  RooRealVar vh("#nu_{r}","high-end exponential"      ,defaults[2],0,10);    vh.setError(0.1);
  RooRealVar s0("#sigma_{0}","peak width value"       ,defaults[3],0,50);    s0.setError(0.1);
  RooRealVar sl("#sigma_{L}","low-end width value"    ,defaults[4],0,20);    sl.setError(0.1);
  RooRealVar sh("#sigma_{H}","high-end width value"   ,defaults[5],20,80);   sh.setError(0.1);
  ModGaus background("background","Modified Gaussian",m,m0,vl,vh,s0,sl,sh);
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
  RooRealVar norm_s("norm_s","N_{s}",0, -100*lumi, lumi*100); norm_s.setError(1);
  RooRealVar norm_b("norm_b","N_{b}",(int)lumi*4000,1000,250000);
  const RooArgList components(signal,background);
  const RooArgList coeffs(norm_s,norm_b);
  RooAddPdf pdf("pdf","f_{s+b}",components,coeffs);
  // Generate & load data
  MakeTree(sigmu,lumi);
  TString inpath = TString::Format("MC/%d/combined-MC-mu%d-%.0fifb.root",year,sigmu,lumi);
  cout << "INPATH: " << inpath << endl;
  TFile *infile = TFile::Open(inpath);
  TTree *tree = (TTree*)infile->Get("tree");
  TH1D *mllg = new TH1D("mllg","mllg",80,100,180);
  tree->Draw("llphoton_m>>mllg","weight");
  RooRealVar weight("weight","weight",1,-10,10);
  RooDataSet data("data","data",tree,RooArgSet(m,weight),"1","weight");
  // Warning, progress, info
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  // Fit data with S+B model
//   pdf.fitTo(data,Extended());
  pdf.fitTo(data, Save(kTRUE), Offset(1),// Minimizer("Minuit2", "migrad"), 
//           Strategy(2), Optimize(1), 
          NumCPU(4), 
          Verbose(0),PrintLevel(-1),
          SumW2Error(1)
          );
  background.fitTo(data,Extended());
  FitPlotter(data, background, m, TString::Format("%dMC",year));
//           Minos(kTRUE));
  // Save resulting model to workspace
  RooWorkspace w("w");
  w.import(pdf);
  w.import(data);
  TString wspacePath = TString::Format("workspaces/weightedMC_mu%d_lumi%.0f.root",sigmu,lumi);
  w.writeToFile(wspacePath);
  return;
}

void WeightedSBFit() {
  WeightedFit(2017,0,42);
  WeightedFit(2016,0,36);
//   for(int i = 0; i < 11; i++)
//     WeightedFit(i);
//   WeightedFit(0,42);
//   WeightedFit(1,36);
//   WeightedFit(2,36);
//   WeightedFit(3,36);
//   WeightedFit(4,36);
//   WeightedFit(5,36);
//   WeightedFit(6,36);
//   WeightedFit(7,36);
//   WeightedFit(8,36);
//   WeightedFit(9,36);
//   WeightedFit(10,36);
//   WeightedFit(50);
//   WeightedFit(100);
//   WeightedFit(150);
//   WeightedFit(200);
//   WeightedFit(250);
//   WeightedFit(300);
  return;
}

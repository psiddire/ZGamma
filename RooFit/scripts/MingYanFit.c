#include "src/ModGaus.cxx"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "scripts/FitPlotter.c"
using namespace RooFit;
void MingYanFit() {
  // Declare background fit pdf and variables
  double defaults[6] = {110,2,2,4,5,50};
//   double defaults[6] = {113.07,1.97,1.67,12.44,7.72,49.2};
//   double defaults[6] = {113.11,1.87,0.71,12.88,7.79,49.1};
  RooRealVar m("llphoton_m","ll#gamma mass [GeV]",100,180);
  RooRealVar m0("m_{0}","mass peak value"             ,defaults[0],105,145);  //m0.setError(0.1);
  RooRealVar vl("#nu_{L}","low-end exponential"       ,defaults[1],  0.1,10); //vl.setError(0.1);
  RooRealVar vr("#Delta#nu","change in exponential"   ,defaults[2],  0.1,10); //vr.setError(0.1);
  RooRealVar s0("#sigma_{0}","peak width value"       ,defaults[3],  1, 20);  //s0.setError(0.1);
  RooRealVar sl("#sigma_{L}","low-end width value"    ,defaults[4],  1, 20);  //sl.setError(0.1);
  RooRealVar sh("#sigma_{H}","high-end width value"   ,defaults[5], 20, 80);  //sh.setError(0.1);
  ModGaus background("background","Modified Gaussian",m,m0,vl,vr,s0,sl,sh);
  // Load data
//   TFile *infile = TFile::Open("/Users/adorsett/Desktop/CERNbox/Zgamma/MingYan-MC-mllg.root");
  TString tag("untag3");
  TString inpath = "/Users/adorsett/Desktop/CERNbox/Zgamma/draw_pico/RooFit/MingYanFiles/";
  inpath += "bkghist_ele_mu_"+tag+".root";
  TFile *infile = TFile::Open(inpath);
  TH1D *hmc = (TH1D*)infile->Get("hbkg");
//   hmc->Rebin();
//   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooDataHist mc("mc","mc",RooArgSet(m),hmc);
  RooFitResult *fit;
//   TVirtualFitter::SetMaxIterations(100000);
  fit = background.fitTo(mc,
                         //Verbose(0),PrintLevel(-1),
                         Save(kTRUE), 
                         SumW2Error(kFALSE),
//                          Offset(2),
                         Minimizer("Minuit2","migrad"),
                         NumCPU(4)); 
  TString label = "MingYan_"+tag;
  FitPlotter(&mc, &background, m, 80,label);
  cout << endl << endl;
  cout << "Background only 80 bin fit status = " << fit->status() << endl;
  cout << "covariance quality = " << fit->covQual() << endl;
  cout << endl << endl;
  // Save resulting model to workspace
  RooWorkspace w("w");
  w.import(background);
  w.import(mc);
  TString wspacePath = TString::Format("workspaces/MingYan-MC.root");
  w.writeToFile(wspacePath);
  return;
}


#include "src/ModGaus.cxx"
#include "src/GausCB.cxx"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "FitPlotter.c"
using namespace RooFit;
void MingYanFit() {
  // Declare background fit pdf and variables
  double defaults[6] = {110,2,2,4,5,50};
//   double defaults[6] = {113.07, 1.97,1.67,12.44,7.72,49.2};
  RooRealVar m("llphoton_m","ll#gamma mass [GeV]",100,180);
  RooRealVar m0("m_{0}","mass peak value"             ,defaults[0],105,145);//m0.setError(0.1);
  RooRealVar vl("#nu_{L}","low-end exponential"       ,defaults[1],0.1,3);  //vl.setError(0.1);
  RooRealVar vr("#nu_{r}","high-end exponential"      ,defaults[2],0.1,3);  //vr.setError(0.1);
  RooRealVar s0("#sigma_{0}","peak width value"       ,defaults[3],2,15);   //s0.setError(0.1);
  RooRealVar sl("#sigma_{L}","low-end width value"    ,defaults[4],2,13);   //sl.setError(0.1);
  RooRealVar sh("#sigma_{H}","high-end width value"   ,defaults[5],40,60);  //sh.setError(0.1);
  ModGaus background("background","Modified Gaussian",m,m0,vl,vr,s0,sl,sh);
//   m0.setConstant(); vl.setConstant(); vr.setConstant(); s0.setConstant(); sl.setConstant(); sh.setConstant();
  // Declare signal fit pdf and variables (then set as constant)
  RooRealVar mu      ("#mu"          ,"mean"          ,124.76,12,126);
  RooRealVar gaus_std("#sigma_{gaus}","Gaussian stdev",2   ,1.0,50);
  RooRealVar cb_std  ("#sigma_{cb}"  ,"CB stdev"      ,1.52,0.1,100);
  RooRealVar cb_alpha("#alpha_{cb}"  ,"CB alpha"      ,1.22,1  ,10);
  RooRealVar cb_n    ("n_{cb}"       ,"CB n"          ,2.49,1  ,40);
  RooRealVar cb_frac ("f_{cb}"       ,"CB fraction"   ,0.78,0  ,1);
  mu.setConstant(); gaus_std.setConstant(); cb_std.setConstant();
  cb_alpha.setConstant(); cb_n.setConstant(); cb_frac.setConstant();
  GausCB signal("signal","Gaussian+Crystal Ball",m,mu,gaus_std,cb_std,cb_alpha,cb_n,cb_frac);
  // Combine signal and background pdfs
  RooRealVar norm_s("norm_s","N_{s}",-1000,1000);
  RooRealVar norm_b("norm_b","N_{b}",0,500000);
  const RooArgList components(signal,background);
  const RooArgList coeffs(norm_s,norm_b);
  RooAddPdf pdf("pdf","f_{s+b}",components,coeffs);
  // Load data
  TFile *infile = TFile::Open("/Users/adorsett/Desktop/CERNbox/Zgamma/MingYan-MC-mllg.root");
  TH1D *hmc = (TH1D*)infile->Get("MC");
//   for(int i = 0; i < hmc->GetNbinsX(); i++) {
//     hmc->SetBinContent(i, (int)hmc->GetBinContent(i));
//     hmc->SetBinError(i, sqrt(hmc->GetBinContent(i)));
//     cout << hmc->GetBinContent(i) << " +- " << hmc->GetBinError(i) << endl;
//     }
    
//   hmc->Rebin();
//   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooDataHist mc("mc","mc",RooArgSet(m),hmc);
  // Fit data with S+B model
//   pdf.fitTo(mc, Save(kTRUE), Offset(1),// Minimizer("Minuit2", "migrad"),
//             Strategy(2), Optimize(1), NumCPU(4));
//   pdf.fitTo(mc,RooFit::Extended());
//   FitPlotter(&mc, &pdf, m, 80,"MingYan-MC");
  RooFitResult *fit;
//   TVirtualFitter::SetMaxIterations(100000);
  fit = background.fitTo(mc,
                         //Verbose(0),PrintLevel(-1),
                         Save(kTRUE), 
                         SumW2Error(0),
                         Minimizer("Minuit2","migrad"),
                         NumCPU(4)); 
  FitPlotter(&mc, &background, m, 80,"MingYan-MC_B-only");
  cout << endl << endl;
  cout << "Background only 80 bin fit status = " << fit->status() << endl;
  cout << "covariance quality = " << fit->covQual() << endl;
  cout << endl << endl;
//   hmc->Rebin();
//   fit = background.fitTo(mc, Verbose(0),PrintLevel(-1), Save(kTRUE), SumW2Error(kFALSE));
//   cout << endl << endl;
//   cout << "Background only 40 bin fit status = " << fit->status() << endl;
//   cout << endl << endl;
//   FitPlotter(&mc, &background, m, 40,"MingYan-MC_B-only");
//   pdf.fitTo(mc,RooFit::Extended());
//   FitPlotter(&mc, &pdf, m, 40,"MingYan-MC");
//   background.fitTo(mc,Extended());
//   background.fitTo(mc, Save(kTRUE), Offset(1),// Minimizer("Minuit2", "migrad"),
//                    Strategy(2), Optimize(1), NumCPU(4));
//   // Save resulting model to workspace
//   RooWorkspace w("w");
//   w.import(pdf);
//   w.import(background);
//   w.import(mc);
//   TString wspacePath = TString::Format("workspaces/MingYan-MC.root");
//   w.writeToFile(wspacePath);
  return;
}


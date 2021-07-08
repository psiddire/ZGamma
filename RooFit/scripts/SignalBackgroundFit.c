#include "src/ModGaus.cxx"
#include "src/GausCB.cxx"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
using namespace RooFit;
void SignalBackgroundFit() {
  int sigmu(0);
  double lumi(1);
  // Declare background fit pdf and variables
  RooRealVar m("llphoton_m","ll#gamma mass [GeV]",125  ,100,180);
  RooRealVar m0("m_{0}","mass peak value"        ,112.8,100,180);
  RooRealVar vl("v_{L}","low-end exponential"    ,3.2  ,0  ,10);
  RooRealVar vh("v_{H}","high-end exponential"   ,4.2  ,0  ,10);
  RooRealVar s0("s_{0}","peak width value"       ,7.3  ,0  ,20);
  RooRealVar sl("s_{L}","low-end width value"    ,9.7  ,0  ,20);
  RooRealVar sh("s_{H}","high-end width value"   ,55.5 ,30 ,70);
  ModGaus background("background","Modified Gaussian",m,m0,vl,vh,s0,sl,sh);
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
  RooRealVar norm_s("norm_s","N_{s}",-50,50);
  RooRealVar norm_b("norm_b","N_{b}",0,15000);
  const RooArgList components(signal,background);
  const RooArgList coeffs(norm_s,norm_b);
  RooAddPdf pdf("pdf","f_{s+b}",components,coeffs);
  // Load data
  TChain chain("tree");
  TString infile = TString::Format("pseudodata/mu%d_%.0fifb.root",sigmu,lumi);
  chain.Add(infile);
  RooDataSet data("data","data",RooArgSet(m),Import(chain));
  // Fit data with S+B model
  pdf.fitTo(data,RooFit::Extended());
  // Plot data with fit
  TCanvas *can = new TCanvas("can","canvas",1200,800);
  can->Divide(1,2);
  can->cd(1);
  RooPlot *plot = m.frame();
  data.plotOn(plot,Binning(80));
  pdf.plotOn(plot,RooFit::Components("background"),RooFit::LineColor(kGreen));
  pdf.plotOn(plot,RooFit::LineColor(kRed));
  pdf.paramOn(plot,&data,"",2,"NELU",0.6,0.9,0.90);
  plot->Draw();
  can->cd(2);
  RooHist *pull = plot->pullHist();
  TLine line(0,0,1,0);
  line.SetNDC(false);
  line.SetLineStyle(2);
  line.SetLineColor(kBlack);
  line.SetLineWidth(2);
  pull->SetMarkerStyle(21);
  pull->GetXaxis()->SetRangeUser(100,180);
  pull->SetTitle("");
  pull->Draw("AP");
  line.Draw("same");
  TString outfile = TString::Format("plots/sb-fit_mu%d.pdf",sigmu);
  can->SaveAs(outfile);
  // Save resulting model to workspace
  RooWorkspace w("w");
  w.import(pdf);
  w.import(data);
  TString wspacePath = TString::Format("workspaces/mu%d_lumi%.0f.root",sigmu,lumi);
  w.writeToFile(wspacePath);
  return;
}


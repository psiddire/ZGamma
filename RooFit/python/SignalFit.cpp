#include "GausCB.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
using namespace RooFit;
// void SignalFit() {
int main() {
  RooRealVar m       ("llphoton_m"   ,"ll#gamma mass [GeV]" ,125 ,110,140);
  RooRealVar mu      ("#mu"          ,"mean"          ,125 ,124,126);
  RooRealVar gaus_std("#sigma_{gaus}","Gaussian stdev", 10 ,1.0,20);
  RooRealVar cb_std  ("#sigma_{cb}"  ,"CB stdev"      ,1.0 ,0.1,10);
  RooRealVar cb_alpha("#alpha_{cb}"  ,"CB alpha"      ,0.7 ,0 ,2);
  RooRealVar cb_n    ("n_{cb}"       ,"CB n"          ,9  ,1  ,20);
  RooRealVar cb_frac ("f_{cb}"       ,"CB fraction"   ,0.9 ,0,1);
  GausCB pdf("pdf","Gaussian+Crystal Ball",m,mu,gaus_std,cb_std,cb_alpha,cb_n,cb_frac);
  // Load data
  TString inpath = "MC/signal-shape.root";
  TFile *infile = TFile::Open(inpath);
  TTree *tree = (TTree*)infile->Get("tree");
  RooRealVar weight("weight","weight",1,-10,10);
  RooDataSet sig("sig","sig",tree,RooArgSet(m,weight),"1","weight");
  pdf.fitTo(sig);
  TCanvas *can = new TCanvas("can","canvas",800,800);
  can->cd();
  RooPlot *plot = m.frame();
  sig.plotOn(plot,Binning(30));
  pdf.plotOn(plot);
  pdf.paramOn(plot,&sig,"",3,"NELU",0.6,0.9,0.90);
  plot->Draw();
  can->SaveAs("plots/MCsig-fit.pdf");
  // Save resulting model to workspace
  RooWorkspace w("w");
  w.import(pdf);
  w.import(sig);
  TString wspacePath = "workspaces/weightedMC_signal.root";
  w.writeToFile(wspacePath);
  return;
}


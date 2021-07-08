#include "ModGaus.h"
#include "ModGaus.cxx"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
using namespace RooFit;
void MCFit() {
  RooRealVar m("llphoton_m","ll#gamma mass",125,100,180);
  RooRealVar m0("m_{0}","mass peak value",115,100,180);
  RooRealVar vl("v_{L}","low-end exponential",1,0,20);
  RooRealVar vh("v_{H}","high-end exponential",1,0,20);
  RooRealVar s0("s_{0}","peak width value",1,0,100);
  RooRealVar sl("s_{L}","low-end width value",2,0,100);
  RooRealVar sh("s_{H}","high-end width value",10,0,100);
//   RooRealVar n("n","Normalization constant",100,0,10000);
  ModGaus pdf("pdf","Modified Gaussian",m,m0,vl,vh,s0,sl,sh);
  TChain chain("tree");
  chain.Add("pseudodata/mu0_36ifb.root");
  RooDataSet mc("mc","mc",RooArgSet(m),Import(chain));
//   RooClassFactory::makePdf("ModGaus", "m,m0,vl,vh,s0,sl,sh,n");
//   RooClassFactory::makePdf("GausCB", "m,mean,gaus_std,cb_std,cb_alpha,cb_n,cb_frac");
  pdf.fitTo(mc);
  TCanvas *can = new TCanvas("can","canvas",1200,800);
  can->cd();
  RooPlot *plot = m.frame();
  mc.plotOn(plot);
  pdf.plotOn(plot);
  pdf.paramOn(plot,&mc,"",2,"NELU",0.6,0.9,0.90);
  plot->Draw();
  can->SaveAs("plots/mc-fit.pdf");
  return;
}


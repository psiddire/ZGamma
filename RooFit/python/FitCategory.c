#include <iostream>
#include <stdlib.h>
#include "ModGaus.cxx"
#include "ModGaus11.cxx"
#include "ModThreeGaus.cxx"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "scripts/PlotFit.c"
using namespace RooFit;

void FitCategory(TString category = "untag4",TString fitType = "ModGaus",
                 int nbin = 80, int xlow = 100, int xhigh = 180, 
                 bool useSumW2 = true) {
  TString outpath = TString::Format("plots/CategoryFits/"+fitType+"_%d-%d_%dbins",xlow,xhigh,nbin);
  if(!useSumW2) outpath += "_Poisson";
  outpath += "/"+category;
  gSystem->mkdir(outpath,kTRUE);
  gSystem->RedirectOutput(outpath+"/log.txt","w");
  RooAbsPdf* background;
  RooRealVar m("llphoton_m","ll#gamma mass [GeV]",xlow,xhigh);
  // Declare background fit pdf and variables
  Double_t defaults[6] = {120,2,0,10,10,50};
  RooRealVar m0("m_{0}","mass peak value [GeV]" ,defaults[0],105,145);  m0.setError(1.0); 
  RooRealVar vl("#nu_{L}","low-end power"       ,defaults[1],  0,  5);  vl.setError(0.1);
  RooRealVar vr("#Delta#nu","power range"       ,defaults[2], -5,  5);  vr.setError(0.1);
  RooRealVar s0("#sigma_{0}","peak width"       ,defaults[3],  1, 40);  s0.setError(1.0);
  RooRealVar sl("#sigma_{L}","low-end width"    ,defaults[4],  1, 40);  sl.setError(1.0);
  RooRealVar sh("#sigma_{H}","high-end width"   ,defaults[5], 20,100);  sh.setError(1.0);
  ModGaus* bM = new ModGaus("background","G_{M}",m,m0,vl,vr,s0,sl,sh);
  
  RooRealVar s1("#sigma_{1}","change in width"        ,40, 0, 100); s1.setError(1);
  ModGaus11* b11 = new ModGaus11("background","G_{1,1}",m,m0,vl,vr,s0,s1);
  if(fitType == "ModGaus11")
    background = b11;
  else background = bM;
  // Load data
  TString inpath = "/Users/adorsett/Desktop/CERNbox/Zgamma/draw_pico/RooFit/CatFiles/";
  if(category == "total") inpath += "mllg.root";
  else inpath += "bkghist_ele_mu_"+category+".root";
  TFile *infile = TFile::Open(inpath);

  TH1D *hmc;
  if(category == "total") hmc = (TH1D*)infile->Get("MC");
  else                    hmc = (TH1D*)infile->Get("hbkg");
  if(nbin == 40) hmc->Rebin();
  if(!useSumW2)
    for(int ib = 0; ib < nbin; ib++)
      if(hmc->GetBinContent(ib) < 0) hmc->SetBinContent(ib,0);

  RooDataHist mc(category,category,m,hmc);
  RooFitResult *fitResult;
  if(useSumW2) fitResult = background->chi2FitTo(mc, Verbose(true), PrintLevel(3),
                                                 Save(kTRUE), Range(xlow, xhigh), 
                                                 Minimizer("Minuit2","migrad"));
  else         fitResult = background->fitTo(    mc, Verbose(true), PrintLevel(3),
                                                 Save(kTRUE), Range(xlow, xhigh), Offset(true),
                                                 Minimizer("Minuit2","migrad"));
  PlotFit(fitResult, &mc, background, m, nbin,xlow,xhigh,outpath,true);
  gSystem->RedirectOutput(0);
  cout << "Category: " << category << "   " << "Fit type: " + fitType;
  if(useSumW2)    cout << " SumW2 errors";
  else            cout << " Poisson errors";
  cout << endl;
  cout << "Status, covQual = " << fitResult->status() << "," << fitResult->covQual() << endl;
  cout << endl;
  return;
}


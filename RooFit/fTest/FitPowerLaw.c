#include <iostream>
#include <stdlib.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
using namespace RooFit;

void FitPowerLaw(int order = 3) {

  TString category = "Inclusive";
  TString formula;
  TString fitType = "PowerLaw";
  int nbin = 80;
  int xlow = 100;
  int xhigh = 180;
  bool useSumW2 = true;

  TString outpath = "plots/";

  gSystem->mkdir(outpath,kTRUE);
  gSystem->RedirectOutput(outpath+"/log.txt","w");

  double sigma_pow, sigma_lpow, sigma_hpow;
  double turnon_pow, turnon_lpow, turnon_hpow;
  double par1, par3, par5;
  double par1_h, par3_h, par5_h;
  double par1_l, par3_l, par5_l;
  double coeff1, coeff3, coeff5;
  double coeff1_h, coeff3_h, coeff5_h;
  double coeff1_l, coeff3_l, coeff5_l;

  turnon_pow = 108;    turnon_lpow = 107.; turnon_hpow = 109.;                                                                                                      
  sigma_pow = 6.;      sigma_lpow = 3.;    sigma_hpow = 8.;                                                                                                         
  par1 = -9.;          par1_l = -15.;      par1_h = -5.;                                                                                                        
  coeff1 = 0.04;       coeff1_l = 0.;      coeff1_h = 1.;                                                                                                       
  par3 = -7.;          par3_l = -10.;      par3_h = -2.;                                                                                                        
  coeff3 = 4.0951e-13; coeff3_l = 0.;      coeff3_h = 1.;                                                                                                       
  par5 = -4.5113;      par5_l = -8;        par5_h = 2.;                                                                                                         
  coeff5 = 4.0951e-13; coeff5_l = 0.;      coeff5_h = 1.;                                                                                                       

  RooRealVar *m = new RooRealVar("CMS_llg_Mass","ll#gamma mass [GeV]",xlow,xhigh);
  RooRealVar *mean = new RooRealVar("Pow_mean","Pow_mean",0.);
  RooRealVar *sigma = new RooRealVar("Pow_sigma","Pow_sigma",sigma_pow,sigma_lpow,sigma_hpow);
  RooRealVar *turnon = new RooRealVar("Pow_turnon","Pow_turnon",turnon_pow,turnon_lpow,turnon_hpow);

  RooRealVar *p1 = new RooRealVar("Pow_p1","Pow_p1",par1,par1_l,par1_h);
  RooRealVar *p3 = new RooRealVar("Pow_p2","Pow_p2",par3,par3_l,par3_h);
  RooRealVar *p5 = new RooRealVar("Pow_p3","Pow_p3",par5,par5_l,par5_h);
  RooRealVar *cp1 = new RooRealVar("Pow_cp1","Pow_cp1",coeff1,coeff1_l,coeff1_h);
  RooRealVar *cp3 = new RooRealVar("Pow_cp3","Pow_cp3",coeff3,coeff3_l,coeff3_h);
  RooRealVar *cp5 = new RooRealVar("Pow_cp5","Pow_cp5",coeff5,coeff5_l,coeff5_h);

  if(order == 4){
    formula = "1e-20+(@0 > @1)*(@3*TMath::Power(@0,@2))";
  }
  else if(order == 5){
    formula = "1e-20+(@0 > @1)*(@3*TMath::Power(@0,@2)+@5*TMath::Power(@0,@4))";
  }
  else if(order == 6){
    formula = "1e-20+(@0 > @1)*(@3*TMath::Power(@0,@2)+@5*TMath::Power(@0,@4)+@7*TMath::Power(@0,@6))";
  }

  RooGenericPdf *step = new RooGenericPdf("Pow_step","Pow_step", formula, RooArgList(*m,*turnon,*p1,*cp1,*p3,*cp3,*p5,*cp5));
  RooGaussModel *gau = new RooGaussModel("Pow_gau","Pow_gau",*m,*mean,*sigma);
  RooFFTConvPdf *gauxpow = new RooFFTConvPdf("Pow_gauxpow","Pow_gauxpow",*m,*step,*gau);

  RooAbsPdf *background = gauxpow;

  // Load data
  TString inpath = "dataws.root";
  TFile *inFile = TFile::Open(inpath);
  RooWorkspace *inWS;
  RooAbsData *data;

  inWS = (RooWorkspace*)inFile->Get("CMS_llg_workspace");
  data = (RooAbsData*)inWS->data("Data_13TeV_inclusive");

  RooFitResult *fitResult;

  fitResult = background->fitTo(*data, Verbose(true), PrintLevel(3), Save(kTRUE), Range(xlow, xhigh), Offset(true), Minimizer("Minuit2","migrad"));

  gSystem->RedirectOutput(0);
  cout << "Category: " << category << "   " << "Fit type: " + fitType;
  cout << endl;
  cout << "Status, covQual = " << fitResult->status() << "," << fitResult->covQual() << endl;
  cout << endl;

  TCanvas *can = new TCanvas("can","canvas",800,800);

  RooPlot *plot = m->frame(Range(xlow,xhigh));
  data->plotOn(plot,Binning(nbin));
  background->plotOn(plot,RooFit::LineColor(kRed),RooFit::LineWidth(2));
  TString chi2Label = TString::Format("#chi^{2}_{R} = %.2f", plot->chiSquare(fitResult->floatParsFinal().getSize()));
  background->paramOn(plot,RooFit::Layout(0.54,0.96,0.89),RooFit::Format("NE",FixedPrecision(1)));
  plot->Draw();

  can->SaveAs(Form("plots/PowerLaw%d.pdf", order));

  return;
}

#include <iostream>
#include <stdlib.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
using namespace RooFit;

void FitLaurent(int order = 3) {

  TString category = "Inclusive";
  TString formula;
  TString fitType = "Laurent";
  int nbin = 80;
  int xlow = 100;
  int xhigh = 180;
  bool useSumW2 = true;

  TString outpath = "plots/";

  gSystem->mkdir(outpath,kTRUE);
  gSystem->RedirectOutput(outpath+"/log.txt","w");

  double sigma_gau, sigma_lgau, sigma_hgau;
  double turnon_gau, turnon_lgau, turnon_hgau;
  double coeff1, coeff2, coeff3, coeff4, coeff5, coeff6;
  double coeff1_h, coeff2_h, coeff3_h, coeff4_h, coeff5_h, coeff6_h;
  double coeff1_l, coeff2_l, coeff3_l, coeff4_l, coeff5_l, coeff6_l;

  turnon_gau = 106;  turnon_lgau = 102.; turnon_hgau = 110.;
  sigma_gau = 6.;    sigma_lgau = 2.;    sigma_hgau = 10.;
  coeff1 = 0.70; coeff1_l = 0.0;  coeff1_h = 1.0;
  coeff2 = 0.34; coeff2_l = 0.0;  coeff2_h = 1.0;
  coeff3 = 0.0038; coeff3_l = 0.0;  coeff3_h = 1.0;
  coeff4 = 0.09; coeff4_l = 0.0;  coeff4_h = 1.0;
  coeff5 = 0.1; coeff5_l = 0.0;  coeff5_h = 1.0;
  coeff6 = 0.1; coeff6_l = 0.0;  coeff6_h = 1.0;

  RooRealVar *m = new RooRealVar("CMS_llg_Mass","ll#gamma mass [GeV]",xlow,xhigh);
  RooRealVar *mean = new RooRealVar("Laurent_mean","Laurent_mean",0.);
  RooRealVar *sigma = new RooRealVar("Laurent_sigma","Laurent_sigma",sigma_gau,sigma_lgau,sigma_hgau);
  RooRealVar *turnon = new RooRealVar("Laurent_turnon","Laurent_turnon",turnon_gau,turnon_lgau,turnon_hgau);

  RooRealVar *cp1 = new RooRealVar("Laurent_cp1","Laurent_cp1",coeff1,coeff1_l,coeff1_h);
  RooRealVar *cp2 = new RooRealVar("Laurent_cp2","Laurent_cp2",coeff2,coeff2_l,coeff2_h);
  RooRealVar *cp3 = new RooRealVar("Laurent_cp3","Laurent_cp3",coeff3,coeff3_l,coeff3_h);
  RooRealVar *cp4 = new RooRealVar("Laurent_cp4","Laurent_cp4",coeff4,coeff4_l,coeff4_h);
  RooRealVar *cp5 = new RooRealVar("Laurent_cp5","Laurent_cp5",coeff5,coeff5_l,coeff5_h);
  RooRealVar *cp6 = new RooRealVar("Laurent_cp6","Laurent_cp6",coeff6,coeff6_l,coeff6_h);

  if(order == 4){
    formula = "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3)+@5*(@0)^(-6))";
  }
  else if(order == 5){
    formula = "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3)+@5*(@0)^(-6)+@6*(@0)^(-2))";
  }
  else if(order == 6){
    formula = "1e-20+(@0 > @1)*(@2*(@0)^(-4)+@3*(@0)^(-5)+@4*(@0)^(-3)+@5*(@0)^(-6)+@6*(@0)^(-2)+@7*(@0)^(-7))";
  }

  RooGenericPdf *step = new RooGenericPdf("Laurent_step","Laurent_step", formula, RooArgList(*m,*turnon,*cp1,*cp2,*cp3,*cp4,*cp5,*cp6));
  RooGaussModel *gau = new RooGaussModel("Laurent_gau","Laurent_gau",*m,*mean,*sigma);
  RooFFTConvPdf *gauxlau = new RooFFTConvPdf("Laurent_gauxlau","Laurent_gauxlau",*m,*step,*gau);

  RooAbsPdf *background = gauxlau;

  // Load data
  TString inpath = "dataws.root";
  TFile *inFile = TFile::Open(inpath);
  RooWorkspace *inWS;
  // RooDataSet *data;
  RooAbsData *data;

  // TH1D *hmc;
  // hmc = (TH1D*)inFile->Get("h_llg");
  inWS = (RooWorkspace*)inFile->Get("CMS_llg_workspace");
  // data = (RooDataSet*)inWS->data("Data_13TeV_inclusive");
  data = (RooAbsData*)inWS->data("Data_13TeV_inclusive");

  // RooDataHist *mc = new RooDataHist(category,category,*m,hmc);
  RooFitResult *fitResult;

  fitResult = background->fitTo(*data, Verbose(true), PrintLevel(3), Save(kTRUE), Range(xlow, xhigh), Offset(true), Minimizer("Minuit2","migrad"));

  //PlotFit(fitResult, data, background, *m, nbin, xlow, xhigh, outpath, true);
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

  background->paramOn(plot,data,chi2Label,2,"NELU",0.6,0.95,0.9);

  plot->Draw();

  can->SaveAs(Form("plots/Laurent%d.pdf", order));

  return;
}

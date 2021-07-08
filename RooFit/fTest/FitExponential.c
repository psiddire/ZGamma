#include <iostream>
#include <stdlib.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
using namespace RooFit;

void FitExponential(int order = 3) {

  TString category = "Inclusive";
  TString fitType = "Exponential";
  TString outpath = "plots/";
  int nbin = 80;
  int xlow = 100;
  int xhigh = 180;
  bool useSumW2 = true;

  gSystem->mkdir(outpath,kTRUE);
  gSystem->RedirectOutput(outpath+"/log.txt","w");

  // Load data
  TString inpath = "dataws.root";
  TFile *inFile = TFile::Open(inpath);
  RooWorkspace *inWS;
  RooAbsData *data;

  inWS = (RooWorkspace*)inFile->Get("CMS_llg_workspace");
  data = (RooAbsData*)inWS->data("Data_13TeV_inclusive");

  if(order == 4){
    double sigma_exp, sigma_lexp, sigma_hexp;
    double turnon_exp, turnon_lexp, turnon_hexp;
    double par1, par1_h, par1_l;
    double coeff1, coeff1_h, coeff1_l;

    sigma_exp = 4.;      sigma_lexp = 2.;    sigma_hexp = 6.;
    turnon_exp = 107.;   turnon_lexp = 105.; turnon_hexp = 109.;
    par1 = -0.026;        par1_l = -0.5;      par1_h = 0.;
    coeff1 = 0.6;        coeff1_l = 0.2;      coeff1_h = 1.0;

    RooRealVar *m = new RooRealVar("CMS_llg_Mass","ll#gamma mass [GeV]",xlow,xhigh);
    RooRealVar *mean = new RooRealVar("Exp_mean","Exp_mean",0.);
    RooRealVar *sigma = new RooRealVar("Exp_sigma","Exp_sigma",sigma_exp,sigma_lexp,sigma_hexp);
    RooRealVar *turnon = new RooRealVar("Exp_turnon","Exp_turnon",turnon_exp,turnon_lexp,turnon_hexp);

    RooRealVar *p1 = new RooRealVar("Exp_p1","Exp_p1",par1,par1_l,par1_h);
    RooRealVar *cp1 = new RooRealVar("Exp_cp1","Exp_cp1",coeff1,coeff1_l,coeff1_h);

    RooGenericPdf *step = new RooGenericPdf("Exp_step","Exp_step", "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2))", RooArgList(*m,*turnon,*p1,*cp1));
    RooGaussModel *gau = new RooGaussModel("Exp_gau","Exp_gau",*m,*mean,*sigma);
    RooFFTConvPdf *gauxexp = new RooFFTConvPdf("Exp_gauxexp","Exp_gauxexp",*m,*step,*gau);

    RooAbsPdf *background = gauxexp;
    RooFitResult *fitResult;
    fitResult = background->fitTo(*data, Verbose(true), PrintLevel(3), Save(kTRUE), Range(xlow, xhigh), Offset(true), Minimizer("Minuit2","migrad"));

    gSystem->RedirectOutput(0);
    cout << "Category: " << category << "   " << "Fit type: " + fitType;
    cout << "Status, covQual = " << fitResult->status() << "," << fitResult->covQual() << endl;

    TCanvas *can = new TCanvas("can","canvas",800,800);

    RooPlot *plot = m->frame(Range(xlow,xhigh));
    data->plotOn(plot,Binning(nbin));
    background->plotOn(plot,RooFit::LineColor(kRed),RooFit::LineWidth(2));
    TString chi2Label = TString::Format("#chi^{2}_{R} = %.2f", plot->chiSquare(fitResult->floatParsFinal().getSize()));
    background->paramOn(plot,RooFit::Layout(0.54,0.96,0.89),RooFit::Format("NE",FixedPrecision(1)));
    plot->Draw();

    can->SaveAs(Form("plots/Exponential%d.pdf", order));
  }
  else if(order == 5){
    double sigma_exp, sigma_lexp, sigma_hexp;
    double turnon_exp, turnon_lexp, turnon_hexp;
    double par1, par1_h, par1_l;
    double coeff1, coeff1_h, coeff1_l;
    double par3, par3_h, par3_l;
    double coeff3, coeff3_h, coeff3_l;

    sigma_exp = 4.;      sigma_lexp = 2.;    sigma_hexp = 6.;
    turnon_exp = 107.;   turnon_lexp = 105.; turnon_hexp = 109.;
    par1 = -0.026;        par1_l = -0.5;      par1_h = 0.;
    coeff1 = 0.6;        coeff1_l = 0.2;      coeff1_h = 1.0;
    par3 = -0.25;        par3_l = -0.5;      par3_h = 0.;
    coeff3 = 0.5;     coeff3_l = 0.;      coeff3_h = 1.;

    RooRealVar *m = new RooRealVar("CMS_llg_Mass","ll#gamma mass [GeV]",xlow,xhigh);
    RooRealVar *mean = new RooRealVar("Exp_mean","Exp_mean",0.);
    RooRealVar *sigma = new RooRealVar("Exp_sigma","Exp_sigma",sigma_exp,sigma_lexp,sigma_hexp);
    RooRealVar *turnon = new RooRealVar("Exp_turnon","Exp_turnon",turnon_exp,turnon_lexp,turnon_hexp);

    RooRealVar *p1 = new RooRealVar("Exp_p1","Exp_p1",par1,par1_l,par1_h);
    RooRealVar *cp1 = new RooRealVar("Exp_cp1","Exp_cp1",coeff1,coeff1_l,coeff1_h);
    RooRealVar *p3 = new RooRealVar("Exp_p3","Exp_p3",par3,par3_l,par3_h);
    RooRealVar *cp3 = new RooRealVar("Exp_cp3","Exp_cp3",coeff3,coeff3_l,coeff3_h);

    RooGenericPdf *step = new RooGenericPdf("Exp_step","Exp_step", "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2)+@5*TMath::Exp(@0*@4))", RooArgList(*m,*turnon,*p1,*cp1,*p3,*cp3));
    RooGaussModel *gau = new RooGaussModel("Exp_gau","Exp_gau",*m,*mean,*sigma);
    RooFFTConvPdf *gauxexp = new RooFFTConvPdf("Exp_gauxexp","Exp_gauxexp",*m,*step,*gau);

    RooAbsPdf *background = gauxexp;
    RooFitResult *fitResult;
    fitResult = background->fitTo(*data, Verbose(true), PrintLevel(3), Save(kTRUE), Range(xlow, xhigh), Offset(true), Minimizer("Minuit2","migrad"));

    gSystem->RedirectOutput(0);
    cout << "Category: " << category << "   " << "Fit type: " + fitType;
    cout << "Status, covQual = " << fitResult->status() << "," << fitResult->covQual() << endl;

    TCanvas *can = new TCanvas("can","canvas",800,800);

    RooPlot *plot = m->frame(Range(xlow,xhigh));
    data->plotOn(plot,Binning(nbin));
    background->plotOn(plot,RooFit::LineColor(kRed),RooFit::LineWidth(2));
    TString chi2Label = TString::Format("#chi^{2}_{R} = %.2f", plot->chiSquare(fitResult->floatParsFinal().getSize()));
    background->paramOn(plot,RooFit::Layout(0.54,0.96,0.89),RooFit::Format("NE",FixedPrecision(1)));
    plot->Draw();

    can->SaveAs(Form("plots/Exponential%d.pdf", order));
  }
  else if(order == 6){
    double sigma_exp, sigma_lexp, sigma_hexp;
    double turnon_exp, turnon_lexp, turnon_hexp;
    double par1, par1_h, par1_l;
    double coeff1, coeff1_h, coeff1_l;
    double par3, par3_h, par3_l;
    double coeff3, coeff3_h, coeff3_l;
    double par5, par5_h, par5_l;
    double coeff5, coeff5_h, coeff5_l;

    sigma_exp = 4.;      sigma_lexp = 2.;    sigma_hexp = 6.;
    turnon_exp = 107.;   turnon_lexp = 105.; turnon_hexp = 109.;
    par1 = -0.026;        par1_l = -0.5;      par1_h = 0.;
    coeff1 = 0.6;        coeff1_l = 0.2;      coeff1_h = 1.0;
    par3 = -0.25;        par3_l = -0.5;      par3_h = 0.;
    coeff3 = 0.5;     coeff3_l = 0.;      coeff3_h = 1.;
    par5 = -0.25;        par5_l = -0.5;      par5_h = 0.;
    coeff5 = 0.5;      coeff5_l = 0.;      coeff5_h = 1.;

    RooRealVar *m = new RooRealVar("CMS_llg_Mass","ll#gamma mass [GeV]",xlow,xhigh);
    RooRealVar *mean = new RooRealVar("Exp_mean","Exp_mean",0.);
    RooRealVar *sigma = new RooRealVar("Exp_sigma","Exp_sigma",sigma_exp,sigma_lexp,sigma_hexp);
    RooRealVar *turnon = new RooRealVar("Exp_turnon","Exp_turnon",turnon_exp,turnon_lexp,turnon_hexp);

    RooRealVar *p1 = new RooRealVar("Exp_p1","Exp_p1",par1,par1_l,par1_h);
    RooRealVar *cp1 = new RooRealVar("Exp_cp1","Exp_cp1",coeff1,coeff1_l,coeff1_h);
    RooRealVar *p3 = new RooRealVar("Exp_p2","Exp_p2",par3,par3_l,par3_h);
    RooRealVar *cp3 = new RooRealVar("Exp_cp3","Exp_cp3",coeff3,coeff3_l,coeff3_h);
    RooRealVar *p5 = new RooRealVar("Exp_p3","Exp_p3",par5,par5_l,par5_h);
    RooRealVar *cp5 = new RooRealVar("Exp_cp5","Exp_cp5",coeff5,coeff5_l,coeff5_h);

    RooGenericPdf *step = new RooGenericPdf("Exp_step","Exp_step", "1e-20+(@0 > @1)*(@3*TMath::Exp(@0*@2)+@5*TMath::Exp(@0*@4)+@7*TMath::Exp(@0*@6))", RooArgList(*m,*turnon,*p1,*cp1,*p3,*cp3,*p5,*cp5));
    RooGaussModel *gau = new RooGaussModel("Exp_gau","Exp_gau",*m,*mean,*sigma);
    RooFFTConvPdf *gauxexp = new RooFFTConvPdf("Exp_gauxexp","Exp_gauxexp",*m,*step,*gau);

    RooAbsPdf *background = gauxexp;
    RooFitResult *fitResult;
    fitResult = background->fitTo(*data, Verbose(true), PrintLevel(3), Save(kTRUE), Range(xlow, xhigh), Offset(true), Minimizer("Minuit2","migrad"));

    gSystem->RedirectOutput(0);
    cout << "Category: " << category << "   " << "Fit type: " + fitType;
    cout << "Status, covQual = " << fitResult->status() << "," << fitResult->covQual() << endl;

    TCanvas *can = new TCanvas("can","canvas",800,800);

    RooPlot *plot = m->frame(Range(xlow,xhigh));
    data->plotOn(plot,Binning(nbin));
    background->plotOn(plot,RooFit::LineColor(kRed),RooFit::LineWidth(2));
    TString chi2Label = TString::Format("#chi^{2}_{R} = %.2f", plot->chiSquare(fitResult->floatParsFinal().getSize()));
    background->paramOn(plot,RooFit::Layout(0.54,0.96,0.89),RooFit::Format("NE",FixedPrecision(1)));
    plot->Draw();

    can->SaveAs(Form("plots/Exponential%d.pdf", order));
  }
  return;
}

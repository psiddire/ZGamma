#include <iostream>
using namespace RooFit;
void GENZfit(int emu)
{
  //set file name
  TString type;
  if(fabs(emu)==11)  type = "ele";
  if(fabs(emu)==13)  type = "mu";
  char fname[500];
  sprintf(fname,"HZg_ggF_125GeV_ext1_M125_13TeV_powheg2_pythia8_%s.txt",type.Data());
  ofstream outfile(fname);
  //import the data
  TFile *f = TFile::Open(Form("../../outfile/GEN/MC_wcut_%s_HZg_ggF_125GeV.root",type.Data()));
  TTree *t = (TTree*)f->Get("tMC");
	
  //build GenZ model
  RooRealVar mcmll("mcmll" ,"mcmll", 40, 120, "GeV");
  RooRealVar eormu("eormu" ,"eormu", 10, 14, "GeV");
  RooDataSet Zmass("Zmass", "Zmass", RooArgSet(mcmll), Import(*t));
  RooRealVar meanCB("meanCB","meanCB", 8.64761e+01, 80,100);
  RooRealVar sigmaCB("sigmaCB","sigmaCB", 4.99069e+00, 4., 10.);
  RooRealVar alphaCB("alphaCB","alphaCB",5.51008e-01, 0.1,1.2);
  RooRealVar nCB("nCB","nCB", 4.99983e+01, 0, 100);
  RooRealVar meanGauss1("meanGauss1","meanGauss1",70., 40, 110);
  RooRealVar sigmaGauss1("sigmaGauss1","sigmaGauss1",8.08747e+00, 0, 20);
  RooRealVar f1("f1","f1", 0.00001,1.0);
  RooRealVar meanGauss2("meanGauss2","meanGauss2",9.10697e+01, 80, 100);
  RooRealVar sigmaGauss2("sigmaGauss2","sigmaGauss2",8.03251e-01, 0.1, 5.1);
  RooRealVar f2("f2","f2",7.40547e-01, 0.0001, 1.0);
  RooRealVar meanGauss3("meanGauss3","meanGauss3",9.05542e+01, 80, 100);
  RooRealVar sigmaGauss3("sigmaGauss3","sigmaGauss3", 1.86069e+00, 0, 10);
  RooRealVar f3("f3","f3",7.38652e-01,0.0001,1.0);
  RooCBShape* singleCB = new RooCBShape("singleCB", "", mcmll, meanCB, sigmaCB, alphaCB, nCB);
  RooGaussian* gaussShape1 = new RooGaussian("gaussShape1", "", mcmll, meanGauss1, sigmaGauss1);
  RooAddPdf* CBplusGauss = new RooAddPdf("CBplusGauss", "", *singleCB, *gaussShape1, f1);
  RooGaussian* gaussShape2 = new RooGaussian("gaussShape2", "", mcmll, meanGauss2, sigmaGauss2);
  RooAddPdf* CBplusGaussplusGauss = new RooAddPdf("CBplusGaussplusGauss", "", *CBplusGauss, *gaussShape2, f2);
  RooGaussian* gaussShape3 = new RooGaussian("gaussShape3", "", mcmll, meanGauss3, sigmaGauss3);
  RooAddPdf* CBplusGaussplusGaussplusGauss = new RooAddPdf("CBplusGaussplusGaussplusGauss", "", *CBplusGaussplusGauss, *gaussShape3, f3);
  RooFitResult* GENZ = CBplusGaussplusGaussplusGauss->fitTo(Zmass, Save(kTRUE));
  RooPlot* xframe4 = mcmll.frame(40,120) ;
  Zmass.plotOn(xframe4,Binning(80), RooFit::Name("Zmass")) ;
  CBplusGaussplusGaussplusGauss->plotOn(xframe4,RooFit::Name("CBplusGaussplusGaussplusGauss"),LineColor(TColor::GetColor("#d9594c")));
  CBplusGaussplusGaussplusGauss->plotOn(xframe4,RooFit::Components("CBplusGaussplusGauss"),LineStyle(kDashed),LineColor(TColor::GetColor("#ce8d66")));
  CBplusGaussplusGaussplusGauss->plotOn(xframe4,RooFit::Components("gaussShape3"),LineStyle(kDashed),LineColor(TColor::GetColor("#b7b868")));
  CBplusGaussplusGaussplusGauss->paramOn(xframe4,Layout(0.6,0.96,0.95));
  xframe4->Draw();

  RooWorkspace *w  = new RooWorkspace("w","");
  w->import(Zmass);
  w->import(*CBplusGaussplusGaussplusGauss);
  w->import(*GENZ);
  w->Print() ;	
  // w->WritetoFile("ele_root");
}

#include <iostream>
#include <stdlib.h>
#include "ModGaus01.cxx"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TCanvas.h"

#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1.h"
 

using namespace RooFit;

void FitCategoryBlind()
{

TString category = "Inclusive";

int nbin = 80;
int xlow = 100;
int xhigh = 180;
bool useSumW2 = true;

TString fitType = "ModGaus01";

RooRealVar m("llphoton_m","ll#gamma mass [GeV]",xlow,xhigh);
// Declare background fit pdf and variables
Double_t defaults[6] = {120,2,0,10,10,50};
RooRealVar m0("m_{0}","mass peak value [GeV]" ,defaults[0],105,145);  m0.setError(1.0); 
RooRealVar vl("#nu_{L}","low-end power"       ,defaults[1],  0,  5);  vl.setError(0.1);
RooRealVar vr("#Delta#nu","power range"       ,defaults[2], -5,  5);  vr.setError(0.1);
RooRealVar s0("#sigma_{0}","peak width"       ,defaults[3],  1, 40);  s0.setError(1.0);
RooRealVar sl("#sigma_{L}","low-end width"    ,defaults[4],  1, 40);  sl.setError(1.0);
RooRealVar sh("#sigma_{H}","high-end width"   ,defaults[5], 20,100);  sh.setError(1.0);
RooRealVar s1("#sigma_{1}", "change in width", 40, 0, 100); s1.setError(1);

ModGaus01* background = new ModGaus01("background", "G_{0,1}", m, m0, vl, s0, s1);

// Load data
TString inpath = "Data.root";
TFile *infile = TFile::Open(inpath);

TTree *tree = (TTree*)infile->Get("opttree");

RooDataSet data(category, category, RooArgSet(m), Import(*tree));

m.setRange("full", 100, 180);
m.setRange("left", 100, 120);
m.setRange("right", 130, 180);

auto blindedData = data.reduce(CutRange("left,right"));

background->fitTo(*blindedData);

TCanvas *canvas=new TCanvas("canvas","canvas",800,600);

RooPlot* plotFrameWithNormRange = m.frame(RooFit::Title("Right: All slices have common normalisation"));

blindedData->plotOn(plotFrameWithNormRange);
background->plotOn(plotFrameWithNormRange, LineColor(kBlue),  Range("left"),  RooFit::NormRange("left,right"));
background->plotOn(plotFrameWithNormRange, LineColor(kGreen), Range("right"), RooFit::NormRange("left,right"));
background->plotOn(plotFrameWithNormRange, LineColor(kRed),   Range("full"),  RooFit::NormRange("left,right"), LineStyle(10));

plotFrameWithNormRange->Draw();
 
canvas->Draw();

}

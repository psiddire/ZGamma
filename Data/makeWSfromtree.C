#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooConstVar.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "TFile.h"
#include <string>

void makeWSfromtree() {
  const int ncats = 1;
  RooRealVar mass("CMS_llg_Mass", "ll#gamma mass [GeV]", 100, 180);
  RooRealVar lumi("IntLumi", "Integrated luminosity", 41.5, "fb^{-1}");
  RooRealVar sqrts("SqrtS","Center of Mass Energy", 13, "TeV");
  RooDataSet * dataset[ncats] = {};
  string catname[ncats] = {"inclusive"};
  for (int i = 0; i < ncats; i++){
    string temp = "Data_13TeV_" + catname[i];
    char * histname = new char [temp.length()+1];
    strcpy (histname, temp.c_str());
    dataset[i] =  new RooDataSet(histname, histname, RooArgList(mass));
  }
  TFile *file = new TFile("/afs/hep.wisc.edu/home/psiddire/CMSSW_10_2_13/src/ZGamma/Data/Tree/Data.root");
  TTree *tree = (TTree*)file->Get("opttree");
  double llphoton_m;
  tree->SetBranchAddress("llphoton_m", &llphoton_m);
  int entries = tree->GetEntries();
  for (int i = 0; i < entries; i++){
    tree->GetEntry(i);
    mass.setVal(llphoton_m);
    dataset[0]->add(RooArgList(mass));
  }
  RooWorkspace *w = new RooWorkspace("CMS_llg_workspace", "CMS_llg_workspace");
  w->import(lumi);
  w->import(sqrts);
  for (int i = 0; i < ncats; i++){
    w->import(*dataset[i]);
  }
  w->writeToFile("dataws.root");
  w->Print();
}

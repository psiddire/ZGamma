#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
void PlotCovariances(RooFitResult *fit, TString outpath);
void PlotLikelihood(RooFitResult *fit, RooAbsPdf *model, RooAbsData *data, TString outpath);
using namespace RooFit;
double FitPlotter(RooFitResult *fit, RooAbsData *data, RooAbsPdf *pdf, 
                  RooRealVar m, int nbins, TString outfile, 
                  bool doCovariance = false, bool reduceRange = false) {
  gErrorIgnoreLevel = kWarning;
  double xlow(0.07), xhigh(0.95);
  double ylow(0.05), yhigh(0.90);
  double hratio(0.25);
  TCanvas *can = new TCanvas("can","canvas",850,800);
  can->Divide(1,2);
  TVirtualPad *top = can->cd(1);
  TVirtualPad *bot = can->cd(2);
  // Margin: L, R, B, T
  // Pad xl, yl, xu, yu
  top->SetMargin(xlow,1-xhigh,0.0 ,1-yhigh);
  bot->SetMargin(xlow,1-xhigh,0.30,0.0);
  top->SetPad(0.0,hratio,1.0,1.0   );
  bot->SetPad(0.0,   0.0,1.0,hratio);
  bot->SetGridy(1);
  top->cd();
  double fitMax(180);
  if(reduceRange) fitMax = 160;
  RooPlot *plot = m.frame(Range(100,fitMax));
//   data->plotOn(plot,Binning(nbins), Range(100,fitMax));
  cout << outfile << endl;
  if(outfile.Contains("NoSumW2")) { 
    data->plotOn(plot,Binning(nbins),DataError(RooAbsData::Poisson));
    cout << "NoSumW2" << endl;
  }
  else data->plotOn(plot,Binning(nbins));
  pdf->plotOn(plot,RooFit::LineColor(kRed),RooFit::LineWidth(2));
  TString chi2Label = TString::Format("#chi^{2}_{R} = %.2f",plot->chiSquare());
  pdf->paramOn(plot,data,chi2Label,2,"NELU",0.6,xhigh,yhigh);
  plot->GetYaxis()->SetTitleOffset(0.9);
  plot->GetYaxis()->SetLabelSize(0.03);
  plot->GetXaxis()->SetLimits(100,fitMax);
  TString funcType = pdf->GetTitle();
  TString category = outfile.Tokenize("_")->At(1)->GetName();
//   plot->SetTitle(TString::Format("%s fit with %s (status = %s, covQual = %s)",
//                                  category, funcType, status, covQual));
  plot->SetTitle(category+" fit with "+funcType+
                 " (status = "+to_string(fit->status())+
                 ", covQual = "+to_string(fit->covQual())+")");
  plot->Draw("same");
  bot->cd();
  RooHist *pull = plot->pullHist();
  TH1D* hpull = new TH1D("hpull","m_{ll#gamma} pull distribution",8,-4,4);
  Double_t ax[nbins], ay[nbins];
  double chisquared(0);
  for(int i = 0; i < nbins; i++) {
    pull->GetPoint(i,ax[i],ay[i]);
    hpull->Fill(ay[i]);
    chisquared += fabs(ay[i]);
  }
  chisquared = chisquared/(nbins-1);
  pull->SetMarkerStyle(21);
  pull->GetXaxis()->SetLimits(100,fitMax);
  pull->SetTitle(TString::Format(";m_{ll#gamma} [GeV]; Pull(%.3f)",
                                   hpull->GetMean()));
  pull->GetXaxis()->SetTitleSize(0.1);
  pull->GetYaxis()->SetTitleSize(0.1);
  pull->GetXaxis()->SetLabelSize(0.08);
  pull->GetYaxis()->SetRangeUser(-4.2,4.2);
  pull->GetYaxis()->SetLabelSize(0.08);
  pull->GetYaxis()->CenterTitle();
  pull->GetYaxis()->SetTickLength(0);
  pull->GetYaxis()->SetTitleOffset(0.2);
  pull->Draw("AP same");
  if(nbins == 40) outfile.ReplaceAll("Fits/","Fits/nb40/");
  gSystem->mkdir(outfile);
  can->SaveAs(outfile+".pdf");
  TCanvas *can2 = new TCanvas("can2","pull canvas",800,800);
  can2->cd();
  hpull->SetLineColor(kBlack);
  hpull->SetFillColor(kAzure);
  hpull->Draw("B");
  can2->SaveAs(outfile+"/pulls.pdf");
  can2->Close();
  can->Close();
  hpull->Delete();
  if(doCovariance) PlotCovariances(fit, outfile);
  PlotLikelihood(fit,pdf,data,outfile);
  return chisquared;
}

void PlotCovariances(RooFitResult *fit, TString outpath){
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;
  RooArgList pars = fit->floatParsFinal();
  TCanvas *can = new TCanvas("can","canvas",900,800); can->cd();
  for(size_t ip1(0); ip1 < pars.getSize(); ip1++) {
    for(size_t ip2(ip1+1); ip2 < pars.getSize(); ip2++) {
      RooRealVar *p1 = (RooRealVar*)pars.find(pars.at(ip1)->namePtr()->GetName());
      RooRealVar *p2 = (RooRealVar*)pars.find(pars.at(ip2)->namePtr()->GetName());
      RooPlot *plot = new RooPlot(*p1,*p2, p1->getVal()-1.1*p1->getError(),
                                           p1->getVal()+1.1*p1->getError(),
                                           p2->getVal()-1.1*p2->getError(),
                                           p2->getVal()+1.1*p2->getError());
      fit->plotOn(plot,*p1, *p2, "ME12ABHV");
      plot->GetYaxis()->SetLabelSize(0.025);
      plot->GetXaxis()->SetLabelSize(0.025);
      plot->Draw();
      TString outfile = outpath+"/Covariance_p"+to_string(ip1)+"-p"+to_string(ip2)+".pdf";
      can->SaveAs(outfile);
      can->Clear();
    }
  }
  TH2* cor = fit->correlationHist();
  cor->GetYaxis()->SetLabelSize(0.1);
  cor->GetXaxis()->SetLabelSize(0.1);
  cor->Draw("colz");
  can->SaveAs(outpath+"/Correlation.pdf");
  can->Close();
  cor->Delete();
  return;
}

void PlotLikelihood(RooFitResult *fit, RooAbsPdf *model, RooAbsData *data, TString outpath) {
//   gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;
  RooArgList pars = fit->floatParsFinal();
  TCanvas *can = new TCanvas("can","canvas",1000,800); can->cd();
  for(size_t ip(0); ip < pars.getSize(); ip++) {
    RooRealVar *par = (RooRealVar*)pars.find(pars.at(ip)->namePtr()->GetName());
    RooAbsReal *nll = model->createNLL(*data);
    TString title = "NLL scan for "+string(par->GetName());
    RooPlot *frame = par->frame(Bins(10),Range(par->getVal()-3*par->getError(),par->getVal()+3*par->getError()), 
                              Title(title));
    nll->plotOn(frame,ShiftToZero());
    frame->Draw();
    TString outfile = outpath+"/"+title+".pdf";
    can->SaveAs(outfile);
    can->Clear();
  }
  return;
}
  





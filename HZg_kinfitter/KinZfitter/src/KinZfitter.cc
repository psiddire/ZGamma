#undef AnalysisClass_cxx
/*************************************************************************
 *  Authors:   Tongguang Cheng
 *  Removed the 2nd Z info for HZg analysis- Ming-Yan Lee
 *************************************************************************/
#ifndef KinZfitter_cpp
#define KinZfitter_cpp
// using namespace RooFit;
/// KinFitter header
#include "ZGamma/HZg_kinfitter/KinZfitter/interface/KinZfitter.h"
#include "ZGamma/HZg_kinfitter/HelperFunction/interface/HelperFunction.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RooWorkspace.h"
#include "RooProduct.h"
#include "RooProdPdf.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "time.h"
///----------------------------------------------------------------------------------------------
/// KinZfitter::KinZfitter - constructor/
///----------------------------------------------------------------------------------------------

KinZfitter::KinZfitter(bool isData)
{    
  // modified with HZg
  PDFName_ = "HZg_ggF_125GeV_ext1_M125_13TeV_powheg2_pythia8";

  debug_ = false;

  if(debug_) std::cout << "KinZfitter. The debug flag is ON with "<<PDFName_<< std::endl;
	
  // / Initialise HelperFunction
  helperFunc_ = new HelperFunction();
  isCorrPTerr_ = 0; 
  isData_ = isData; 

}

KinZfitter::~KinZfitter()
{
  delete helperFunc_;
}

void KinZfitter::Setup(TreeReader &data, vector<int> &lepID, vector<int> &fsrID, vector<double> &errpT, int eormu, bool isData,  std::map<unsigned int, TLorentzVector> selectedLeptons, std::map<unsigned int, TLorentzVector> selectedFsrPhotons){

  // reset everything for each event
  idsZ1_.clear(); /* idsZ2_.clear();      */ 
  idsFsrZ1_.clear(); /* idsFsrZ2_.clear(); */
	
  p4sZ1_.clear();p4sZ1ph_.clear(); /*  p4sZ2_.clear(); p4sZ2ph_.clear(); */
  p4sZ1REFIT_.clear();  p4sZ1phREFIT_.clear(); /* p4sZ2REFIT_.clear();p4sZ2phREFIT_.clear(); */
     
  pTerrsZ1_.clear();  pTerrsZ1ph_.clear(); /* pTerrsZ2_.clear(); pTerrsZ2ph_.clear(); */
  pTerrsZ1REFIT_.clear();  pTerrsZ1phREFIT_.clear(); /* pTerrsZ2REFIT_.clear();pTerrsZ2phREFIT_.clear(); */

  
  initZs( data, lepID, fsrID, errpT, eormu, isData,  selectedLeptons, selectedFsrPhotons);
  fs_=""; 
  if(eormu== 0 ) fs_="ele";
  else  fs_="mu";   
  if(debug_) cout<<"fs is "<<fs_<<endl;    
  if(debug_) cout<<"list ids"<<endl;
  if(debug_) cout<<"IDs[0] "<<idsZ1_[0]<<" IDs[1] "<<idsZ1_[1]/* <<" IDs[2] "<<idsZ2_[0]<<" IDs[3] "<<idsZ2_[1] */<<endl;


  /////////////
  string paramZ1_dummy = "KinZfitter/KinZfitter/ParamZ1/dummy.txt";
    
  TString paramZ1 = TString( paramZ1_dummy.substr(0,paramZ1_dummy.length() - 9));

  paramZ1+=PDFName_;
  paramZ1+="_";
  paramZ1+=+fs_;
  paramZ1+=".txt";

  if(debug_) cout<<"paramZ1 in "<<paramZ1<<endl;
  
  std::ifstream input(paramZ1);
  std::string line;
  while (!input.eof() && std::getline(input,line))
    {
      std::istringstream iss(line);
      string p; double val;
      if(iss >> p >> val) {
	if(p=="sg")  { sgVal_ = val;}
	if(p=="a" )  { aVal_ = val;}
	if(p=="n" )  { nVal_ = val;}
	if(p=="f")   { fVal_ = val;}

	if(p=="mean" )  { meanVal_ = val;}
	if(p=="sigma" )  { sigmaVal_ = val;}
	if(p=="f1")   { f1Val_ = val;}

	if(p=="meanCB")   { meanCB_ = val;}
	if(p=="sigmaCB")   { sigmaCB_ = val;}
	if(p=="alphaCB")   { alphaCB_ = val;}
	if(p=="nCB")   { nCB_ = val;}
	if(p=="meanGauss1")   { meanGauss1_ = val;}
	if(p=="sigmaGauss1")   { sigmaGauss1_ = val;}
	if(p=="f1")   { f1_ = val;}
	if(p=="meanGauss2")   { meanGauss2_ = val; }
	if(p=="sigmaGauss2")   { sigmaGauss2_ = val;}
	if(p=="f2")   { f2_ = val;}
	if(p=="meanGauss3")   { meanGauss3_ = val;}
	if(p=="sigmaGauss3")   { sigmaGauss3_ = val;}
	if(p=="f3")   { f3_ = val;}

      }
    }
}

///----------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------

void KinZfitter::initZs( TreeReader &data, vector<int> &lepID, vector<int> &fsrID, vector<double> &errpT,int eormu, bool isData, std::map<unsigned int, TLorentzVector>  selectedLeptons, std::map<unsigned int, TLorentzVector> selectedFsrPhotons){
  if(debug_) cout<<"init leptons"<<endl;


  for(unsigned int il = 0; il<selectedLeptons.size(); il++)
    {
      double pTerr = 0; TLorentzVector p4;
      p4 = selectedLeptons[il];
      pTerr = helperFunc_->pterr(data, lepID[il],-99, eormu, isData); //put lepton err first
      idsZ1_.push_back(eormu);
      pTerrsZ1_.push_back(pTerr);
      p4sZ1_.push_back(p4);            
    }
  errpT = pTerrsZ1_ ;//!!!!!!!REMEMBER TO PUT BACK
  
  if(debug_) cout<<"init fsr photons"<<endl;
  if(selectedFsrPhotons.size()>0)
    {
      for(unsigned int ifsr = 0; ifsr<2; ifsr++)
	{
	  TLorentzVector p4 = selectedFsrPhotons[ifsr];
	  if(selectedFsrPhotons[ifsr].Pt()==0) continue;
	  if(fsrID[ifsr]==-99) continue;
	  if(debug_) cout<<"ifsr "<<ifsr<<endl;
	  double pTerr = 0;
	  pTerr = helperFunc_->pterr(data, lepID[ifsr], fsrID[ifsr], eormu, isData); 
	  // pTerr = 0.1;

	  if(debug_) cout<<" pt err is "<<pTerr<<endl;

	  cout<<pTerr<<endl;
	  pTerrsZ1ph_.push_back(pTerr);
	  p4sZ1ph_.push_back(p4);
	  idsFsrZ1_.push_back(idsZ1_[ifsr]);
	}
      if(debug_) cout<<"p4sZ1ph_ "<<p4sZ1ph_.size()<<endl;
    }  
}

void KinZfitter::SetZResult(double l1, double l2, double lph1, double lph2)
{

  if(debug_) cout<<"start set Z result"<<endl;

  // pT scale after refitting w.r.t. reco pT
  lZ1_l1_ = l1; lZ1_l2_ = l2;

  if(debug_) cout<<"l1 "<<l1<<" l2 "<<l2<<endl;

  lZ1_ph1_ = lph1; lZ1_ph2_ = lph2;

  TLorentzVector Z1_1 = p4sZ1_[0]; TLorentzVector Z1_2 = p4sZ1_[1];

  TLorentzVector Z1_1_True(0,0,0,0);
  Z1_1_True.SetPtEtaPhiM(lZ1_l1_*Z1_1.Pt(),Z1_1.Eta(),Z1_1.Phi(),Z1_1.M());
  TLorentzVector Z1_2_True(0,0,0,0);
  Z1_2_True.SetPtEtaPhiM(lZ1_l2_*Z1_2.Pt(),Z1_2.Eta(),Z1_2.Phi(),Z1_2.M());

  p4sZ1REFIT_.push_back(Z1_1_True); p4sZ1REFIT_.push_back(Z1_2_True);

  for(unsigned int ifsr1 = 0; ifsr1 < p4sZ1ph_.size(); ifsr1++){

    TLorentzVector Z1ph = p4sZ1ph_[ifsr1];
    TLorentzVector Z1phTrue(0,0,0,0);
  
    double l = 1.0;
    if(ifsr1==0) l = lZ1_ph1_; if(ifsr1==1) l = lZ1_ph2_;
  
    Z1phTrue.SetPtEtaPhiM(l*Z1ph.Pt(),Z1ph.Eta(),Z1ph.Phi(),Z1ph.M());
  
    p4sZ1phREFIT_.push_back(Z1phTrue);
    pTerrsZ1phREFIT_.push_back(pTerrsZ1ph_[ifsr1]);

  }

  if(debug_) cout<<"end set Z1 result"<<endl;
}

double KinZfitter::GetRefitMZ1()
{

  vector<TLorentzVector> p4s = GetRefitP4s();

  TLorentzVector pZ1(0,0,0,0);

  pZ1 = p4s[0] + p4s[1];

  return pZ1.M();

}

double KinZfitter::GetMZ1()
{

  vector<TLorentzVector> p4s = GetP4s();

  TLorentzVector pZ1(0,0,0,0);

  pZ1 = p4s[0] + p4s[1];

  return pZ1.M();

}

double KinZfitter::GetMZ1Err()
{

  vector<TLorentzVector> p4s;
  vector<double> pTErrs;

  p4s.push_back(p4sZ1_[0]);p4s.push_back(p4sZ1_[1]);
  pTErrs.push_back(pTerrsZ1_[0]); pTErrs.push_back(pTerrsZ1_[1]);
  for(unsigned int ifsr1 = 0; ifsr1<p4sZ1ph_.size(); ifsr1++){

    p4s.push_back(p4sZ1ph_[ifsr1]);
    pTErrs.push_back(pTerrsZ1ph_[ifsr1]);

  }
  if(debug_) cout<<"Get MZ1 Err"<<endl;
  return helperFunc_->masserror(p4s,pTErrs);
  // return 0.5;
}


vector<TLorentzVector> KinZfitter::GetRefitP4s()
{

  TLorentzVector Z1_1 = p4sZ1REFIT_[0]; TLorentzVector Z1_2 = p4sZ1REFIT_[1];

  // fsr photons

  for(unsigned int ifsr1 = 0; ifsr1<p4sZ1phREFIT_.size(); ifsr1++){

    int id_fsr1 = idsFsrZ1_[ifsr1];
    TLorentzVector Z1ph = p4sZ1phREFIT_[ifsr1];

    if(id_fsr1==idsZ1_[0]) Z1_1 = Z1_1 + Z1ph;
    if(id_fsr1==idsZ1_[1]) Z1_2 = Z1_2 + Z1ph;

  }

  vector<TLorentzVector> p4s;
  p4s.push_back(Z1_1); p4s.push_back(Z1_2);
  return p4s;

}

vector<TLorentzVector> KinZfitter::GetP4s()
{

  TLorentzVector Z1_1 = p4sZ1_[0]; TLorentzVector Z1_2 = p4sZ1_[1];

  // fsr photons

  for(unsigned int ifsr1 = 0; ifsr1<p4sZ1ph_.size(); ifsr1++){

    int id_fsr1 = idsFsrZ1_[ifsr1];
    TLorentzVector Z1ph = p4sZ1ph_[ifsr1];
    if(id_fsr1==idsZ1_[0]) Z1_1 = Z1_1 + Z1ph;
    if(id_fsr1==idsZ1_[1]) Z1_2 = Z1_2 + Z1ph;

  }

  vector<TLorentzVector> p4s;
  p4s.push_back(Z1_1); 
  p4s.push_back(Z1_2);
  return p4s;

}

void KinZfitter::KinRefitZ()
{
  double l1,l2,lph1,lph2;

  l1 = 1.0; l2 = 1.0; lph1 = 1.0; lph2 = 1.0;

  SetFitInput(fitInput1, p4sZ1_, pTerrsZ1_, p4sZ1ph_, pTerrsZ1ph_);
  if (debug_) cout<<"done input"<<endl;
  Driver(fitInput1, fitOutput1);
  if (debug_) cout<<"done driver"<<endl;
  SetFitOutput(fitInput1, fitOutput1, l1, l2, lph1, lph2, pTerrsZ1REFIT_, pTerrsZ1phREFIT_, covMatrixZ1_);

  if(debug_) cout<<"l1 "<<l1<<"; l2 "<<l2<<" lph1 "<<lph1<<" lph2 "<<lph2<<endl;

  SetZResult(l1, l2, lph1, lph2/* , l3, l4, lph3, lph4 */);

  if(debug_) cout<<"Z refit done"<<endl;
}

void  KinZfitter::Driver(KinZfitter::FitInput &input, KinZfitter::FitOutput &output) {

  MakeModel(input, output);

}


void  KinZfitter::SetFitInput(KinZfitter::FitInput &input, 
                              vector<TLorentzVector> ZLep, vector<double> ZLepErr,
                              vector<TLorentzVector> ZGamma, vector<double> ZGammaErr) {

  TLorentzVector lep1 = ZLep[0]; TLorentzVector lep2 = ZLep[1];

  input.pTRECO1_lep = lep1.Pt(); input.pTRECO2_lep = lep2.Pt();
  input.pTErr1_lep = ZLepErr[0]; input.pTErr2_lep = ZLepErr[1];
  input.theta1_lep = lep1.Theta(); input.theta2_lep = lep2.Theta();
  input.phi1_lep = lep1.Phi(); input.phi2_lep = lep2.Phi();
  input.m1 = lep1.M(); input.m2 = lep2.M();

  input.nFsr = 0;
  TLorentzVector nullFourVector(0, 0, 0, 0);
  TLorentzVector Gamma1, Gamma2;
  Gamma1 = nullFourVector; Gamma2 = nullFourVector;

  input.pTRECO1_gamma = Gamma1.Pt(); input.pTErr1_gamma = 0;
  input.theta1_gamma = Gamma1.Theta(); input.phi1_gamma = Gamma1.Phi();
  input.pTRECO2_gamma = Gamma2.Pt(); input.pTErr2_gamma = 0;
  input.theta2_gamma = Gamma2.Theta(); input.phi2_gamma = Gamma1.Phi();


  if (int(ZGamma.size()) >= 1) {

    input.nFsr = 1;
    TLorentzVector gamma1 = ZGamma[0];
    input.pTRECO1_gamma = gamma1.Pt(); input.pTErr1_gamma = ZGammaErr[0];
    input.theta1_gamma = gamma1.Theta(); input.phi1_gamma = gamma1.Phi();

  }

  if (int(ZGamma.size()) == 2) {

    input.nFsr = 2;
    TLorentzVector gamma2 = ZGamma[1];
    input.pTRECO2_gamma = gamma2.Pt(); input.pTErr2_gamma = ZGammaErr[1];
    input.theta2_gamma = gamma2.Theta(); input.phi2_gamma = gamma2.Phi();

  }

  //      /*if (debug_)*/ cout << "nFsr: " << input.nFsr << endl;
}


void KinZfitter::SetFitOutput(KinZfitter::FitInput &input, KinZfitter::FitOutput &output,
                              double &l1, double &l2, double &lph1, double &lph2, 
                              vector<double> &pTerrsREFIT_lep, vector<double> &pTerrsREFIT_gamma,
                              TMatrixDSym &covMatrixZ) {

  l1 = output.pT1_lep/input.pTRECO1_lep;
  l2 = output.pT2_lep/input.pTRECO2_lep;
  pTerrsREFIT_lep.push_back(output.pTErr1_lep);
  pTerrsREFIT_lep.push_back(output.pTErr2_lep);

  if (debug_) {

    cout << "lep1 pt before: " << input.pTRECO1_lep << ", lep1 pt after: " << output.pT1_lep << endl;
    cout << "lep2 pt before: " << input.pTRECO2_lep << ", lep2 pt after: " << output.pT2_lep << endl;

  }

  if (input.nFsr >= 1) {

    lph1 = 1;
    pTerrsREFIT_gamma.push_back(input.pTErr1_gamma);

  }

  if (input.nFsr == 2) {

    lph2 = 1;
    pTerrsREFIT_gamma.push_back(input.pTErr2_gamma);

  }

  int size = output.covMatrixZ.GetNcols();
  covMatrixZ.ResizeTo(size,size);
  covMatrixZ = output.covMatrixZ;

}


void KinZfitter::MakeModel(/*RooWorkspace &w,*/ KinZfitter::FitInput &input, KinZfitter::FitOutput &output) {

  //lep
  RooRealVar pTRECO1_lep("pTRECO1_lep", "pTRECO1_lep", input.pTRECO1_lep, 5, 500);
  RooRealVar pTRECO2_lep("pTRECO2_lep", "pTRECO2_lep", input.pTRECO2_lep, 5, 500);
  RooRealVar pTMean1_lep("pTMean1_lep", "pTMean1_lep", 
			 input.pTRECO1_lep, max(5.0, input.pTRECO1_lep-2*input.pTErr1_lep), input.pTRECO1_lep+2*input.pTErr1_lep);
  RooRealVar pTMean2_lep("pTMean2_lep", "pTMean2_lep", 
			 input.pTRECO2_lep, max(5.0, input.pTRECO2_lep-2*input.pTErr2_lep), input.pTRECO2_lep+2*input.pTErr2_lep);
  RooRealVar pTSigma1_lep("pTSigma1_lep", "pTSigma1_lep", input.pTErr1_lep);
  RooRealVar pTSigma2_lep("pTSigma2_lep", "pTSigma2_lep", input.pTErr2_lep);
  RooRealVar theta1_lep("theta1_lep", "theta1_lep", input.theta1_lep);
  RooRealVar theta2_lep("theta2_lep", "theta2_lep", input.theta2_lep);
  RooRealVar phi1_lep("phi1_lep", "phi1_lep", input.phi1_lep);
  RooRealVar phi2_lep("phi2_lep", "phi2_lep", input.phi2_lep);
  RooRealVar m1("m1", "m1", input.m1);
  RooRealVar m2("m2", "m2", input.m2);

	
  //gamma
  RooRealVar pTRECO1_gamma("pTRECO1_gamma", "pTRECO1_gamma", input.pTRECO1_gamma, 5, 500);
  RooRealVar pTRECO2_gamma("pTRECO2_gamma", "pTRECO2_gamma", input.pTRECO2_gamma, 5, 500);
  RooRealVar pTMean1_gamma("pTMean1_gamma", "pTMean1_gamma", 
			   input.pTRECO1_gamma, -0.5, 0.5);
  RooRealVar pTMean2_gamma("pTMean2_gamma", "pTMean2_gamma", 
			   input.pTRECO2_gamma, -0.5, 0.5);
  RooRealVar pTSigma1_gamma("pTSigma1_gamma", "pTSigma1_gamma", input.pTErr1_gamma);
  RooRealVar pTSigma2_gamma("pTSigma2_gamma", "pTSigma2_gamma", input.pTErr2_gamma);
  RooRealVar theta1_gamma("theta1_gamma", "theta1_gamma", input.theta1_gamma);
  RooRealVar theta2_gamma("theta2_gamma", "theta2_gamma", input.theta2_gamma);
  RooRealVar phi1_gamma("phi1_gamma", "phi1_gamma", input.phi1_gamma);
  RooRealVar phi2_gamma("phi2_gamma", "phi2_gamma", input.phi2_gamma);

  //gauss
  RooGaussian gauss1_lep("gauss1_lep", "gauss1_lep", pTRECO1_lep, pTMean1_lep, pTSigma1_lep);
  RooGaussian gauss2_lep("gauss2_lep", "gauss2_lep", pTRECO2_lep, pTMean2_lep, pTSigma2_lep);
  RooGaussian gauss1_gamma("gauss1_gamma", "gauss1_gamma", pTRECO1_gamma, pTMean1_gamma, pTSigma1_gamma);
  RooGaussian gauss2_gamma("gauss2_gamma", "gauss2_gamma", pTRECO2_gamma, pTMean2_gamma, pTSigma2_gamma);


  TString makeE_lep = "TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1)))+@2*@2)";
  RooFormulaVar E1_lep("E1_lep", makeE_lep, RooArgList(pTMean1_lep, theta1_lep, m1));  //w.import(E1_lep);
  RooFormulaVar E2_lep("E2_lep", makeE_lep, RooArgList(pTMean2_lep, theta2_lep, m2));  //w.import(E2_lep);

  TString makeE_gamma = "TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1))))";
  RooFormulaVar E1_gamma("E1_gamma", makeE_gamma, RooArgList(pTMean1_gamma, theta1_gamma));  //w.import(E1_gamma);
  RooFormulaVar E2_gamma("E2_gamma", makeE_gamma, RooArgList(pTMean2_gamma, theta2_gamma));  //w.import(E2_gamma);

  //dotProduct 3d
  TString dotProduct_3d = "@0*@1*( ((TMath::Cos(@2))*(TMath::Cos(@3)))/((TMath::Sin(@2))*(TMath::Sin(@3)))+(TMath::Cos(@4-@5)))";
  RooFormulaVar p1v3D2("p1v3D2", dotProduct_3d, RooArgList(pTMean1_lep, pTMean2_lep, theta1_lep, theta2_lep, phi1_lep, phi2_lep));
  RooFormulaVar p1v3Dph1("p1v3Dph1", dotProduct_3d, RooArgList(pTMean1_lep, pTMean1_gamma, theta1_lep, theta1_gamma, phi1_lep, phi1_gamma));
  RooFormulaVar p2v3Dph1("p2v3Dph1", dotProduct_3d, RooArgList(pTMean2_lep, pTMean1_gamma, theta2_lep, theta1_gamma, phi2_lep, phi1_gamma));
  RooFormulaVar p1v3Dph2("p1v3Dph2", dotProduct_3d, RooArgList(pTMean1_lep, pTMean2_gamma, theta1_lep, theta2_gamma, phi1_lep, phi2_gamma));
  RooFormulaVar p2v3Dph2("p2v3Dph2", dotProduct_3d, RooArgList(pTMean2_lep, pTMean2_gamma, theta2_lep, theta2_gamma, phi2_lep, phi2_gamma));
  RooFormulaVar ph1v3Dph2("ph1v3Dph2", dotProduct_3d, RooArgList(pTMean1_gamma, pTMean2_gamma, theta1_gamma, theta2_gamma, phi1_gamma, phi2_gamma));

  TString dotProduct_4d = "@0*@1-@2";
  RooFormulaVar p1D2("p1D2", dotProduct_4d, RooArgList(E1_lep, E2_lep, p1v3D2));  //w.import(p1D2);
  RooFormulaVar p1Dph1("p1Dph1", dotProduct_4d, RooArgList(E1_lep, E1_gamma, p1v3Dph1));//  w.import(p1Dph1);
  RooFormulaVar p2Dph1("p2Dph1", dotProduct_4d, RooArgList(E2_lep, E1_gamma, p2v3Dph1)); // w.import(p2Dph1);
  RooFormulaVar p1Dph2("p1Dph2", dotProduct_4d, RooArgList(E1_lep, E2_gamma, p1v3Dph2));  //w.import(p1Dph2);
  RooFormulaVar p2Dph2("p2Dph2", dotProduct_4d, RooArgList(E2_lep, E2_gamma, p2v3Dph2));  //w.import(p2Dph2);
  RooFormulaVar ph1Dph2("ph1Dph2", dotProduct_4d, RooArgList(E1_gamma, E2_gamma, ph1v3Dph2)); // w.import(ph1Dph2);

  RooFormulaVar* mZ;

  mZ = new RooFormulaVar("mZ", "TMath::Sqrt(2*@0+@1*@1+@2*@2)", RooArgList(p1D2, m1, m2));

  //true shape
  RooRealVar meanCB("meanCB","",meanCB_);
  RooRealVar sigmaCB("sigmaCB","",sigmaCB_);
  RooRealVar alphaCB("alphaCB","",alphaCB_);
  RooRealVar nCB("nCB","",nCB_);
  RooRealVar meanGauss1("meanGauss1","",meanGauss1_);
  RooRealVar sigmaGauss1("sigmaGauss1","",sigmaGauss1_);
  RooRealVar f1("f1","",f1_);
  RooRealVar meanGauss2("meanGauss2","",meanGauss2_);
  RooRealVar sigmaGauss2("sigmaGauss2","",sigmaGauss2_);
  RooRealVar f2("f2","",f2_);
  RooRealVar meanGauss3("meanGauss3","",meanGauss3_);
  RooRealVar sigmaGauss3("sigmaGauss3","",sigmaGauss3_);
  RooRealVar f3("f3","",f3_);

  RooCBShape* singleCB = new RooCBShape("singleCB", "", *mZ, meanCB, sigmaCB, alphaCB, nCB);
  RooGaussian* gaussShape1 = new RooGaussian("gaussShape1", "", *mZ, meanGauss1, sigmaGauss1);
  RooAddPdf* CBplusGauss = new RooAddPdf("CBplusGauss", "", *singleCB, *gaussShape1, f1);
  RooGaussian* gaussShape2 = new RooGaussian("gaussShape2", "", *mZ, meanGauss2, sigmaGauss2);
  RooAddPdf* CBplusGaussplusGauss = new RooAddPdf("CBplusGaussplusGauss", "", *CBplusGauss, *gaussShape2, f2);
  RooGaussian* gaussShape3 = new RooGaussian("gaussShape3", "", *mZ, meanGauss3, sigmaGauss3);
  RooAddPdf* CBplusGaussplusGaussplusGauss = new RooAddPdf("CBplusGaussplusGaussplusGauss", "", *CBplusGaussplusGauss, *gaussShape3, f3);

  RooProdPdf *model = new RooProdPdf("model","",RooArgList(gauss1_lep, gauss2_lep,*CBplusGaussplusGaussplusGauss));

  //make fit
  RooArgSet *rastmp;
  rastmp = new RooArgSet(pTRECO1_lep, pTRECO2_lep);

  RooDataSet* pTs = new RooDataSet("pTs","pTs", *rastmp);
  pTs->add(*rastmp);

  RooFitResult* r;
  if (mass4lRECO_ > 140) {
    // r = RelBW->fitTo(*pTs,RooFit::Save(),RooFit::PrintLevel(-1),RooFit::PrintEvalErrors(-1));
    r = model->fitTo(*pTs,RooFit::Save(),RooFit::Constrain(*mZ),RooFit::PrintLevel(-1));

  } else {
    // r = PDFRelBWxCBxgauss->fitTo(*pTs,RooFit::Save(),RooFit::PrintLevel(-1),RooFit::PrintEvalErrors(-1));
    r = model->fitTo(*pTs,RooFit::Save(),RooFit::Constrain(*mZ),RooFit::PrintLevel(-1));

  }
  //save fit result
  const TMatrixDSym& covMatrix = r->covarianceMatrix();
  const RooArgList& finalPars = r->floatParsFinal();

  for (int i=0 ; i<finalPars.getSize(); i++){
 
    TString name = TString(((RooRealVar*)finalPars.at(i))->GetName());
    if(debug_) cout<<"name list of RooRealVar for covariance matrix "<<name<<endl;

  }

  int size = covMatrix.GetNcols();
  output.covMatrixZ.ResizeTo(size,size);
  output.covMatrixZ = covMatrix;
    
  output.pT1_lep = pTMean1_lep.getVal();
  output.pT2_lep = pTMean2_lep.getVal();
  output.pTErr1_lep = pTMean1_lep.getError();
  output.pTErr2_lep = pTMean2_lep.getError();

  if (r) 
    delete r;
  delete rastmp;
  delete pTs;
  delete mZ;
  delete singleCB;
  delete gaussShape1;
  delete CBplusGauss;
  delete gaussShape2;
  delete CBplusGaussplusGauss;
  delete gaussShape3;
  delete CBplusGaussplusGaussplusGauss;
  delete model;

}


#endif

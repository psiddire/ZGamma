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
//#include "DataFormats/Math/interface/deltaR.h"
#include "RooWorkspace.h"
#include "RooProduct.h"
#include "RooProdPdf.h"
#include "RooFormulaVar.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "time.h"
///----------------------------------------------------------------------------------------------
/// KinZfitter::KinZfitter - constructor/
///----------------------------------------------------------------------------------------------

KinZfitter::KinZfitter(bool isData)
{    
  // modified with HZg
  PDFName_ = "HZg_ggF_125GeV_ext1_M125_13TeV_powheg2_pythia8_new";

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
	//old version HZZ
	/*if(p=="sg")  { sgVal_ = val;}
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
	if(p=="f3")   { f3_ = val;}*/
	
	if(p=="alpha")aVal_ = val;                                                                      
	if(p=="mean")meanVal_ = val; 
        if(p=="n" )  nVal_ = val;                                                                       
        if(p=="sigma" )  sgVal_ = val;                                                                 
        if(p=="bwgamma")  bwsigVal_ = val;
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
      pTerr = helperFunc_->pterr(data, lepID[il], fsrID[il], eormu, isData); //put lepton err first
     idsZ1_.push_back(eormu);
      pTerrsZ1_.push_back(pTerr);
      p4sZ1_.push_back(p4);            
      /* else{

	 idsZ2_.push_back(pdgId);
	 pTerrsZ2_.push_back(pTerr);
	 p4sZ2_.push_back(p4);
	 } */
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
	  pTerrsZ1ph_.push_back(pTerr);
	  p4sZ1ph_.push_back(p4);
	  idsFsrZ1_.push_back(idsZ1_[ifsr]);
	}
      if(debug_) cout<<"p4sZ1ph_ "<<p4sZ1ph_.size()<<endl;
    }  

}

void KinZfitter::SetZResult(double l1, double l2, double lph1, double lph2
			    /*,  double l3, double l4, double lph3, double lph4 */)
{

  if(debug_) cout<<"start set Z result"<<endl;

  // pT scale after refitting w.r.t. reco pT

  lZ1_l1_ = l1; lZ1_l2_ = l2;
  /* lZ2_l1_ = l3; lZ2_l2_ = l4; */

  if(debug_) cout<<"l1 "<<l1<<" l2 "<<l2<<endl;
  /* if(debug_) cout<<"l3 "<<l3<<" l4 "<<l4<<endl; */

  lZ1_ph1_ = lph1; lZ1_ph2_ = lph2;
  /* lZ2_ph1_ = lph1; lZ2_ph2_ = lph2; */

  TLorentzVector Z1_1 = p4sZ1_[0]; TLorentzVector Z1_2 = p4sZ1_[1];
  /* TLorentzVector Z2_1 = p4sZ2_[0]; TLorentzVector Z2_2 = p4sZ2_[1]; */

  TLorentzVector Z1_1_True(0,0,0,0);
  Z1_1_True.SetPtEtaPhiM(lZ1_l1_*Z1_1.Pt(),Z1_1.Eta(),Z1_1.Phi(),Z1_1.M());
  TLorentzVector Z1_2_True(0,0,0,0);
  Z1_2_True.SetPtEtaPhiM(lZ1_l2_*Z1_2.Pt(),Z1_2.Eta(),Z1_2.Phi(),Z1_2.M());
  /*   TLorentzVector Z2_1_True(0,0,0,0);
       Z2_1_True.SetPtEtaPhiM(lZ2_l1_*Z2_1.Pt(),Z2_1.Eta(),Z2_1.Phi(),Z2_1.M());
       TLorentzVector Z2_2_True(0,0,0,0);
       Z2_2_True.SetPtEtaPhiM(lZ2_l2_*Z2_2.Pt(),Z2_2.Eta(),Z2_2.Phi(),Z2_2.M()); */

  p4sZ1REFIT_.push_back(Z1_1_True); p4sZ1REFIT_.push_back(Z1_2_True);

  /* p4sZ2REFIT_.push_back(Z2_1_True); p4sZ2REFIT_.push_back(Z2_2_True); */

  for(unsigned int ifsr1 = 0; ifsr1 < p4sZ1ph_.size(); ifsr1++){

    TLorentzVector Z1ph = p4sZ1ph_[ifsr1];
    TLorentzVector Z1phTrue(0,0,0,0);
  
    double l = 1.0;
    if(ifsr1==0) l = lZ1_ph1_; if(ifsr1==1) l = lZ1_ph2_;
  
    Z1phTrue.SetPtEtaPhiM(l*Z1ph.Pt(),Z1ph.Eta(),Z1ph.Phi(),Z1ph.M());
  
    p4sZ1phREFIT_.push_back(Z1phTrue);
    pTerrsZ1phREFIT_.push_back(pTerrsZ1ph_[ifsr1]);

  }

  // since it is Z1 refit result, Z2 kinematics keep at what it is 
  //  p4sZ2REFIT_.push_back(p4sZ2_[0]); p4sZ2REFIT_.push_back(p4sZ2_[1]);
  //  pTerrsZ2REFIT_.push_back(pTerrsZ2_[0]); pTerrsZ2REFIT_.push_back(pTerrsZ2_[1]);

  /* for(unsigned int ifsr2 = 0; ifsr2 < p4sZ2ph_.size(); ifsr2++){

     TLorentzVector Z2ph = p4sZ2ph_[ifsr2];
     TLorentzVector Z2phTrue(0,0,0,0);

     double l = 1.0;
     if(ifsr2==0) l = lZ2_ph1_; if(ifsr2==1) l = lZ2_ph2_;

     Z2phTrue.SetPtEtaPhiM(l*Z2ph.Pt(),Z2ph.Eta(),Z2ph.Phi(),Z2ph.M());

     p4sZ2phREFIT_.push_back(Z2phTrue);
     pTerrsZ2phREFIT_.push_back(pTerrsZ2ph_[ifsr2]);
     } */

  if(debug_) cout<<"end set Z1 result"<<endl;
}

double KinZfitter::GetM4l()
{

  vector<TLorentzVector> p4s = GetP4s();

    TLorentzVector pH(0,0,0,0);
  for(unsigned int i = 0; i< p4s.size(); i++){
    pH = pH + p4s[i];
  }

  return pH.M();

}


double KinZfitter::GetRefitM4l()//removed for HZg
{

  vector<TLorentzVector> p4s = GetRefitP4s();

  TLorentzVector pH(0,0,0,0); 
  for(unsigned int i = 0; i< p4s.size(); i++){
    pH = pH + p4s[i];
  }

  return pH.M();

}

double KinZfitter::GetRefitMZ1()
{

  vector<TLorentzVector> p4s = GetRefitP4s();

  TLorentzVector pZ1(0,0,0,0);

  pZ1 = p4s[0] + p4s[1];

  return pZ1.M();

}
/* double KinZfitter::GetRefitMZ2()// removed for Zg
   {

   vector<TLorentzVector> p4s = GetRefitP4s();

   TLorentzVector pZ2(0,0,0,0);

   pZ2 = p4s[2] + p4s[3];

   return pZ2.M();

   } */

double KinZfitter::GetMZ1()
{

  vector<TLorentzVector> p4s = GetP4s();

  TLorentzVector pZ1(0,0,0,0);

  pZ1 = p4s[0] + p4s[1];

  return pZ1.M();

}
/* double KinZfitter::GetMZ2()//removed for HZg
   {

   vector<TLorentzVector> p4s = GetP4s();

   TLorentzVector pZ2(0,0,0,0);

   pZ2 = p4s[2] + p4s[3];

   return pZ2.M();

   } */


double KinZfitter::GetRefitM4lErr()//removed
{

  vector<TLorentzVector> p4s;
  vector<double> pTErrs;

  p4s.push_back(p4sZ1REFIT_[0]);p4s.push_back(p4sZ1REFIT_[1]);
  // p4s.push_back(p4sZ2REFIT_[0]);p4s.push_back(p4sZ2REFIT_[1]);

  // patch when MINIUT FAILS
  if(pTerrsZ1REFIT_[0]==0||pTerrsZ1REFIT_[1]==0)
    return GetM4lErr();

  pTErrs.push_back(pTerrsZ1REFIT_[0]); pTErrs.push_back(pTerrsZ1REFIT_[1]);
  // pTErrs.push_back(pTerrsZ2REFIT_[0]); pTErrs.push_back(pTerrsZ2REFIT_[1]);

  for(unsigned int ifsr1 = 0; ifsr1<p4sZ1phREFIT_.size(); ifsr1++){

    p4s.push_back(p4sZ1phREFIT_[ifsr1]);
    pTErrs.push_back(pTerrsZ1phREFIT_[ifsr1]);

  }

  // for(unsigned int ifsr2 = 0; ifsr2<p4sZ2phREFIT_.size(); ifsr2++){

  // p4s.push_back(p4sZ2phREFIT_[ifsr2]);
  // pTErrs.push_back(pTerrsZ2phREFIT_[ifsr2]);

  // }
  return helperFunc_->masserror(p4s,pTErrs);
  // return 0.1;
}
/*
  double KinZfitter::GetRefitM4lErrFullCov()
  {


  vector<TLorentzVector> Lp4s = GetRefitP4s();
  vector<TLorentzVector> p4s;
  vector<double> pTErrs;

  p4s.push_back(p4sZ1REFIT_[0]);p4s.push_back(p4sZ1REFIT_[1]);
  pTErrs.push_back(pTerrsZ1REFIT_[0]); pTErrs.push_back(pTerrsZ1REFIT_[1]);
  // patch when MINUIT FAILS
  if(pTerrsZ1REFIT_[0]==0||pTerrsZ1REFIT_[1]==0)

  return GetM4lErr();


  if(p4sZ1phREFIT_.size()>=1){
  p4s.push_back(p4sZ1phREFIT_[0]); pTErrs.push_back(pTerrsZ1phREFIT_[0]);
  }
  if(p4sZ1phREFIT_.size()==2){
  p4s.push_back(p4sZ1phREFIT_[1]); pTErrs.push_back(pTerrsZ1phREFIT_[1]);
  }

  p4s.push_back(p4sZ2REFIT_[0]);p4s.push_back(p4sZ2REFIT_[1]);
  pTErrs.push_back(pTerrsZ2REFIT_[0]); pTErrs.push_back(pTerrsZ2REFIT_[1]);

  if(p4sZ2phREFIT_.size()>=1){
  p4s.push_back(p4sZ2phREFIT_[0]); pTErrs.push_back(pTerrsZ2phREFIT_[0]);
  }
  if(p4sZ2phREFIT_.size()==2){
  p4s.push_back(p4sZ2phREFIT_[1]); pTErrs.push_back(pTerrsZ2phREFIT_[1]);
  }


  double errorUncorr = helperFunc_->masserror(p4s,pTErrs);

  vector<double> pTErrs1; vector<double> pTErrs2;
  vector<double> pTErrs3; vector<double> pTErrs4; 

  for(unsigned int i = 0; i<pTErrs.size(); i++){
  if(i==0) pTErrs1.push_back(pTErrs[i]);
  else pTErrs1.push_back(0.0);
  }
  for(unsigned int i = 0; i<pTErrs.size(); i++){
  if(i==1) pTErrs2.push_back(pTErrs[i]);
  else pTErrs2.push_back(0.0);
  }
  for(unsigned int i = 0; i<pTErrs.size(); i++){
  if(i==2) pTErrs3.push_back(pTErrs[i]);
  else pTErrs3.push_back(0.0);
  }
  for(unsigned int i = 0; i<pTErrs.size(); i++){
  if(i==3) pTErrs4.push_back(pTErrs[i]);
  else pTErrs4.push_back(0.0);
  } 

  double error1 = helperFunc_->masserror(p4s,pTErrs1);
  double error2 = helperFunc_->masserror(p4s,pTErrs2);
  double error3 = helperFunc_->masserror(p4s,pTErrs3);
  double error4 = helperFunc_->masserror(p4s,pTErrs4); 

  
  double errorph1 = 0.0; double errorph2 = 0.0;
  if(p4sZ2phREFIT_.size()>=1){

  vector<double> pTErrsph1;
  for(unsigned int i = 0; i<pTErrs.size(); i++){
  if(i==2) pTErrsph1.push_back(pTErrs[i]);
  else pTErrsph1.push_back(0.0);
  }

  errorph1 = helperFunc_->masserror(p4s,pTErrsph1);

  }

  if(p4sZ2phREFIT_.size()>=2){

  vector<double> pTErrsph2;
  for(unsigned int i = 0; i<pTErrs.size(); i++){
  if(i==3) pTErrsph2.push_back(pTErrs[i]);
  else pTErrsph2.push_back(0.0);
  }

  errorph2 = helperFunc_->masserror(p4s,pTErrsph2);

  }
  

  if(debug_) cout<<"error1 "<<error1<<" error2 "<<error2<<endl;
  // if(debug_) cout<<"error3 "<<error3<<" error4 "<<error4<<endl;

  // covariance matrix
  double delta12 = error1*error2*covMatrixZ1_(0,1)/sqrt(covMatrixZ1_(0,0)*covMatrixZ1_(1,1));
  // double delta34=0.0;
  if (mass4lRECO_ > cutoff_) {
  delta34 = error3*error4*covMatrixZ2_(0,1)/sqrt(covMatrixZ2_(0,0)*covMatrixZ2_(1,1));
  }
  if (debug_) cout<<"delta34 "<<delta34<<endl; 

  double delta1ph1 = 0.0; double delta1ph2 = 0.0;
  double delta2ph1 = 0.0; double delta2ph2 = 0.0;
  double deltaph1ph2 = 0.0;
  if(p4sZ1phREFIT_.size()>=1){

  delta1ph1 = error1*errorph1*covMatrixZ1_(0,2)/sqrt(covMatrixZ1_(0,0)*covMatrixZ1_(2,2));
  delta2ph1 = error2*errorph1*covMatrixZ1_(1,2)/sqrt(covMatrixZ1_(1,1)*covMatrixZ1_(2,2));
  }

  if(p4sZ1phREFIT_.size()>=2){
  delta1ph2 = error1*errorph2*covMatrixZ1_(0,3)/sqrt(covMatrixZ1_(0,0)*covMatrixZ1_(3,3));
  delta2ph2 = error2*errorph2*covMatrixZ1_(1,3)/sqrt(covMatrixZ1_(1,1)*covMatrixZ1_(3,3));
  delta1ph2 = errorph1*errorph2*covMatrixZ1_(2,3)/sqrt(covMatrixZ1_(2,2)*covMatrixZ1_(3,3));
  }
  

  //double correlation = delta12;//+delta1ph1+delta1ph2+delta2ph1+delta2ph2+deltaph1ph2;
  double correlation = delta12 +delta34 ;
  double err = sqrt(errorUncorr*errorUncorr+correlation);

  return err;
  
  }*/

double KinZfitter::GetM4lErr()
{
  
  vector<TLorentzVector> p4s;
  vector<double> pTErrs;
  
  p4s.push_back(p4sZ1_[0]);p4s.push_back(p4sZ1_[1]);
  // p4s.push_back(p4sZ2_[0]);p4s.push_back(p4sZ2_[1]);
  
  pTErrs.push_back(pTerrsZ1_[0]); pTErrs.push_back(pTerrsZ1_[1]);
  // pTErrs.push_back(pTerrsZ2_[0]); pTErrs.push_back(pTerrsZ2_[1]);
  
  for(unsigned int ifsr1 = 0; ifsr1<p4sZ1ph_.size(); ifsr1++){
      
    p4s.push_back(p4sZ1ph_[ifsr1]);
    pTErrs.push_back(pTerrsZ1ph_[ifsr1]);
  
  }
  
  // for(unsigned int ifsr2 = 0; ifsr2<p4sZ2ph_.size(); ifsr2++){
      
  // p4s.push_back(p4sZ2ph_[ifsr2]);
  // pTErrs.push_back(pTerrsZ2ph_[ifsr2]);
  
  // }
  
  // return helperFunc_->masserror(p4s,pTErrs);
  return 0.5;
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
  // TLorentzVector Z2_1 = p4sZ2REFIT_[0]; TLorentzVector Z2_2 = p4sZ2REFIT_[1];

  /// fsr photons

  for(unsigned int ifsr1 = 0; ifsr1<p4sZ1phREFIT_.size(); ifsr1++){

    int id_fsr1 = idsFsrZ1_[ifsr1];
    TLorentzVector Z1ph = p4sZ1phREFIT_[ifsr1];

    if(id_fsr1==idsZ1_[0]) Z1_1 = Z1_1 + Z1ph;
    if(id_fsr1==idsZ1_[1]) Z1_2 = Z1_2 + Z1ph;

  }

  /*    for(unsigned int ifsr2 = 0; ifsr2<p4sZ2phREFIT_.size(); ifsr2++){

        int id_fsr2 = idsFsrZ2_[ifsr2];
        TLorentzVector Z2ph = p4sZ2phREFIT_[ifsr2];

        if(id_fsr2==idsZ2_[0]) Z2_1 = Z2_1 + Z2ph;
        if(id_fsr2==idsZ2_[1]) Z2_2 = Z2_2 + Z2ph;

	}
  */
  vector<TLorentzVector> p4s;
  p4s.push_back(Z1_1); p4s.push_back(Z1_2);
  /* p4s.push_back(Z2_1); p4s.push_back(Z2_2); */

  return p4s;

}

vector<TLorentzVector> KinZfitter::GetP4s()
{

  TLorentzVector Z1_1 = p4sZ1_[0]; TLorentzVector Z1_2 = p4sZ1_[1];

  /// fsr photons

  for(unsigned int ifsr1 = 0; ifsr1<p4sZ1ph_.size(); ifsr1++){

    int id_fsr1 = idsFsrZ1_[ifsr1];
    TLorentzVector Z1ph = p4sZ1ph_[ifsr1];
    if(id_fsr1==idsZ1_[0]) Z1_1 = Z1_1 + Z1ph;
    if(id_fsr1==idsZ1_[1]) Z1_2 = Z1_2 + Z1ph;

  }

  /*  for(unsigned int ifsr2 = 0; ifsr2<p4sZ2ph_.size(); ifsr2++){

      int id_fsr2 = idsFsrZ2_[ifsr2];
      TLorentzVector Z2ph = p4sZ2ph_[ifsr2];

      if(id_fsr2==idsZ2_[0]) Z2_1 = Z2_1 + Z2ph;
      if(id_fsr2==idsZ2_[1]) Z2_2 = Z2_2 + Z2ph;

      } */

  vector<TLorentzVector> p4s;
  p4s.push_back(Z1_1); 
  p4s.push_back(Z1_2);
  /* p4s.push_back(Z2_1); 
     p4s.push_back(Z2_2); */
  return p4s;

}

void KinZfitter::KinRefitZ()
{
  double l1,l2,lph1,lph2;
  /* double l3,l4,lph3,lph4; */

  l1 = 1.0; l2 = 1.0; lph1 = 1.0; lph2 = 1.0;
  /* l3 = 1.0; l4 = 1.0; lph3 = 1.0; lph4 = 1.0; */

  // bool fourEfourMu = IsFourEFourMu(idsZ1_, idsZ2_);

  // mass4lRECO_ = GetM4l();

  //cout << mass4lRECO << endl;

  // if (mass4lRECO_ <= cutoff_) {//fit Z1

  SetFitInput(fitInput1, p4sZ1_, pTerrsZ1_, p4sZ1ph_, pTerrsZ1ph_);
  if (debug_) cout<<"done input"<<endl;
  Driver(fitInput1, fitOutput1);
  if (debug_) cout<<"done driver"<<endl;
  SetFitOutput(fitInput1, fitOutput1, l1, l2, lph1, lph2, pTerrsZ1REFIT_, pTerrsZ1phREFIT_, covMatrixZ1_);

  // pTerrsZ2REFIT_.push_back(pTerrsZ2_[0]); pTerrsZ2REFIT_.push_back(pTerrsZ2_[1]);

  // } else {//fit two Zs

  // if (fourEfourMu) {//4e,4mu, do reshuffle

  // RepairZ1Z2(p4sZ1_, pTerrsZ1_, p4sZ1ph_, pTerrsZ1ph_, p4sZ2_, pTerrsZ2_, p4sZ2ph_, pTerrsZ2ph_, idsZ1_, idsZ2_);

  // }
       
  // SetFitInput(fitInput1, p4sZ1_, pTerrsZ1_, p4sZ1ph_, pTerrsZ1ph_);
  // Driver(fitInput1, fitOutput1);
  // SetFitOutput(fitInput1, fitOutput1, l1, l2, lph1, lph2, pTerrsZ1REFIT_, pTerrsZ1phREFIT_, covMatrixZ1_);

  /*   SetFitInput(fitInput2, p4sZ2_, pTerrsZ2_, p4sZ2ph_, pTerrsZ2ph_);
       Driver(fitInput2, fitOutput2);
       SetFitOutput(fitInput2, fitOutput2, l3, l4, lph3, lph4, pTerrsZ2REFIT_, pTerrsZ2phREFIT_, covMatrixZ2_);
  */
  // }

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

	if (debug_) cout << "nFsr: " << input.nFsr << endl;
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
 //gauss
  RooGaussian gauss1_lep("gauss1_lep", "gauss1_lep", pTRECO1_lep, pTMean1_lep, pTSigma1_lep);
  RooGaussian gauss2_lep("gauss2_lep", "gauss2_lep", pTRECO2_lep, pTMean2_lep, pTSigma2_lep);

  TString makeE_lep = "TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1)))+@2*@2)";
  RooFormulaVar E1_lep("E1_lep", makeE_lep, RooArgList(pTMean1_lep, theta1_lep, m1));  //w.import(E1_lep);
  RooFormulaVar E2_lep("E2_lep", makeE_lep, RooArgList(pTMean2_lep, theta2_lep, m2));  //w.import(E2_lep);
  //dotProduct 3d
  TString dotProduct_3d = "@0*@1*( ((TMath::Cos(@2))*(TMath::Cos(@3)))/((TMath::Sin(@2))*(TMath::Sin(@3)))+(TMath::Cos(@4-@5)))";
  RooFormulaVar p1v3D2("p1v3D2", dotProduct_3d, RooArgList(pTMean1_lep, pTMean2_lep, theta1_lep, theta2_lep, phi1_lep, phi2_lep));
  TString dotProduct_4d = "@0*@1-@2";
  RooFormulaVar p1D2("p1D2", dotProduct_4d, RooArgList(E1_lep, E2_lep, p1v3D2));  //w.import(p1D2);
  
  //gamma ---treat seperately, if no photon--> don't use it
  TString makeE_gamma = "TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1))))";
  if(input.pTErr1_gamma==0.)input.pTErr1_gamma = 0.001;
  if(input.pTErr2_gamma==0.)input.pTErr2_gamma = 0.001;
  
	RooRealVar pTRECO1_gamma("pTRECO1_gamma", "pTRECO1_gamma", input.pTRECO1_gamma, 5, 500);
	RooRealVar pTMean1_gamma("pTMean1_gamma", "pTMean1_gamma", 
	input.pTRECO1_gamma, max(0.5, input.pTRECO1_gamma-2*input.pTErr1_gamma), max(0.5+0.001,input.pTRECO1_gamma+2*input.pTErr1_gamma)); 
	RooRealVar pTSigma1_gamma("pTSigma1_gamma", "pTSigma1_gamma", input.pTErr1_gamma);
	RooRealVar theta1_gamma("theta1_gamma", "theta1_gamma", input.theta1_gamma);
	RooRealVar phi1_gamma("phi1_gamma", "phi1_gamma", input.phi1_gamma);
	RooGaussian gauss1_gamma("gauss1_gamma", "gauss1_gamma", pTRECO1_gamma, pTMean1_gamma, pTSigma1_gamma);
	RooFormulaVar E1_gamma("E1_gamma", makeE_gamma, RooArgList(pTMean1_gamma, theta1_gamma));  //w.import(E1_gamma);
	RooFormulaVar p1v3Dph1("p1v3Dph1", dotProduct_3d, RooArgList(pTMean1_lep, pTMean1_gamma, theta1_lep, theta1_gamma, phi1_lep, phi1_gamma));
	RooFormulaVar p2v3Dph1("p2v3Dph1", dotProduct_3d, RooArgList(pTMean2_lep, pTMean1_gamma, theta2_lep, theta1_gamma, phi2_lep, phi1_gamma));
	RooFormulaVar p1Dph1("p1Dph1", dotProduct_4d, RooArgList(E1_lep, E1_gamma, p1v3Dph1));//  w.import(p1Dph1);
	RooFormulaVar p2Dph1("p2Dph1", dotProduct_4d, RooArgList(E2_lep, E1_gamma, p2v3Dph1)); // w.import(p2Dph1);
	RooRealVar pTRECO2_gamma("pTRECO2_gamma", "pTRECO2_gamma", input.pTRECO2_gamma, 5, 500);
	RooRealVar pTMean2_gamma("pTMean2_gamma", "pTMean2_gamma", 
	input.pTRECO2_gamma, max(0.5, input.pTRECO2_gamma-2*input.pTErr2_gamma), max(0.5+0.001,input.pTRECO2_gamma+2*input.pTErr2_gamma)); 
	RooRealVar pTSigma2_gamma("pTSigma2_gamma", "pTSigma2_gamma", input.pTErr2_gamma);
	RooRealVar theta2_gamma("theta2_gamma", "theta2_gamma", input.theta2_gamma);
	RooRealVar phi2_gamma("phi2_gamma", "phi2_gamma", input.phi2_gamma);
	RooGaussian gauss2_gamma("gauss2_gamma", "gauss2_gamma", pTRECO2_gamma, pTMean2_gamma, pTSigma2_gamma);
	RooFormulaVar E2_gamma("E2_gamma", makeE_gamma, RooArgList(pTMean2_gamma, theta2_gamma));  //w.import(E2_gamma);
	RooFormulaVar p1v3Dph2("p1v3Dph2", dotProduct_3d, RooArgList(pTMean1_lep, pTMean2_gamma, theta1_lep, theta2_gamma, phi1_lep, phi2_gamma));
	RooFormulaVar p2v3Dph2("p2v3Dph2", dotProduct_3d, RooArgList(pTMean2_lep, pTMean2_gamma, theta2_lep, theta2_gamma, phi2_lep, phi2_gamma));
	RooFormulaVar p1Dph2("p1Dph2", dotProduct_4d, RooArgList(E1_lep, E2_gamma, p1v3Dph2));  //w.import(p1Dph2);
	RooFormulaVar p2Dph2("p2Dph2", dotProduct_4d, RooArgList(E2_lep, E2_gamma, p2v3Dph2));  //w.import(p2Dph2);
	RooFormulaVar ph1v3Dph2("ph1v3Dph2", dotProduct_3d, RooArgList(pTMean1_gamma, pTMean2_gamma, theta1_gamma, theta2_gamma, phi1_gamma, phi2_gamma));
	RooFormulaVar ph1Dph2("ph1Dph2", dotProduct_4d, RooArgList(E1_gamma, E2_gamma, ph1v3Dph2)); // w.import(ph1Dph2);


  RooFormulaVar* mZ;

  //mZ
  if (input.nFsr == 1) {
    mZ = new RooFormulaVar("mZ", "TMath::Sqrt(2*@0+2*@1+2*@2+@3*@3+@4*@4)", RooArgList(p1D2, p1Dph1, p2Dph1, m1, m2));
    } 
    else if (input.nFsr == 2) {
    mZ = new RooFormulaVar("mZ", "TMath::Sqrt(2*@0+2*@1+2*@2+2*@3+2*@4+2*@5+@6*@6+@7*@7)", RooArgList(p1D2,p1Dph1,p2Dph1,p1Dph2,p2Dph2,ph1Dph2, m1, m2));
    }
    else {
  mZ = new RooFormulaVar("mZ", "TMath::Sqrt(2*@0+@1*@1+@2*@2)", RooArgList(p1D2, m1, m2));
  }

  //true shape
  /*RooRealVar meanCB("meanCB","",meanCB_);
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
  
  RooProdPdf *model = new RooProdPdf("model","",RooArgList(gauss1_lep, gauss2_lep,*CBplusGaussplusGaussplusGauss));*/
  RooAbsReal *Z =mZ;
  RooRealVar *Zmass = new RooRealVar("Zmass","Zmass",Z->getValV(),-1000,10000);
  RooRealVar bwMean("bwMean", "m_{Z^{0}}", 91.187);
  RooRealVar bwGamma("bwGamma", "#Gamma", bwsigVal_);
  RooRealVar sg("sg", "sg", sgVal_);
  RooRealVar a("a", "a", aVal_);
  RooRealVar n("n", "n", nVal_);
  RooRealVar mean("mean","mean",meanVal_);
  RooCBShape CB("CB","CB",*mZ,mean,sg,a,n);
  RooGenericPdf RelBW("RelBW","1/( pow(Zmass*Zmass-bwMean*bwMean,2)+pow(Zmass,4)*pow(bwGamma/bwMean,2) )",RooArgSet(*mZ,bwMean,bwGamma) );
  RooFFTConvPdf *RelBWxCB = new RooFFTConvPdf("RelBWxCB","RelBWxCB", *Z, RelBW,CB);
  RooProdPdf *model = new RooProdPdf("model","",RooArgList(gauss1_lep, gauss2_lep,*RelBWxCB));


  //make fit
  RooArgSet *rastmp;
  rastmp = new RooArgSet(pTRECO1_lep, pTRECO2_lep);
  
    if(input.nFsr == 1) {
    rastmp = new RooArgSet(pTRECO1_lep, pTRECO2_lep, pTRECO1_gamma);
    }

    if(input.nFsr == 2) {
    rastmp = new RooArgSet(pTRECO1_lep, pTRECO2_lep, pTRECO1_gamma, pTRECO2_gamma);
    }
    RooDataSet* pTs = new RooDataSet("pTs","pTs", *rastmp);
  pTs->add(*rastmp);
  //pTs->addColumn(*Zmass);
  //RooRealVar *plotmZ = (RooRealVar*)pTs->addColumn(*Z);
  RooRealVar *plotmZ = Zmass;
  RooPlot *xframe = plotmZ->frame(50.,120.);
   //model->plotOn(xframe,RooFit::Name("model"),LineColor(kBlue));
   RelBWxCB->plotOn(xframe,RooFit::Name("RelBWxCB"));
  /*CBplusGaussplusGaussplusGauss->plotOn(xframe,RooFit::Name("CBplusGaussplusGaussplusGauss"),LineStyle(kDashed),LineColor(kRed));*/
   xframe->Draw(); 
   ///gPad->Print(Form("new_%.2f.pdf",input.pTRECO2_lep));
  //RooRealVar *plotmZ = (RooRealVar*)pTs->addColumn(*mZ);
   RooFitResult* r;
   RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
   RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
   RooMsgService::instance().setSilentMode(true);
   
   r = model->fitTo(*pTs,RooFit::Save(),RooFit::Constrain(*mZ),RooFit::PrintLevel(-1));

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
    if (input.nFsr >= 1) 
	{
		output.pT1_gamma = pTMean1_gamma.getVal();
		output.pTErr1_gamma = pTMean1_gamma.getError();
    }

    if (input.nFsr == 2) 
	{
		output.pT2_gamma = pTMean2_gamma.getVal();
		output.pTErr2_gamma = pTMean2_gamma.getError();
    }

  if (r) delete r;
  //delete Z;
  delete Zmass; 
  delete rastmp;
  delete pTs;
  delete mZ;
  /*delete singleCB;
  delete gaussShape1;
  delete CBplusGauss;
  delete gaussShape2;
  delete CBplusGaussplusGauss;
  delete gaussShape3;
  delete CBplusGaussplusGaussplusGauss;*/
  delete RelBWxCB;
  delete model;  
}

/*bool KinZfitter::IsFourEFourMu(vector<int> &Z1id, vector<int> &Z2id) {
bool flag = false;
if (abs(Z1id[0]) == abs(Z2id[0])) flag = true;
return flag;
}*/
/* 
   void  KinZfitter::RepairZ1Z2(vector<TLorentzVector> &Z1Lep, vector<double> &Z1LepErr,
   vector<TLorentzVector> &Z1Gamma, vector<double> &Z1GammaErr,
   vector<TLorentzVector> &Z2Lep, vector<double> &Z2LepErr,
   vector<TLorentzVector> &Z2Gamma, vector<double> &Z2GammaErr,
   vector<int> &Z1id, vector<int> &Z2id) {

   typedef pair<int, TLorentzVector> Lep;
   typedef pair<Lep, Lep> Z;
   typedef pair<double, double> ZLepErr;

   Lep lep1, lep2, lep3, lep4;
   lep1 = make_pair(Z1id[0], Z1Lep[0]);
   lep2 = make_pair(Z1id[1], Z1Lep[1]);
   lep3 = make_pair(Z2id[0], Z2Lep[0]);
   lep4 = make_pair(Z2id[1], Z2Lep[1]);

   Z Z1_cfg1, Z2_cfg1, Z1_cfg2, Z2_cfg2;
   Z1_cfg1 = make_pair(lep1, lep2); 
   Z2_cfg1 = make_pair(lep3, lep4);

   ZLepErr Z1LepErr_cfg1, Z2LepErr_cfg1, Z1LepErr_cfg2, Z2LepErr_cfg2;
   Z1LepErr_cfg1 = make_pair(Z1LepErr[0], Z1LepErr[1]);
   Z2LepErr_cfg1 = make_pair(Z2LepErr[0], Z2LepErr[1]);

   //      Z1_cfg2 = make_pair(lep1, (lep1.first + lep3.first == 0) ? lep3 : lep4);
   //      Z2_cfg2 = make_pair(lep2, (lep2.first + lep4.first == 0) ? lep4 : lep3);
   if (lep1.first + lep3.first == 0) {

   Z1_cfg2 = make_pair(lep1, lep3);
   Z1LepErr_cfg2 = make_pair(Z1LepErr[0], Z2LepErr[0]);
   Z2_cfg2 = make_pair(lep2, lep4);
   Z2LepErr_cfg2 = make_pair(Z1LepErr[1], Z2LepErr[1]);

   } else { 

   Z1_cfg2 = make_pair(lep1, lep4);
   Z1LepErr_cfg2 = make_pair(Z1LepErr[0], Z2LepErr[1]);
   Z2_cfg2 = make_pair(lep2, lep3);
   Z2LepErr_cfg2 = make_pair(Z1LepErr[1], Z2LepErr[0]);

   }

   //       ZLepErr Z1LepErr_cfg1, Z2LepErr_cfg1, Z1LepErr_cfg2, Z2LepErr_cfg2;
   // Z1LepErr_cfg1 = make_pair(Z1LepErr[0], Z1LepErr[1]);
   // Z2LepErr_cfg1 = make_pair(Z2LepErr[0], Z2LepErr[1]);

   // Z1LepErr_cfg2 = make_pair(Z1LepErr[0], (lep1.first + lep3.first == 0) ? Z2LepErr[0] : Z2LepErr[1]);
   // Z2LepErr_cfg2 = make_pair(Z1LepErr[1], (lep2.first + lep4.first == 0) ? Z2LepErr[1] : Z2LepErr[0]);

   double massZ1_cfg1 = (Z1_cfg1.first.second + Z1_cfg1.second.second).M();
   double massZ2_cfg1 = (Z2_cfg1.first.second + Z2_cfg1.second.second).M();
   double massZ1_cfg2 = (Z1_cfg2.first.second + Z1_cfg2.second.second).M();
   double massZ2_cfg2 = (Z2_cfg2.first.second + Z2_cfg2.second.second).M();

   //      double massZDiff_cfg1 = abs(massZ1_cfg1-massZ2_cfg1);
   //      double massZDiff_cfg2 = abs(massZ1_cfg2-massZ2_cfg2);
   double massZDiff_cfg1 = abs(massZ1_cfg1-91.2) + abs(massZ2_cfg1-91.2);
   double massZDiff_cfg2 = abs(massZ1_cfg2-91.2) + abs(massZ2_cfg2-91.2);

   if (debug_) cout << "massZdiff_cfg1: " << massZDiff_cfg1 << ", massZdiff_cfg2: " << massZDiff_cfg2 << endl;
   if (debug_) cout << "Z1lep2Pt_cfg2: "  << Z1_cfg2.second.second.Pt() << ", Z2lep2Pt_cfg2: " << Z2_cfg2.second.second.Pt() << endl;

   if (massZDiff_cfg1 > massZDiff_cfg2) {

   Z1id[0] = Z1_cfg2.first.first;
   Z1Lep[0] = Z1_cfg2.first.second; 
   Z1LepErr[0] = Z1LepErr_cfg2.first; 

   Z1id[1] = Z1_cfg2.second.first;
   Z1Lep[1] = Z1_cfg2.second.second; 
   Z1LepErr[1] = Z1LepErr_cfg2.second;

   Z2id[0] = Z2_cfg2.first.first;
   Z2Lep[0] = Z2_cfg2.first.second; 
   Z2LepErr[0] = Z2LepErr_cfg2.first;

   Z2id[1] = Z2_cfg2.second.first;
   Z2Lep[1] = Z2_cfg2.second.second;
   Z2LepErr[1] = Z2LepErr_cfg2.second;

   }
   }
*/
//old version
/* 
   int KinZfitter::PerZ1Likelihood(double & l1, double & l2, double & lph1, double & lph2)
   {

   l1= 1.0; l2 = 1.0;
   lph1 = 1.0; lph2 = 1.0;

   if(debug_) cout<<"start Z refit"<<endl;

   TLorentzVector Z1_1 = p4sZ1_[0]; TLorentzVector Z1_2 = p4sZ1_[1];

   double RECOpT1 = Z1_1.Pt(); double RECOpT2 = Z1_2.Pt();
   double pTerrZ1_1 = pTerrsZ1_[0]; double pTerrZ1_2 = pTerrsZ1_[1];

   if(debug_)cout<<"pT1 "<<RECOpT1<<" pTerrZ1_1 "<<pTerrZ1_1<<endl;
   if(debug_)cout<<"pT2 "<<RECOpT2<<" pTerrZ1_2 "<<pTerrZ1_2<<endl;

   //////////////

   TLorentzVector Z1_ph1, Z1_ph2;
   double pTerrZ1_ph1, pTerrZ1_ph2;
   double RECOpTph1, RECOpTph2;

   TLorentzVector nullFourVector(0, 0, 0, 0);
   Z1_ph1=nullFourVector; Z1_ph2=nullFourVector;
   RECOpTph1 = 0; RECOpTph2 = 0;
   pTerrZ1_ph1 = 0; pTerrZ1_ph2 = 0;

   if(p4sZ1ph_.size()>=1){

   Z1_ph1 = p4sZ1ph_[0]; pTerrZ1_ph1 = pTerrsZ1ph_[0];
   RECOpTph1 = Z1_ph1.Pt();
   if(debug_) cout<<"put in Z1 fsr photon 1 pT "<<RECOpTph1<<" pT err "<<pTerrZ1_ph1<<endl; 
   }
   if(p4sZ1ph_.size()==2){
   //if(debug_) cout<<"put in Z1 fsr photon 2"<<endl;
   Z1_ph2 = p4sZ1ph_[1]; pTerrZ1_ph2 = pTerrsZ1ph_[1];
   RECOpTph2 = Z1_ph2.Pt();     
   }

   RooRealVar* pT1RECO = new RooRealVar("pT1RECO","pT1RECO", RECOpT1, 5, 500);
   RooRealVar* pT2RECO = new RooRealVar("pT2RECO","pT2RECO", RECOpT2, 5, 500);
   
   double RECOpT1min = max(5.0, RECOpT1-2*pTerrZ1_1);
   double RECOpT2min = max(5.0, RECOpT2-2*pTerrZ1_2);

   RooRealVar* pTph1RECO = new RooRealVar("pTph1RECO","pTph1RECO", RECOpTph1, 5, 500);
   RooRealVar* pTph2RECO = new RooRealVar("pTph2RECO","pTph2RECO", RECOpTph2, 5, 500);

   double RECOpTph1min = max(0.5, RECOpTph1-2*pTerrZ1_ph1);
   double RECOpTph2min = max(0.5, RECOpTph2-2*pTerrZ1_ph2);

   // observables pT1,2,ph1,ph2
   RooRealVar* pT1 = new RooRealVar("pT1", "pT1FIT", RECOpT1, RECOpT1min, RECOpT1+2*pTerrZ1_1 );
   RooRealVar* pT2 = new RooRealVar("pT2", "pT2FIT", RECOpT2, RECOpT2min, RECOpT2+2*pTerrZ1_2 );

   RooRealVar* m1 = new RooRealVar("m1","m1", Z1_1.M());
   RooRealVar* m2 = new RooRealVar("m2","m2", Z1_2.M());

   if(debug_) cout<<"m1 "<<m1->getVal()<<" m2 "<<m2->getVal()<<endl;

   double Vtheta1, Vphi1, Vtheta2, Vphi2;
   Vtheta1 = (Z1_1).Theta(); Vtheta2 = (Z1_2).Theta();
   Vphi1 = (Z1_1).Phi(); Vphi2 = (Z1_2).Phi();

   RooRealVar* theta1 = new RooRealVar("theta1","theta1",Vtheta1);
   RooRealVar* phi1   = new RooRealVar("phi1","phi1",Vphi1);
   RooRealVar* theta2 = new RooRealVar("theta2","theta2",Vtheta2);
   RooRealVar* phi2   = new RooRealVar("phi2","phi2",Vphi2);

   // dot product to calculate (p1+p2+ph1+ph2).M()
   RooFormulaVar E1("E1","TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1)))+@2*@2)",
   RooArgList(*pT1,*theta1,*m1));
   RooFormulaVar E2("E2","TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1)))+@2*@2)",
   RooArgList(*pT2,*theta2,*m2));
   if(debug_) cout<<"E1 "<<E1.getVal()<<"; E2 "<<E2.getVal()<<endl;

   /////

   RooRealVar* pTph1 = new RooRealVar("pTph1", "pTph1FIT", RECOpTph1, RECOpTph1min, RECOpTph1+2*pTerrZ1_ph1 );
   RooRealVar* pTph2 = new RooRealVar("pTph2", "pTph2FIT", RECOpTph2, RECOpTph2min, RECOpTph2+2*pTerrZ1_ph2 );

   double Vthetaph1, Vphiph1, Vthetaph2, Vphiph2;
   Vthetaph1 = (Z1_ph1).Theta(); Vthetaph2 = (Z1_ph2).Theta();
   Vphiph1 = (Z1_ph1).Phi(); Vphiph2 = (Z1_ph2).Phi();

   RooRealVar* thetaph1 = new RooRealVar("thetaph1","thetaph1",Vthetaph1);
   RooRealVar* phiph1   = new RooRealVar("phiph1","phiph1",Vphiph1);
   RooRealVar* thetaph2 = new RooRealVar("thetaph2","thetaph2",Vthetaph2);
   RooRealVar* phiph2   = new RooRealVar("phiph2","phi2",Vphiph2);

   RooFormulaVar Eph1("Eph1","TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1))))", 
   RooArgList(*pTph1,*thetaph1));
   RooFormulaVar Eph2("Eph2","TMath::Sqrt((@0*@0)/((TMath::Sin(@1))*(TMath::Sin(@1))))", 
   RooArgList(*pTph2,*thetaph2));

   //// dot products of 4-vectors

   // 3-vector DOT
   RooFormulaVar* p1v3D2 = new RooFormulaVar("p1v3D2",
   "@0*@1*( ((TMath::Cos(@2))*(TMath::Cos(@3)))/((TMath::Sin(@2))*(TMath::Sin(@3)))+(TMath::Cos(@4-@5)))",
   RooArgList(*pT1,*pT2,*theta1,*theta2,*phi1,*phi2));    
   if(debug_) cout<<"p1 DOT p2 is "<<p1v3D2->getVal()<<endl;
   // 4-vector DOT metric 1 -1 -1 -1
   RooFormulaVar p1D2("p1D2","@0*@1-@2",RooArgList(E1,E2,*p1v3D2));

   //lep DOT fsrPhoton1

   // 3-vector DOT
   RooFormulaVar* p1v3Dph1 = new RooFormulaVar("p1v3Dph1",
   "@0*@1*( (TMath::Cos(@2)*TMath::Cos(@3))/(TMath::Sin(@2)*TMath::Sin(@3))+TMath::Cos(@4-@5))",
   RooArgList(*pT1,*pTph1,*theta1,*thetaph1,*phi1,*phiph1));

   // 4-vector DOT metric 1 -1 -1 -1
   RooFormulaVar p1Dph1("p1Dph1","@0*@1-@2",RooArgList(E1,Eph1,*p1v3Dph1));

   // 3-vector DOT
   RooFormulaVar* p2v3Dph1 = new RooFormulaVar("p2v3Dph1",
   "@0*@1*( (TMath::Cos(@2)*TMath::Cos(@3))/(TMath::Sin(@2)*TMath::Sin(@3))+TMath::Cos(@4-@5))",
   RooArgList(*pT2,*pTph1,*theta2,*thetaph1,*phi2,*phiph1));
   // 4-vector DOT metric 1 -1 -1 -1
   RooFormulaVar p2Dph1("p2Dph1","@0*@1-@2",RooArgList(E2,Eph1,*p2v3Dph1));

   // lep DOT fsrPhoton2 

   // 3-vector DOT
   RooFormulaVar* p1v3Dph2 = new RooFormulaVar("p1v3Dph2",
   "@0*@1*( (TMath::Cos(@2)*TMath::Cos(@3))/(TMath::Sin(@2)*TMath::Sin(@3))+TMath::Cos(@4-@5))",
   RooArgList(*pT1,*pTph2,*theta1,*thetaph2,*phi1,*phiph2));

   // 4-vector DOT metric 1 -1 -1 -1
   RooFormulaVar p1Dph2("p1Dph2","@0*@1-@2",RooArgList(E1,Eph2,*p1v3Dph2));

   // 3-vector DOT
   RooFormulaVar* p2v3Dph2 = new RooFormulaVar("p2v3Dph2",
   "@0*@1*( (TMath::Cos(@2)*TMath::Cos(@3))/(TMath::Sin(@2)*TMath::Sin(@3))+TMath::Cos(@4-@5))",
   RooArgList(*pT2,*pTph2,*theta2,*thetaph2,*phi2,*phiph2));
   // 4-vector DOT metric 1 -1 -1 -1
   RooFormulaVar p2Dph2("p2Dph2","@0*@1-@2",RooArgList(E2,Eph2,*p2v3Dph2));

   // fsrPhoton1 DOT fsrPhoton2

   // 3-vector DOT
   RooFormulaVar* ph1v3Dph2 = new RooFormulaVar("ph1v3Dph2",
   "@0*@1*( (TMath::Cos(@2)*TMath::Cos(@3))/(TMath::Sin(@2)*TMath::Sin(@3))+TMath::Cos(@4-@5))",
   RooArgList(*pTph1,*pTph2,*thetaph1,*thetaph2,*phiph1,*phiph2));    
   // 4-vector DOT metric 1 -1 -1 -1
   RooFormulaVar ph1Dph2("ph1Dph2","@0*@1-@2",RooArgList(Eph1,Eph2,*ph1v3Dph2));

   // mZ1

   RooFormulaVar* mZ1;
   if(p4sZ1ph_.size()==1)
   mZ1 = new RooFormulaVar("mZ1","TMath::Sqrt(2*@0+2*@1+2*@2+@3*@3+@4*@4)",
   RooArgList(p1D2, p1Dph1, p2Dph1, *m1,*m2));
   else if(p4sZ1ph_.size()==2)
   mZ1 = new RooFormulaVar("mZ1","TMath::Sqrt(2*@0+2*@1+2*@2+2*@3+2*@4+2*@5+@6*@6+@7*@7)",
   RooArgList(p1D2,p1Dph1,p2Dph1,p1Dph2,p2Dph2,ph1Dph2, *m1,*m2));
   else
   mZ1 = new RooFormulaVar("mZ1","TMath::Sqrt(2*@0+@1*@1+@2*@2)",RooArgList(p1D2,*m1,*m2));

   if(debug_) cout<<"mZ1 is "<<mZ1->getVal()<<endl;

   // pTerrs, 1,2,ph1,ph2
   RooRealVar sigmaZ1_1("sigmaZ1_1", "sigmaZ1_1", pTerrZ1_1);
   RooRealVar sigmaZ1_2("sigmaZ1_2", "sigmaZ1_2", pTerrZ1_2);

   RooRealVar sigmaZ1_ph1("sigmaZ1_ph1", "sigmaZ1_ph1", pTerrZ1_ph1);
   RooRealVar sigmaZ1_ph2("sigmaZ1_ph2", "sigmaZ1_ph2", pTerrZ1_ph2);

   // resolution for decay products
   RooGaussian gauss1("gauss1","gaussian PDF", *pT1RECO, *pT1, sigmaZ1_1);
   RooGaussian gauss2("gauss2","gaussian PDF", *pT2RECO, *pT2, sigmaZ1_2);

   RooGaussian gaussph1("gaussph1","gaussian PDF", *pTph1RECO, *pTph1, sigmaZ1_ph1);
   RooGaussian gaussph2("gaussph2","gaussian PDF", *pTph2RECO, *pTph2, sigmaZ1_ph2);

   RooRealVar bwMean("bwMean", "m_{Z^{0}}", 91.187);
   RooRealVar bwGamma("bwGamma", "#Gamma", 2.5);

   RooRealVar sg("sg", "sg", sgVal_);
   RooRealVar a("a", "a", aVal_);
   RooRealVar n("n", "n", nVal_);

   RooCBShape CB("CB","CB",*mZ1,bwMean,sg,a,n);
   RooRealVar f("f","f", fVal_);

   RooRealVar mean("mean","mean",meanVal_);
   RooRealVar sigma("sigma","sigma",sigmaVal_);
   RooRealVar f1("f1","f1",f1Val_);

   RooGenericPdf RelBW("RelBW","1/( pow(mZ1*mZ1-bwMean*bwMean,2)+pow(mZ1,4)*pow(bwGamma/bwMean,2) )", RooArgSet(*mZ1,bwMean,bwGamma) );

   RooAddPdf RelBWxCB("RelBWxCB","RelBWxCB", RelBW, CB, f);
   RooGaussian gauss("gauss","gauss",*mZ1,mean,sigma);
   RooAddPdf RelBWxCBxgauss("RelBWxCBxgauss","RelBWxCBxgauss", RelBWxCB, gauss, f1);

   RooProdPdf *PDFRelBWxCBxgauss;

   if(p4sZ1ph_.size()==1)    
   PDFRelBWxCBxgauss = new RooProdPdf("PDFRelBWxCBxgauss","PDFRelBWxCBxgauss", 
   RooArgList(gauss1, gauss2, gaussph1, RelBWxCBxgauss) );
   else if(p4sZ1ph_.size()==2)
   PDFRelBWxCBxgauss = new RooProdPdf("PDFRelBWxCBxgauss","PDFRelBWxCBxgauss", 
   RooArgList(gauss1, gauss2, gaussph1, gaussph2, RelBWxCBxgauss) );
   else
   PDFRelBWxCBxgauss = new RooProdPdf("PDFRelBWxCBxgauss","PDFRelBWxCBxgauss", 
   RooArgList(gauss1, gauss2, RelBWxCBxgauss) );

   // observable set
   RooArgSet *rastmp;
   if(p4sZ1ph_.size()==1)
   rastmp = new RooArgSet(*pT1RECO,*pT2RECO,*pTph1RECO);
   else if(p4sZ1ph_.size()>=2)
   rastmp = new RooArgSet(*pT1RECO,*pT2RECO,*pTph1RECO,*pTph2RECO);
   else
   rastmp = new RooArgSet(*pT1RECO,*pT2RECO);

   RooDataSet* pTs = new RooDataSet("pTs","pTs", *rastmp);
   pTs->add(*rastmp); 

   //RooAbsReal* nll;
   //nll = PDFRelBWxCBxgauss->createNLL(*pTs);
   //RooMinuit(*nll).migrad();

   //RooFitResult* r = PDFRelBWxCBxgauss->fitTo(*pTs,RooFit::Save(),RooFit::PrintLevel(-1),RooFit::PrintEvalErrors(-1));
   RooFitResult* r = PDFRelBWxCBxgauss->fitTo(*pTs,RooFit::Save(),RooFit::Constrain(*mZ1),RooFit::PrintLevel(-1));
   const TMatrixDSym& covMatrix = r->covarianceMatrix();
   
   const RooArgList& finalPars = r->floatParsFinal();
   for (int i=0 ; i<finalPars.getSize(); i++){
   TString name = TString(((RooRealVar*)finalPars.at(i))->GetName());

   if(debug_) cout<<"name list of RooRealVar for covariance matrix "<<name<<endl;

   }

   int size = covMatrix.GetNcols();
   //TMatrixDSym covMatrixTest_(size);
   covMatrixZ1_.ResizeTo(size,size);
   covMatrixZ1_ = covMatrix;   

   if(debug_) cout<<"save the covariance matrix"<<endl;
    
   l1 = pT1->getVal()/RECOpT1; l2 = pT2->getVal()/RECOpT2;
   double pTerrZ1REFIT1 = pT1->getError(); double pTerrZ1REFIT2 = pT2->getError();

   pTerrsZ1REFIT_.push_back(pTerrZ1REFIT1);
   pTerrsZ1REFIT_.push_back(pTerrZ1REFIT2);

   if(p4sZ1ph_.size()>=1){

   if(debug_) cout<<"set refit result for Z1 fsr photon 1"<<endl;

   lph1 = pTph1->getVal()/RECOpTph1;
   double pTerrZ1phREFIT1 = pTph1->getError();
   if(debug_) cout<<"scale "<<lph1<<" pterr "<<pTerrZ1phREFIT1<<endl;  
   
   pTerrsZ1phREFIT_.push_back(pTerrZ1phREFIT1);

   } 
   if(p4sZ1ph_.size()==2){

   lph2 = pTph2->getVal()/RECOpTph2;
   double pTerrZ1phREFIT2 = pTph2->getError();
   pTerrsZ1phREFIT_.push_back(pTerrZ1phREFIT2);

   }

   //delete nll;
   if (r) delete r;
   delete mZ1;
   delete pT1; delete pT2; delete pTph1; delete pTph2;
   delete pT1RECO; delete pT2RECO; delete pTph1RECO; delete pTph2RECO;
   delete ph1v3Dph2; delete p1v3Dph1; delete p2v3Dph1; delete p1v3Dph2; delete p2v3Dph2;
   delete PDFRelBWxCBxgauss;
   delete pTs;
   delete rastmp;
   delete m1; delete m2;
   delete theta1; delete theta2; delete phi1; delete phi2;
   delete thetaph1; delete thetaph2; delete phiph1; delete phiph2;
   delete p1v3D2;

   if(debug_) cout<<"end Z1 refit"<<endl;

   return 0;

   } */

/*
  TMatrixDSym GetRefitZ1BigCov(){

  


  }
*/
     
/*
  RooFormulaVar KinZfitter::p1DOTp2(RooRealVar pT1, RooRealVar theta1, RooRealVar phi1, RooRealVar m1, TString index1, RooRealVar pT2, RooRealVar theta2, RooRealVar phi2, RooRealVar m2, TString index2)
  {

  RooFormulaVar cosDeltaPhi("cosDeltaPhi"+index1+index2,"TMath::Cos(@0-@1)",RooArgList(phi1,phi2));

  RooFormulaVar cosTheta1("cosTheta1"+index1,"TMath::Cos(@0)",RooArgList(theta1));
  RooFormulaVar sinTheta1("sinTheta1"+index1,"TMath::Sin(@0)",RooArgList(theta1));
  RooFormulaVar cosTheta2("cosTheta2"+index2,"TMath::Cos(@0)",RooArgList(theta2));
  RooFormulaVar sinTheta2("sinTheta2"+index2,"TMath::Sin(@0)",RooArgList(theta2));

  RooFormulaVar E1("E"+index1,"TMath::Sqrt((@0*@0)/(@1*@1)+@2*@2)", RooArgList(pT1,sinTheta1,m1));
  RooFormulaVar E2("E"+index2,"TMath::Sqrt((@0*@0)/(@1*@1)+@2*@2)", RooArgList(pT2,sinTheta2,m2));
  // 3-vector DOT
  RooFormulaVar p1Dp2("p"+index1+"Dp"+index2,"@0*@1*(@2+(@3*@4)/(@5*@6))",RooArgList(pT1,pT2,cosDeltaPhi,cosTheta1,cosTheta2,sinTheta1,sinTheta2));
  // 4-vector DOT metric 1 -1 -1 -1
  RooFormulaVar p1Dv4p2("p"+index1+"DOTp"+index2,"@0*@1-@2",RooArgList(E1,E2,p1Dp2));


  return p1Dv4p2;
  }
*/

#endif

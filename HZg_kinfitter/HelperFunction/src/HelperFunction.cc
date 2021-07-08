
// -*- C++ -*-
//
// Package:     Subsystem/Package
// Class  :     HelperFunction
// 
// Implementation:
//     [Notes on implementation]
//
// Original Author:  Tongguang Cheng
//         Created:  Mon, 21 Dec 2015 12:47:33 GMT
//

#ifndef HelperFunction_CC
#define HelperFunction_CC
// system include files

// user include files
#include "KinZfitter/HelperFunction/interface/HelperFunction.h"

// fileinPath
#include "FWCore/ParameterSet/interface/FileInPath.h"


HelperFunction::HelperFunction()
{

  //declarations
  debug_ = 1;

  // MORIOND 17
  TString s_corr_e_1 = TString("KinZfitter/HelperFunction/hists/DYJetsToLL_M-50_m2eLUT_m2e_1.root");
  TString s_corr_e_2 = TString("KinZfitter/HelperFunction/hists/DYJetsToLL_M-50_m2eLUT_m2e_2.root");
  TString s_corr_e_3 = TString("KinZfitter/HelperFunction/hists/DYJetsToLL_M-50_m2eLUT_m2e_3.root");
  TString s_corr_mu = TString("KinZfitter/HelperFunction/hists/DYJetsToLL_M-50_m2muLUT_m2mu.root");

  f_corr_e_1 = boost::shared_ptr<TFile>( new TFile(s_corr_e_1)); 
  f_corr_e_2 = boost::shared_ptr<TFile>( new TFile(s_corr_e_2)); 
  f_corr_e_3 = boost::shared_ptr<TFile>( new TFile(s_corr_e_3)); 
  f_corr_mu = boost::shared_ptr<TFile>( new TFile(s_corr_mu));
        
  el_corr_1 = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(f_corr_e_1->Get("2e")->Clone() )) );
  el_corr_2 = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(f_corr_e_2->Get("2e")->Clone() )) );
  el_corr_3 = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(f_corr_e_3->Get("2e")->Clone() )) );

  mu_corr = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(f_corr_mu->Get("2mu")->Clone() )) );
        
  x_elpTaxis_1 = el_corr_1->GetXaxis(); y_eletaaxis_1 = el_corr_1->GetYaxis();
  maxPtEl_1 = x_elpTaxis_1->GetXmax(); minPtEl_1 = x_elpTaxis_1->GetXmin();

  x_elpTaxis_2 = el_corr_2->GetXaxis(); y_eletaaxis_2 = el_corr_2->GetYaxis();
  maxPtEl_2 = x_elpTaxis_2->GetXmax(); minPtEl_2 = x_elpTaxis_2->GetXmin();

  x_elpTaxis_3 = el_corr_3->GetXaxis(); y_eletaaxis_3 = el_corr_3->GetYaxis();
  maxPtEl_3 = x_elpTaxis_3->GetXmax(); minPtEl_3 = x_elpTaxis_3->GetXmin();

  x_mupTaxis = mu_corr->GetXaxis(); y_muetaaxis = mu_corr->GetYaxis();
  maxPtMu = x_mupTaxis->GetXmax(); minPtMu = x_mupTaxis->GetXmin();

}


HelperFunction::~HelperFunction(){}

double HelperFunction:: masserrorFullCov(std::vector<TLorentzVector> p4s, TMatrixDSym covMatrix)
{
  int ndim = 3*p4s.size();
  //  if(debug_) cout<<""<<endl;

  TMatrixD jacobian(1,ndim);

  double e = 0; double mass = 0;
  double px = 0; double py = 0; double pz = 0;
  for (unsigned int ip = 0; ip < p4s.size(); ip++) {
         
    e = e + p4s[ip].E();
    px = px + p4s[ip].Px();
    py = py + p4s[ip].Py();
    pz = pz + p4s[ip].Pz();
  }

  mass = TMath::Sqrt(e*e-px*px-py*py-pz*pz);

  for (unsigned int i = 0, o = 0; i < p4s.size(); i++, o += 3) {

    double pxi = p4s[i].Px();
    double pyi = p4s[i].Py();
    double pzi = p4s[i].Pz();
    double ei = p4s[i].E();

    jacobian(0, o+0) = (e*(pxi/ei) - px)/mass;
    jacobian(0, o+1) = (e*(pyi/ei) - py)/mass;
    jacobian(0, o+2) = (e*(pzi/ei) - pz)/mass;
  }

  TMatrixDSym massCov = covMatrix.Similarity(jacobian);

  double dm2 = massCov(0,0);
  return (dm2 > 0 ? std::sqrt(dm2) : 0.0);
}


double HelperFunction::masserror( std::vector<TLorentzVector> Lep, std::vector<double> pterr){
  // if(Lep.size()!= pterr.size()!=4) {std::cout<<" Lepsize="<<Lep.size()<<", "<<pterr.size()<<std::endl;}
  TLorentzVector compositeParticle ;
  for(unsigned int i=0; i<Lep.size(); i++){
    compositeParticle+=Lep[i];
  }
  double mass  =  compositeParticle.M();

  double masserr = 0;

  for(unsigned int i=0; i<Lep.size(); i++){
    TLorentzVector variedLep; // = Lep[i];

    variedLep.SetPtEtaPhiM(Lep[i].Pt()+ pterr[i], Lep[i].Eta(), Lep[i].Phi(), Lep[i].M());
    TLorentzVector compositeParticleVariation ;
    for(unsigned int j=0; j<Lep.size(); j++){
      if(i!=j)compositeParticleVariation+=Lep[j];
      else compositeParticleVariation+=variedLep;
    }

    masserr += (compositeParticleVariation.M()-mass)*(compositeParticleVariation.M()-mass);
  }

  return sqrt(masserr);
}


double HelperFunction::pterr(TreeReader &data, Int_t lepID, Int_t fsrID, Int_t channel,  bool isData){
  // reco::GsfElectron *gsf; reco::Muon *mu;
  // reco::PFCandidate *pf;
  // pat::Muon *patmu;
  double pterrLep = 0.0;
  float *elePt = data.GetPtrFloat("elePt");
  float *muPt = data.GetPtrFloat("muPt");
  float *eleEta = data.GetPtrFloat("eleEta");
  float *muEta = data.GetPtrFloat("muEta");
  Int_t *eleEcalDrivenSeed = data.GetPtrInt("eleEcalDrivenSeed");
  float *pfPhoEtErr = data.GetPtrFloat("pfPhoEtErr");
  
  
  if (channel == 0) {
  double pT_e = elePt[lepID];
  double eta_e = eleEta[lepID];
    pterrLep = correlepterr(data,lepID);
    if(eleEcalDrivenSeed[lepID]==1){		//ecal driven
      if (fabs(eta_e) < 1) {
	if (pterrLep/pT_e < 0.03) { // hardcode 1
	  int xbin = x_elpTaxis_1->FindBin(pT_e);
	  int ybin = y_eletaaxis_1->FindBin(fabs(eta_e));
	  if(pT_e > minPtEl_1 && pT_e < maxPtEl_1 ){
	    pterrLep*=el_corr_1->GetBinContent(xbin,ybin);
	  } else {
	    pterrLep*=1.0;
	  }
	} else {pterrLep*=1.03;} // hardcode 2        
      } else if (fabs(eta_e) > 1 && fabs(eta_e) < 2.5) {
	if (pterrLep/pT_e < 0.07) { // hardcode 3
	  int xbin = x_elpTaxis_2->FindBin(pT_e);
	  int ybin = y_eletaaxis_2->FindBin(fabs(eta_e));
	  if(pT_e > minPtEl_2 && pT_e < maxPtEl_2 ){
	    pterrLep*=el_corr_2->GetBinContent(xbin,ybin);
	  } else {
	    pterrLep*=1.0;
	  }
	} else {pterrLep*=0.74;} // hardcode 4
      } // 1 < |eta| < 2.5      
    }
    else //tracker driven
      {
	int xbin = x_elpTaxis_3->FindBin(pT_e);
	int ybin = y_eletaaxis_3->FindBin(fabs(eta_e));
	if(pT_e > minPtEl_3 && pT_e < maxPtEl_3 ){
	  pterrLep*=el_corr_3->GetBinContent(xbin,ybin);
	} else {
	  pterrLep*=1.0;
	}
      }
  }
  else 
    {
      float* muBestTrkPtError = data.GetPtrFloat("muBestTrkPtError");
      pterrLep= muBestTrkPtError[lepID];
      if(debug_)cout<<"reco pt err is "<<pterrLep<<endl;

      int xbin = x_mupTaxis->FindBin(muPt[lepID]);
      int ybin = y_muetaaxis->FindBin(fabs(muEta[lepID]));
      if(muPt[lepID]>minPtMu && muPt[lepID]<maxPtMu ){
	pterrLep*=mu_corr->GetBinContent(xbin,ybin);
      } else {
	pterrLep*=1.0;
      }
	  if(fsrID!=-99)//FSR constribution
	  {
		  pterrLep = pfPhoEtErr[fsrID];
	  }
    }

  return pterrLep;

}
double HelperFunction::correlepterr(TreeReader &data, Int_t lepID)
{
	float *ecalEnergy = data.GetPtrFloat("eleCalibEn");
	float *eleSCEta = data.GetPtrFloat("eleSCEta");
	float *elePt = data.GetPtrFloat("elePt");
	float *eleEta = data.GetPtrFloat("eleEta");
	float *elePhi = data.GetPtrFloat("elePhi");
	Int_t *eleEcalDrivenSeed = data.GetPtrInt("eleEcalDrivenSeed");
	float *elePtError = data.GetPtrFloat("elePtError");
	double pterr = 0.;
	if(eleEcalDrivenSeed[lepID]==1)
	  {
	    pterr = elePtError[lepID];
	  }
	else {
	TLorentzVector ele;
	ele.SetPtEtaPhiM(elePt[lepID],eleEta[lepID],elePhi[lepID],0.511*0.001);
	float err2 = 0.0;
	if(fabs(eleSCEta[lepID])< 1.479)
	{
		 err2 += (5.24e-02*5.24e-02)/ecalEnergy[lepID];
		 err2 += (2.01e-01*2.01e-01)/(ecalEnergy[lepID]*ecalEnergy[lepID]);
		 err2 += 1.00e-02*1.00e-02;
	}
	else 
	{
		err2 += (1.46e-01*1.46e-01)/ecalEnergy[lepID];
		err2 += (9.21e-01*9.21e-01)/(ecalEnergy[lepID]*ecalEnergy[lepID]);
		err2 += 1.94e-03*1.94e-03;	
	}
	float perr = ecalEnergy[lepID]*sqrt(err2);
	pterr = perr*elePt[lepID]/ele.P();
	}
	return pterr;
}
#endif

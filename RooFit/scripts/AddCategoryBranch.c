
void AddCategoryBranch(double sigmu = 0) {
  TChain chain("tree");
  chain.Add("/Users/adorsett/Desktop/CERNbox/pico/NanoAODv5/zgamma_channelIslandsv3/2016/mc/merged_zgmc_llg/*_DYJets*");
  chain.Add("/Users/adorsett/Desktop/CERNbox/pico/NanoAODv5/zgamma_channelIslandsv3/2016/mc/merged_zgmc_llg/*_ZGToLLG*");
  chain.Add("/Users/adorsett/Desktop/CERNbox/pico/NanoAODv5/zgamma_channelIslandsv3/2016/mc/merged_zgmc_llg/*_LLAJJ*");
  if(sigmu > 0) {
    chain.Add("/Users/adorsett/Desktop/CERNbox/pico/NanoAODv5/zgamma_channelIslandsv3/2016/HToZG/merged_zgmc_llg/*GluGluH*");
    chain.Add("/Users/adorsett/Desktop/CERNbox/pico/NanoAODv5/zgamma_channelIslandsv3/2016/HToZG/merged_zgmc_llg/*VBFH*");
  }
  chain.SetBranchStatus("*",1);
  float vbf_mva, weight;
  int nel, nmu, type;
  vector<float> *el_pt(nullptr),       *mu_pt(nullptr),        *llphoton_m(nullptr), *ll_m(nullptr);
  vector<float> *photon_pt(nullptr),   *photon_eta(nullptr),   *photon_phi(nullptr), *photon_drmin(nullptr);
  vector<float> *llphoton_pt(nullptr), *llphoton_eta(nullptr), *llphoton_phi(nullptr);
  vector<float> *ll_pt(nullptr),       *ll_eta(nullptr),       *ll_phi(nullptr);
  vector<int>   *ll_lepid(nullptr),    *ll_i1(nullptr),        *ll_i2(nullptr);
  chain.SetBranchAddress("weight",      &weight);
  chain.SetBranchAddress("type",        &type);
  chain.SetBranchAddress("nel",         &nel);
  chain.SetBranchAddress("nmu",         &nmu);
  chain.SetBranchAddress("vbf_mva",     &vbf_mva);
  chain.SetBranchAddress("el_pt",       &el_pt);
  chain.SetBranchAddress("mu_pt",       &mu_pt);
  chain.SetBranchAddress("photon_pt",   &photon_pt);
  chain.SetBranchAddress("photon_eta",  &photon_eta);
  chain.SetBranchAddress("photon_phi",  &photon_phi);
  chain.SetBranchAddress("photon_drmin",&photon_drmin);
  chain.SetBranchAddress("llphoton_m",  &llphoton_m);
  chain.SetBranchAddress("llphoton_pt", &llphoton_pt);
  chain.SetBranchAddress("llphoton_eta",&llphoton_eta);
  chain.SetBranchAddress("llphoton_phi",&llphoton_phi);
  chain.SetBranchAddress("ll_m",        &ll_m);
  chain.SetBranchAddress("ll_pt",       &ll_pt);
  chain.SetBranchAddress("ll_eta",      &ll_eta);
  chain.SetBranchAddress("ll_phi",      &ll_phi);
  chain.SetBranchAddress("ll_lepid",    &ll_lepid);
  chain.SetBranchAddress("ll_i1",       &ll_i1);
  chain.SetBranchAddress("ll_i2",       &ll_i2);
  TH1D *cat1 = new TH1D("cat1","ll#gamma mass spectrum (VBF);            m_{ll#gamma} [GeV]; Events",80,100,180);
  TH1D *cat2 = new TH1D("cat2","ll#gamma mass spectrum (High rel p_{T,#gamma}); m_{ll#gamma} [GeV]; Events",80,100,180);
  TH1D *cat3 = new TH1D("cat3","ll#gamma mass spectrum (High p_{Tt} e);  m_{ll#gamma} [GeV]; Events",80,100,180);
  TH1D *cat4 = new TH1D("cat4","ll#gamma mass spectrum (Low p_{Tt} e);   m_{ll#gamma} [GeV]; Events",80,100,180);
  TH1D *cat5 = new TH1D("cat5","ll#gamma mass spectrum (High p_{Tt} #mu);m_{ll#gamma} [GeV]; Events",80,100,180);
  TH1D *cat6 = new TH1D("cat6","ll#gamma mass spectrum (Low p_{Tt} #mu); m_{ll#gamma} [GeV]; Events",80,100,180);
  for(int i = 0; i < chain.GetEntries(); i++) {
    chain.GetEntry(i);
    bool el = nel > 1 && ll_lepid->at(0) == 11 && el_pt->at(ll_i1->at(0)) > 25 && el_pt->at(ll_i2->at(0)) > 15;
    bool mu = nmu > 1 && ll_lepid->at(0) == 13 && mu_pt->at(ll_i1->at(0)) > 20 && mu_pt->at(ll_i2->at(0)) > 10;
    bool others = llphoton_m->at(0)+ll_m->at(0) >= 185 && 
                  photon_pt->at(0)/llphoton_m->at(0) >= 15./110 &&
                  photon_drmin->at(0) > 0.4;
    bool baseline = (el || mu) && others;
    if(!baseline) continue;
    TVector3 g, h, z;
    g.SetPtEtaPhi(photon_pt->at(0),photon_eta->at(0),photon_phi->at(0));
    h.SetPtEtaPhi(llphoton_pt->at(0),llphoton_eta->at(0),llphoton_phi->at(0));
    z.SetPtEtaPhi(ll_pt->at(0),ll_eta->at(0),ll_phi->at(0));
    g.SetZ(0); h.SetZ(0); z.SetZ(0);
    double pTt = h.Cross((z-g).Unit()).Mag();
    double relpT = photon_pt->at(0)/llphoton_m->at(0);
    double llgm = llphoton_m->at(0);
    if(type >= 200000 && type <= 206000) weight = weight*sigmu;
    if(vbf_mva > -0.9)      cat1->Fill(llgm,weight);
    else if(relpT > 0.4)    cat2->Fill(llgm,weight);
    else if(pTt > 40 && el) cat3->Fill(llgm,weight);
    else if(pTt < 40 && el) cat4->Fill(llgm,weight);
    else if(pTt > 40 && mu) cat5->Fill(llgm,weight);
    else if(pTt < 40 && mu) cat6->Fill(llgm,weight);
  }
  TFile out("Fits/categories.root","RECREATE");
  cat1->Write();
  cat2->Write();
  cat3->Write();
  cat4->Write();
  cat5->Write();
  cat6->Write();
  out.Close();
  return;
}

void FitGausCats() {
//   MakeCats();
  TFile *in = TFile::Open("Fits/categories.root");
  TH1D *cat1 = (TH1D*)in->Get("cat1");
  TH1D *cat2 = (TH1D*)in->Get("cat2");
  TH1D *cat3 = (TH1D*)in->Get("cat3");
  TH1D *cat4 = (TH1D*)in->Get("cat4");
  TH1D *cat5 = (TH1D*)in->Get("cat5");
  TH1D *cat6 = (TH1D*)in->Get("cat6");
  GausFit(cat1,4);
  GausFit(cat2,4);
  GausFit(cat3,4);
  GausFit(cat4,1);
  GausFit(cat5,4);
  GausFit(cat6,1);
  return;
}


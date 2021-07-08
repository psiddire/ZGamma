void MakeCategories(double lumi, int mu) {
  TString inpath = TString::Format("pseudodata/categories_mu%d.root",mu);
  TFile *infile = TFile::Open(inpath);
  TH1D *cat1 = (TH1D*)infile->Get("cat1");
  TH1D *cat2 = (TH1D*)infile->Get("cat2");
  TH1D *cat3 = (TH1D*)infile->Get("cat3");
  TH1D *cat4 = (TH1D*)infile->Get("cat4");
  TH1D *cat5 = (TH1D*)infile->Get("cat5");
  TH1D *cat6 = (TH1D*)infile->Get("cat6");
  vector<TH1D*> cats = {cat1,cat2,cat3,cat4,cat5,cat6};
  TString outname = TString::Format("pseudodata/category_mu%d_%.0fifb.root",mu,lumi);
  TFile *out = new TFile(outname,"RECREATE");
  TTree *tree = new TTree("tree","tree");
  double llphoton_m;
  int category;
  tree->Branch("llphoton_m",&llphoton_m);
  tree->Branch("category",&category);
  for(size_t ic(0); ic < cats.size(); ic++) 
    for(int i = 0; i < lumi*cats.at(ic)->Integral(); i++) {
      llphoton_m = cats.at(ic)->GetRandom();
      category = ic + 1;
      tree->Fill();
    }
  tree->Write();
  out->Close();
  return;
}

void MakePseudodata(double lumi, int mu, bool includeBackground = true){
  TChain chain("tree");
  chain.Add("/Users/adorsett/Desktop/CERNbox/pico/NanoAODv5/zgamma_channelIslandsv3/2016/mc/merged_zgmc_zgsr/*_DYJets*");
  chain.Add("/Users/adorsett/Desktop/CERNbox/pico/NanoAODv5/zgamma_channelIslandsv3/2016/mc/merged_zgmc_zgsr/*_ZGToLLG*");
  chain.Add("/Users/adorsett/Desktop/CERNbox/pico/NanoAODv5/zgamma_channelIslandsv3/2016/mc/merged_zgmc_zgsr/*_LLAJJ*");
  chain.Add("/Users/adorsett/Desktop/CERNbox/pico/NanoAODv5/zgamma_channelIslandsv3/2016/HToZG/merged_zgmc_zgsr/*");
  TTree tree("tree","tree");
  double llphoton_m;
  tree.Branch("llphoton_m",&llphoton_m);
  TH1D *sig = new TH1D("sig","ll#gamma mass spectrum; m_{ll#gamma}; Events",80,100,180);
  chain.SetBranchStatus("*",1);
  TString draw_cut = TString::Format("%f*weight*(ll_lepid[0]==13)*(%d*(type>=200000 && type <= 207000) + %d*(type < 200000))",
                                      lumi,mu,includeBackground);
  chain.Draw("llphoton_m[0]>>sig",draw_cut);
  for(int i = 0; i < sig->Integral(); i++) {
    llphoton_m = sig->GetRandom();
    tree.Fill();
  }
  TString outname = TString::Format("pseudodata/mu%d_%.0fifb.root",mu,lumi);
  if(!includeBackground) outname = TString::Format("pseudodata/signal_mu%d_%.0fifb.root",mu,lumi);
  TFile out(outname,"RECREATE");
  tree.Write();
  out.Close();
  return;
}

void GeneratePseudodata() {
  MakeCategories(1,0);
  MakeCategories(1,100);
  MakeCategories(36,0);
  MakeCategories(36,100);
//   MakePseudodata(1,100);
//   MakePseudodata(36,0);
//   MakePseudodata(36,10);
//   MakePseudodata(36,100);
//   MakePseudodata(36,100,false);
  return;
}


void MakeTree(int year = 2016, int sigmu = 0, double lumi = 36) {
  if(year == 2017 && lumi == 36) lumi = 42;
  TChain chain("tree");
  TString base = "/Users/adorsett/Desktop/CERNbox/pico/NanoAODv5/zgamma_channelIslandsv3/"+to_string(year);
  chain.Add(base+"/mc/merged_zgmc_zgsr/*_DYJets*");
  chain.Add(base+"/mc/merged_zgmc_zgsr/*_ZGToLLG*");
  chain.Add(base+"/mc/merged_zgmc_zgsr/*_LLAJJ*");
  chain.Add(base+"/HToZG/merged_zgmc_zgsr/*");
  chain.SetBranchStatus("*",0);
  chain.SetBranchStatus("llphoton_m",1);
  chain.SetBranchStatus("weight",1);
  chain.SetBranchStatus("type",1);
  vector<float> *llphoton_m(nullptr);
  float weight;
  int type;
  chain.SetBranchAddress("llphoton_m",&llphoton_m);
  chain.SetBranchAddress("weight",&weight);
  chain.SetBranchAddress("type",&type);
  TString outname = TString::Format("MC/%d/combined-MC-mu%d-%.0fifb.root",year,sigmu,lumi);
  TFile out(outname,"RECREATE");
  TTree tree("tree","tree");
  double llgm;
  tree.Branch("llphoton_m",&llgm);
  tree.Branch("weight",&weight);
  tree.Branch("type",&type);
  for(int i(0); i < chain.GetEntries(); i++) {
    chain.GetEntry(i);
    llgm = llphoton_m->at(0);
    if(type >= 200000 && type <= 207000) weight = weight*sigmu;
    weight = weight*lumi;
    tree.Fill();
  }
  tree.Write();
  out.Close();
  return tree;
}

void CombinePicos() {
  MakeTree(2017,0);
  MakeTree(2017,1);
  MakeTree(2017,2);
  MakeTree(2017,3);
  MakeTree(2017,4);
  MakeTree(2017,5);
  MakeTree(2017,6);
  MakeTree(2017,7);
  MakeTree(2017,8);
  MakeTree(2017,9);
  MakeTree(2017,10);
  return;
}

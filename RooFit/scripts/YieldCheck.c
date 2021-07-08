
void YieldCheck() {
  int year(2017);
  TChain chain("tree");
  TString base = "/Users/adorsett/Desktop/CERNbox/pico/NanoAODv5/zgamma_channelIslandsv3/"+to_string(year);
  chain.Add(base+"/mc/merged_zgmc_zgsr/*_DYJets*");
  chain.Add(base+"/mc/merged_zgmc_zgsr/*_ZGToLLG*");
  chain.Add(base+"/mc/merged_zgmc_zgsr/*_LLAJJ*");
  TH1D *mllg = new TH1D("mllg","llg mass",80,100,180);
  chain.Draw("llphoton_m[0]>>mllg","weight");
  cout << "Background yield/ifb: " << mllg->Integral() << endl;
  TChain schain("tree");
  schain.Add(base+"/HToZG/merged_zgmc_zgsr/*");
  TH1D *smllg = new TH1D("smllg","llg mass",80,100,180);
  schain.Draw("llphoton_m[0]>>smllg","weight");
  cout << "Signal yield/ifb: " << smllg->Integral() << endl;
  return;
}

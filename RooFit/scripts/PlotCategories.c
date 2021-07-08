
void PlotCategory(int icat = 0) {
  vector<TString> catfiles = {"bkghist_ele_mu_dijet1.root",
                              "bkghist_ele_mu_dijet2.root",
                              "bkghist_ele_mu_dijet3.root",
                              "bkghist_ele_mu_untag1.root",
                              "bkghist_ele_mu_untag2.root",
                              "bkghist_ele_mu_untag3.root",
                              "bkghist_ele_mu_untag4.root"};
  TString inpath = "MingYanFiles/"+catfiles.at(icat);
  TFile *infile = TFile::Open(inpath);
  TH1D *hmc = (TH1D*)infile->Get("hbkg");
  TCanvas *can = new TCanvas("can","canvas",800,800);
  TVirtualPad *pad = can->cd();
  pad->SetMargin(0.12,0.05,0.10,0.10);
  pad->SetPad(0,0,1,1);
  pad->cd();
  TString title("Background in ");
  if(icat < 3) title += TString::Format("dijet%d category",icat+1);
  else title += TString::Format("untagged%d category",icat-2);
  title += "; m_{ll#gamma} [GeV]; Events/GeV";
  hmc->SetTitle(title);
  hmc->GetYaxis()->SetLabelSize(0.03);
  hmc->Draw();
  can->SaveAs(TString::Format("MingYanFiles/plots/cat%d.pdf",icat));
  return;
}

void PlotCategories() {
  for(int i = 0; i < 7; i++) PlotCategory(i);
  return;
}


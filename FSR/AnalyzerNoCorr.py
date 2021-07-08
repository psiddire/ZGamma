import ROOT
import argparse
import WeightLumi

parser = argparse.ArgumentParser()
parser.add_argument('input_path')
parser.add_argument('output_path')
args = parser.parse_args()

ROOT.ROOT.EnableImplicitMT()

f = ROOT.TFile(args.input_path)
t = f.Get("tree")

h_ll = ROOT.TH1F("h_ll", "ll_mass", 100, 50, 150)
h_llg = ROOT.TH1F("h_llg", "llg_mass", 80, 100, 180)

weight = WeightLumi.weight_lumi(args.input_path)

for i in range(0, t.GetEntries()):

    if i%1000==0:
        print i
    t.GetEntry(i)
    if not (t.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 or t.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8):
        continue
    if t.nll < 1 or t.nphoton < 1:
        continue
    if all(dr < 0.4 for dr in t.photon_drmin):
        continue
    if t.ll_charge.size()==1 and t.ll_charge[0]!=0:
        continue
    if sum(map(bool, t.mu_id)) < 2 or t.nmu < 2:
        continue
    if t.nfsrphoton < 1 or t.nfsrphoton > 1:
        continue

    idll = []
    for z in range(t.nll):
        x = t.ll_i1[z]
        y = t.ll_i2[z]
        if (t.mu_charge[x]*t.mu_charge[y] == -1 and t.mu_pt[x] > 20 and t.mu_pt[y] > 10 and t.mu_id[x] and t.mu_id[y] \
            and abs(t.mu_dz[x]) < 1.0 and abs(t.mu_dxy[x]) < 0.5 and abs(t.mu_dz[y]) < 1.0 and abs(t.mu_dxy[y]) < 0.5 \
            and t.mu_reliso[x] < 0.25 and t.mu_reliso[y] < 0.25):
            idll.append(z)
            break
    if len(idll) == 0:
        continue
    massZ = t.ll_m[idll[0]]

    idllg = []
    iph = 0
    ph = ROOT.TLorentzVector()
    for z in range(t.nllphoton):
        if t.llphoton_ill[z]!=idll[0]:
            continue
        if (bool(t.photon_id[iph]) and t.photon_drmin[iph] > 0.4 and t.photon_pt[iph] > 15):
            ph.SetPtEtaPhiM(t.photon_pt[iph], t.photon_eta[iph], t.photon_phi[iph], 0)
            idllg.append(z)
            break
        iph = iph + 1
    if len(idllg) == 0:
        continue
    massH = t.llphoton_m[idllg[0]]

    if (massZ + massH) < 185:
        continue
    if (ph.E()/massH) < 15/110:
        continue
    h_ll.Fill(massZ, weight)
    h_llg.Fill(massH, weight)

f1 = ROOT.TFile(args.output_path, "RECREATE")
h_ll.Write()
h_llg.Write()
f1.Close

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
    if not t.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL:
        continue
    if t.nll < 1 or t.nphoton < 1:
        continue
    if all(dr < 0.4 for dr in t.photon_drmin):
        continue
    if t.ll_charge.size()==1 and t.ll_charge[0]!=0:
        continue
    if sum(map(bool, t.el_id)) < 2:
        continue
        
    idll = []
    for z in range(t.nll):
        x = t.ll_i1[z]
        y = t.ll_i2[z]
        if (t.el_charge[x]*t.el_charge[y] == -1 and \
            t.el_pt[x] > 25 and t.el_dz[x] < 0.01 and t.el_dxy[x] < 0.005 and bool(t.el_id[x]) and \
            t.el_pt[y] > 15 and t.el_dz[y] < 0.01 and t.el_dxy[y] < 0.005 and bool(t.el_id[y])):
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

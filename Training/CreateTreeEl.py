import ROOT
from array import array
import argparse
import WeightLumi
from math import sqrt, pi, exp, cos
import os
from RooFunctorFromWS import FunctorFromMVA

var_d_star = ['photon_mva', 'dR', 'pt_mass', 'cosTheta', 'costheta1', 'photon_res', 'photon_rapidity', 'el1_rapidity', 'el2_rapidity']
xml_name = os.path.join(os.getcwd(), "dataset/weights/TMVAClassification_BDT.weights.xml")
functor = FunctorFromMVA('BDT method', xml_name, *var_d_star)

def var_d(photon_mva, dR, pt_mass, cosTheta, costheta1, photon_res, photon_rapidity, el1_rapidity, el2_rapidity):
    return {'photon_mva' : photon_mva, 'dR' : dR, 'pt_mass' : pt_mass, 'cosTheta' : cosTheta, 'costheta1' : costheta1, 'photon_res' : photon_res, 'photon_rapidity' : photon_rapidity, 'el1_rapidity' : el1_rapidity, 'el2_rapidity' : el2_rapidity}

parser = argparse.ArgumentParser()
parser.add_argument('input_path')
parser.add_argument('output_path')
args = parser.parse_args()

ROOT.ROOT.EnableImplicitMT()

f = ROOT.TFile(args.input_path)
t = f.Get("tree")

weightL = WeightLumi.weight_lumi(args.input_path)

f1 = ROOT.TFile(args.output_path, "RECREATE")
tree1 = ROOT.TTree("opttree", "opttree")
mva = array('d', [0])
tree1.Branch('mva', mva, 'mva/D')
llphoton_m = array('d', [0])
tree1.Branch('llphoton_m', llphoton_m, 'llphoton_m/D')
ll_m = array('d', [0])
tree1.Branch('ll_m', ll_m, 'll_m/D')
weight = array('d', [0])
tree1.Branch('weight', weight, 'weight/D')

for i in range(0, t.GetEntries()):
    if i%1000==0:
       print i
    t.GetEntry(i)
    if not t.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL:
        continue
    if t.nll < 1 or t.nphoton < 1 or t.nel < 2 or t.nmu > 0:
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
        if (t.el_charge[x]*t.el_charge[y] == -1 and t.el_pt[x] > 25 and t.el_pt[y] > 15 and t.el_id[x] and t.el_id[y] \
            and abs(t.el_dz[x]) < 1 and abs(t.el_dxy[x]) < 0.5 and abs(t.el_dz[y]) < 1 and abs(t.el_dxy[y]) < 0.5):
            idll.append(z)
            el1_rapidity = t.el_eta[x]
            el2_rapidity = t.el_eta[y]
            break
    if len(idll) == 0:
        continue
    massZ = t.ll_m[idll[0]]

    if t.ll_dr[idll[0]] < 0.4:
        continue

    idllg = []
    iph = 0
    ph = ROOT.TLorentzVector()
    for z in range(t.nllphoton):
        if t.llphoton_ill[z]!=idll[0]:
            continue
        if (t.photon_drmin[iph] > 0.4 and t.photon_pt[iph] > 15 and \
            bool(bool(abs(t.photon_eta[iph]) < 1.4442 and t.photon_idmva[iph] > -0.4) or \
                 bool(1.566 < abs(t.photon_eta[iph]) < 2.5 and t.photon_idmva[iph] > -0.58))):
            ph.SetPtEtaPhiM(t.photon_pt[iph], t.photon_eta[iph], t.photon_phi[iph], 0)
            photon_mva = t.photon_idmva[iph]
            photon_rapidity = t.photon_eta[iph]
            photon_res = t.photon_pterr[iph]/t.photon_pt[iph]
            idllg.append(z)
            break
        iph = iph + 1
    if len(idllg) == 0:
        continue
    massH = t.llphoton_m[idllg[0]]
    ptH = t.llphoton_pt[idllg[0]]

    if massZ < 50:
        continue
    if (massZ + massH) < 185:
        continue
    if (ph.E()/massH) < 15/110:
        continue
    if massH < 100 or massH > 180:
        continue

    pt_mass = ptH/massH
    cosTheta = t.llphoton_cosTheta[idllg[0]]
    costheta1 = t.llphoton_costheta1[idllg[0]]
    dR = t.llphoton_dr[idllg[0]]

    mva[0] = functor(**var_d(photon_mva, dR, pt_mass, cosTheta, costheta1, photon_res, photon_rapidity, el1_rapidity, el2_rapidity))
    llphoton_m[0] = massH
    ll_m[0] = massZ

    if "Double" in args.input_path:
        weight[0] = weightL
    else:
        weight[0] = weightL*t.weight

    tree1.Fill()

f1.Write()

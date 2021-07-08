import ROOT
import rootnotes
c1 = rootnotes.default_canvas()

# Load file and print
f = ROOT.TFile("testWS.root")
f.ls()

# Load Workspace
ws = f.Get("CMS_emu_workspace")
ws.Print()

# Variable and Data
mass = ROOT.RooRealVar("CMS_emu_Mass","CMS_emu_Mass", 100, 180)
data = ws.data("Data_13TeV_vbf")

# Plotting
xframe = mass.frame()
data.plotOn(xframe)
xframe.Draw()
c1

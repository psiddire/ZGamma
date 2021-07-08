import ROOT
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('input_file')
args = parser.parse_args()

f = ROOT.TFile(args.input_file)
t = f.Get("opttree")
print t.Print()

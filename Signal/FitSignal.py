import ROOT
ROOT.TGaxis.SetMaxDigits(3)

def PlotFit(fit, data, pdf, m, nbins, xlow, xhigh, outfile, axis_name, doCovariance = False):
    gErrorIgnoreLevel = ROOT.kWarning
    # Declare canvas with 2 pads
    can = ROOT.TCanvas("can", "canvas", 800, 800)
    can.cd()
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1)
    pad1.Draw()
    pad1.cd()
    pad1.SetFillColor(0)
    pad1.SetBorderMode(0)
    pad1.SetBorderSize(10)
    pad1.SetTickx(1)
    pad1.SetTicky(1)
    pad1.SetLeftMargin(0.14)
    pad1.SetRightMargin(0.05)
    pad1.SetTopMargin(0.122)
    pad1.SetBottomMargin(0.026)
    pad1.SetFrameFillStyle(0)
    pad1.SetFrameLineStyle(0)
    pad1.SetFrameLineWidth(3)
    pad1.SetFrameBorderMode(0)
    pad1.SetFrameBorderSize(10)

    plot = m.frame(ROOT.RooFit.Range(xlow, xhigh))
    # Use Poisson errors if specified
    # Otherwise, will use errors from data histogram
    if "Poisson" in outfile:
        data.plotOn(plot, ROOT.RooFit.Binning(nbins), DataError(ROOT.Poisson))
    else:
        data.plotOn(plot, ROOT.RooFit.Binning(nbins))
    pdf.plotOn(plot, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineWidth(2))
    # Add reduced chi-squared to parameter legend by making it the title
    chi2Label = "#chi^{2}_{R} = " + str(plot.chiSquare(fit.floatParsFinal().getSize()))
    pdf.paramOn(plot, data, chi2Label, 2, "NELU", 0.6, 0.85, 0.85)

    plot.GetXaxis().SetTitle("")
    plot.GetXaxis().SetTitleSize(0)
    plot.GetXaxis().SetLabelSize(0)
    plot.GetXaxis().SetNdivisions(505)
    plot.GetYaxis().SetLabelFont(42)
    plot.GetYaxis().SetLabelOffset(0.01)
    plot.GetYaxis().SetLabelSize(0.06)
    plot.GetYaxis().SetTitleSize(0.07)
    plot.GetYaxis().SetTitleOffset(0.8)
    plot.SetMarkerStyle(20)
    plot.SetMarkerSize(1)
    plot.SetLineWidth(1)
    plot.GetXaxis().SetLimits(xlow, xhigh)
    # Retrieve function type and fit category
    funcType = pdf.GetTitle()
    category = data.GetTitle()
    plot.SetTitle(category+" fit with "+funcType+" (status = "+str(fit.status())+", covQual = "+str(fit.covQual())+")")
    plot.SetTitleSize(0.15)
    plot.Draw("same")
    pad1.RedrawAxis()

    can.cd()
    pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.3)
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(0.14)
    pad2.SetRightMargin(0.05)
    pad2.SetTickx(1)
    pad2.SetTicky(1)
    pad2.SetFrameLineWidth(3)
    pad2.SetGridy()
    pad2.Draw()
    pad2.cd()

    pull = plot.pullHist()
    # Set plot options and draw pulls
    pull.SetMarkerStyle(20)
    pull.GetXaxis().SetLimits(xlow, xhigh)
    pull.SetTitle(";" + axis_name + " [GeV]; Pull")
    pull.GetXaxis().SetLabelSize(0.08)
    pull.GetYaxis().SetLabelSize(0.08)
    pull.GetXaxis().SetNdivisions(505)
    pull.GetYaxis().SetNdivisions(5)
    pull.GetXaxis().SetTitleSize(0.16)
    pull.GetYaxis().SetTitleSize(0.16)
    pull.GetYaxis().SetTitleOffset(0.35)
    pull.GetXaxis().SetTitleOffset(1.04)
    pull.GetXaxis().SetLabelSize(0.11)
    pull.GetYaxis().SetLabelSize(0.11)
    pull.GetXaxis().SetTitleFont(42)
    pull.GetYaxis().SetTitleFont(42)
    pull.Draw("AP same")

    can.cd()
    pad1.Draw()

    ROOT.gPad.RedrawAxis()

    # Make output directory and save the plot
    can.SaveAs(outfile+".pdf");
    return

def makepdf(hm, outpath="signal"):
    mass = ROOT.RooRealVar("CMS_zgamma_Mass", "m_{ee#gamma}", 110, 140, "GeV")
    a1_dcb = ROOT.RooRealVar("a1", "a1", 1., 0.5, 2)
    a2_dcb = ROOT.RooRealVar("a2", "a2", 1., 0.5, 2)
    dm_dcb = ROOT.RooRealVar("dm", "dm", -0.1, -.5, 0.)
    mean_dcb = ROOT.RooFormulaVar("mean", "mean", "(125+@0)", ROOT.RooArgList(dm_dcb))
    n1_dcb = ROOT.RooRealVar("n1", "n1", 30., 10., 50.)
    n2_dcb = ROOT.RooRealVar("n2", "n2", 30., 10., 50.)
    sigma_dcb = ROOT.RooRealVar("sigma", "sigma", 2, 0.1, 2.5)
    pdf = ROOT.RooDoubleCBFast("pdf", "pdf", mass, mean_dcb, sigma_dcb, a1_dcb, n1_dcb, a2_dcb, n2_dcb)
    dh = ROOT.RooDataHist("data", "data", ROOT.RooArgList(mass), ROOT.RooFit.Import(hm))
    fitResult = pdf.fitTo(dh, ROOT.RooFit.Minimizer("Minuit2", "minimize"), \
                          ROOT.RooFit.Save(1), ROOT.RooFit.Range(110, 140), \
                          ROOT.RooFit.SumW2Error(ROOT.kTRUE))
    nbin = 30
    xlow = 110
    xhigh = 140
    print a1_dcb.getVal()
    print a2_dcb.getVal()
    print mean_dcb.getVal()
    print n1_dcb.getVal()
    print n2_dcb.getVal()
    print sigma_dcb.getVal()
    print fitResult.minNll()
    PlotFit(fitResult, dh, pdf, mass, nbin, xlow, xhigh, outpath, "m_{ll#gamma}", True)
    return

def makepdfZ(hm, outpath="signal"):
    mass = ROOT.RooRealVar("CMS_z_Mass", "m_{ee}", 70, 110, "GeV")
    a1_dcb = ROOT.RooRealVar("a1", "a1", 1., 0.5, 2)
    a2_dcb = ROOT.RooRealVar("a2", "a2", 1., 0.5, 2)
    dm_dcb = ROOT.RooRealVar("dm", "dm", -0.1, -.5, 0.)
    mean_dcb = ROOT.RooFormulaVar("mean", "mean", "(91+@0)", ROOT.RooArgList(dm_dcb))
    n1_dcb = ROOT.RooRealVar("n1", "n1", 30., 10., 50.)
    n2_dcb = ROOT.RooRealVar("n2", "n2", 30., 10., 50.)
    sigma_dcb = ROOT.RooRealVar("sigma", "sigma", 2, 0.1, 2.5)
    pdf = ROOT.RooDoubleCBFast("pdf", "pdf", mass, mean_dcb, sigma_dcb, a1_dcb, n1_dcb, a2_dcb, n2_dcb)
    dh = ROOT.RooDataHist("data", "data", ROOT.RooArgList(mass), ROOT.RooFit.Import(hm))
    fitResult = pdf.fitTo(dh, ROOT.RooFit.Minimizer("Minuit2", "minimize"), \
                          ROOT.RooFit.Save(1), ROOT.RooFit.Range(70, 110), \
                          ROOT.RooFit.SumW2Error(ROOT.kTRUE))
    nbin = 40
    xlow = 70
    xhigh = 110
    print a1_dcb.getVal()
    print a2_dcb.getVal()
    print mean_dcb.getVal()
    print n1_dcb.getVal()
    print n2_dcb.getVal()
    print sigma_dcb.getVal()
    print fitResult.minNll()
    PlotFit(fitResult, dh, pdf, mass, nbin, xlow, xhigh, outpath, "m_{ll}", True)
    return


def makeSumOfGaussians(hm):

    mass = ROOT.RooRealVar("CMS_zgamma_Mass", "m_{ee#gamma}", 110, 140, "GeV")

    dh = ROOT.RooDataHist("data", "data", ROOT.RooArgList(mass), ROOT.RooFit.Import(hm))

    mh = 125.0

    pdfs = ROOT.RooArgList()

    # massHypothesis + deltaM
    expr = "%f + @0" % mh

    dmuvars0 = ROOT.RooRealVar("dmu_0", "delta mu 0", 0, -2.5, +2.5)
    meanVar0 = ROOT.RooFormulaVar("mu_0", "mean Gaussian 0", expr, ROOT.RooArgList(dmuvars0))
    sigmaVar0 = ROOT.RooRealVar("sigma_0", "sigma Gaussian 0", 1, 0.01, 10)
    pdf0 = ROOT.RooGaussian("sigpdf_0", "Gaussian 0", mass, meanVar0, sigmaVar0)

    dmuvars1 = ROOT.RooRealVar("dmu_1", "delta mu 1", 1, -2.5, +2.5)
    meanVar1 = ROOT.RooFormulaVar("mu_1", "mean Gaussian 1", expr, ROOT.RooArgList(dmuvars1))
    sigmaVar1 = ROOT.RooRealVar("sigma_1", "sigma Gaussian 1", 2, 0.01, 10)
    pdf1 = ROOT.RooGaussian("sigpdf_1", "Gaussian 1", mass, meanVar0, sigmaVar1)

    mix = ROOT.RooRealVar("frac", "fraction variable for Gaussian sum", 0.5, 0, 1)
    pdf = ROOT.RooAddPdf("sigpdf_2", "sigpdf_2", ROOT.RooArgList(pdf0, pdf1), \
                         ROOT.RooArgList(mix), ROOT.kTRUE)

    fitResult = pdf.fitTo(dh, ROOT.RooFit.Minimizer("Minuit2", "minimize"),
                          ROOT.RooFit.Save(1), ROOT.RooFit.Range(110, 140),
                          ROOT.RooFit.SumW2Error(ROOT.kTRUE))
    nbin = 30
    xlow = 110
    xhigh = 140
    outpath = "signal_gauss"
    FitCategory.PlotFit(fitResult, dh, pdf, mass, nbin, xlow, xhigh, outpath, True)

    return

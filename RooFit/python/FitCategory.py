import ROOT
ROOT.gSystem.Load('ModGaus_cxx.so')
ROOT.gSystem.Load('ModGaus01_cxx.so')
ROOT.gSystem.Load('ModGaus10_cxx.so')
ROOT.gSystem.Load('ModGaus11_cxx.so')
ROOT.gSystem.Load('ModGaus21_cxx.so')
ROOT.TGaxis.SetMaxDigits(3)

def PlotFit(fit, data, pdf, m, nbins, xlow, xhigh, outfile, doCovariance = False):
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
    pull.SetTitle(";m_{ll#gamma} [GeV]; Pull")
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
    ROOT.gSystem.mkdir(outfile);
    can.SaveAs(outfile+".pdf");
    return


def PlotCovariances(fit, outpath):
    ROOT.gStyle.SetOptStat(0)
    gErrorIgnoreLevel = ROOT.kWarning
    pars = fit.floatParsFinal()
    can = ROOT.TCanvas("can", "canvas", 800, 800)
    can.cd()
    for ip1 in range(pars.getSize()):
        for ip2 in range(ip1+1, pars.getSize()):
            p1 = pars.find(pars.at(ip1).namePtr().GetName())
            p2 = pars.find(pars.at(ip2).namePtr().GetName())
            plot = ROOT.RooPlot(p1, p2, p1.getVal()-1.1*p1.getError(),
                                p1.getVal()+1.1*p1.getError(),
                                p2.getVal()-1.1*p2.getError(),
                                p2.getVal()+1.1*p2.getError())
            fit.plotOn(plot, p1, p2, "ME12ABHV")
            plot.GetYaxis().SetLabelSize(0.025)
            plot.GetXaxis().SetLabelSize(0.025)
            plot.Draw()
            outfile = outpath + "/Covariance_p" + str(ip1) + "-p" + str(ip2) + ".pdf"
            can.SaveAs(outfile)
            can.Clear()
    cor = fit.correlationHist()
    cor.GetYaxis().SetLabelSize(0.1)
    cor.GetXaxis().SetLabelSize(0.1)
    cor.Draw("colz")
    can.SaveAs(outpath+"/Correlation.pdf")
    can.Close()
    cor.Delete()
    return


def PlotLikelihood(fit, model, data, outpath):
    gErrorIgnoreLevel = ROOT.kWarning
    pars = fit.floatParsFinal()
    can = ROOT.TCanvas("can", "canvas", 800, 800)
    can.cd()
    for ip1 in range(pars.getSize()):
        par = pars.find(pars.at(ip1).namePtr().GetName())
        nll = model.createNLL(data)
        title = "NLL-scan-for-"+ str(ip1)
        frame = par.frame(ROOT.RooFit.Bins(10), ROOT.RooFit.Range(par.getVal()-3*par.getError(),par.getVal()+3*par.getError()), ROOT.RooFit.Title(title))
        nll.plotOn(frame, ROOT.RooFit.ShiftToZero())
        frame.Draw()
        outfile = outpath+"/"+title+".pdf"
        can.SaveAs(outfile)
        can.Clear()
    return


def FitCategory(category = "untag4", fitType = "ModGaus", nbin = 80, xlow = 100, xhigh = 180, useSumW2 = True, defaults = []):
    outpath = "plots/CategoryFits/" + fitType + str(xlow) + "-" + str(xhigh) + "_" + str(nbin)
    if not useSumW2:
        outpath = outpath + "_Poisson"
    outpath = outpath + "/"+category
    ROOT.gSystem.mkdir(outpath, True)
    ROOT.gSystem.RedirectOutput(outpath+"/log.txt", "w")

    if defaults == []:
        defaults = [120, 2.5, 0, 10, 10, 50, 0, 0, 50]

    m = ROOT.RooRealVar("llphoton_m", "ll#gamma mass [GeV]", xlow, xhigh)
    m0 = ROOT.RooRealVar("m_{0}", "mass peak value [GeV]", defaults[0], 105, 145)
    m0.setError(1.0)
    vl = ROOT.RooRealVar("#nu_{L}", "low-end power", defaults[1], 0, 5)
    vl.setError(0.1)
    vr = ROOT.RooRealVar("#Delta#nu", "power range", defaults[2], -5, 5)
    vr.setError(0.1)
    s0 = ROOT.RooRealVar("#sigma_{0}", "peak width", defaults[3], 1, 40)
    s0.setError(1.0)
    sl = ROOT.RooRealVar("#sigma_{L}", "low-end width", defaults[4], 1, 40)
    sl.setError(1.0)
    sh = ROOT.RooRealVar("#sigma_{H}", "high-end width", defaults[5], 20, 100)
    sh.setError(1.0)
    s1 = ROOT.RooRealVar("#sigma_{1}", "1st order width", defaults[6], 1, 40)
    s1.setError(1.0)
    s2 = ROOT.RooRealVar("#sigma_{2}", "2nd order width", defaults[7], -50, 50)
    s2.setError(1.0)
    v0 = ROOT.RooRealVar("#nu_{0}", "0th order power", 2.5, 0, 5)
    v0.setError(0.1)
    v1 = ROOT.RooRealVar("#nu_{1}", "1st order power", 2.5, 0, 5)
    v1.setError(0.1)

    if fitType == "ModGaus21":
        background = ROOT.ModGaus21("background", "G_{2, 1}", m, m0, vl, vr, s0, s1, s2)
    elif fitType == "ModGaus11":
        background = ROOT.ModGaus11("background", "G_{1, 1}", m, m0, vl, vr, s0, s1)
    elif fitType == "ModGaus10":
        background = ROOT.ModGaus10("background", "G_{1, 0}", m, m0, vl, vr, s0)
    elif fitType == "ModGaus01":
        background = ROOT.ModGaus01("background", "G_{0, 1}", m, m0, vl, s0, s1)
    else:
        background = ROOT.ModGaus("background", "G_{M}", m, m0, vl, vr, s0, sl, sh)

    inpath = "CatFiles/"
    inpath = inpath + "bkghist_ele_mu_"+category+".root"
    infile = ROOT.TFile(inpath)

    hmc = infile.Get("hbkg")

    if nbin==40:
        hmc.Rebin()

    if not useSumW2:
        for ib in range(nbin):
            if(hmc.GetBinContent(ib) < 0):
                hmc.SetBinContent(ib, 0)

    mc = ROOT.RooDataHist(category, category, ROOT.RooArgList(m), ROOT.RooFit.Import(hmc))

    fitResult = background.fitTo(mc, ROOT.RooFit.Minimizer("Minuit2", "Migrad"), ROOT.RooFit.Save(1), ROOT.RooFit.Range(xlow, xhigh), ROOT.RooFit.SumW2Error(ROOT.kTRUE))

    PlotFit(fitResult, mc, background, m, nbin, xlow, xhigh, outpath, True)

    PlotCovariances(fitResult, outpath)
    
    PlotLikelihood(fitResult, background, mc, outpath)

    return

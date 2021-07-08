def weight_lumi(fName):
    if "GluGluHToZG" in fName:
        entries = 400000
        xsec = 48.58
    elif "VBFHToZG" in fName:
        entries = 200000
        xsec = 3.782
    elif "WminusH" in fName:
        entries = 299276
        xsec = 0.533
    elif "WplusH" in fName:
        entries = 299978
        xsec = 0.840
    elif "ZH" in fName:
        entries = 297389
        xsec = 0.884
    #elif "ttH" in fName:
    #    entries = 200000
    #    xsec = ?
    else:
        raise ValueError('Not a valid input file')
    equi_lumi = entries/xsec
    lumi = 41500
    weight = lumi/equi_lumi
    return weight 

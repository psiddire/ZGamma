def weight_lumi(fName):
    if "DYJetsToLL" in fName:
        entries = 102486448
        xsec = 6025.0
    elif "ZGToLLG" in fName:
        entries = 29885702
        xsec = 55.48
    elif "TTGJets" in fName:
        entries = 3534208
        xsec = 4.078
    elif "WW" in fName:
        entries = 15883000
        xsec = 118.7
    elif "WZ" in fName:
        entries = 7898000
        xsec = 51.11
    elif "ZZ" in fName:
        entries = 2708000
        xsec = 16.91
    elif "TTJets" in fName:
        entries = 248504443
        xsec = 831.76
    elif "WJets" in fName:
        entries = 81254459
        xsec = 61526.7
    elif "GluGluHToZG" in fName:
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
    elif "Double" in fName:
        return 1.0
    else:
        raise ValueError('Not a valid input file')
    equi_lumi = entries/xsec
    lumi = 41500
    weight = lumi/equi_lumi
    return weight

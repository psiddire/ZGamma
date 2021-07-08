def weight_lumi(fName):
    if "DYJetsToLL" in fName:
        entries = 102486448
        xsec = 6025
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
    else:
        raise ValueError('Not a valid input file')
    equi_lumi = entries/xsec
    lumi = 41500
    weight = lumi/equi_lumi
    return weight 

import ROOT as r
import fedrarootlogon
import numpy as np 
from matplotlib import pyplot as plt
import time

IDBRICK = 333
S0, SL = 1, 2
off_name = str(IDBRICK) + "_S" + str(S0) + "_S" + str(SL) + "_offsets.root"

fileIn = r.TFile(off_name, "READ")
t_merge = fileIn.Get("t_merge")

### CALCULATING OFFSETS
offName2 = str(IDBRICK) + "_S" + str(S0) + "_S" + str(SL) + "_offsets_fit.root"
outFile2 = r.TFile(offName2, "RECREATE")
off_tup = r.TNtuple("merge_offsets", "", "Xoff:Yoff:TXoff:TYoff")

X_offset, Y_offset, TX_offset, TY_offset = 0, 0, 0, 0

# Create Histograms
t_merge.Draw("DXp>>h_DXp(600, -300, 300")
h_DXp = r.gDirectory.Get("h_DXp")
t_merge.Draw("DYp>>h_DYp(600, -300, 300")
h_DYp = r.gDirectory.Get("h_DYp")
t_merge.Draw("DTX>>h_DTX(1000, -.1, .1")
h_DTX = r.gDirectory.Get("h_DTX")
t_merge.Draw("DTY>>h_DTY(1000, -.1, .1")
h_DTY = r.gDirectory.Get("h_DTY")

# DX
bin_max = h_DXp.GetMaximumBin()
binc = h_DXp.GetXaxis().GetBinCenter(bin_max)
g1 = r.TF1("g1", "gaus(0)", binc-50, binc+50)
g1.SetParameter(0, 100)
g1.SetParameter(1, binc)
g1.SetParameter(2, 5)
h_DXp.Fit("g1", "L","", binc-10, binc+10)
fit_prob = g1.GetProb()
mean = g1.GetParameter(1)
if (fit_prob>0.0001):
    print(" Used Fit Value " + str(mean) + " instead of max value " + str(binc))
    X_offset = float(mean)
else:
    print(" WARNING: low fit probability " + str(fit_prob) + " using binc " + str(binc))
    X_offset = float(binc)
#c_DXp.SaveAs("merge_off/c_DXp.C")

# DY
bin_max = h_DYp.GetMaximumBin()
binc = h_DYp.GetXaxis().GetBinCenter(bin_max)
g1 = r.TF1("g1", "gaus(0)", binc-50, binc+50)
g1.SetParameter(0, 100)
g1.SetParameter(1, binc)
g1.SetParameter(2, 5)
h_DYp.Fit("g1", "L","", binc-10, binc+10)
fit_prob = g1.GetProb()
mean = g1.GetParameter(1)
if (fit_prob>0.0001):
    print(" Used Fit Value " + str(mean) + " instead of max value " + str(binc))
    Y_offset = float(mean)
else:
    print(" WARNING: low fit probability " + str(fit_prob) + " using binc " + str(binc))
    Y_offset = float(binc)
#c_DYp.SaveAs("merge_off/c_DYp.C")

# DTX
bin_max = h_DTX.GetMaximumBin()
binc = h_DTX.GetXaxis().GetBinCenter(bin_max)
g1 = r.TF1("g1", "gaus(0)", binc-0.02, binc+0.02)
g1.SetParameter(0, 100)
g1.SetParameter(1, binc)
g1.SetParameter(2, .01)
h_DTX.Fit("g1", "L","", binc-0.01, binc+0.01)
fit_prob = g1.GetProb()
mean = g1.GetParameter(1)
if (fit_prob>0.0001):
    print(" Used Fit Value " + str(mean) + " instead of max value " + str(binc))
    TX_offset = float(mean)
else:
    print(" WARNING: low fit probability " + str(fit_prob) + " using binc " + str(binc))
    TX_offset = float(binc)
#c_DTX.SaveAs("merge_off/c_DTX.C")

# DTY
bin_max = h_DTY.GetMaximumBin()
binc = h_DTY.GetXaxis().GetBinCenter(bin_max)
g1 = r.TF1("g1", "gaus(0)", binc-0.01, binc+0.01)
g1.SetParameter(0, 100)
g1.SetParameter(1, binc)
g1.SetParameter(2, .01)
h_DTY.Fit("g1", "L","", binc-0.01, binc+0.01)
fit_prob = g1.GetProb()
mean = g1.GetParameter(1)
if (fit_prob>0.0001):
    print(" Used Fit Value " + str(mean) + " instead of max value " + str(binc))
    TY_offset = float(mean)
else:
    print(" WARNING: low fit probability " + str(fit_prob) + " using binc " + str(binc))
    TY_offset = float(binc)
#c_DTY.SaveAs("merge_off/c_DTY.C")

off_tup.Fill(X_offset, Y_offset, TX_offset, TY_offset)
outFile2.cd()

c_DXp = r.TCanvas()
c_DXp.cd()
h_DXp.SetTitle("DXp;DXp[#mum];Entries")
h_DXp.Draw("COLZ")
c_DXp.Write("c_DXp")

c_DYp = r.TCanvas()
c_DYp.cd()
h_DYp.SetTitle("DYp;DYp[#mum];Entries")
h_DYp.Draw("COLZ")
c_DYp.Write("c_DYp")

c_DTX = r.TCanvas()
c_DTX.cd()
h_DTX.SetTitle("DTX;DTX[mrad];Entries")
h_DTX.Draw("COLZ")
c_DTX.Write("c_DTX")

c_DTY = r.TCanvas()
c_DTY.cd()
h_DTY.SetTitle("DTY;DTY[mrad];Entries")
h_DTY.Draw("COLZ")
c_DTY.Write("c_DTY")

off_tup.Write()
outFile2.Close()
fileIn.Close()

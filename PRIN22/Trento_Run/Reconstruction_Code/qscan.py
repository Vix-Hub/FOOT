# script to have a first look at reconstructed grains 

import ROOT as r 
import fedrarootlogon 

inputFile = r.TFile("dm_grains.dm.root", "READ")
Vdmr = inputFile.Get("Vdmr")

c = r.TCanvas()
c.Divide(3,3)

# XYZ Distributions

c.cd(1)
Vdmr.Draw("gr.y+y:gr.x+x", "", "colz")

c.cd(2)
Vdmr.Draw("gr.z+z:gr.x+x", "", "colz")

c.cd(3)
Vdmr.Draw("gr.z+z:gr.y+y", "", "colz")

c.cd(4)
Vdmr.Draw("gr.ncl", "")

c.cd(5)
Vdmr.Draw("gr.y+y:gr.x+x", "gr.ncl>=3", "colz")

c.SaveAs("qscan.pdf")


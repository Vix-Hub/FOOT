## script to add Nincoming branch to vtx tree
import ROOT as r 
import fedrarootlogon 
import numpy as np 

VtxFileName = "vertices_improved_fast_3_new_temp_2.root"
VtxFile = r.TFile(VtxFileName, "READ")
vtx = VtxFile.Get("vtx")

OutFileName = "vertices_incoming.root"
OutFile = r.TFile(OutFileName, "RECREATE")

vtx2 = vtx.CloneTree(0)

n_outgoing = np.zeros(1, np.intc)
vtx2.Branch("n_outgoing", n_outgoing, "n_outgoing/I")

for i in range(vtx.GetEntries()):
    vtx.GetEntry(i)
    n_temp = 0
    for j in range(vtx.n):
        if (vtx.incoming[j] == 1):
            n_temp += 1
    n_outgoing[0] = n_temp
    vtx2.Fill()

OutFile.cd()
vtx2.Write("vtx")
OutFile.Close()
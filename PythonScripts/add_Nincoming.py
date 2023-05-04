## script to add Nincoming branch to vtx tree
import ROOT as r 
import fedrarootlogon 
import numpy as np 

VtxFileName = "vertices_improved_fast_3_new_temp_2.root"
VtxFile = r.TFile(VtxFileName, "READ")
vtx = VtxFile.Get("vtx")
vrec = VtxFile.Get("EdbVertexRec")

OutFileName = "vertices_incoming.root"
OutFile = r.TFile(OutFileName, "RECREATE")

vtx2 = vtx.CloneTree(0)

n_outgoing = np.zeros(1, np.intc)
n_S2 = np.zeros(1, np.intc)
vtx2.Branch("n_outgoing", n_outgoing, "n_outgoing/I")
vtx2.Branch("n_S2", n_S2, "n_S2/I")

for i in range(vtx.GetEntries()):
    vtx.GetEntry(i)
    n_temp, n_temp_s2 = 0, 0
    for j in range(vtx.n):
        if (vtx.incoming[j] == 1):
            n_temp += 1
        if (vtx.npl[i] + vtx.plate[i] > 31):
            n_temp_s2 += 1
    n_outgoing[0] = n_temp
    n_S2[0] = n_temp_s2
    vtx2.Fill()

OutFile.cd()
vtx2.Write("vtx")
vrec.Write("EdbVertexRec")
OutFile.Close()
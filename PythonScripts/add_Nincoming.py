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
n_10 = np.zeros(1, np.intc)
n_5 = np.zeros(1, np.intc)
nseg5 = np.zeros(1, np.intc)
nseg10 = np.zeros(1, np.intc)

vtx2.Branch("n_outgoing", n_outgoing, "n_outgoing/I")
vtx2.Branch("n_S2", n_S2, "n_S2/I")
vtx2.Branch("n_10", n_10, "npl10/I")
vtx2.Branch("n_5", n_5, "npl5/I")
vtx2.Branch("nseg_5", nseg5, "nseg5/I")
vtx2.Branch("nseg_10", nseg10, "nseg10/I")

for i in range(vtx.GetEntries()):
    vtx.GetEntry(i)
    vertex = vrec.eVTX.At(i)
    n_temp, n_temp_s2 = 0, 0
    npl5temp, npl10temp = 0, 0
    nseg5temp, nseg10temp = 0, 0
    for j in range(vtx.n):
        if (vtx.incoming[j] == 1):
            n_temp += 1
        if (vtx.npl[j] + vtx.plate[j] > 31):
            n_temp_s2 += 1
        track = vertex.GetTrack(j)
        npl = track.GetSegmentLast().Plate() - track.GetSegmentFirst().Plate()
        nseg = track.N()
        if (npl>=5):
            npl5temp += 1
        if (npl>=10):
            npl10temp += 1
        if (nseg>=5):
            nseg5temp += 1
        if (nseg>=10):
            nseg10temp += 1
    n_outgoing[0] = n_temp
    n_S2[0] = n_temp_s2
    n_5[0] = npl5temp
    n_10[0] = npl10temp
    nseg5[0] = nseg5temp
    nseg10[0] = nseg10temp
    vtx2.Fill()

OutFile.cd()
vtx2.Write("vtx")
vrec.Write("EdbVertexRec")
OutFile.Close()
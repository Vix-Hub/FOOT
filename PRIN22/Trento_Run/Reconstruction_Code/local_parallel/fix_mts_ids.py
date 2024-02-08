# script to add branch with unique ID to each mt
import ROOT as r 
import numpy as np

MergedFile = r.TFile("mt.merged.all.root", "READ")
mtracks = MergedFile.Get("mtracks")

MergedFile2 = r.TFile("mt.merged.all.new.root", "RECREATE")
mtracks2 = mtracks.CloneTree(0)

id2 = np.zeros(1, dtype=np.intc)

mtracks2.Branch("id2", id2, "id2/I")
counter = 0
for i in range(mtracks.GetEntries()):
    mtracks.GetEntry(i)
    counter += 1
    id2[0] = counter
    mtracks2.Fill()

MergedFile2.cd()
mtracks2.Write("mtracks")
MergedFile2.Close()
import ROOT as r 
import fedrarootlogon 
import numpy as np 
import time

inputName = "vertices_improved_fast_3_new_2.root"
inputFile = r.TFile(inputName, "READ")
inputVtx = inputFile.Get("vtx")
inputVtx.BuildIndex("vID")
vrec = inputFile.Get("EdbVertexRec")

outputName = "vertices_improved_fast_3_new_2_MCr.root"
outputFile = r.TFile(outputName, "RECREATE")

IDBRICKMC = 2

# CREATE NEW OBJECTS
vrec2 = vrec
outVtx = inputVtx.CloneTree(0)

Z_rec = np.zeros(1000, dtype = np.intc)
#n2 = np.zeros(0, dtype=np.intc)
#t.Branch('mynum', n, 'mynum/I')
outVtx.Branch("Z_flag2", Z_rec, "Z_flag2[n]/I")

# CHARGE INFO
chargeName = "b" + str(IDBRICKMC) + "_vtx.root"
chargeFile = r.TFile(chargeName, "READ")
chargeTree = chargeFile.Get("tracks")

chargeTree.BuildIndex("s[0].ID()", "s[0].Plate()")

# LOOP ON VERTEX TRACKS
n_vertices = vrec2.eVTX.GetEntries()
Z = 0
t0 = time.time()
for ivtx in range(n_vertices):
	vertex = vrec2.eVTX.At(ivtx)
	charges = []
	for itrk in range(vertex.N()):
		track = vertex.GetTrack(itrk)
		skip, pl0, id0, S2_segs = 0, 0, 0, 0
		for iseg in range(track.N()):
			segp = track.GetSegment(iseg)
			segf = track.GetSegmentF(iseg)
			if (segp.Plate() >= 31 and skip == 0):
				pl0, id0 = segp.Plate(), segp.ID()
				skip = 1
		if (skip>0): #if track gets to S2
			chargeTree.GetEntryWithIndex(int(id0), int(pl0)) #find track in S2
			Z = int(chargeTree.s[0].W()-70)
			if (chargeTree.in_vtx == 1):
				charges.append(Z)
			for iseg in range(track.N()):
				track.GetSegment(iseg).SetFlag(Z)
				track.GetSegmentF(iseg).SetFlag(Z)
			track.SetFlag(Z) # saving Z to EdbVertexRec segments and tracks
		else:
			charges.append(0)
	#print(charges)
	inputVtx.GetEntryWithIndex(vertex.ID())
	for j, charge in enumerate(charges):
		Z_rec[j] = np.intc(charge)
	#print(Z_rec)
	outVtx.Fill() #saving new Z in Ttree
	if (ivtx%1000 == 0):
		print("Completed " + str(round(100*ivtx/n_vertices, 2)) + " %")
		print(" Time: " + str(round(time.time()-t0, 2)) + " s")

outputFile.cd()
outVtx.Write("vtx")
vrec2.Write()
outputFile.Close()
inputFile.Close()

'''
#### CHECK
print("")
print("--- Checking --- ")
inFile2 = r.TFile(outputName, "READ")
vtxrec = inFile2.Get("EdbVertexRec")

n_vertices = vtxrec.eVTX.GetEntries()
flags = []
flags2 = []
for ivtx in range(n_vertices):
	vertex = vtxrec.eVTX.At(ivtx)
	for itrk in range(vertex.N()):
		track = vertex.GetTrack(itrk)
		flags.append(track.GetSegment(0).Flag())
		flags2.append(track.Flag())

from matplotlib import pyplot as plt
plt.hist(flags)
plt.show()

plt.hist(flags2)
plt.show()
'''

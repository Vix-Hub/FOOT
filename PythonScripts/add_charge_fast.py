import ROOT as r 
import fedrarootlogon 
import numpy as np 
import time

#script to add charge branch to vtx tree 
IDBRICK = 2
MC = 0
if (IDBRICK<10):
	MC = 1

inputName = "vertices_incoming.root"
inputFile = r.TFile(inputName, "READ")
inputVtx = inputFile.Get("vtx")
inputVtx.BuildIndex("vID")
vrec = inputFile.Get("EdbVertexRec")

outputName = "vertices_final_Z.root"
outputFile = r.TFile(outputName, "RECREATE")

# CREATE NEW OBJECTS
vrec2 = vrec
outVtx = inputVtx.CloneTree(0)

Z_rec = np.zeros(1000, dtype = np.intc)
Z_sum_out = np.zeros(1, dtype=np.intc)
outVtx.Branch("Z_flag2", Z_rec, "Z_flag2[n]/I")
outVtx.Branch("Z_sum", Z_sum_out, "Z_sum/I")

# CHARGE INFO
chargeName = "b" + str(IDBRICK).zfill(6) + ".2.0.0.trk.root"
if (not MC):
	chargeName = "b" + str(IDBRICK).zfill(6) + ".2.0.0.trk_updateV.root"
chargeFile = r.TFile(chargeName, "READ")
chargeTree = chargeFile.Get("tracks")

charge_couples = []
charges_S2 = []
for track in chargeTree:
	charge_couples.append((track.s[0].ID(), track.s[0].Plate()))
	if (MC):
		charges_S2.append(track.s[0].W()-70)
	else:
		charges_S2.append(track.s[0].Flag())

print(" Saved Charge info from track file \n")

# LOOP ON VERTEX TRACKS
n_vertices = vrec2.eVTX.GetEntries()
Z = 0
t0 = time.time()
for ivtx in range(n_vertices):
	vertex = vrec2.eVTX.At(ivtx)
	zsum = 0
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
			position = charge_couples.index((int(id0), int(pl0))) #find track in S2
			Z = int(charges_S2[position])
			charges.append(Z)
			zsum += Z
			for iseg in range(track.N()):
				track.GetSegment(iseg).SetFlag(Z)
				track.GetSegmentF(iseg).SetFlag(Z)
			track.SetFlag(Z) # saving Z to EdbVertexRec segments and tracks
		else:
			charges.append(-2)
	#print(charges)
	inputVtx.GetEntryWithIndex(vertex.ID())
	for j, charge in enumerate(charges):
		Z_rec[j] = np.intc(charge)
	Z_sum_out[0] = zsum
	outVtx.Fill() #saving new Z in Ttree
	if (ivtx%1000 == 0):
		print("Completed " + str(round(100*ivtx/n_vertices, 2)) + " %")
		print(" Time: " + str(round(time.time()-t0, 2)) + " s")

outputFile.cd()
outVtx.Write("vtx")
vrec2.Write()
outputFile.Close()
inputFile.Close()


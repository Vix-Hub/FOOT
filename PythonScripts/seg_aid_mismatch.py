import ROOT as r 
import fedrarootlogon
import numpy as np 

VtxFile = r.TFile("vertices_improved_fast_3_new_temp_2.root", "READ")
vrec = VtxFile.Get("EdbVertexRec")
n_vertices = vrec.eVTX.GetEntries()

TrkFile = r.TFile("b000002.0.1.2.trk.root", "READ")
tracks = TrkFile.Get("tracks")

Vtx_Trks = []
Vtx_Trks_couples = []
mismatch, match0, tot = 0, 0, 0

for i in range(n_vertices):
    vertex = vrec.eVTX.At(i)
    for j in range(vertex.N()):
        track = vertex.GetTrack(j)
        Vtx_Trks.append((track.GetSegmentFirst().Aid(1), track.GetSegmentLast().Aid(1)))
        Vtx_Trks_couples.append((track.GetSegmentFirst().Plate(), track.GetSegmentFirst().ID()))

for track in tracks:
    tot += 1
    if (track.s[0]<30):
        continue
    if ((track.s[0].Plate(), track.s[0].ID()) in Vtx_Trks_couples):
        index = Vtx_Trks_couples.index((track.s[0].Plate(), track.s[0].ID()))
        if ( (track.s[0].Aid(1), track.s[track.nseg-1].Aid(1)) != Vtx_Trks[index] ):
            mismatch += 1 
        else:
            match0 += 1
    if (tot%10000 == 0):
        print(" Completed " + str(100.*tot/tracks.GetEntries()) + " % ")

print(" Found Mismatch in segments' Aid(1) for " + str(mismatch) + " Tracks, out of " + str(match0+mismatch) + " tracks connected to vertices ")

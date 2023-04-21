import ROOT as r 
import fedrarootlogon 
import numpy as np 

VtxFile1 = r.TFile("vertices_merged.root", "READ")
VtxFile2 = r.TFile("vertices_AnaFake.root", "READ")

vrec1 = VtxFile1.Get("EdbVertexRec")
n1 = vrec1.eVTX.GetEntries()

vrec2 = VtxFile2.Get("EdbVertexRec")
n2 = vrec2.eVTX.GetEntries()

long_tracks, initial_vIDs = [], []
kept_tracks = []

initial_N, init_X, init_Y, init_Z = [], [], [], []
vids, nlongs = [], []

outFile = r.TFile("deleted_vertices.root", "RECREATE")
tup = r.TNtuple("t_del", "info about vertices deleted in FAST=1 and that were in vertices_merged", "vID:N:vx:vy:vz:Nlong")

for i in range(n1):
    v1 = vrec1.eVTX.At(i)
    nlong = 0
    for j in range(v1.N()):
        track = v1.GetTrack(j)
        if (track.GetSegmentFirst().Plate() < 30 and track.GetSegmentLast().Plate()>30):
            long_tracks.append((track.GetSegmentFirst().Plate(), track.GetSegmentFirst().ID()))
            initial_vIDs.append(v1.ID())
            initial_N.append(v1.N())
            init_X.append(v1.X())
            init_Y.append(v1.Y())
            init_Z.append(v1.Z())
            nlong += 1
    vids.append(v1.ID())
    nlongs.append(nlong)

# find long tracks removed 

for i in range(n2):
    v2 = vrec2.eVTX.At(i)
    for j in range(v2.N()):
        track = v2.GetTrack(j)
        if (track.GetSegmentFirst().Plate() < 30 and track.GetSegmentLast().Plate()>30):
            kept_tracks.append((track.GetSegmentFirst().Plate(), track.GetSegmentFirst().ID()))

appended_vids = []

for j, long_track in enumerate(long_tracks):
    if (not long_track in kept_tracks):
        if (not (initial_vIDs[j] in appended_vids)):
            tup.Fill(initial_vIDs[j], initial_N[j], init_X[j], init_Y[j], init_Z[j], nlongs[vids.index(initial_vIDs[j])])
            appended_vids.append(initial_vIDs[j])

outFile.cd()
tup.Write()
outFile.Close()
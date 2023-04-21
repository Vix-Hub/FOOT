## checks which tracks that were not present in VertexFile1 are present in VertexFile2 
## some cuts are applied to the tracks (example: select only tracks merged between S1 and S2 for proveV)

import ROOT as r 
import fedrarootlogon
import numpy as np 
from debug_libraries import FillVertexMergedTracks

VertexFile1 = r.TFile("vertices_merged.root", "READ")
vrec1 = VertexFile1.Get("EdbVertexRec")

VertexFile2 = r.TFile("vertices_improved_fast_3_new_temp_2.root", "READ")
vrec2 = VertexFile2.Get("EdbVertexRec")

n_vertices1 = vrec1.eVTX.GetEntries()
n_vertices2 = vrec2.eVTX.GetEntries()

vtx_tracks_1 = FillVertexMergedTracks(vrec1)
vtx_tracks_2 = FillVertexMergedTracks(vrec2)

outFile = r.TFile("recovered_tracks.root", "RECREATE")
outTup = r.TNtuple("tup", "Info if track is present in both vertex files in find_recovered_tracks.py", "s0id:s0plate:is_rec")

for track in vtx_tracks_1:
    recovered = 0
    if (track in vtx_tracks_2):
        recovered = 1
    outTup.Fill(track[1], track[0], recovered)

outFile.cd()
outTup.Write()
outFile.Close()
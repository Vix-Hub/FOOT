import ROOT as r 
import fedrarootlogon
import numpy as np 

# finds merged tracks between S0 and SL not connected to vertices and saves info

IDBRICK = 2
S0, SL = 1, 2
S0_END_PLATE = 30

start_string = "{:06}".format(IDBRICK)
TrackFile = r.TFile( "b" + start_string + ".0." + str(S0) + "." + str(SL) + ".trk_merged_2.root", "READ")
tracks = TrackFile.Get("tracks")

VertexFile = r.TFile("vertices.root", "READ")
vrec = VertexFile.Get("EdbVertexRec")

outfile = r.TFile("missing_tracks_vertices.root", "RECREATE")
tup = r.TNtuple("t_missing", "Info Concerning Merged Tracks not connected to any vertices in vertices.root ", "nseg:npl:s0plate:s0id:s0X:s0Y:s0Z:s0TX:s0TY:s0_MCEvt:s0_MCtrk:s0flag:Z_first:Z_last")
aus_tup = r.TNtuple("t_missing2", "more info", "sLX:sLY:sLZ:sLTX:sLTY")

vertex_trk_couples = [] 

n_vertices = vrec.eVTX.GetEntries()

for i in range(n_vertices):
    vertex = vrec.eVTX.At(i)
    for j in range(vertex.N()):
        track = vertex.GetTrack(j)
        for iseg in range(vertex.N()):
            vertex_trk_couples.append((track.GetSegment(iseg).Plate(), track.GetSegment(iseg).ID()))

for j, track in enumerate(tracks): 
    if (track.s[0].Plate()>S0_END_PLATE or track.s[track.nseg-1].Plate()<S0_END_PLATE):
        continue
    elif (track.s[0].Plate()<S0_END_PLATE and track.s[track.nseg-1].Plate()>S0_END_PLATE):
        current_couple = (track.s[0].Plate(), track.s[0].ID())
        if (not current_couple in vertex_trk_couples):
            tup.Fill(track.nseg, track.npl, track.s[0].Plate(), track.sf[0].ID(), track.sf[0].X(), track.sf[0].Y(), track.sf[0].Z(), track.sf[0].TX(), track.sf[0].TY(), track.s[0].MCEvt(), track.s[0].MCTrack(), tracks.s[0].Flag(), tracks.s[0].W()-70, tracks.s[tracks.nseg-1].W()-70)
            aus_tup.Fill(track.s[track.nseg-1].X(), track.s[track.nseg-1].Y(), track.s[track.nseg-1].Z(), track.s[track.nseg-1].TX(), track.s[track.nseg-1].TY())
    if (j%10000==0):
        print(" Completed " + str(100.*j/tracks.GetEntries()) + " % ")

outfile.cd()
tup.Write()
aus_tup.Write()
outfile.Close()

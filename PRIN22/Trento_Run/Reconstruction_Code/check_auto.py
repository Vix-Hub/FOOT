import ROOT as r 
import fedrarootlogon 

# script to test auto analysis of NITs

vtxfile = r.TFile("manual.root", "READ")
mtracks = vtxfile.Get("mtracks")

inter_candfile = r.TFile("interaction_cand.root", "READ")
tup2 = inter_candfile.Get("tup2")

outfile = r.TFile("auto.root", "RECREATE")
tup = r.TNtuple("tup", "Imp of true vtx couples", "imp:d1:d2:id1:id2")

vids = []

for track in mtracks:
    if (track.vID<0 or track.vID in vids):
        continue
    vids.append(track.vID)

vertex_tracks = []

for vID in vids:
    track_list = []
    for track in mtracks:
        if track.vID==vID:
            track_list.append(track.id)
    vertex_tracks.append(track_list)

print(" Collected info about vertex tracks, now checking estimated IPs ")

for j, tracks in enumerate(vertex_tracks):
    for couple in tup2:
        if (not (couple.id1 in tracks) or not(couple.id2 in tracks)):
            continue
        tup.Fill(couple.imp, couple.d1, couple.d2, couple.id1, couple.id2)
    if (j%10==0):
        print(" Completed ", 100.*j/len(vertex_tracks), " %")

outfile.cd()
tup.Write()
outfile.Close()
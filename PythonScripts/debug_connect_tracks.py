import ROOT as r 
import fedrarootlogon
import numpy as np 

merge_file = r.TFile("2_S1_S2_offsets.root", "READ")
merge_best = merge_file.Get("merge_best")

s0_idS1s, s0_plateS1s, s0_ids, s0_plates, bs, b_backs = [], [], [], [], [], []
s0_couples_S1 = []
total_S2_couples = []

for entry in merge_best:
    s0_idS1s.append(entry.s0_idS1)
    s0_plateS1s.append(entry.s0_plateS1)
    s0_couples_S1.append((entry.s0_plateS1, entry.s0_idS1))
    s0_ids.append(entry.s0_id)
    s0_plates.append(entry.s0_plate)
    total_S2_couples.append((entry.s0_plate, entry.s0_id))

    bs.append(entry.b)
    b_backs.append(entry.b_back)

count = 0
count_good = 0
lost_tracks = []

for couple in s0_couples_S1:
    if (s0_couples_S1.count(couple)>1):
        count += 1
        indices = [i for i, x in enumerate(s0_couples_S1) if x == couple]
        ips, ipsb = [], []
        s0_couplesS2 = []
        good = 0
        for pos in indices:
            ips.append(bs[pos])
            ipsb.append(b_backs[pos])
            s0_couplesS2.append((s0_plates[pos], s0_ids[pos]))
            if (bs[pos]<100 and b_backs[pos]<100):
                good += 1
        if (len(indices)>1):
            lost_tracks.append(s0_couplesS2)
        if (good == len(indices)):
            count_good += len(indices) - 1
            print(s0_couplesS2)
             
print("Shared best candidates: " + str(count))
print(" Shared best candidates all with good b, b_back (tracks that end with no candidate) " + str(count_good))


## Now check which tracks are not connected 
merged_file = r.TFile("b000002.0.1.2.trk_merged.root", "READ")
tracks = merged_file.Get("tracks")

unconnected_couples = []

for track in tracks:
    if (track.nseg>2 and track.s[0].Plate()>30 and track.s[0].Plate()<40): #not connected with same cuts in merge_offsets
        position = total_S2_couples.index(((track.s[0].Plate(), track.s[0].ID())))
        if ( (bs[position]<100 and b_backs[position]<100)):
            unconnected_couples.append((track.s[0].Plate(), track.s[0].ID()))

print(" Printing tracks to use to debug connect tracks ")

for unconnected in unconnected_couples:
    if (not (unconnected in lost_tracks)):
        print(unconnected)
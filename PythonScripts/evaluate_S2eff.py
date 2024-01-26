import ROOT as r 
import fedrarootlogon
import numpy as np

TrackFile = r.TFile("b000111.2.0.0.trk.root", "READ")
tracks = TrackFile.Get("tracks")

outfile = r.TFile("S2_eff.root", "RECREATE")
tup = r.TNtuple("tup", "S2 efficiency check", "exp_segs:segs:plate")

r0_plates, r0_expected, r0_total = [], [], []
p0 = 31
for i in range(31, 67):
    if ((i-p0)%4==0):
        r0_expected.append(0)
        r0_total.append(0)
        r0_plates.append(i)
counter = 0
for track in tracks:
    counter+=1
    if track.nseg<5:
        continue
    for plate in r0_plates:
        index = r0_plates.index(plate)
        if (track.s[0].Plate()<= plate and track.s[track.nseg-1].Plate()>=plate):
            r0_expected[index] += 1
            for seg in track.s:
                if (seg.Plate()==plate):
                    r0_total[index] += 1
                    break
    if (counter%10000==0):
        print(" Completed " + str(100.*counter/tracks.GetEntries()) + " %")

efficiency = []
for i in range(len(r0_plates)):
    tup.Fill(r0_expected[i], r0_total[i], r0_plates[i])
    if (i!=0 and i!=len(r0_plates)-1):
        efficiency.append(r0_total[i]/r0_expected[i])


print(" Average Efficiency (excluding first and last): " + str(np.mean(efficiency)))

outfile.cd()
tup.Write()
outfile.Close()
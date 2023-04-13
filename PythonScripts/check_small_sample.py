import ROOT as r 
import fedrarootlogon
import numpy as np

## Simple Script to check How many True MC events are in track file from small sample
## usage: python check_small_sample.py in the same directory of the tracks file

file = r.TFile("b000003.2.0.0.trk.root", "READ")
tracks = file.Get("tracks")

true_mc_events = []

for track in tracks:
    mcevt = track.s[0].MCEvt()
    true_mc_events.append(mcevt)

print(" Number of MC True Events: " + str(len(np.unique(true_mc_events))))
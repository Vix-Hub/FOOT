import ROOT as r 
import fedrarootlogon 
import numpy as np 

OutFile = r.TFile("events.root", "RECREATE")
out_tree = r.TTree("events", "Events Tree")

s0plates = np.zeros(1, np.intc)
s0ids = np.zeros(1000, np.intc)
vIDs = np.zeros(1000, np.intc)
N = np.zeros(1, np.intc)
MCevts = np.zeros(1, np.intc)
trids = np.zeros(1, np.intc)
Xs = np.zeros(1, np.single)
Ys = np.zeros(1, np.single)
Zs = np.zeros(1, np.single)
TXs = np.zeros(1, np.single)
TYs = np.zeros(1, np.single)

out_tree.Branch("s0plate", s0plates, "s0plate/I")
out_tree.Branch("N", N, "N/I") #number of vertices conected to same MC event 
out_tree.Branch("s0id", s0ids, "s0id/I") #first segment info
out_tree.Branch("vIDs", vIDs, "vIDs[N]/I")
out_tree.Branch("MCEvt", MCevts, "MCEvt/I")
out_tree.Branch("s0X", Xs, "s0X/F")
out_tree.Branch("s0Y", Ys, "s0Y/F")
out_tree.Branch("s0Z", Zs, "s0Z/F")
out_tree.Branch("s0TX", TXs, "s0TX/F")
out_tree.Branch("s0TY", TYs, "s0TY/F")



infile = open("events2.txt", "r")
lines = infile.readlines()

for k, line in enumerate(lines):
    if (k == 0):
        continue
    values = line.split()
    
    trids[0] = int(values[0])
    s0plates[0] = int(values[1])
    s0ids[0] = int(values[2])
    Xs[0] = float(values[3])
    Ys[0] = float(values[4])
    Zs[0] = float(values[5])
    TXs[0] = float(values[6])
    TYs[0] = float(values[7])
    MCevts[0] = int(values[8])

    N[0] = len(values) - 9
    for j in range(N[0]):
        vIDs[j] = int(values[9+j])
    out_tree.Fill()


out_tree.Write()
OutFile.Close()
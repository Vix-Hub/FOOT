import ROOT as r 
import fedrarootlogon
import sys
sys.path.append("..")  # Add the parent directory to the sys.path
from functions import * 
import numpy as np
import pickle

# confronto con tup->Scan("d1:id1:id2:d2:start1:start2", "imp<10 && TMath::Min(d1,d2)<10 && d1<30 && d2<30 ")

CheckFile = r.TFile("mt.merged2.root", "READ")
mtracks = CheckFile.Get("mtracks")
mtracks.BuildIndex("id")

CheckFile2 = r.TFile("interaction_cand.root", "READ")
tup = CheckFile2.Get("tup")

with open("manual_list.pkl", "rb") as file:
    manual_list = pickle.load(file)

all_list, all_counter = [], []

EMU_START, EMU_END = 110, 130 #approx
plastic = 0

outFile = r.TFile("manual.root", "RECREATE")
outTup = r.TNtuple("tup", "Properties of Vtx Tracks ", "len:theta:phi:id")
outTup2 = r.TNtuple("tup2", " Properties of True Vtx Candidates", "imp:d1:d2:id1:id2")

mtracks_info, mtracks_ids = [], []
estimated_pvs = []  #this should be changed for better accuracy

# Manual Check on Tracks
for event_number, event in enumerate(manual_list):
    test_segs = []
    for trackID in event:
        trackID = int(trackID)
        mtracks.GetEntryWithIndex(trackID)
        all_list.append(trackID)
        all_counter.append(event_number)

        outTup.Fill(mtracks.len, mtracks.theta, mtracks.phi, trackID)

        seg0 = r.EdbSegP()
        seg0.SetX(mtracks.grx[0])
        seg0.SetY(mtracks.gry[0])
        seg0.SetZ(mtracks.grz[0])

        tx0, ty0 = Convert_to_TX_TY(mtracks.theta, mtracks.phi)
        seg0.SetTX(tx0)
        seg0.SetTY(ty0)
        seg0.SetID(trackID)
        test_segs.append(seg0)


        id, length, theta, phi, Ngr = mtracks.id, mtracks.len, mtracks.theta, mtracks.phi, mtracks.Ngr  #eventually add score
        grx, gry, grz = [], [], []

        for i in range(Ngr):
            grx.append(mtracks.grx[i])
            gry.append(mtracks.gry[i])
            grz.append(mtracks.grz[i])

        mtracks_info.append([id, length, theta, phi, Ngr, grx, gry, grz])
        mtracks_ids.append(id)

    imp, pv, parallel = CheckImpactN(test_segs[0], test_segs[1], 10)
    estimated_pvs.append(pv)
    print("pv ", pv)
    if (pv[2]<EMU_START):
        plastic +=1

print(" Total interactions in plastic ", plastic)

print(" now Creating Vtx Tree ")

VtxTree = r.TTree("vtx", "vertices")

multi = np.zeros(1, dtype=np.intc)
vIDs = np.zeros(1, dtype=np.intc)
trk_ID = np.zeros(1000, dtype=np.intc)
lengths, thetas, phis, Ngrs = np.zeros(1000, dtype=np.float32), np.zeros(1000, dtype=np.float32), np.zeros(1000, dtype=np.float32), np.zeros(1000, dtype=np.intc)
vx, vy = np.zeros(1, dtype=np.float32), np.zeros(1, dtype=np.float32)
vz = np.zeros(1, dtype=np.float32)

VtxTree.Branch("n", multi, "n/I")
VtxTree.Branch("vID", vIDs, "vID/I")
# Mt properties
VtxTree.Branch("IDTrack", trk_ID, "IDTrack[n]/I")
VtxTree.Branch("len", lengths, "len[n]/F")
VtxTree.Branch("theta", thetas, "theta[n]/F")
VtxTree.Branch("phi", phis, "phi[n]/F")
VtxTree.Branch("Ngr", Ngrs, "Ngr[n]/I")
VtxTree.Branch("vx", vx, "vx/F")
VtxTree.Branch("vy", vy, "vy/F")
VtxTree.Branch("vz", vz, "vz/F")

for iev, event in enumerate(manual_list):
    event = list(event)
    multi[0] = int(len(event))
    vIDs[0] = VtxTree.GetEntries()

    for j, track_ID in enumerate(event):
        position = mtracks_ids.index(track_ID)
        info = mtracks_info[position]

        trk_ID[j] = track_ID

        lengths[j] = info[1]
        thetas[j] = info[2]
        phis[j] = info[3]
        Ngrs[j] = info[4]

    vx[0] = estimated_pvs[iev][0]
    vy[0] = estimated_pvs[iev][1]
    vz[0] = estimated_pvs[iev][2]
        
    VtxTree.Fill()


print(" Adding info to tracks tuple ")

mtracks2 = mtracks.CloneTree(0)
trk_vID = np.zeros(1, dtype=np.intc)

mtracks2.Branch("vID", trk_vID, "vID/I")

for i in range(mtracks.GetEntries()):
    mtracks.GetEntry(i)
    if (mtracks.id in all_list):
        trk_vID[0] = all_counter[all_list.index(mtracks.id)]
    else:
        trk_vID[0] = -99
    mtracks2.Fill()


MIN_THETA, MAX_THETA = 0.001, 1.5
MIN_LEN = 60
MAX_LEN_OVER_NGR = 3 

print(" Looking for n=1 candidates with cuts: ")
print(" MAX Theta: ", MAX_THETA, " Min Theta ", MIN_THETA, " MIN LEN ", MIN_LEN, " MAX LEN OVER NGR ", MAX_LEN_OVER_NGR)

candidates_ids = []
for entry in mtracks:
    if (entry.len > MIN_LEN and entry.theta>MIN_THETA and entry.theta < MAX_THETA and entry.len/entry.Ngr < MAX_LEN_OVER_NGR):
        if (not (entry.id in mtracks_ids)):
            if (entry.grz[0]>=EMU_START and entry.grz[0]<EMU_END):
                candidates_ids.append(entry.id)

print(" Found N = " + str(len(candidates_ids)) + " n=1 candidates ")
print(candidates_ids)


outFile.cd()
outTup.Write()
outTup2.Write()
VtxTree.Write()
mtracks2.Write()
outFile.Close()

import ROOT as r 
import fedrarootlogon
import sys
sys.path.append("..")  # Add the parent directory to the sys.path
from functions import * 
import numpy as np

# confronto con tup->Scan("d1:id1:id2:d2:start1:start2", "imp<10 && TMath::Min(d1,d2)<10 && d1<30 && d2<30 ")

CheckFile = r.TFile("mt.merged2.root", "READ")
mtracks = CheckFile.Get("mtracks")
mtracks.BuildIndex("id")

CheckFile2 = r.TFile("interaction_cand.root", "READ")
tup = CheckFile2.Get("tup")

manual_list =  [[1653, 3660], [3764, 3762], [2735, 2676, 2686],  [4839, 4841], [2666, 5540], [8209, 8205], [16252, 16995], [17723, 18926], [17919, 17915, 17922],
  [17728, 18897, 18783], [20075, 20391], [21064, 21056, 21066], [21518, 21712], [22653, 22700], [23097, 23093, 23095], 
  [26776, 29243], [31690, 31770, 33407], [23901, 23897], [31174, 31175, 31182], [31886, 33268, 33287], [34821, 34820],
  [36753, 36726, 36730],  [40512, 40473], [42781, 42784, 42787], [43337, 44740], [43472, 44618], [43552, 41866], [44892, 44871],
  [47102, 48517], [47770, 47803], [48264, 49661], [49726, 50381], [50665, 49453], [50691, 51383], [51336, 51338], [51702, 52181], 
  [52202, 52179], [52975, 54333, 54338], [53206, 54062], [53975, 53982], [54130, 54131], [54363, 54707], [54827, 55873], [58272, 58285], [60788, 61049], [61027, 61022, 61028], [61200, 61199, 61204], [61400, 61404], [61488, 61844], [61717, 61586], [62242, 62482], [62510, 63499], [62702, 63288], [63814, 63812], [64600, 65221], [67079, 67078], [67717, 67719, 67724], [67774, 67777], [68145, 69125], [68751, 69641, 69642], [68816, 69585], [70771, 71069], [73877, 73878], [74061, 75211, 75212], [74894, 74895, 75533], [74915, 74914], 
  [80697, 80695], [80075, 80078], [82373, 82370, 82374], [82419, 82421], [84296, 84297], [84495, 84403, 84496, 84500], [85217, 85643], [86363, 86364], [88046, 89104], [88879, 88245], [89072, 89169], [89077, 89076], [89439, 89660], [90033, 90037], [90279, 90274], 
  [91596, 91595], [92481, 92480], [94196, 94198], [97566, 97560], [98626, 98623], [103351, 103651], [110346, 110347]]



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

# script to produce plots to analyze interaction candidates from interaction_cand.root 

import ROOT as r 
import fedrarootlogon
import numpy as np 
import pickle 

DO_PLOTS = 1

# Selection Cuts
MAX_IMP, MAX_D1D2, MIN_D1D2 = 5, 20, 5

r.gROOT.SetBatch(True)

# Ausiliary function
def WriteCut(interaction_ids):
    start = " id == "
    for entry in interaction_ids:
        start += str(entry) + " || id == "
    start = start + " - 999 "
    return start

def FindXY_Limits(ids_list, tree):
    x, y = [], []
    for entry in tree:
        if (entry.id in ids_list):
            x.append(entry.grx[0])
            y.append(entry.gry[0])

    return min(x), max(x), min(y), max(y)

def GetIds_with_Limits(tree, min_x, max_x, min_y, max_y):
    ids = []
    for entry in tree:
        if (entry.grx[0]>min_x and entry.grx[0]<max_x and entry.gry[0]>min_y and entry.gry[0]<max_y):
            ids.append(entry.id)
    return ids
    


IntFile = r.TFile("interaction_cand.root", "READ")
tup2 = IntFile.Get("tup2")

TrackFile = r.TFile("mt.merged2.root", "READ")
mtracks = TrackFile.Get("mtracks")
mtracks.SetMarkerStyle(6)

outFile = r.TFile("interaction_cand_plots.root", "RECREATE")
vtxTree = r.TTree("vtx", "")

vID = np.zeros(1, dtype=np.intc)
n = np.zeros(1, dtype=np.intc)
trk_ID = np.zeros(1000, dtype=np.intc)

vtxTree.Branch("vID", vID, "vID/I")
vtxTree.Branch("n", n, "n/I")
vtxTree.Branch("IDTrack", trk_ID, "IDTrack[n]/I")

interaction_couples = []
interaction_couples_info = []
for entry in tup2:
    if (entry.imp >= MAX_IMP or entry.d1 > MAX_D1D2 or entry.d2 > MAX_D1D2 or min([entry.d1, entry.d2]) > MIN_D1D2):
        continue
    interaction_couples.append((entry.id1, entry.id2))
    interaction_couples_info.append((entry.imp, entry.d1, entry.d2))

interaction_candidates = []
for i in range(len(interaction_couples)):
    couple = interaction_couples[i]
    id1 = couple[0]
    final_candidate = list(couple)
    for j in range(i+1, len(interaction_couples)):
        couple2 = interaction_couples[j]
        id2, id3 = couple2[0], couple2[1]
        if (id1 == id2 and (not (id3 in final_candidate))):
            final_candidate.append(id3)
        if (id1 == id3 and (not (id2 in final_candidate))):
            final_candidate.append(id2)
    interaction_candidates.append(final_candidate)

# remove repeated candidates
final_interaction_candidates = []
print(" Removing Duplicates ")
for i in range(len(interaction_candidates)):
    candidate1 = interaction_candidates[i]
    shared = 0
    longest = 1
    for j in range(len(interaction_candidates)):
        if (i==j):
            continue
        candidate2 = interaction_candidates[j]

        set1 = set(candidate1)
        set2 = set(candidate2)
        common_elements = set1.intersection(set2)

        if (candidate1==[2662.0, 2677.0, 2679.0, 2677.0]):
            print(" debug ", set1, set2, common_elements)

        if (len(common_elements)>0):
            shared = 1
            lengths = [len(candidate1), len(candidate2)]
            maximum_index = lengths.index(max(lengths))
            if (maximum_index != 0):
                longest = 0

    if (shared == 0):
        final_interaction_candidates.append(candidate1)
    elif (shared == 1 and longest == 1):
        final_interaction_candidates.append(candidate1)

    if (i%1000==0):
        print(" Completed " + str(100.*i / len(interaction_candidates)) + " %")

print(" Interaction Candidates ", final_interaction_candidates)

with open("final_candidates.pkl", "wb") as file:
    pickle.dump(final_interaction_candidates, file)


if (DO_PLOTS):
    outFile = r.TFile("interaction_cand_plots.root", "RECREATE")
    vtxTree = r.TTree("vtx", "")

    vID = np.zeros(1, dtype=np.intc)
    n = np.zeros(1, dtype=np.intc)
    trk_ID = np.zeros(1000, dtype=np.intc)

    vtxTree.Branch("vID", vID, "vID/I")
    vtxTree.Branch("n", n, "n/I")
    vtxTree.Branch("IDTrack", trk_ID, "IDTrack[n]/I")

    counter = 0
    for candidate in final_interaction_candidates:
        n[0] = len(candidate)
        vID[0] = counter
        for j, id in enumerate(candidate):
            trk_ID[j] = id
        vtxTree.Fill()

        counter += 1

    N = len(final_interaction_candidates)

    texts = []
    for i2 in range(N):
            temp = r.TPaveText(0.6, 0.6, 0.8, 0.8)
            temp.AddText(str(i2))
            texts.append(temp)
    N_canvas = int(N/9)
    canvas_list = []
    canvas_list2 = []
    counter = 0

    for i in range(N_canvas+1):
        c = r.TCanvas()
        c.Divide(3,3)

        c2 = r.TCanvas()
        c2.Divide(3,3)
        counters = []
        for j1 in range(1, 10):
            c.cd(j1)
            if (counter==len(final_interaction_candidates)):
                break
            cut = WriteCut(final_interaction_candidates[counter])
            counters.append(counter)
            mtracks.Draw("grz:gry:grx", cut)
            idx = (i)*9 + j1-1
            texts[idx].Draw()

            x_min, x_max, y_min, y_max = FindXY_Limits(final_interaction_candidates[counter], mtracks)
            ids = GetIds_with_Limits(mtracks, x_min-50, x_max+50, y_min-50, y_max+50)
            cut2 = WriteCut(ids)
            c2.cd(j1)
            mtracks.Draw("grz:gry:grx", cut2)

            counter += 1

            temp_title = ""
        for number in counters:
            temp_title += " " + str(number) + " "
        c.SetTitle(temp_title)
        c.Update()
        canvas_list.append(c)
        canvas_list2.append(c2)
        print(" Drawing canvas # " + str(i) + " out of " + str(N_canvas))

        if (i==1):
            break

    outFile.cd()
    for j, canvas in enumerate(canvas_list):
        canvas_list[j].Write("c"+str(j))
        canvas_list2[j].Write("k"+str(j))

    vtxTree.Write()
    outFile.Close()

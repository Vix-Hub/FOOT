import ROOT as r 
import fedrarootlogon
import numpy as np 

import sys
sys.path.append("..")  # Add the parent directory to the sys.path

from functions import *

def isin(element, test_elements, assume_unique=False):
    "..."
    element = np.asarray(element)
    return np.in1d(element, test_elements, assume_unique=assume_unique).reshape(element.shape)

def Evaluate_d(pv, segment):
    return np.sqrt( (pv[0]-segment.X())**2 + (pv[1]-segment.Y())**2 + (pv[2]-segment.Z())**2 )

# script to find close MTs and estimate the number of vertices in a sample
# modified version to use coordinates of first grains
# added check with last grains (backscattering)
# to add no sovrapposizione in coordinate tra le tracce + taglio su d0 e d1
# + rimuovere tracce che hanno troppi candidati (>=4)

DEBUG_TRKID1, DEBUG_TRKID2 = 2541, 2797

MT_File = r.TFile("mt.merged2.root", "READ")
mtracks = MT_File.Get("mtracks")
MAX_LEN_OVER_NGR = 3. 
THRESHOLD = 20
MIN_NGR = 25
DZMAX = 10

#output
outFile = r.TFile("interaction_cand.root", "RECREATE")
outTree = r.TTree("intree", "Views with interaction candidates") 

viewID = np.zeros(1, dtype=np.intc)
Ncands = np.zeros(1, dtype=np.intc)

cand_x = np.zeros(10000, dtype=np.float32)
cand_y = np.zeros(10000, dtype=np.float32)
cand_d = np.zeros(10000, dtype=np.float32)

outTree.Branch("vID", viewID, "vID/I")
outTree.Branch("ncand", Ncands, "Ncand/I")
outTree.Branch("cx", cand_x, "cx[Ncand]/F")
outTree.Branch("cy", cand_y, "cy[Ncand]/F")
outTree.Branch("d", cand_d, "d[Ncand]/F")

outTuple = r.TNtuple("tup", "Preliminary Checks", "imp:id1:id2:start1:start2:d1:d2")
outTuple2 = r.TNtuple("tup2", "Minimum IP ", "imp:id1:id2:d1:d2")


# store first grains XY coordinates of good MTs
selected_first_xyz = []
all_first_xyz = []
selected_views, views = [], []

for mt in mtracks:
    if (mt.is_dup==1 or mt.to_break==1):
        continue
    all_first_xyz.append([mt.hdid, mt.grx[0], mt.gry[0], mt.grz[0], mt.theta, mt.phi, mt.id, mt.grx[mt.Ngr-1], mt.gry[mt.Ngr-1], mt.grz[mt.Ngr-1]])
    views.append(mt.hdid)
    if (mt.Ngr<MIN_NGR):
        continue
    elif (mt.len/mt.Ngr>MAX_LEN_OVER_NGR and mt.Ngr<30):
        continue
    selected_first_xyz.append([mt.hdid, mt.grx[0], mt.gry[0], mt.grz[0], mt.theta, mt.phi, mt.id, mt.grx[mt.Ngr-1], mt.gry[mt.Ngr-1], mt.grz[mt.Ngr-1]])
    views.append(mt.hdid)

print(" Found ", len(selected_first_xyz), " MTs ")
# look for close tracks
counter, totalviews, totalcountperview = 0, 0, 0
list_of_views = []
list_of_counts = []
list_of_distances = []

for selected_track in selected_first_xyz:
    grx0, gry0, grz0 = selected_track[1], selected_track[2], selected_track[3]
    id0, theta0, phi0 = selected_track[6], selected_track[4], selected_track[5]
    grxN0, gryN0, grzN0 = selected_track[7], selected_track[8], selected_track[9]

    if (id0==DEBUG_TRKID1):
        print(" Here for debug trk ", DEBUG_TRKID1)

    candidates = []
    current_counter = 0
    for track in all_first_xyz:
        if (abs(grx0-track[1])<400 and abs(gry0-track[2])<300):
            #if max(grx0, track[1]) >= min(grxN0, track[7]) and max(gry0, track[2]) >= min(gryN0, track[8]):
                candidates.append(track)
    bs = []

    #create EdbSegP Object
    seg0 = r.EdbSegP()
    seg0.SetX(grx0)
    seg0.SetY(gry0)
    seg0.SetZ(grz0)

    tx0, ty0 = Convert_to_TX_TY(theta0, phi0)
    seg0.SetTX(tx0)
    seg0.SetTY(ty0)
    seg0.SetID(id0)

    segN0 = r.EdbSegP()
    segN0.SetX(grxN0)
    segN0.SetY(gryN0)
    segN0.SetZ(grzN0)

    segN0.SetTX(tx0)
    segN0.SetTY(ty0)
    segN0.SetID(id0)

    for candidate in candidates:
        if (candidate[1] == selected_track[1]):
            continue
        grx1, gry1, grz1 = candidate[1], candidate[2], candidate[3]
        id1, theta1, phi1 = candidate[6], candidate[4], candidate[5]
        grxN1, gryN1, grzN1 = candidate[7], candidate[8], candidate[9]

        seg1 = r.EdbSegP()
        seg1.SetX(grx1)
        seg1.SetY(gry1)
        seg1.SetZ(grz1)
        tx1, ty1 = Convert_to_TX_TY(theta1, phi1)
        seg1.SetTX(tx1)
        seg1.SetTY(ty1)
        seg1.SetID(id1)

        segN1 = r.EdbSegP()
        segN1.SetX(grxN1)
        segN1.SetY(gryN1)
        segN1.SetZ(grzN1)
        segN1.SetTX(tx1)
        segN1.SetTY(ty1)
        segN1.SetID(id1)

        d = np.sqrt((grz1-grz0)**2 + (gry1-gry0)**2 + (grx1-grx0)**2)

        # check 4 combinations (first first, first last, last first, last last)
        imp1, imp2, imp3, imp4 = 0, 0, 0, 0
        distances = []
        imps = []

        if (CheckProb(seg0, seg1, DZMAX, 1, 0.1, 1, 1)):
            imp1, pv, parallel = CheckImpactN(seg0, seg1, 10)
            d0, d1 = Evaluate_d(pv, seg0), Evaluate_d(pv, seg1)
            if (imp1<500):
                outTuple.Fill(imp1, seg0.ID(), seg1.ID(), 1, 1, d0, d1)
                distances.append((d0,d1))
                imps.append(imp1)

        if (CheckProb(seg0, segN1, DZMAX, 1, 0.1, 1, 0)):
            imp2, pv, parallel = CheckImpactN(seg0, segN1, 10)
            d0, d1 = Evaluate_d(pv, seg0), Evaluate_d(pv, segN1)
            if (imp2<500):
                outTuple.Fill(imp2, seg0.ID(), seg1.ID(), 1, 0, d0, d1)
                distances.append((d0,d1))
                imps.append(imp2)

        if (CheckProb(segN0, seg1, DZMAX, 1, 0.1, 0, 1)):
            imp3, pv, parallel = CheckImpactN(segN0, seg1, 10)
            d0, d1 = Evaluate_d(pv, segN0), Evaluate_d(pv, seg1)
            if (imp3<500):
                outTuple.Fill(imp3, seg0.ID(), seg1.ID(), 0, 1, d0, d1)
                distances.append((d0,d1))
                imps.append(imp3)

        if (CheckProb(segN0, segN1, DZMAX, 1, 0.1, 0, 0)):
            imp4, pv, parallel = CheckImpactN(segN0, segN1, 10)
            d0, d1 = Evaluate_d(pv, segN0), Evaluate_d(pv, segN1)
            if (imp4<500):
                outTuple.Fill(imp4, seg0.ID(), seg1.ID(), 0, 0, d0, d1)
                distances.append((d0,d1))
                imps.append(imp4)

        if (len(imps)>0):
            outTuple2.Fill(min(imps), seg0.ID(), seg1.ID(), distances[imps.index(min(imps))][0], distances[imps.index(min(imps))][1])
        if (seg0.ID()==DEBUG_TRKID1 and seg1.ID()==DEBUG_TRKID2):
            print(" debugging ", pv, " segs ", seg0.X(), seg0.Y(), seg0.Z())
            print( seg1.X(), seg1.Y(), seg1.Z())
            print( " d1 ", d)
            print(" d0 ", np.sqrt((pv[0]-seg0.X())**2 + (pv[1]-seg0.Y())**2 + (pv[2]-seg0.Z())**2))
        
        if (d<THRESHOLD):
            bs.append(d)
            counter += 1

print("Found", counter, " Possible interactions")

outFile.cd()
outTree.Write()
outTuple.Write()
outTuple2.Write()
outFile.Close()


"""
for view in np.unique(views):
    first_grain_xy = []
    view_counter = 0
    totalcountperview = 0
    for entry in first_xy:
        if (entry[0] == view):
            first_grain_xy.append([entry[1], entry[2]])  #recover all MT grains from a given view
    ds, xy_couples = [], []
    for i in range(len(first_grain_xy)):
        x0, y0 = first_grain_xy[i][0], first_grain_xy[i][1]
        for j in range(i+1, len(first_grain_xy)):
            x1, y1 = first_grain_xy[j][0], first_grain_xy[j][1]
            d = np.sqrt((x1-x0)**2 + (y1-y0)**2 )
            if (d<THRESHOLD):
                counter += 1
                totalcountperview += 1
                xy_couples.append([x0,y0,x1,y1])
                ds.append(d)

    if (len(ds)>0):
       
        viewID[0] = view
        ds2, xy_couples2 = [], []
        for i in range(len(ds)):
            if (ds[i]<1e-30):
                continue
            ds2.append(ds[i])
            xy_couples2.append(xy_couples[i])
        Ncands[0] = len(ds2)
        for i in range(len(ds2)):
            cand_x[i] = (xy_couples2[i][0] + xy_couples2[i][2])/2
            cand_y[i] = (xy_couples2[i][1] + xy_couples2[i][3])/2
            cand_d[i] = ds[i]
            print(ds[i])
        outTree.Fill()


print("Found", counter, " Possible interactions")

outFile.cd()
outTree.Write()
outFile.Close()
"""

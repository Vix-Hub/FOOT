# script to produce plots to analyze interaction candidates from interaction_cand.root with added cuts

import ROOT as r 
import numpy as np 
import pickle 
import sys
sys.path.append("..")  # Add the parent directory to the sys.path
from functions import *
from scipy.optimize import minimize

#remember to check if one is using mtracks.id2 (parallel) or mtracks.id (serial)
DO_PLOTS = 1
OPTIONAL_CANVASES = 0 # adds canvases with ALL mts near the vertex candidates but slowes down the script
PARALLEL = 1 #id branch issue, PARALLEL = 1 -> use id

# Selection Cuts
MAX_IMP, MAX_D1D2, MIN_D1D2 = 5, 30, 5  #5, 20, 5
NGR_MIN = 5
MIN_NGR_OVER_LEN, MAX_NGR_OVER_LEN = 0.2, 1 #1
MAX_THETA = 1. #1
MIN_R2_XYZ = 0.75 #0.75
MIN_DELTA = 3 #microns
MAX_IP_OVER_N = 3
MIN_DTHETA_COUPLE = 0. 


r.gROOT.SetBatch(True)

id_string = "id"
if (PARALLEL==0):
    id_string = "id2"

# Ausiliary functions
def WriteCut(interaction_ids):
    start = " " + id_string + " == "
    for entry in interaction_ids:
        start += str(entry) + " || " + id_string + " == "
    start = start + " - 999 "
    return start

def FindXY_Limits(ids_list, tree):
    x, y = [], []
    for entry in tree:
        current_id = entry.id
        if (PARALLEL==0):
            current_id = entry.id2
        if (current_id in ids_list):
            x.append(entry.grx[0])
            y.append(entry.gry[0])

    return min(x), max(x), min(y), max(y)

def GetIds_with_Limits(tree, min_x, max_x, min_y, max_y):
    ids = []
    for entry in tree:
        if (entry.grx[0]>min_x and entry.grx[0]<max_x and entry.gry[0]>min_y and entry.gry[0]<max_y):
            if (PARALLEL==0):
                ids.append(entry.id2)
            else:
                ids.append(entry.id)
    return ids
    
# Polynomial Regression
def polyfit(x, y, degree):
    results = {}

    coeffs = np.polyfit(x, y, degree)

     # Polynomial Coefficients
    results['polynomial'] = coeffs.tolist()

    # r-squared
    p = np.poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)                         # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    results['determination'] = ssreg / sstot

    return results

def Evaluate_IP(pv, track):

    ips = []
    for grain_index in range(len(track)):
        dz = -track[grain_index][2]+pv[2] #using first grain
        tx, ty = Convert_to_TX_TY(track[0][3], track[0][4])
        #tx, ty = track[grain_index][5], track[grain_index][6]
        x = track[grain_index][0] + dz*tx
        y = track[grain_index][1] + dz*ty
        ip = np.sqrt( (pv[0]-x)**2 + (pv[1]-y)**2 )
        ips.append(ip)
    return min(ips)

def Evaluate_d(pv, track):
    ds = []
    for grain_index in range(len(track)):
        ds.append( np.sqrt( (pv[0]-track[grain_index][0])**2 + (pv[1]-track[grain_index][1])**2 + (pv[2]-track[grain_index][2])**2 ) )
    return min(ds)

def residual(params, data):
    # Define the residual function to minimize (sum of impact parameters)
    # params: vertex coordinates
    # data: list of arrays, each array contains 3D points of a track
    # track_lengths: list of integers, representing the number of points in each track

    # Reshape params for easier access
    vertex_coords = params.reshape((-1, 3))

    # Initialize the total residual
    total_residual = 0.0

    # Loop over tracks
    for track in data:
        IP = Evaluate_IP(vertex_coords[0], track)
        # Sum the impact parameters for the current track
        total_residual += IP

    return total_residual

def fit_vertex(initial_guess, data):
    # Optimize the residual function to find optimal vertex coordinates
    #print(data)
    result = minimize(residual, initial_guess, args=(data,), method='L-BFGS-B')

    # Extract optimal vertex coordinates
    optimal_vertex_coords = result.x.reshape((-1, 3))

    # Calculate quality of fit (sum of impact parameters)
    fit_quality = result.fun

    return optimal_vertex_coords, fit_quality

IntFile = r.TFile("interaction_cand.root", "READ")
tup2 = IntFile.Get("tup2")

TrackFile = r.TFile("mt.merged2.root", "READ")
mtracks = TrackFile.Get("mtracks")
mtracks.SetMarkerStyle(6)
mtracks.BuildIndex(id_string)

outFile = r.TFile("interaction_cand_plots_cuts.root", "RECREATE")

positions_tup = r.TNtuple("pos_tup", "Candidates positions", "grx:gry:grz:n_cand")
positions_tup.SetMarkerStyle(20)
positions_tup.SetMarkerColor(r.kRed)

check_tup = r.TNtuple("debug_tup", "", "R2_xy:R2_xz:R2_yz:id:theta_new")

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

        #if (candidate1==[2662.0, 2677.0, 2679.0, 2677.0]):
        #    print(" debug ", set1, set2, common_elements)

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

    #if (i%1000==0):
    #    print(" Completed " + str(100.*i / len(interaction_candidates)) + " %")

#print(" Interaction Candidates ", final_interaction_candidates)

final_interaction_candidates_cuts = []
#additional cuts (track quality + vertex quality)
print(" Additional cuts ")
removed_track, removed_vtx = 0, 0
added_ids = []
for icand, cand in enumerate(final_interaction_candidates):
    bad_tracks = 0
    long_tracks = 0
    tracks_data = []
    test_thetas = []
    for entry in cand: #track quality
        if (entry in added_ids):
            keep = 0
            continue
        else:
            added_ids.append(entry)

        fit_xs, fit_ys, fit_zs = [], [], []
        current_track = []
        mtracks.GetEntryWithIndex(int(entry))
        for i in range(mtracks.Ngr):
            fit_xs.append(mtracks.grx[i])
            fit_ys.append(mtracks.gry[i])
            fit_zs.append(mtracks.grz[i])
            current_track.append([mtracks.grx[i], mtracks.gry[i], mtracks.grz[i], mtracks.theta, mtracks.phi])

        if (len(fit_xs)==0 or len(fit_ys)==0 or len(fit_zs)==0):
            continue

        tracks_data.append(np.array(current_track))

        res_xz = polyfit(fit_xs, fit_zs, 1) #fitting projections
        res_yz = polyfit(fit_ys, fit_zs, 1)
        res_xy = polyfit(fit_xs, fit_ys, 1)

        R2_xz = res_xz["determination"]
        R2_yz = res_yz["determination"]
        R2_xy = res_xy["determination"]

        fit_xs, fit_ys, fit_zs = sorted(fit_xs), sorted(fit_ys), sorted(fit_zs)
        DX, DY, DZ = np.abs(fit_xs[-1]-fit_xs[0]), np.abs(fit_ys[-1]-fit_ys[0]), np.abs(fit_zs[-1]-fit_zs[0])
        if (DX<MIN_DELTA): #avoid using fit results if interval is too small
            R2_xy, R2_xz = 1, 1 
        if (DY<MIN_DELTA):
            R2_xy, R2_yz = 1, 1
        if (DZ<MIN_DELTA):
            R2_xz, R2_yz = 1, 1

        Ngr_test, len_test = mtracks.Ngr, mtracks.len

        #evaluate theta again (new ref system)
        theta_test = np.pi/2 - np.arctan(np.sqrt(DX*DX/(DZ*DZ) + DY*DY/(DZ*DZ)))
        test_thetas.append(theta_test)

        check_tup.Fill(R2_xy, R2_xz, R2_yz, entry, theta_test)

        if (Ngr_test/len_test <= MIN_NGR_OVER_LEN or Ngr_test/len_test >= MAX_NGR_OVER_LEN 
            or theta_test >= MAX_THETA or R2_xz <= MIN_R2_XYZ or R2_yz <= MIN_R2_XYZ or R2_xy <= MIN_R2_XYZ
            or Ngr_test<NGR_MIN): #cut on Ngr maybe too strong
                bad_tracks += 1 

        if (Ngr_test>NGR_MIN):
            long_tracks += 1

    if (len(cand) - bad_tracks < 2 or long_tracks == 0): #if there are less than 2 good tracks or only short tracks
        removed_track += 1
        continue

    # check for n=2 vertices (one broken track seen as a vertex), cut not tested properly
    """
    if (len(cand)==2):
        dtheta = np.abs(test_thetas[0] - test_thetas[1])
        if (dtheta<= MIN_DTHETA_COUPLE):
            continue
    """
    #overall vertex quality
    if (len(tracks_data)==0):
        continue
            
    initial_guess = [0, 0, 0]
    qualities, positions = [], []
    for track in tracks_data:
        for index in [0, len(track)-1]:
            initial_guess = [track[index][0], track[index][1], track[index][2]]
            test, testq = fit_vertex(initial_guess, tracks_data)
            qualities.append(testq)
            positions.append(test)

    testq = min(qualities)
    test = positions[qualities.index(testq)][0] #use for vertex position
    if (testq/len(cand) >= MAX_IP_OVER_N):
        removed_vtx += 1
        continue
    #print(" Estimated vertex position for candidate # " + str(len(final_interaction_candidates_cuts)) + " " + str(test))
    
    positions_tup.Fill(test[0], test[1], test[2], len(final_interaction_candidates_cuts))
    final_interaction_candidates_cuts.append(cand)

    if (icand%50==0):
        print(" Completed " + str(100.*icand/len(final_interaction_candidates)) + " %")

print(" Removed Track Cut ", removed_track)
print(" Removed Vertex Cut ", removed_vtx)

with open("final_candidates_cuts.pkl", "wb") as file:
    pickle.dump(final_interaction_candidates_cuts, file)

print(len(final_interaction_candidates))
print(len(final_interaction_candidates_cuts))

if (DO_PLOTS):
    
    N = len(final_interaction_candidates_cuts)

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
            if (counter==len(final_interaction_candidates_cuts)):
                break
            cut = WriteCut(final_interaction_candidates_cuts[counter])
            #print(cut)
            
            counters.append(counter)
            mtracks.Draw("grz:gry:grx", cut)
            idx = (i)*9 + j1-1
            texts[idx].Draw()
            positions_tup.Draw("grz:gry:grx", "n_cand=="+str(counter), "same")

            if (OPTIONAL_CANVASES):
                x_min, x_max, y_min, y_max = FindXY_Limits(final_interaction_candidates_cuts[counter], mtracks)
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

        #if (i==1):
        #    break

    outFile.cd()
    for j, canvas in enumerate(canvas_list):
        canvas_list[j].Write("c"+str(j))
        canvas_list2[j].Write("k"+str(j))
    positions_tup.Write()
    check_tup.Write()
    outFile.Close()

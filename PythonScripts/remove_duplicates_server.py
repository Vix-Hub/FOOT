# -*- coding: utf-8 -*-
import ROOT as r 
import fedrarootlogon 
import numpy as np 
from ROOT import EdbPattern
import time

#DEBUG_TRKID = 94853

# INPUT FILE, OUTPUT FILE
InFileName = "vertices_improved_fast_5_twosteps.root"
InFile = r.TFile(InFileName, "READ")
vrec = InFile.Get("EdbVertexRec")
old_vtx = InFile.Get("vtx")

print("vrec starting with " + str(vrec.eVTX.GetEntries()))

OutFileName = "vertices_improved_fast_5.root"
OutFile = r.TFile(OutFileName, "RECREATE")

dproc = r.EdbDataProc()
gAli = dproc.PVR()
scancond = r.EdbScanCond()
gAli.SetScanCond(scancond)
vrec.SetPVRec(gAli)

# NEW EDBVERTEX REC
vertexrec = r.EdbVertexRec()
vertexrec.SetPVRec(gAli)
#vertexrec.eDZmax=3000.
#vertexrec.eProbMin=0.01
#vertexrec.eImpMax=15.
vertexrec.eUseMom=False
vertexrec.eUseSegPar=True
vertexrec.eQualityMode=0

# FIND DUPLICATED TRACKS TO REMOVE FROM INPUT FILE
n_vertices = vrec.eVTX.GetEntries()
track_list, track_vertex_list = [], []

for i in range(n_vertices):
    vertex = vrec.eVTX.At(i)
    for j in range(vertex.N()):
        track = vertex.GetTrack(j)
        track_list.append(track.Track())
        track_vertex_list.append((vertex.ID(), track.Track(), vertex.GetVTa(j).Imp(), j))

# some vertices contain the same track twice -> remove them
for i in range(n_vertices):
    vertex = vrec.eVTX.At(i)
    vtx_trk_list, pos_list = [], []
    for j in range(vertex.N()):
        track = vertex.GetTrack(j)
        vtx_trk_list.append(track.Track())
        pos_list.append(j)
    for k, ttrack in enumerate(vtx_trk_list):
        if (vtx_trk_list.count(ttrack)>1):
            vta_to_remove = vertex.GetVTa(pos_list[k])
            vertex.RemoveVTA(vta_to_remove)
            break

to_remove = []

for ivtx in range(n_vertices):
    vertex = vrec.eVTX.At(ivtx)
    old_vID = vertex.ID()
    for itrk in range(vertex.N()):
        track = vertex.GetTrack(itrk)
        if (track_list.count(track.Track())>1): # se la traccia è collegata a più vertici
            vIDs, bs, positions = [], [], []
            for info in track_vertex_list:
                if (info[1] == track.Track()):
                    vIDs.append(info[0])
                    bs.append(info[2])
                    positions.append(info[3])
            idx = bs.index(min(bs))
            best_vID = vIDs[idx]  # vertice che si collega meglio

            if (best_vID == vertex.ID()):
                continue
            else:
                to_remove.append((vertex.ID(), itrk))
    if (ivtx%(1000)==0):
        print(" Completed " + str(100.*ivtx/n_vertices) + " %")

to_remove2 = []
done_vids = []

for entry in to_remove:
    vID = entry[0]
    if (vID in done_vids):
        continue
    done_vids.append(vID)
    itrk_remove = entry[1]
    final_list = [vID, itrk_remove]
    for entry2 in to_remove:
        if (entry[0]==entry2[0] and entry!=entry2):                
            final_list.append(entry2[1])
    to_remove2.append(final_list)
    
#print(" to remove2 ")
#print(to_remove2)
added_vids = []

for entry in to_remove2:
    old_vID = entry[0]
    vertex = vrec.eVTX.At(old_vID)
    itrks_to_remove = []
    for i in range(len(entry)):
        if (i!=0):
            itrks_to_remove.append(entry[i])

    if (len(itrks_to_remove)==1):
        vta_to_remove = vertex.GetVTa(itrks_to_remove[0]) #only remove one track to avoid seg fault
        vertex.RemoveVTA(vta_to_remove)
    else:
        for j2, duplicated_track in enumerate(itrks_to_remove):
            vta_to_remove = vertex.GetVTa(duplicated_track-j2) #need to change the index every time I remove a track
            vertex.RemoveVTA(vta_to_remove)

    vertexrec.AddVertex(vertex)
    vertex.SetID(old_vID)
    if (old_vID==6100):
            print("Here Added 6100 (I)")
    added_vids.append(old_vID)

# Add remaining vertices (that did not share tracks from the beginning)
for i in range(n_vertices):
    vertex = vrec.eVTX.At(i)
    old_vID = vertex.ID()
    if (not (vertex.ID() in added_vids)):
        vertexrec.AddVertex(vertex)
        vertex.SetID(old_vID)
        added_vids.append(old_vID)

print(" Now I have " + str(vertexrec.eVTX.GetEntries()) + " vertices ")

# CHECK FOR DUPLICATES
track_list2 = []
still_dup, dup_list = 0, []
dup_id = 0
for i in range(vertexrec.eVTX.GetEntries()):
    vertex = vertexrec.eVTX.At(i)
    for j in range(vertex.N()):
        track = vertex.GetTrack(j)
        track_list2.append(track.Track())
for el in track_list2:
    if (track_list2.count(el)>1):
        print("Still duplicated " + str(el))
        dup_list.append(el)
        still_dup+=1


print("Number of duplicates after changes: " + str(still_dup))

OutFile.cd()
vertexrec.Write()

# CREATE TREE
FAST = 3
PLMIN, PLMAX = 0, 66
LASTLAYER = [1,30,66,76,83,90,120,140]

print("Creating Vertex Tree with Modified Vertices ")

vx, vy, vz = np.zeros(1, np.single), np.zeros(1, np.single), np.zeros(1, np.single)
maxaperture, probability = np.zeros(1, np.single), np.zeros(1, np.single)

n = np.zeros(1, np.intc)
v_flag, vplate = np.zeros(1, np.intc), np.zeros(1, np.intc)
vID = np.zeros(1, np.intc)

max_dim = 5000
IDTrack, TrackTrack = np.zeros(max_dim, np.intc), np.zeros(max_dim, np.intc)
maxgap, nseg_vec, nseg_S1 = np.zeros(max_dim, np.intc), np.zeros(max_dim, np.intc), np.zeros(max_dim, np.intc)
npl, npl_S1, nholes, nholes_S1 = np.zeros(max_dim, np.intc), np.zeros(max_dim, np.intc), np.zeros(max_dim, np.intc), np.zeros(max_dim, np.intc)
plate = np.zeros(max_dim, np.intc)

X, Y, Z = np.zeros(max_dim, np.single), np.zeros(max_dim, np.single), np.zeros(max_dim, np.single)
TX, TY = np.zeros(max_dim, np.single), np.zeros(max_dim, np.single)
Theta, impactparameter = np.zeros(max_dim, np.single), np.zeros(max_dim, np.single)

incoming, Z_flag = np.zeros(max_dim, np.intc), np.zeros(max_dim, np.intc)
MC_Charge_first, MC_Charge_last, MC_Charge_S2 = np.zeros(max_dim, np.intc), np.zeros(max_dim, np.intc), np.zeros(max_dim, np.intc)
MC_evID_last, MC_evID_first, MC_trackID_first = np.zeros(max_dim, np.intc), np.zeros(max_dim, np.intc), np.zeros(max_dim, np.intc)
MC_trackID_last = np.zeros(max_dim, np.intc)
MC_IDPart, MC_nseg = np.zeros(max_dim, np.intc), np.zeros(max_dim, np.intc)
MC_mother_first, MC_mother_last = np.zeros(max_dim, np.intc), np.zeros(max_dim, np.intc)
MC_pdgcode_first, MC_pdgcode_last = np.zeros(max_dim, np.intc), np.zeros(max_dim, np.intc)
MC_firstplate, MC_lastplate = np.zeros(max_dim, np.intc), np.zeros(max_dim, np.intc)

out_tree = r.TTree("vtx", " Vertices Tree")

out_tree.Branch("vID", vID, "vID/I")
out_tree.Branch("vx", vx, "vx/F")
out_tree.Branch("vy", vy, "vy/F")
out_tree.Branch("vz", vz, "vz/F")
out_tree.Branch("vplate", vplate, "vplate/I")
out_tree.Branch("v_flag", v_flag, "v_flag/I")
out_tree.Branch("maxaperture", maxaperture, "maxaperture/F")
out_tree.Branch("probability", probability, "probability/F")
out_tree.Branch("n", n, "n/I")
out_tree.Branch("IDTrack", IDTrack, "IDTrack[n]/I")
out_tree.Branch("TrackTrack", TrackTrack, "TrackTrack[n]/I")
out_tree.Branch("X", X, "X[n]/F")
out_tree.Branch("Y", Y, "Y[n]/F")
out_tree.Branch("Z", Z, "Z[n]/F")
out_tree.Branch("TX", TX, "TX[n]/F")
out_tree.Branch("TY", TY, "TY[n]/F")
out_tree.Branch("Theta", Theta, "Theta[n]/F")
out_tree.Branch("nseg", nseg_vec, "nseg_vec[n]/I")
out_tree.Branch("npl", npl, "npl[n]/I")
out_tree.Branch("nholes", nholes, "nholes[n]/I")
out_tree.Branch("nseg_S1", nseg_S1, "nseg_S1[n]/I")
out_tree.Branch("nholes_S1", nholes_S1, "nholes_S1[n]/I")
out_tree.Branch("plate", plate, "plate[n]/I")
out_tree.Branch("maxgap", maxgap, "maxgap[n]/I")
out_tree.Branch("incoming", incoming, "incoming[n]/I")
out_tree.Branch("impactparameter", impactparameter, "impactparameter[n]/F")
out_tree.Branch("Z_flag", Z_flag, "Z_flag[n]/I")
out_tree.Branch("MC_Charge_first", MC_Charge_first, "MC_Charge_first[n]/I")

  
out_tree.Branch("plate", plate, "plate[n]/I")
out_tree.Branch("maxgap", maxgap, "maxgap[n]/I")
out_tree.Branch("incoming", incoming, "incoming[n]/I")
out_tree.Branch("impactparameter", impactparameter, "impactparameter[n]/F")
out_tree.Branch("Z_flag", Z_flag, "Z_flag[n]/I")
out_tree.Branch("MC_Charge_first", MC_Charge_first, "MC_Charge_first[n]/I")
out_tree.Branch("MC_Charge_last", MC_Charge_last, "MC_Charge_last[n]/I")
out_tree.Branch("MC_Charge_S2", MC_Charge_S2, "MC_Charge_S2[n]/I")
out_tree.Branch("MC_evID_first", MC_evID_first, "MC_evID_first[n]/I")
out_tree.Branch("MC_evID_last", MC_evID_last, "MC_evID_last[n]/I")
out_tree.Branch("MC_trackID_first", MC_trackID_first, "MC_trackID_first[n]/I")
out_tree.Branch("MC_trackID_last", MC_trackID_last, "MC_trackID_last[n]/I")
out_tree.Branch("MC_mother_first", MC_mother_first, "MC_mother_first[n]/I")
out_tree.Branch("MC_mother_last", MC_mother_last, "MC_mother_last[n]/I")
out_tree.Branch("MC_pdgcode_first", MC_pdgcode_first, "MC_pdgcode_first[n]/I")
out_tree.Branch("MC_pdgcode_last", MC_pdgcode_last, "MC_pdgcode_last[n]/I")
out_tree.Branch("MC_firstplate", MC_firstplate, "MC_firstplate[n]/I")
out_tree.Branch("MC_lastplate", MC_lastplate, "MC_lastplate[n]/I")

old_vtx.BuildIndex("vID")

for ivtx in range(vertexrec.eVTX.GetEntries()):
    vertex = vertexrec.eVTX.At(ivtx)
    #if (vertex.ID() == 0 and FAST != 1):
    #    vID[0] = ivtx + 900000
    #else:
    #    vID[0] = ivtx
    vID[0] = vertex.ID()
    #vertex.SetID(int(vID[0]))
    vx[0] = vertex.VX()
    vy[0] = vertex.VY()
    vz[0] = vertex.VZ()
    n[0] = vertex.N()
    old_vtx.GetEntryWithIndex(int(vertex.ID()))
    vplate[0] = int(old_vtx.vplate)
    vz[0] = float(old_vtx.vz)
    #if (vplate[0]<0):
    #    old_vtx.GetEntry(vID[0])
    #    vplate[0] = int(old_vtx.vplate)
    #    vz[0] = float(old_vtx.vz)
    v_flag[0] = vertex.Flag()

    maxaperture[0] = vertex.MaxAperture()
    probability[0] = vertex.V().prob()
    for itrk in range(vertex.N()):
        track = vertex.GetTrack(itrk)

        IDTrack[itrk] = track.Track()
        TrackTrack[itrk] = track.Track()
        nseg_vec[itrk] = track.N()
        npl[itrk] = track.Npl()
        zpos = vertex.GetVTa(itrk).Zpos()
        incoming[itrk] = zpos
        X[itrk] = track.GetSegment(0).X()
        Y[itrk] = track.GetSegment(0).Y()
        Z[itrk] = track.GetSegment(0).Z()
        TX[itrk] = track.TX()
        TY[itrk] = track.TY()
        Theta[itrk] = track.Theta()
        plate[itrk] = track.GetSegment(0).Plate()
        nholes[itrk] = track.N0()
        maxgap[itrk] = track.CheckMaxGap()

        nseg_S1[itrk] = 0
        npl_S1[itrk] = 0
        nholes_S1[itrk] = 0

        if plate[itrk] < 31:
            for iseg in range(track.N()):
                if track.GetSegment(iseg).Plate() <= LASTLAYER[1]:
                    nseg_S1[itrk] += 1
            npl_S1[itrk] = LASTLAYER[1] - plate[itrk] + 1
            nholes_S1[itrk] = npl_S1[itrk] - nseg_S1[itrk]

        impactparameter[itrk] = vertex.GetVTa(itrk).Imp()

        Z_flag[itrk] = track.Flag()
        for idx in range(track.N()):
            if track.GetSegment(idx).Plate() > LASTLAYER[1] + 1:
                MC_Charge_S2[itrk] = int(track.GetSegment(idx).W()) - 70
                continue

        MC_evID_first[itrk] = track.GetSegmentFirst().MCEvt()
        MC_evID_last[itrk] = track.GetSegmentLast().MCEvt()
        MC_trackID_first[itrk] = track.GetSegmentFirst().MCTrack()
        MC_trackID_last[itrk] = track.GetSegmentLast().MCTrack()
        MC_pdgcode_first[itrk] = track.GetSegmentFirst().Aid(0)
        MC_pdgcode_last[itrk] = track.GetSegmentLast().Aid(0)
        MC_Charge_first[itrk] = int(track.GetSegmentFirst().W()) - 70
        MC_Charge_last[itrk] = int(track.GetSegmentLast().W()) - 70
        MC_mother_first[itrk] = track.GetSegmentFirst().Aid(1)
        MC_mother_last[itrk] = track.GetSegmentLast().Aid(1)
        MC_firstplate[itrk] = track.GetSegmentFirst().Vid(0)
        MC_lastplate[itrk] = track.GetSegmentLast().Vid(0)

    out_tree.Fill()

OutFile.cd()
out_tree.Write()
    


OutFile.Close()
InFile.Close()

import ROOT as r 
import fedrarootlogon
import numpy as np 
import time
import sys
sys.path.append("..")  # Add the parent directory to the sys.path
from functions import *
import time

# script to improve post processing after linking (from merge_mts_new.py)
# completely remove duplicates 
# break tracks if needed? for now, simply ignore them
# add check on GRY and GRZ

EVERBOSE = -100
DEBUG_TRK_ID = 2541

MTFile = r.TFile("mt.merged.root", "READ")
mtracks = MTFile.Get("mtracks")

MTFile2 = r.TFile("mt.merged2.root", "RECREATE")
mtracks2 = mtracks.CloneTree(0)

grains_list = [] #used to check if grains are shared by check on X coordinates
entries, ids = [], []

for i in range(mtracks.GetEntries()):
    mtracks.GetEntry(i)
    if (mtracks.is_dup==1 or mtracks.to_break==1):
        continue
    current_list = []
    for igrain in range(mtracks.Ngr):
        current_list.append((mtracks.grx[igrain], mtracks.gry[igrain]))
    grains_list.append(current_list)
    entries.append(i)
    ids.append(mtracks.id)

t0 = time.time()

for i in range(len(entries)):

    avoid_saving = 0
    grx1 = grains_list[i]
    id = ids[i]

    if (id==DEBUG_TRK_ID):
        print(" here debugging ")

    for j in range(len(entries)):
        if (i!=j):
            grx2 = grains_list[j]

            if (abs(grx1[0][0]-grx2[0][0])>400 or abs(grx1[0][1]-grx2[0][1])>300):
                continue

            id2 = ids[j]
            common_elements = set(grx1) & set(grx2)

            if common_elements:
                # If there are common elements, select the longest list
                longest_list = max([grx1, grx2], key=len)
                if (longest_list!=grx1):
                    avoid_saving = 1

            if (id==DEBUG_TRK_ID and common_elements):
                print(" here debugging ", common_elements, avoid_saving, " other id", id2)
    if (i%1000 == 0):
        print(" Completed " + str(100.*i/len(entries)) + " % " + str(time.time()-t0) + " s ")
    
    if (avoid_saving==0):
        mtracks.GetEntry(entries[i])
        mtracks2.Fill()
        #print(i)

MTFile2.cd()
mtracks2.Write()
MTFile2.Close() 

import ROOT as r 
#import fedrarootlogon
import numpy as np 
import time
import sys
sys.path.append("..")  # Add the parent directory to the sys.path
from functions import *
import time

XY_LIMIT = 100 #max absolute difference

# script to improve post processing after linking (from merge_mts_new.py)
# completely remove duplicates 
# break tracks if needed? for now, simply ignore them
# add check on GRY and GRZ
MTFile = r.TFile("mt.merged.all.root", "READ")
mtracks = MTFile.Get("mtracks")


# Removing duplicates
MTFile2 = r.TFile("mt.merged2.root", "RECREATE")

outTree2 = mtracks.CloneTree(0)

is_dup = np.zeros(1, dtype=np.intc)
to_break = np.zeros(1, dtype=np.float32)
outTree2.Branch("is_dup", is_dup, "is_dup/I")
outTree2.Branch("to_break", to_break, "to_break/F")
unique_entries = {}

print(" Removing duplicates and marking unwanted links ")

for entry in range(mtracks.GetEntries()):
	mtracks.GetEntry(entry)

	entry_values = []
	# Convert the values of the entry to a tuple to use as a key in the dictionary
	for coord_x in mtracks.grx:
		entry_values.append(coord_x)
	entry_values = tuple(sorted(entry_values))

	# Check if the entry is already in the dictionary
	if entry_values in unique_entries:
		# Entry is a duplicate
		is_dup[0] = 1
	else:
		# Entry is unique
		is_dup[0] = 0
		# Add the entry to the dictionary
		unique_entries[entry_values] = None

	if (entry==7552 or entry==6821):
		print(" Entry ", entry, " dup ", is_dup[0], " entry values ", entry_values)

	break_index = -1
	grains_x, grains_y, grains_z = mtracks.grx, mtracks.gry, mtracks.grz
	sorted_grains_z = np.sort(grains_z)

	grains_x = list(grains_x)
	grains_z = list(grains_z)

	for i in range(len(grains_x)-1):
		idx = grains_z.index(sorted_grains_z[i])
		next_idx = grains_z.index(sorted_grains_z[i+1])
		dx = abs(grains_x[idx]-grains_x[next_idx])
		dy = abs(grains_y[idx]-grains_y[next_idx])
		dz = abs(sorted_grains_z[i]-sorted_grains_z[i+1])
		d = np.sqrt(dx**2 + dy**2 + dz**2)
		if (dx>XY_LIMIT or dy>XY_LIMIT or dz>XY_LIMIT or d>50):
			break_index = i
			if (mtracks.id == 2677):
				print(" check ", dx, dy, dz, break_index, i, grains_x[i], grains_x[i+1])
			break

	
	if (break_index < 0):
		to_break[0] = 0 #save mt normally
	else:
		to_break[0] = 1. # mark first unwanted link
		if (mtracks.id == 330):
			print(" check ", dx, dy, dz, break_index, to_break)

	# Fill the modified TTree with the current entry
	#hdid[0] = iview

	outTree2.Fill()


EVERBOSE = -100
DEBUG_TRK_ID = 2541
 
mtracks2 = outTree2.CloneTree(0)  #not really needed to split the steps, optimization needed!

grains_list = [] #used to check if grains are shared by check on X coordinates
entries, ids = [], []

for i in range(outTree2.GetEntries()):
    outTree2.GetEntry(i)
    if (outTree2.is_dup==1 or outTree2.to_break==1):
        continue
    current_list = []
    for igrain in range(outTree2.Ngr):
        current_list.append((outTree2.grx[igrain], outTree2.gry[igrain]))
    grains_list.append(current_list)
    entries.append(i)
    ids.append(outTree2.id)

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

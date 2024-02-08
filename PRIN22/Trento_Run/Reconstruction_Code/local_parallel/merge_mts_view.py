import ROOT as r 
#import fedrarootlogon
import numpy as np 
import time
import sys
#sys.path.append("..")  # Add the parent directory to the sys.path
from functions import *
import pickle
import argparse

parser = argparse.ArgumentParser(description="Process a View")
parser.add_argument("viewID", type=int, help="view ID to process")
args = parser.parse_args()
print(args.viewID)
SINGLE_VIEW = args.viewID
N_VIEWS = 200

EVERBOSE = -100
ANGLE_THRESHOLD = 0.15 #tested
IP_THRESHOLD = 80
MIN_IP_THRESHOLD = 50 #at least one should be lower than this
XY_LIMIT = 100 #max absolute difference
Z_LIMIT = 100
DEBUG_TRK_ID = -435

ViewFile = r.TFile("dm_tracks2.dm.root", "READ")
Vdmr = ViewFile.Get("Vdmr")
# organize grains in an easy-to-access way
gr_imt = Vdmr.GetLeaf("gr.imt")
gr_x = Vdmr.GetLeaf("gr.x")
gr_y = Vdmr.GetLeaf("gr.y")
gr_z = Vdmr.GetLeaf("gr.z")
xleaf = Vdmr.GetLeaf("x")
yleaf = Vdmr.GetLeaf("y")
zleaf = Vdmr.GetLeaf("z")

score_treshold = -999
n_grains_min = 3
N = Vdmr.GetEntries()

outFile = r.TFile("mt.merged.temp."+str(SINGLE_VIEW)+".root", "RECREATE")
outTree = r.TTree("mtracks", "Merged MTs")

Ngrs = np.zeros(1, dtype=np.intc)
scores = np.zeros(1, dtype=np.float32)
thetas = np.zeros(1, dtype=np.float32)
phis = np.zeros(1, dtype=np.float32)
lens = np.zeros(1, dtype=np.float32)
ids = np.zeros(1, dtype=np.intc)
hdid = np.zeros(1, dtype=np.intc)

grxs = np.zeros(2000, dtype=np.float32)
grys = np.zeros(2000, dtype=np.float32)
grzs = np.zeros(2000, dtype=np.float32)

outTree.Branch("Ngr", Ngrs, "Ngr/I")
outTree.Branch("score", scores, "score/F")
outTree.Branch("theta", thetas, "theta/F")
outTree.Branch("phi", phis, "phi/F")
outTree.Branch("len", lens, "len/F")
outTree.Branch("id", ids, "id/I")

outTree.Branch("grx", grxs, "grx[Ngr]/F")
outTree.Branch("gry", grys, "gry[Ngr]/F")
outTree.Branch("grz", grzs, "grz[Ngr]/F")
outTree.Branch("hdid", hdid, "hdid/I")


with open("view_info.pkl", "rb") as file:
	view_data = pickle.load(file)
	
views_x, views_y, actual_views = view_data["views_x"], view_data["views_y"], view_data["actual_views"]

t0 = time.time()
total_used_links = []

for k1, iview in enumerate(actual_views):

	if (SINGLE_VIEW>=0 and (iview<SINGLE_VIEW or iview>=SINGLE_VIEW+N_VIEWS or iview>=Vdmr.GetEntries())): #option for debugging
		continue

	Vdmr.GetEntry(iview)
	hdflag = Vdmr.GetLeaf("hd.flag").GetValue(iview)
	if (hdflag !=0):
		continue

	# select closest views based on XY values
	current_view_x = views_x[k1]
	current_view_y = views_y[k1]

	close_views = [iview]
	for k, iview2 in enumerate(actual_views):
		if(iview==iview2):
			continue
		candidate_view_x = views_x[k]
		candidate_view_y = views_y[k]
		if (abs(candidate_view_x-current_view_x)<=400 and abs(candidate_view_y-current_view_y)<=300 ):
			close_views.append(iview2)

	if (EVERBOSE==100):
		print(" Closest Views ", close_views)

	microtracks_grains_coordinates, microtracks_ids, microtracks_angles = [], [], []
	microtracks_coordinates = []
	N_total_microtracks = 0
	total_microtracks_grains_coordinates = []

	current_view_microtracks_grains_coordinates, current_view_microtracks_ids = [], [] #used to store info about view #iview
	current_view_microtracks_angles, current_view_microtracks_lens = [], []
	current_N_microtracks = 0

	for selected_view in close_views:
		Vdmr.GetEntry(selected_view)

		# check again to avoid using preliminary views
		hdflag = Vdmr.GetLeaf("hd.flag").GetValue(selected_view) 
		if (hdflag !=0):
			continue

		N_microtracks = len(Vdmr.mt)
		N_total_microtracks += N_microtracks

		if (selected_view==iview):
			current_N_microtracks = N_microtracks
		
		total_microtracks_grains_coordinates.append([]) #this list is organized as [ [list of MTs from close view 1],  [ same for close view 2 ], ...]
		current_close_view_microtracks_grains_coordinates, current_close_view_microtracks_ids = [], [] #updated for each view

		for itrk in range(N_microtracks): #first save MT info like ID, theta, phi, X, Y, Z
			microtracks_ids.append(int(Vdmr.GetLeaf("mt.id").GetValue(itrk)))
			microtracks_angles.append((Vdmr.GetLeaf("mt.theta").GetValue(itrk), Vdmr.GetLeaf("mt.phi").GetValue(itrk)))
			microtracks_coordinates.append((Vdmr.GetLeaf("mt.x").GetValue(itrk) + xleaf.GetValue(itrk), Vdmr.GetLeaf("mt.y").GetValue(itrk) + yleaf.GetValue(itrk), Vdmr.GetLeaf("mt.z").GetValue(itrk) + zleaf.GetValue(itrk), selected_view, itrk))
		
			if (selected_view==iview): #here do the same but only for current view (iview)
				current_view_microtracks_grains_coordinates.append([])
				current_view_microtracks_angles.append(microtracks_angles[-1])
				current_view_microtracks_ids.append(microtracks_ids[-1])
				current_view_microtracks_lens.append((Vdmr.GetLeaf("mt.len").GetValue(itrk), Vdmr.GetLeaf("mt.ngr").GetValue(itrk)))

			current_close_view_microtracks_grains_coordinates.append([])
			current_close_view_microtracks_ids.append(microtracks_ids[-1])

		# find the grains corresponding to the MTs
		for igr in range(len(Vdmr.gr)):
			imt = gr_imt.GetValue(igr)
			x, y, z = gr_x.GetValue(igr), gr_y.GetValue(igr), gr_z.GetValue(igr)
			x += xleaf.GetValue(igr)
			y += yleaf.GetValue(igr)
			z += zleaf.GetValue(igr)
			if(imt>-1):
				position = current_close_view_microtracks_ids.index(int(imt))
				current_close_view_microtracks_grains_coordinates[position].append((x, y, z))
				if (selected_view == iview):
					current_view_microtracks_grains_coordinates[position].append((x, y, z))

		# store the grains info of all selected views in this list 
		for idx in range(N_microtracks):	
			total_microtracks_grains_coordinates[-1].append([current_close_view_microtracks_ids[idx], current_close_view_microtracks_grains_coordinates[idx]])  # [ [id1, grains1], ...]

		if (EVERBOSE==100):
			print(" For View # " + str(iview) + " there are " + str(N_microtracks) + " MTs ")
			if (N_microtracks>0):
				print(" First MT: " + str(microtracks_ids[0]))

	if (EVERBOSE==100):
		for j, entry in enumerate(total_microtracks_grains_coordinates):
			print(" close view ", close_views[j])
			print(" total grains ", entry)

	new_microtrack_grains_coordinates = []
	to_merge = []

	# loop to calculate IP -> find links (= couples of mergeable MTs) like [(view1, itrk1), (view2, itrk2)]

	for i1 in range(N_total_microtracks):
		theta1, phi1 = microtracks_angles[i1][0], microtracks_angles[i1][1]
		tx1, ty1 = Convert_to_TX_TY(theta1, phi1)
		mt1_x, mt1_y, mt1_z = microtracks_coordinates[i1][0], microtracks_coordinates[i1][1], microtracks_coordinates[i1][2]
		view1 = microtracks_coordinates[i1][3]
		itrk1 = microtracks_coordinates[i1][4]

		#still to be tested
		#grains1 = FindGrains(total_microtracks_grains_coordinates, view1, itrk1, close_views)
		#intx1, inty1, intz1 = GetGrainsIntervals(grains1)

		for i2 in range(i1+1, N_total_microtracks):
			theta2, phi2 = microtracks_angles[i2][0], microtracks_angles[i2][1]
			view2 = microtracks_coordinates[i2][3]
			itrk2 = microtracks_coordinates[i2][4]
			tx2, ty2 = Convert_to_TX_TY(theta2, phi2)

			#grains2 = FindGrains(total_microtracks_grains_coordinates, view2, itrk2, close_views)
			#intx2, inty2, intz2 = GetGrainsIntervals(grains2)
			#cond1, cond2, cond3 = TestIntervals(intx1, intx2), TestIntervals(inty1, inty2), TestIntervals(intz1, intz2)

			mt2_x, mt2_y, mt2_z = microtracks_coordinates[i2][0], microtracks_coordinates[i2][1], microtracks_coordinates[i2][2]

			dz = mt2_z - mt1_z
			dx = mt2_x - dz*tx1 - mt1_x
			dy = mt2_y - dz*ty1 - mt1_y
			b_for = np.sqrt(dx**2 + dy**2)

			dz = mt1_z - mt2_z
			dx = mt1_x - dz*tx2 - mt2_x
			dy = mt1_y - dz*ty2 - mt2_y
			b_back = np.sqrt(dx**2 + dy**2)

			if (EVERBOSE == 100):
				print(" Calculated IP using MT coodinates between ", (view1, itrk1, view2, itrk2), " IP: ", b_for, " IP_back: ", b_back)
				print( " dz ", dz, " dx, dy ", dx, dy, " ")
				print(" ty1, ty2 ", ty1, ty2, " tx1, tx2 ", tx1, tx2)
				print(" theta1, phi1 ", theta1, phi1, " theta2, phi2 ", theta2, phi2)
				print(" x1, x2, y1, y2 ", mt1_x, mt2_x, mt1_y, mt2_y)
				#print(" cond1, cond2, cond3 ", cond1, cond2, cond3, "\n")

			# actual check on IP
			if (b_for<IP_THRESHOLD and b_back<IP_THRESHOLD and abs(theta1-theta2)<ANGLE_THRESHOLD 
	   			and abs(phi1-phi2)<ANGLE_THRESHOLD and abs(mt1_x-mt2_x)<XY_LIMIT 
	   			and abs(mt2_y-mt1_y)<XY_LIMIT and abs(dz)<Z_LIMIT and min([b_for, b_back])<MIN_IP_THRESHOLD):
				to_merge.append((view1, itrk1, view2, itrk2))

	if(EVERBOSE==100):
		print("TO MERGE ", to_merge)

	# at first we consider all the links, even if there will be duplicates
	to_merge_mts_with_dup = [] 
	for i in range(len(to_merge)):
		tracks1 = extract_tracks_from_couples(to_merge[i])
		for j in range(1, len(to_merge)):
			if (i==j):
				continue
			tracks2 =  extract_tracks_from_couples(to_merge[j])
			duplicates = find_duplicates(tracks1, tracks2)
			if (len(duplicates)>0): # check for shared tracks
				tracks3 = merge_and_remove_duplicates(tracks1, tracks2) # if tracks are shared merge links
				tracks1 = tracks3
		#save merged link 
		to_merge_mts_with_dup.append(tracks1)

	to_merge_mts = eliminate_equivalent_lists(to_merge_mts_with_dup) #removing duplicates

	if (EVERBOSE==100):
		print(" merged mts with duplicates ")
		print(to_merge_mts_with_dup)
		print( "to merge mts without dup ")
		print(to_merge_mts)

	#actually merge grains using the MTs' Z as order
	for i in range(len(to_merge_mts)):
		to_merge_links = to_merge_mts[i]
		
		if (EVERBOSE==100):
			print(" to merge ", to_merge_links)

		tracks = extract_tracks_from_couples(to_merge_links)
		for el in tracks:
			if (not (el in total_used_links)):
				total_used_links.append(el)

		merged_microtracks_grains_coordinates = []
		merged_grains_list = []
		for mtrack_identifier in to_merge_links:
			view_index = close_views.index(mtrack_identifier[0])

			if (EVERBOSE==100):
				print(" len ", len(total_microtracks_grains_coordinates), " view index ", view_index, " close views ", close_views, " view ", mtrack_identifier[0])
				print(" is uguale ? ", total_microtracks_grains_coordinates[0]==total_microtracks_grains_coordinates[1])
			
			for mtinfo in total_microtracks_grains_coordinates[view_index]:
				if (mtinfo[0]==mtrack_identifier[1] and len(mtinfo[1])>0 and close_views[view_index]==mtrack_identifier[0]):
					for entry in mtinfo[1]:
						merged_grains_list.append(entry)
						
		sorted_grains_list = sort_tuples_by_z(merged_grains_list) #sort list by Z

		if (len(sorted_grains_list)==0):
			continue

		# fit should still be properly tested
		ty, tx, L = line_fit_and_properties(sorted_grains_list) # evaluate new XYZ, theta, phi, len
		#theta, phi = Convert_to_Theta_Phi(np.deg2rad(tx), np.deg2rad(ty))

		merged_ngr = len(sorted_grains_list)
		DX, DY, DZ = sorted_grains_list[merged_ngr-1][0]-sorted_grains_list[0][0], sorted_grains_list[merged_ngr-1][1]-sorted_grains_list[0][1], sorted_grains_list[merged_ngr-1][2]-sorted_grains_list[0][2]
		if (DX==0 or DZ==0):
			continue

		theta, phi = np.arctan( np.sqrt( (DX/DZ)**2 + (DY/DZ)**2) ), np.arctan( (DY/DX) )

		thetas[0] = theta 
		phis[0] = phi
		lens[0] = L
		ids[0] = outTree.GetEntries()
		Ngrs[0] = len(sorted_grains_list)
		hdid[0] = iview

		if len(sorted_grains_list) > 2000:
			continue


		for j, grain in enumerate(sorted_grains_list):
			grxs[j] = sorted_grains_list[j][0]
			grys[j] = sorted_grains_list[j][1]
			grzs[j] = sorted_grains_list[j][2]

		if (ids[0]==DEBUG_TRK_ID):
			print(" MERGED here id ", DEBUG_TRK_ID, " for iview ", iview)
			print(" sorted grains ", sorted_grains_list)

		if (EVERBOSE==100):
			print(" MERGED - Here added ID ", ids[0])
		outTree.Fill()

	total_used_mts = [] #ausiliary list to avoid adding the same MTs more than once (may still need testing)
	for link in total_used_links:
		for mytuple in link:
			if (not (mytuple in total_used_mts)):
				total_used_mts.append(mytuple)

	# add remaining MTs
	for itrk in range(current_N_microtracks):
		mt_info = (iview,current_view_microtracks_ids[itrk])

		if (EVERBOSE==100):
			print("mt_info ", mt_info)

		if (not (mt_info in total_used_mts)):
			total_used_mts.append(mt_info)

			if (EVERBOSE==100):
				print(" OK not used ")

			thetas[0] = current_view_microtracks_angles[itrk][0] 
			phis[0] = current_view_microtracks_angles[itrk][1]
			lens[0] = current_view_microtracks_lens[itrk][0]
			ids[0] = outTree.GetEntries()
			Ngrs[0] = current_view_microtracks_lens[itrk][1]

			for j, grain in enumerate(current_view_microtracks_grains_coordinates[itrk]):
				grxs[j] = current_view_microtracks_grains_coordinates[itrk][j][0]
				grys[j] = current_view_microtracks_grains_coordinates[itrk][j][1]
				grzs[j] = current_view_microtracks_grains_coordinates[itrk][j][2]

			if (EVERBOSE==100):
				print(" NOT MERGED - Here added ID ", ids[0])

			if (ids[0]==DEBUG_TRK_ID):
				print(" here id ", DEBUG_TRK_ID, " for iview ", iview)

			outTree.Fill()

print(" Completed process with viewID ", SINGLE_VIEW)

outFile.cd()
outTree.Write()
outFile.Close()

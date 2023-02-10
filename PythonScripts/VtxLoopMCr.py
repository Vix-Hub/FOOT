
import ROOT as r
import fedrarootlogon
import numpy as np
import time

IDBRICK = 2
charge_name = "b00000" + str(IDBRICK) + ".2.0.0.trk.root"
charge_file = r.TFile(charge_name, "READ")

tracks = charge_file.Get("tracks")

f = open('vtx_tracks_MC.txt', 'r')
lines = f.readlines()
check_list, check_list2 = [], []
check_ulteriore, vids = [], []
for line in lines:
	#check_list.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2]), float(line.split()[3]), int(line.split()[4])])  #X,Y,TX,TY,Plate_ID, seg_ID, vz, S2_segs, vid
	check_list.append([float(line.split()[4]), float(line.split()[5])]) # Plate_ID, Seg_ID
	check_ulteriore.append( float(line.split()[7]) )  #S2_segs
	check_list2.append(float(line.split()[6])) #vz
	vids.append(float(line.split()[8]))
out_file2 = r.TFile("b" + str(IDBRICK) + "_vtx.root", "RECREATE")

in_vtx = np.zeros(1, dtype = np.intc)
v_id = np.zeros(1, dtype = np.intc)
vz = np.zeros(1, dtype = np.single)

output_tree = tracks.CloneTree(0)
output_tree.Branch("in_vtx", in_vtx, "in_vtx/I")  #tracce con match in S2
output_tree.Branch("vz", vz, "vz/F")
output_tree.Branch("v_ID", v_id, "v_ID/I")

tracks.BuildIndex("s[0].ID()", "s[0].Plate()")

n_trk = 0

start = time.time()
match = 0
n_check, not_found = 0, 0

check_tuple = r.TNtuple("check_t", "", "s_ID:s_Plate:S2_segs:v_ID")

check_max = 10

for i in range(len(check_list)):
	
	check1 = check_list[i]
	if (tracks.GetEntryWithIndex(int(check1[1]), int(check1[0]))>0):
		in_vtx[0] = 1
		vz[0] = check_list2[i]
		match = match + 1
		v_id[0] = vids[i]
		output_tree.Fill()
		if (n_check < check_max):
			print("Number of segments in S2 file: " + str(tracks.nseg))
			print("Number of segments in S2 from vtx file: " + str(check_ulteriore[i]))
			print("v_ID = " + str(vids[i]))

			print("Seg IDs: " + str(check1[1]) + " " + str(tracks.s[0].ID()))
			print("Seg Plates: " + str(check1[0]) + " " + str( tracks.s[0].Plate() )) 
			n_check = n_check + 1
	elif (tracks.GetEntryWithIndex(int(check1[1]), int(check1[0])) == -1):
		in_vtx[0] = 0
		vz[0] = -50
		v_id[0] = -500
		if (not_found < check_max):
			print("NOT FOUND: seg_ID = " + str(check1[1]) + ", seg_Plate = " + str(check1[0]) + ", S2_segs = " + str(check_ulteriore[i])  + ", V_ID = " + str(vids[i]) )
			not_found += 1
		check_tuple.Fill(check1[1], check1[0], check_ulteriore[i], vids[i])
		output_tree.Fill()

	if ( (i)%1000 == 0 ):
		print(" Completed " + str(round(100*i/len(check_list), 2)) + " %")

check_tuple.Write("check_t")
output_tree.Write("tracks")
out_file2.Close()
charge_file.Close()

print("Total Number of Tracks Connected to Vertices with Segments in S2 = " + str(len(check_list)))
print(" Total Number of Matches with tracks in S2 = " + str(match) + " (" + str(100*match/len(check_list)) + " %)" )
print(" Time: " + str(time.time()-start) + " seconds ")

import ROOT as r 
import fedrarootlogon 
import pickle

# script to save to lists the XY coordinates of each View and the entries of the scanning views 

ViewFile = r.TFile("dm_tracks2.dm.root", "READ")
Vdmr = ViewFile.Get("Vdmr")

views_x, views_y = [], [] # find XY coordinates of all views so that later one can find views close to the selected one
actual_views = []
N = Vdmr.GetEntries()
for iview in range(N):
	Vdmr.GetEntry(iview)
	hdflag = Vdmr.GetLeaf("hd.flag").GetValue(iview)
	if (hdflag==0):
		views_x.append(Vdmr.GetLeaf("x").GetValue(iview))
		views_y.append(Vdmr.GetLeaf("y").GetValue(iview))
		actual_views.append(iview)
		
dataset = {}
dataset["views_x"] = views_x
dataset["views_y"] = views_y
dataset["actual_views"] = actual_views

with open("view_info.pkl", "wb") as file:
	pickle.dump(dataset, file, protocol=2)
	
print(" Len Actual Views ", len(actual_views))

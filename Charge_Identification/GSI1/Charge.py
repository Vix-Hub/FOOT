import ROOT as r 
import fedrarootlogon 
import uproot 
import numpy as np 
import sys 

sys.path.insert(0, "/home/baronunix/Scripts/")
from Clustering_Cosmici_Frammenti import *
import time

brick_id = "GSI1"
track_name = "b000111.2.0.0.trk.root"
file_name = "b11_vol.root"

#TAGLIO SU k0 e VR0_max
k0_min = 7
VR0_max = 25000
#ki_min = 2

charge_name = "output_b111_k0_min" + str(k0_min) + ".root"
Calcolo_Variabili_Volume_New(file_name, track_name) #calcolo VRi_av, ki

#trk_file = uproot.open(file_name)
#tracks = trk_file['tracks']
#track_info = tracks.arrays(['s.eX', 's.eY', 's.eTX', 's.eTY', 's.eScanID.ePlate', 's.eID', 'nseg', 'k0', 'k1', 'k2', 'k3', 'VR0_av', 'VR1_av', 'VR2_av', 'VR3_av', 'tan'])

#TAGLIO COSMICI
a0, b0 = 2600, 0.6
a1, b1 = 3850, 0.7
a2, b2 = (a1+a0)/2, 0.65

#TAGLI FRAMMENTI
frag_cut = "VR0_av >= " + str(a2) + "*(1+ TMath::Exp(" + str(b2) + " *tan*tan"

#TAGLIO VR1 vs VR0
c1 = 4750
c2 = 5000
c0 = (c1+c2)/2

Z_rec = np.zeros(1, np.intc)

t0 = time.time()

file2 = r.TFile(file_name, "READ")
tracks_2 = file2.Get("tracks")

## PCA ANALYSIS
file_name01 = '01.root'
tup_name01 = '01_c'
file_name012 = '012.root'
tup_name012 = '012_c'
file_name123 = '123.root'
tup_name123 = '123_c'
file_name013 = '013.root'
tup_name013 = '013_c'

################################# VP123
principal = r.TPrincipal(3, "ND")

file_pca = r.TFile("123t.root", "RECREATE")
cluster_v1 = tracks_2.CopyTree("VR0_av >= " + str(a2) + "*(1 + TMath::Exp(" + str(b2) + " *tan*tan)) ")
campione_pca = cluster_v1.CopyTree("k1 >=2 && k2>=3 && k3>=2 || k1>=2 && k2>=2 && k3>=3")
campione_pca.Write("pca")

principal = r.TPrincipal(3, "ND")
for track in campione_pca:
	vr0, vr1, vr2 = track.VR1_av, track.VR2_av, track.VR3_av
	vrs = np.zeros(3)
	vrs[0] = vr0
	vrs[1] = vr1
	vrs[2] = vr2
	principal.AddRow(vrs)
	
principal.MakePrincipals()
principal.MakeCode()
r.gInterpreter.ProcessLine('.L pca.C+')
r.gSystem.Load("pca_C.so")

vr123s = []
for track in campione_pca:
	vr2, vr3, vr0 = track.VR1_av, track.VR2_av, track.VR3_av
	vrs = np.zeros(3)
	vrs[0] = vr2
	vrs[1] = vr3
	vrs[2] = vr0
	princ = np.zeros(3)
	principal.X2P(vrs, princ)
	vr123s.append(princ)

vr123, vr123b = [], []
for i in vr123s:
	vr123.append(i[0])
	vr123b.append(i[1])

pca_1 = r.TNtuple("pca_1", "", "VR123:VR123b")

for i in range(len(vr123)):
	pca_1.Fill(vr123[i], vr123b[i]) 

pca_1.Write("V123")

## FIT 123
pca_1.Draw("VR123>>v123s(60, -4., 4.5.)", "", "COLZ")

r.gStyle.SetOptStat(0)
histos = r.gDirectory.Get("v123s")
histos.SetTitle("VP_{123} [GSI1]; VP_{123}; Entries")
fit_func = r.TF1("fit_func", "[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [4], [5]) + [6]*TMath::Gaus(x, [7], [8])", -4, 7)
fit_func.SetParameters(501, -1.5, .5, 110, 0.9, 0.7, 110, 3.5, 0.5)

#ampiezze
fit_func.SetParLimits(0, 400, 520)
fit_func.SetParLimits(3, 100, 200)
fit_func.SetParLimits(6, 100, 200)

#punto medio
fit_func.SetParLimits(1, -1.6, -1.3)
fit_func.SetParLimits(4, 0., 1.2)
fit_func.SetParLimits(7, 2., 4.)

#deviazione_st
fit_func.SetParLimits(2, 0., 2.)
fit_func.SetParLimits(5, 0.2, 2)
fit_func.SetParLimits(8, 0.3, 2.5)
tr = -1.8
histos.Fit("fit_func", "", "", tr, 4.5)
params = fit_func.GetParameters()

prova = make_classification_123_ALL(pca_1, fit_func, file_name123, tup_name123, 3, 2, 3, 4)
file_pca.Close()

##############################################################
### VP012
file_pca = r.TFile("012t.root", "RECREATE")
cluster_v1 = tracks_2.CopyTree("VR0_av >= " + str(a2) + "*(1 + TMath::Exp(" + str(b2) + " *tan*tan)) ")
campione_pca = cluster_v1.CopyTree("k1 >=2 && k2>=3")
campione_pca.Write("pca")

principal = r.TPrincipal(3, "ND")

for track in campione_pca:
	vr0, vr1, vr2 = track.VR0_av, track.VR1_av, track.VR2_av
	vrs = np.zeros(3)
	vrs[0] = vr0
	vrs[1] = vr1
	vrs[2] = vr2
	principal.AddRow(vrs)
	
principal.MakePrincipals()
principal.MakeCode()
r.gInterpreter.ProcessLine('.L pca.C+')
r.gSystem.Load("pca_C.so")

vr123_tris = []
for track in campione_pca:
	vr2, vr3, vr0 = track.VR0_av, track.VR1_av, track.VR2_av
	vrs = np.zeros(3)
	vrs[0] = vr2
	vrs[1] = vr3
	vrs[2] = vr0
	princ = np.zeros(3)
	principal.X2P(vrs, princ)
	vr123_tris.append(princ)

vr123_vr2_check = []
for i in vr123_tris:
	vr123_vr2_check.append(i[0])

pca_3 = r.TNtuple("pca_1", "vb", "VR012")

for i in range(len(vr123_vr2_check)):
	pca_3.Fill(vr123_vr2_check[i])
pca_3.Write()

pca_3.Draw("VR012>>v123s(75, -4., 7.)", "", "COLZ")
histos = r.gDirectory.Get("v123s")
histos.SetTitle("VP_{012} [k_{1}>=1 && k_{2}>=2]; VP_{012}; Entries")
fit_func = r.TF1("fit_func", "[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [4], [5]) + [6]*TMath::Gaus(x, [7], [8])", -4, 7)
fit_func.SetParameters(120, -1.5, .5, 60, 0.9, 0.7, 40, 3., 0.5)

#ampiezze
fit_func.SetParLimits(0, 300, 650)
fit_func.SetParLimits(3, 100, 350)
fit_func.SetParLimits(6, 10, 250)
#punto medio
fit_func.SetParLimits(1, -1.4, -1.)
fit_func.SetParLimits(4, -1, 2.)
fit_func.SetParLimits(7, 2., 4.)
#deviazione_st
fit_func.SetParLimits(2, 0., 2.5)
fit_func.SetParLimits(5, 0.2, 2)
fit_func.SetParLimits(8, 0.3, 2.5)
tr = -1.65

histos.Fit("fit_func", "S", "", tr, 7.)
params = fit_func.GetParameters()

prova = make_classification_012(pca_3, fit_func, 3, -1.65)
prova2 = make_classification_012_X(pca_3, fit_func, file_name012, tup_name012, 3, -1.65, 2, 3, 4)
file_pca.Close()
##############################################################
### VP013
file_pca = r.TFile("013t.root", "RECREATE")
cluster_v1 = tracks_2.CopyTree("VR0_av >= " + str(a2) + "*(1 + TMath::Exp(" + str(b2) + " *tan*tan)) ")
campione_pca = cluster_v1.CopyTree("k1 >=2 && k3>=2")
campione_pca.Write("pca")

principal = r.TPrincipal(3, "ND")

for track in campione_pca:
	vr0, vr1, vr2 = track.VR0_av, track.VR1_av, track.VR3_av
	vrs = np.zeros(3)
	vrs[0] = vr0
	vrs[1] = vr1
	vrs[2] = vr2
	principal.AddRow(vrs)
	
principal.MakePrincipals()
principal.MakeCode()
r.gInterpreter.ProcessLine('.L pca.C+')
r.gSystem.Load("pca_C.so")

vr123_tris = []

for track in campione_pca:
	vr2, vr3, vr0 = track.VR0_av, track.VR1_av, track.VR3_av
	vrs = np.zeros(3)
	vrs[0] = vr2
	vrs[1] = vr3
	vrs[2] = vr0
	princ = np.zeros(3)
	principal.X2P(vrs, princ)
	vr123_tris.append(princ)

vr123_vr2_check = []
for i in vr123_tris:
	vr123_vr2_check.append(i[0])

pca_3 = r.TNtuple("pca_1", "vb", "VR013")

for i in range(len(vr123_vr2_check)):
	pca_3.Fill(vr123_vr2_check[i])
	
pca_3.Write()

pca_3.Draw("VR013>>v123s(75, -4., 7.)", "", "COLZ")
histos = r.gDirectory.Get("v123s")
histos.SetTitle("VR013; VR013; Conteggi")
fit_func = r.TF1("fit_func", "[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [4], [5]) + [6]*TMath::Gaus(x, [7], [8])", -4, 7)
fit_func.SetParameters(120, -1.5, .5, 60, 0.9, 0.7, 40, 3., 0.5)

#ampiezze
fit_func.SetParLimits(0, 300, 650)
fit_func.SetParLimits(3, 100, 350)
fit_func.SetParLimits(6, 10, 250)
#punto medio
fit_func.SetParLimits(1, -2.5, -1.)
fit_func.SetParLimits(4, -1, 2.)
fit_func.SetParLimits(7, 2., 4.)
#deviazione_st
fit_func.SetParLimits(2, 0., 2.5)
fit_func.SetParLimits(5, 0.2, 2)
fit_func.SetParLimits(8, 0.3, 2.5)

tr = -1.653
histos.Fit("fit_func", "S", "", tr, 7.)
params = fit_func.GetParameters()

prova2 = make_classification_013_X(pca_3, fit_func, file_name013, tup_name013, 3, tr, 2, 3, 4)
file_pca.Close()

### CHARGE
## Definizione tree in uscita
file0 = r.TFile(track_name, "READ")
tracks = file0.Get("tracks")

file_out = r.TFile(charge_name, "RECREATE")
Z_rec = np.zeros(1, dtype = np.intc)
output_tree = tracks.CloneTree(0)

#tracks.SetBranchAddress("n0", Z_rec)
output_tree.Branch("Z", Z_rec, "Z/I")
tracks.BuildIndex("trid")

## tree con variabili di volume
nseg0 = np.zeros(1, dtype = np.intc)
nseg1 = np.zeros(1, dtype = np.intc)
nseg2 = np.zeros(1, dtype = np.intc)
nseg3 = np.zeros(1, dtype = np.intc)

volume0 = np.zeros(1, dtype = np.double)
volume1 = np.zeros(1, dtype = np.double)
volume2 = np.zeros(1, dtype = np.double)
volume3 = np.zeros(1, dtype = np.double)

output_tree2 = tracks.CloneTree(0)
output_tree2.Branch("VR0_av", volume0, "VR0_av/D")
output_tree2.Branch("VR1_av", volume1, "VR1_av/D")
output_tree2.Branch("VR2_av", volume2, "VR2_av/D")
output_tree2.Branch("VR3_av", volume3, "VR3_av/D")

output_tree2.Branch("k0", nseg0, "k0/I")
output_tree2.Branch("k1", nseg1, "k1/I")
output_tree2.Branch("k2", nseg2, "k2/I")
output_tree2.Branch("k3", nseg3, "k3/I")

output_tree2.Branch("Z", Z_rec, "Z/I")
tracks_2.BuildIndex("trid", "npl")
i = 0
t0 = time.time()
for track in tracks_2:
	k0, k1, k2, k3 = track.k0, track.k1, track.k2, track.k3
	VR0_av, VR1_av, VR2_av, VR3_av = track.VR0_av, track.VR1_av, track.VR2_av, track.VR3_av
	tan = track.tan

	trid_to_assign, npl_to_assign = track.trid, track.npl
	tracks_2.GetEntryWithIndex(trid_to_assign, npl_to_assign)
	tracks.GetEntry(tracks_2.trid)

	if (k0 >= k0_min and VR0_av <= VR0_max):
		if(VR0_av >= a2*(1+np.exp(b2*tan*tan))):
			if (k1<2 and k2<2 and k3<2):
				Z_rec[0] = 1 #high energy protons
				output_tree.Fill()
				output_tree2.Fill()
			elif(k1>0 and VR1_av<=c0 and k2<2 and k3<2):
				Z_rec[0] = 1 #low energy protons
				output_tree.Fill()
				output_tree2.Fill()
			elif(k1>0 and VR1_av>c0 and k2<2 and k3<2):
				Z_rec[0] = 2 #high energy Z=2
				output_tree.Fill()
				output_tree2.Fill()
			else:
				i = i - 1 #non conto
				continue # -> PCA
		elif(VR0_av < a2*(1+np.exp(b2*tan*tan))):
			Z_rec[0] = -1 #cosmic rays
			output_tree.Fill()
			output_tree2.Fill()
	else:
		Z_rec[0] = -2
		output_tree.Fill()
		output_tree2.Fill()
	i = i+1
	if (i%1000 == 0):
			print("Carica assegnata a " + str(i) + " tracce, su " + str(tracks_2.GetEntries()))
	

file_pca123 = r.TFile(file_name123, "READ")
pca_3 = file_pca123.Get(tup_name123)
file_infopca123 = r.TFile("123t.root", "READ")
info_pca123 = file_infopca123.Get("pca")

pca_3.AddFriend(info_pca123)
tracks_2.BuildIndex("trid", "npl")
count3 = 0

for track_123 in pca_3:
	
	trid_to_assign, npl_to_assign = track_123.trid, track_123.npl
	tracks_2.GetEntryWithIndex(trid_to_assign, npl_to_assign)
	tracks.GetEntry(tracks_2.trid)
	
	check0 = track_123.k2>1 and track_123.k3>1
	check1 = track_123.tan == tracks_2.tan and track_123.k0 == tracks_2.k0 and track_123.VR0_av == tracks_2.VR0_av
	check2 = track_123.k1 == tracks_2.k1 and track_123.VR1_av == tracks_2.VR1_av
	check3 = track_123.k2 == tracks_2.k2 and track_123.VR2_av == tracks_2.VR2_av
	check4 = track_123.k3 == tracks_2.k3 and track_123.VR3_av == tracks_2.VR3_av
	check5 = track_123.nseg == tracks_2.nseg and track_123.npl == tracks_2.npl
	check = check1 and check2 and check3 and check4 and check5 and check0
	
	nseg0[0], nseg1[0], nseg2[0], nseg3[0] = tracks_2.k0, tracks_2.k1, tracks_2.k2, tracks_2.k3
	volume0[0], volume1[0], volume2[0], volume3[0] = tracks_2.VR0_av, tracks_2.VR1_av, tracks_2.VR2_av, tracks_2.VR3_av
	
	if (check):
		Z_rec[0] = int(track_123.Z_c)
		output_tree.Fill()
		output_tree2.Fill()
		count3 = count3 +1
	
	if (count3%1000 == 0):
		print("Carica assegnata a " + str(i+count3) + " tracce, su " + str(tracks_2.GetEntries()))
file_pca123.Close()
file_infopca123.Close()


file_pca012 = r.TFile(file_name012, "READ")
pca_012 = file_pca012.Get(tup_name012)
file_infopca012 = r.TFile("012t.root", "READ")
info_pca012 = file_infopca012.Get("pca")


pca_012.AddFriend(info_pca012)
count4 = 0

for track_012 in pca_012:
	
	countx = 0
	if (track_012.k2 > 1 and track_012.k3<=1):
		trid_to_assign, npl_to_assign = track_012.trid, track_012.npl
		tracks_2.GetEntryWithIndex(trid_to_assign, npl_to_assign)
		tracks.GetEntry(tracks_2.trid)
	
		check1 = track_012.tan == tracks_2.tan and track_012.k0 == tracks_2.k0 and track_012.VR0_av == tracks_2.VR0_av
		check2 = track_012.k1 == tracks_2.k1 and track_012.VR1_av == tracks_2.VR1_av
		check3 = track_012.k2 == tracks_2.k2 and track_012.VR2_av == tracks_2.VR2_av
		check4 = track_012.k3 == tracks_2.k3 and track_012.VR3_av == tracks_2.VR3_av
		check5 = track_012.nseg == tracks_2.nseg and track_012.npl == tracks_2.npl
		check = check1 and check2 and check3 and check4 and check5
		
		nseg0[0], nseg1[0], nseg2[0], nseg3[0] = tracks_2.k0, tracks_2.k1, tracks_2.k2, tracks_2.k3
		volume0[0], volume1[0], volume2[0], volume3[0] = tracks_2.VR0_av, tracks_2.VR1_av, tracks_2.VR2_av, tracks_2.VR3_av
		
		if (check):
			Z_rec[0] = int(track_012.Z_012) 
			output_tree.Fill()
			output_tree2.Fill()
			countx = countx + 1
			count4 = count4 +1
		
		if (count4%100 == 0):
		   print("Carica assegnata a " + str(i+count3+count4) + " tracce, su " + str(tracks_2.GetEntries()))
file_pca012.Close()
file_infopca012.Close()

file_pca013 = r.TFile(file_name013, "READ")
pca_013 = file_pca013.Get(tup_name013)
file_infopca013 = r.TFile("013t.root", "READ")
info_pca013 = file_infopca013.Get("pca")

pca_013.AddFriend(info_pca013)
#tracks_2.BuildIndex("trid", "npl")
count5 = 0

for track_012 in pca_013:
	
	if (track_012.k3 > 1 and track_012.k2<=1):
		trid_to_assign, npl_to_assign = track_012.trid, track_012.npl
		tracks_2.GetEntryWithIndex(trid_to_assign, npl_to_assign)
		tracks.GetEntry(tracks_2.trid)
	
		check1 = track_012.tan == tracks_2.tan and track_012.k0 == tracks_2.k0 and track_012.VR0_av == tracks_2.VR0_av
		check2 = track_012.k1 == tracks_2.k1 and track_012.VR1_av == tracks_2.VR1_av
		check3 = track_012.k2 == tracks_2.k2 and track_012.VR2_av == tracks_2.VR2_av
		check4 = track_012.k3 == tracks_2.k3 and track_012.VR3_av == tracks_2.VR3_av
		check5 = track_012.nseg == tracks_2.nseg and track_012.npl == tracks_2.npl
		check = check1 and check2 and check3 and check4 and check5
		
		nseg0[0], nseg1[0], nseg2[0], nseg3[0] = tracks_2.k0, tracks_2.k1, tracks_2.k2, tracks_2.k3
		volume0[0], volume1[0], volume2[0], volume3[0] = tracks_2.VR0_av, tracks_2.VR1_av, tracks_2.VR2_av, tracks_2.VR3_av
		
		if (check):
			Z_rec[0] = int(track_012.Z_013)
			#print(Z_rec)
			output_tree.Fill()
			output_tree2.Fill()
			count5 = count5 +1
		elif (check and track_012.VR2_av > 10000):
			Z_rec[0] = 11
			#print(Z_rec)
			output_tree.Fill()
			output_tree2.Fill()
			count5 = count5 +1
		#if (count == 0):
			#print("Error")
			#break
		
		if (count5%100 == 0):
		   print("Carica assegnata a " + str(i+count3+count4+count5) + " tracce, su " + str(tracks_2.GetEntries()))
file_pca013.Close()
file_infopca013.Close()

file_out.cd()
#output_tree.Write("tracks")
output_tree2.Write("tracks")   #con info su VRi_av e ki
file_out.Close()

print("Time needed (root): " + str(time.time()-t0))








		

import ROOT as r 
import fedrarootlogon 
import sys
sys.path.insert(0, "/home/baronunix/Scripts/")
from Clustering_Cosmici_Frammenti import *

track_name = "b000333.2.0.0.trk.root"
path = "~/Scripts/GSI3/"
file_name = "~/Scripts/GSI3/b33_vol.root"
file_name_cuts = "~/Scripts/GSI3/b33_vol_cuts.root"
brick_id = "GSI3"
charge_name = "output_b333_k0_min_1.root"
vtx_name = "output_b333_vtx.root"

## SET OPTIONS
ADD_VTX = 0  #includo nella classificazione anche le tracce collegate ai vertici *se non soddisfano normalmente i tagli
CALCULATE_VOLUME = 1
DO_PCA = 1

treename = "tracks_n"
if (ADD_VTX==1):
    treename = "tracks"
    file_name_cuts = "b33_vtx.root"
    charge_name = "output_v2_b333.root"
    file_name = path + vtx_name #prendo dal file con info di volume e vertici

VERBOSE = 0 
if (VERBOSE == 1):
	file_out = r.TFile("Check_Plots.root", "RECREATE")

#CALCOLO VARIABILI VOLUME
print(" --- Calculating Volume Variables --- ")
if (CALCULATE_VOLUME == 1):
    Calcolo_Variabili_Volume_New(file_name, track_name)
print( " --- Done --- ")

#TAGLIO k0 / AGGIUNTA TRACCE VERTICI
k0_min = 1
VR0_max = 25000

file1 = r.TFile(file_name, "READ")
tracks_V = file1.Get(treename)
file2 = r.TFile(file_name_cuts, "RECREATE")

if (ADD_VTX == 0):
    tracks_V2 = tracks_V.CopyTree("k0>=" + str(k0_min) + " && VR0_av<"+str(VR0_max))
    tracks_V2.Write("tracks_cuts")
else:
    tracks_V2 = tracks_V.CopyTree("k0>=" + str(k0_min) + " && VR0_av<"+str(VR0_max) + " || in_vtx == 1")
    tracks_V2.Write("tracks_cuts")
    
file2.Close()
file1.Close()

file2 = r.TFile(file_name_cuts, "READ")
tracks_2 = file2.Get("tracks_cuts")

#TAGLIO COSMICI
a2, b2 = 2400, 0.95  #2400, 0.95 x k0>=4
c2, d2 = 100, .45
frag_cut = "VR0_av >= " + str(a2) + "*(1 + TMath::Exp(" + str(b2) + " *s[0].Theta()*s[0].Theta()))"
nfrag_cut = "VR0_av < " + str(a2) + "*(1 + TMath::Exp(" + str(b2) + " *s[0].Theta()*s[0].Theta()))"

#NOMI FILE PCA
file_name01 = '01.root'
tup_name01 = '01_c'
file_name012 = '012.root'
tup_name012 = '012_c'
file_name123 = '123.root'
tup_name123 = '123_c'
file_name013 = '013.root'
tup_name013 = '013_c'

################################################################ PCA: VP01
print(" --- Calculating PCA variables and fitting --- ")
if (DO_PCA):
    do_PCA_VP01(tracks_2, VERBOSE, a2, b2, brick_id, file_name01, tup_name01)

################################################################ PCA: VP123
    do_PCA_VP123(tracks_2, VERBOSE, a2, b2, brick_id, file_name123, tup_name123)

################################################################ PCA: VP012
    do_PCA_VP012(tracks_2, VERBOSE, a2, b2, brick_id, file_name012, tup_name012)

################################################################ PCA: VP013
    do_PCA_VP013(tracks_2, VERBOSE, a2, b2, brick_id, file_name013, tup_name013)

print(" --- Created PCA files --- ")

######################################################### SCRITTURA CARICA OUTPUT

file_b = r.TFile(track_name, "READ")
tracks = file_b.Get("tracks")
tracks.BuildIndex("trid")

## Definizione tree in uscita

file_out = r.TFile(charge_name, "RECREATE")

Z_rec = np.zeros(1, dtype = np.intc)
in_vtx = np.zeros(1, dtype = np.intc)
vid = np.zeros(1, dtype = np.intc)
vz = np.zeros(1, dtype = np.single)

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
if (ADD_VTX == 1):
    output_tree2.Branch("in_vtx", in_vtx, "in_vtx/I")
    output_tree2.Branch("vID", vid, "vID/I")
    output_tree2.Branch("vz", vz, "vz/F")

out_tree, vtx_info = 0, 0

if (ADD_VTX == 1):
    vtx_file = r.TFile(vtx_name, 'READ')
    out_tree = vtx_file.Get('tracks')
    out_tree.BuildIndex('s[0].Plate()', 's[0].ID()')

    vtx_info_file = r.TFile("vtx_tracks.root", "READ")
    vtx_info = vtx_info_file.Get("vtx_info")
    vtx_info.BuildIndex("s0_id", "s0_plate")

for i in range(tracks.GetEntries()):
    tracks.GetEntry(i)
    k0, k1, k2, k3 = 0, 0, 0, 0
    VR0av, VR1av, VR2av, VR3av = 0, 0, 0, 0
    vr0, vr1, vr2, vr3 = 0, 0, 0, 0
    for s in tracks.s:
        if((s.Plate()-31)%4 - 0 == 0):
            k0+=1
            vr0 += s.Volume()
        if((s.Plate()-31)%4 - 1 == 0):
            k1+=1
            vr1 += s.Volume()   
        if((s.Plate()-31)%4 - 2 == 0):
            k2+=1
            vr2 += s.Volume()
        if((s.Plate()-31)%4 - 3 == 0):
            k3+=1
            vr3 += s.Volume()      
    if(k0!=0):
            VR0av = vr0/k0
    if(k1!=0):
            VR1av = vr1/k1
    if(k2!=0):
            VR2av = vr2/k2
    if(k3!=0):
            VR3av = vr3/k3
            
    nseg0[0], nseg1[0], nseg2[0], nseg3[0] = k0, k1, k2, k3
    volume0[0], volume1[0], volume2[0], volume3[0] = VR0av, VR1av, VR2av, VR3av
    cond = k0<k0_min or VR0av>=VR0_max
    in_vtx0 = 0
    vz0, vid0 = -1, -1
    id0, plate0 = tracks.s[0].ID(), tracks.s[0].Plate()

    if (ADD_VTX==1):  #check if entry is connected to vertex
        if(vtx_info.GetEntryWithIndex(int(id0), int(plate0))>=0):
            in_vtx0 = 1
            vz0 = vtx_info.vz
            vid0 = vtx_info.vid
        cond = cond and (in_vtx0 == 0)  
        in_vtx[0] = in_vtx0
        vz[0], vid[0] = vz0, vid0

    if (cond):
        Z_rec[0] = -1
        output_tree.Fill()
        output_tree2.Fill()
    if (i%10000 ==0):
        print(" Completed " + str(round(100*i/tracks.GetEntries(), 2) ) + " %")

count = 0
g1, g2, g3, g4, g5, g6, g7 = 0, 0, 0, 0, 0, 0, 0
check, check2 = [], []
for i in range(tracks_2.GetEntries()): #tracks_2.GetEntries()
    tracks_2.GetEntry(i)
    tracks.GetEntry(tracks_2.trid)
    if (ADD_VTX==1):
        #print(" ### ENTRY " + str(i))
        #print(" Tracks entry S0Plate S0ID " + str(tracks.s[0].Plate()) + " " + str(tracks.s[0].ID()))
        #print(" Tracks 2 TRID " + str(tracks_2.trid))
        out_tree.GetEntryWithIndex(tracks.s[0].Plate(), tracks.s[0].ID())
        #print(" s0x " + str(out_tree.s[0].X()))
        in_vtx[0] = out_tree.in_vtx
        vz[0], vid[0] = out_tree.vz, out_tree.v_ID
        #print(" invtx " + str(out_tree.in_vtx))
        #print(" ")
        check.append(in_vtx[0])
        if (in_vtx[0] == 1):
            check2.append(in_vtx[0])
    frag_cut = (tracks_2.VR0_av > a2*(1 + np.exp(b2 * tracks_2.s[0].Theta()*tracks_2.s[0].Theta())))
    cond1 = (tracks_2.k2 <=1 and tracks_2.k3<=1)
    cond2 = tracks_2.VR1_av>0
    cond3 = tracks_2.VR1_av>6000
    cond0 = tracks_2.VR1_av == 0
    cond4 = tracks_2.k2 == 0
    cond_vp123 = tracks_2.k1>0 and tracks_2.k2>1 and tracks_2.k3>1
    cond_vp012 = tracks_2.k1>0 and tracks_2.k2>1 and tracks_2.k3<=1
    cond_vp013 = tracks_2.k1>0 and tracks_2.k2<=1 and tracks_2.k3>1
    
    nseg0[0], nseg1[0], nseg2[0], nseg3[0] = tracks_2.k0, tracks_2.k1, tracks_2.k2, tracks_2.k3
    volume0[0], volume1[0], volume2[0], volume3[0] = tracks_2.VR0_av, tracks_2.VR1_av, tracks_2.VR2_av, tracks_2.VR3_av
    
    if ((i%10000)==0):
        print("Carica Assegnata a " + str(round(count,2)) + " tracce, su " + str(tracks_2.GetEntries()))
    if (frag_cut and cond0):
        Z_rec[0] = 1
        output_tree.Fill()
        output_tree2.Fill()
        count = count + 1
        g1 = g1 + 1
    elif (frag_cut and cond1 and cond2 and tracks_2.VR1_av<7500):
        #Z_rec[0] = 2
        #print("Z= "+str(Z_rec))
        #output_tree.Fill()
        #count = count + 1
        g2 = g2 + 1 
        continue
    elif (frag_cut and cond1 and cond2 and tracks_2.VR1_av>7500):
        Z_rec[0] = 11 
        #print("Z= "+str(Z_rec))
        output_tree.Fill()
        output_tree2.Fill()
        count = count + 1
        g3 = g3 +1
    #elif (frag_cut and cond4 and cond2):
        #Z_rec[0] = 1
        #output_tree.Fill()
        #count = count + 1
    elif (not frag_cut):
        Z_rec[0] = -2
        #print("Z=" +str(Z_rec))
        output_tree.Fill()
        output_tree2.Fill()
        count = count + 1
        g4 = g4 + 1
    elif (cond_vp123):
        g5 = g5 + 1
        continue
    elif (cond_vp012):
        g6 = g6 + 1
        continue
    elif (cond_vp013):
        g7 = g7 +1
        continue
    else:
        Z_rec[0] = 10
        output_tree.Fill()
        output_tree2.Fill()
    if (count == 0):
        print("Errore")
    #print("count= " + str(count))
print(g1, g2, g3, g4, g5 ,g6, g7)   
#print(" LEN CHECK = " + str(len(check))) 
#print(" LEN CHECK 2 = " + str(len(check2)))

file_pca01 = r.TFile(file_name01, "READ")
pca_3 = file_pca01.Get(tup_name01)

file_info_pca01 = r.TFile("PCA_01.root", "READ")
info_pca01 = file_info_pca01.Get("pca")

pca_3.AddFriend(info_pca01)
tracks_2.BuildIndex("trid", "npl")
count2 = 0
m1, m2 = 0, 0

for track_123 in pca_3:
    
    trid_to_assign, npl_to_assign = track_123.trid, track_123.npl
    tracks_2.GetEntryWithIndex(trid_to_assign, npl_to_assign)
    tracks.GetEntry(tracks_2.trid)
    if (ADD_VTX==1):
        out_tree.GetEntryWithIndex(tracks.s[0].Plate(), tracks.s[0].ID()) 
        in_vtx[0] = out_tree.in_vtx
        vz[0], vid[0] = out_tree.vz, out_tree.v_ID
        
    
    #check0 = track_123.k2>1 and track_123.k3>1
    check1 = track_123.s[0].Theta() == tracks_2.s[0].Theta() and track_123.k0 == tracks_2.k0 and track_123.VR0_av == tracks_2.VR0_av
    check2 = track_123.k1 == tracks_2.k1 and track_123.VR1_av == tracks_2.VR1_av
    check3 = track_123.k2 == tracks_2.k2 and track_123.VR2_av == tracks_2.VR2_av
    check4 = track_123.k3 == tracks_2.k3 and track_123.VR3_av == tracks_2.VR3_av
    check5 = track_123.nseg == tracks_2.nseg and track_123.npl == tracks_2.npl
    check = check1 and check2 and check3 and check4 and check5
    
    nseg0[0], nseg1[0], nseg2[0], nseg3[0] = tracks_2.k0, tracks_2.k1, tracks_2.k2, tracks_2.k3
    volume0[0], volume1[0], volume2[0], volume3[0] = tracks_2.VR0_av, tracks_2.VR1_av, tracks_2.VR2_av, tracks_2.VR3_av
    
    if (check):
        if (int(track_123.Z_c == 1)):
            Z_rec[0] = 2 #sono invertiti in questa classificazione
            #print(Z_rec)
            output_tree.Fill()
            output_tree2.Fill()
            m2 = m2 + 1
            count2 = count2 +1
        elif(int(track_123.Z_c == 2) ):
            Z_rec[0] = 1  #sono invertiti in questa classificazione
            #print(Z_rec)
            output_tree.Fill()
            output_tree2.Fill()
            m1 = m1 +1 
            count2 = count2 +1
    #if (count == 0):
        #print("Error")
        #break
    
    if (count2%1000 == 0):
        print("Carica assegnata a " + str(count2+count) + " tracce, su " + str(tracks_2.GetEntries()))
file_pca01.Close()
file_info_pca01.Close()
print("100%")
print(count2)
print(m1, m2)

file_pca123 = r.TFile(file_name123, "READ")
pca_3 = file_pca123.Get(tup_name123)

file_info_pca123 = r.TFile("PCA_123.root", "READ")
info_pca123 = file_info_pca123.Get("pca")

pca_3.AddFriend(info_pca123)
tracks_2.BuildIndex("trid", "npl")
count3 = 0

for track_123 in pca_3:
    
    trid_to_assign, npl_to_assign = track_123.trid, track_123.npl
    tracks_2.GetEntryWithIndex(trid_to_assign, npl_to_assign)
    tracks.GetEntry(tracks_2.trid)

    if (ADD_VTX==1):
        out_tree.GetEntryWithIndex(tracks.s[0].Plate(), tracks.s[0].ID()) 
        in_vtx[0] = out_tree.in_vtx
        vz[0], vid[0] = out_tree.vz, out_tree.v_ID
    
    check0 = track_123.k2>1 and track_123.k3>1
    check1 = track_123.s[0].Theta() == tracks_2.s[0].Theta() and track_123.k0 == tracks_2.k0 and track_123.VR0_av == tracks_2.VR0_av
    check2 = track_123.k1 == tracks_2.k1 and track_123.VR1_av == tracks_2.VR1_av
    check3 = track_123.k2 == tracks_2.k2 and track_123.VR2_av == tracks_2.VR2_av
    check4 = track_123.k3 == tracks_2.k3 and track_123.VR3_av == tracks_2.VR3_av
    check5 = track_123.nseg == tracks_2.nseg and track_123.npl == tracks_2.npl
    check = check1 and check2 and check3 and check4 and check5 and check0
    
    nseg0[0], nseg1[0], nseg2[0], nseg3[0] = tracks_2.k0, tracks_2.k1, tracks_2.k2, tracks_2.k3
    volume0[0], volume1[0], volume2[0], volume3[0] = tracks_2.VR0_av, tracks_2.VR1_av, tracks_2.VR2_av, tracks_2.VR3_av
    
    if (check):
        Z_rec[0] = int(track_123.Z_c)
        #print(Z_rec)
        output_tree.Fill()
        output_tree2.Fill()
        count3 = count3 +1
    #if (count == 0):
        #print("Error")
        #break
    
    if (count3%1000 == 0):
        print("Carica assegnata a " + str(count2+count+count3) + " tracce, su " + str(tracks_2.GetEntries()))
file_pca123.Close()
file_info_pca123.Close()
print("100%")
print(count3)


file_pca012 = r.TFile(file_name012, "READ")
pca_012 = file_pca012.Get(tup_name012)

file_info_pca012 = r.TFile("PCA_012.root", "READ")
info_pca012 = file_info_pca012.Get("pca")

pca_012.AddFriend(info_pca012)
#tracks_2.BuildIndex("trid", "npl")
count4 = 0

for track_012 in pca_012:
    
    countx = 0
    if (track_012.k2 > 1 and track_012.k3<=1):
        trid_to_assign, npl_to_assign = track_012.trid, track_012.npl
        tracks_2.GetEntryWithIndex(trid_to_assign, npl_to_assign)
        tracks.GetEntry(tracks_2.trid)

        if (ADD_VTX==1):
            out_tree.GetEntryWithIndex(tracks.s[0].Plate(), tracks.s[0].ID()) 
            in_vtx[0] = out_tree.in_vtx
            vz[0], vid[0] = out_tree.vz, out_tree.v_ID
    
        check1 = track_012.s[0].Theta() == tracks_2.s[0].Theta() and track_012.k0 == tracks_2.k0 and track_012.VR0_av == tracks_2.VR0_av
        check2 = track_012.k1 == tracks_2.k1 and track_012.VR1_av == tracks_2.VR1_av
        check3 = track_012.k2 == tracks_2.k2 and track_012.VR2_av == tracks_2.VR2_av
        check4 = track_012.k3 == tracks_2.k3 and track_012.VR3_av == tracks_2.VR3_av
        check5 = track_012.nseg == tracks_2.nseg and track_012.npl == tracks_2.npl
        check = check1 and check2 and check3 and check4 and check5
        
        nseg0[0], nseg1[0], nseg2[0], nseg3[0] = tracks_2.k0, tracks_2.k1, tracks_2.k2, tracks_2.k3
        volume0[0], volume1[0], volume2[0], volume3[0] = tracks_2.VR0_av, tracks_2.VR1_av, tracks_2.VR2_av, tracks_2.VR3_av
        
        if (check and track_012.VR2_av < 8000):
            Z_rec[0] = int(track_012.Z_012) 
            #print(Z_rec)
            output_tree.Fill()
            output_tree2.Fill()
            countx = countx + 1
            count4 = count4 +1
        elif (check and track_012.VR2_av > 8000):
            Z_rec[0] = 11
            #print(Z_rec)
            output_tree.Fill()
            output_tree2.Fill()
            countx = countx + 1
            count4 = count4 +1
        #if (count == 0):
            #print("Error")
            #break
        
        if (count4%100 == 0):
           print("Carica assegnata a " + str(count2+count+count3+count4) + " tracce, su " + str(tracks_2.GetEntries()))
file_pca012.Close()
file_info_pca012.Close()
print("100%")
print(count4)

file_pca013 = r.TFile(file_name013, "READ")
pca_013 = file_pca013.Get(tup_name013)

file_info_pca013 = r.TFile("PCA_013.root", "READ")
info_pca013 = file_info_pca013.Get("pca")

pca_013.AddFriend(info_pca013)
#tracks_2.BuildIndex("trid", "npl")
count5 = 0

for track_012 in pca_013:
    
    if (track_012.k3 > 1 and track_012.k2<=1):
        trid_to_assign, npl_to_assign = track_012.trid, track_012.npl
        tracks_2.GetEntryWithIndex(trid_to_assign, npl_to_assign)
        tracks.GetEntry(tracks_2.trid)

        if (ADD_VTX==1):
            out_tree.GetEntryWithIndex(tracks.s[0].Plate(), tracks.s[0].ID()) 
            in_vtx[0] = out_tree.in_vtx
            vz[0], vid[0] = out_tree.vz, out_tree.v_ID
    
        check1 = track_012.s[0].Theta() == tracks_2.s[0].Theta() and track_012.k0 == tracks_2.k0 and track_012.VR0_av == tracks_2.VR0_av
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

file_pca013.Close()
file_info_pca013.Close()
print("100%")
print(count5)

file_out.cd()
output_tree2.Write("tracks")   #con info su VRi_av e ki
file_out.Close()




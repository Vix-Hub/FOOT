
import ROOT as r
import fedrarootlogon 
from matplotlib import pyplot as plt
import matplotlib
import awkward as ak
import uproot
import numpy as np
from numpy import where
from sklearn.datasets import make_classification


def Bhattacharyya_Distance(modello):
    mu1, mu2 = modello.means_[0], modello.means_[1]
    Sigma1, Sigma2 = modello.covariances_[0], modello.covariances_[1]
    Sigma = (Sigma1 + Sigma2)/2
    D =1/8* np.matmul( np.matmul((mu1-mu2).T, np.linalg.inv(Sigma)), (mu1-mu2)  )+ 0.5*np.log( np.linalg.det(Sigma)/ np.sqrt(np.linalg.det(Sigma1)*np.linalg.det(Sigma2)) )
    return D


def Calcolo_Variabili_Volume_New(file_name, trk_name):
    
    track_file = r.TFile(trk_name, "READ")
    tracks = track_file.Get("tracks") 
    
    VR0, VR1, VR2, VR3 = [], [], [], []

    k0s, k1s, k2s, k3s = [], [], [], [] 
    tan = []
    
    file = r.TFile(file_name, "RECREATE") 
    
    vr0n = np.zeros(1, dtype = np.float32)
    vr1n = np.zeros(1, dtype = np.float32)
    vr2n = np.zeros(1, dtype = np.float32)
    vr3n = np.zeros(1, dtype = np.float32)
    tann = np.zeros(1, dtype = np.float32)
    k0n = np.zeros(1, dtype = np.intc)
    k1n = np.zeros(1, dtype = np.intc)
    k2n = np.zeros(1, dtype = np.intc)
    k3n = np.zeros(1, dtype = np.intc)
    
    output_tree = tracks.CloneTree(0)
    output_tree.Branch("VR0_av", vr0n, "VR0_av/F")
    output_tree.Branch("VR1_av", vr1n, "VR1_av/F")
    output_tree.Branch("VR2_av", vr2n, "VR2_av/F")
    output_tree.Branch("VR3_av", vr3n, "VR3_av/F")
    output_tree.Branch("k0", k0n, "k0/I")
    output_tree.Branch("k1", k1n, "k1/I")
    output_tree.Branch("k2", k2n, "k2/I")
    output_tree.Branch("k3",k3n, "k3/I")
    output_tree.Branch("tan", tann, "tan/F")
    

    for track in tracks:
        plate_0 = 31
        vr0, vr1, vr2, vr3 = 0,0,0,0
        k0, k1, k2, k3 = 0,0,0,0
        for s in track.s:
            if((s.Plate()-plate_0)%4 - 0 == 0):
                vr0 += s.Volume()
                k0+=1
            elif((s.Plate()-plate_0)%4 - 1 == 0):
                vr1 += s.Volume()
                k1+=1
            elif((s.Plate()-plate_0)%4 - 2 == 0):
                vr2 += s.Volume()
                k2+=1
            elif((s.Plate()-plate_0)%4 - 3 == 0):
                vr3+= s.Volume()
                k3+=1
        if(k0!=0):
            VR0.append(vr0/k0)
        else:
            VR0.append(0)
        if(k1!=0):
            VR1.append(vr1/k1)
        else:
            VR1.append(0)
        if(k2!=0):
            VR2.append(vr2/k2)
        else:
            VR2.append(0)
        if(k3!=0):
            VR3.append(vr3/k3)
        else:
            VR3.append(0)
        k0s.append(k0)
        k1s.append(k1)
        k2s.append(k2)
        k3s.append(k3)
        tan.append(np.sqrt(track.s[0].TX()*track.s[0].TX() + track.s[0].TY()*track.s[0].TY() ))
        
    for i in range(len(VR0)):
        tracks.GetEntry(i)
        
        vr0n[0] = VR0[i]
        vr1n[0] = VR1[i]
        vr2n[0] = VR2[i]
        vr3n[0] = VR3[i]
        
        k0n[0] = k0s[i]
        k1n[0] = k1s[i]
        k2n[0] = k2s[i]
        k3n[0] = k3s[i]
        
        tann[0] = tan[i]
        
        output_tree.Fill()
        
    output_tree.Write("tracks_n")
    file.Close()
    track_file.Close()
    
    return 1



def Calcolo_Variabili_Volume1(file_name, trk_name, k0_min = 0, k1_min = 0, k2_min = 0, k3_min = 0):
    
    track_file = r.TFile(trk_name, "READ")
    tracks = track_file.Get("tracks") 
    
    VR0, VR1, VR2, VR3 = [], [], [], []

    k0s, k1s, k2s, k3s = [], [], [], [] 
    tan = []

    for track in tracks:
        plate_0 = 31
        vr0, vr1, vr2, vr3 = 0,0,0,0
        k0, k1, k2, k3 = 0,0,0,0
        for s in track.s:
            if((s.Plate()-plate_0)%4 - 0 == 0):
                vr0 += s.Volume()
                k0+=1
            elif((s.Plate()-plate_0)%4 - 1 == 0):
                vr1 += s.Volume()
                k1+=1
            elif((s.Plate()-plate_0)%4 - 2 == 0):
                vr2 += s.Volume()
                k2+=1
            elif((s.Plate()-plate_0)%4 - 3 == 0):
                vr3+= s.Volume()
                k3+=1
        if(k0!=0):
            VR0.append(vr0/k0)
        else:
            VR0.append(0)
        if(k1!=0):
            VR1.append(vr1/k1)
        else:
            VR1.append(0)
        if(k2!=0):
            VR2.append(vr2/k2)
        else:
            VR2.append(0)
        if(k3!=0):
            VR3.append(vr3/k3)
        else:
            VR3.append(0)
        k0s.append(k0)
        k1s.append(k1)
        k2s.append(k2)
        k3s.append(k3)
        tan.append(np.sqrt(track.s[0].TX()*track.s[0].TX() + track.s[0].TY()*track.s[0].TY() ))
    
    file = r.TFile(file_name, "RECREATE")   
    t2s = r.TNtuple("Tupla", "VR_i","VR0_av:VR1_av:VR2_av:VR3_av:tan:k0:k1:k2:k3")
    t2n = r.TNtuple("Tupla_k0>=" + str(k0_min) + str(k1_min) + str(k2_min) + str(k3_min), "VR_i","VR0_av:VR1_av:VR2_av:VR3_av:tan:k0:k1:k2:k3")

    entries = []
    for i in range(len(VR0)):
        t2s.Fill(VR0[i], VR1[i], VR2[i], VR3[i], tan[i],  k0s[i], k1s[i], k2s[i], k3s[i] )
        if (k0s[i] >= k0_min and k1s[i] >= k1_min and k2s[i] >= k2_min and k3s[i] >= k3_min):
            t2n.Fill(VR0[i], VR1[i], VR2[i], VR3[i], tan[i],  k0s[i], k1s[i], k2s[i], k3s[i])
            entries.append(i)

    s, nseg, npl, n0 = [], [], [], []
    for k in entries:
        tracks.GetEntry(k)
        s.append(tracks.s)
        nseg.append(tracks.nseg)
        npl.append(tracks.npl)
        n0.append(tracks.n0)
        
    seg_info = r.TNtuple("Seg_info", "seg", "nseg:npl:n0")
    for i in range(len(entries)):
        seg_info.Fill(nseg[i], npl[i], n0[i])
        
    s0 = s[0]
    segment_tree = r.TTree("s_tree", "t")
    segment_tree.Branch("s", s0)
    for i in range(len(s)):
        s0 = s[i]
        segment_tree.Fill()
        
    segment_tree.Write()
    
    seg_info.Write()
    
    t2s.Write()
    t2n.Write()
    file.Close()
    
    return 1



def Clustering_GM1(file_name, variables, k0_min = 0, k1_min= 0, k2_min= 0, k3_min = 0, want_plots=0, n_clust=0):
    
    file2 = uproot.open(file_name)
    vrs2 = file2['Tupla_k0>='+str(k0_min) + str(k1_min) + str(k2_min) + str(k3_min)]
    
    y = vrs2.arrays(variables[0])[:]
    x = vrs2.arrays(variables[1])[:]

    x, y = ak.to_numpy(x), ak.to_numpy(y)
    x, y = list(x), list(y)

    x_n, y_n = [], []
    for i in range(len(x)):
        x_n.append(x[i][0])
        y_n.append(y[i][0])

    X, cnt, c0_s = [], [], []
    for i in range(len(x_n)):
        X.append([x_n[i], y_n[i]])
        cnt.append(i)
    
    from sklearn.mixture import GaussianMixture
    model = GaussianMixture(n_components=2, covariance_type = 'full', random_state = 5)
    model.fit(X)
    
    yhat = model.predict(X)
    clusters = np.unique(yhat)

    cluster1_x, cluster1_y = [], [] 
    cluster2_x, cluster2_y = [], [] 

    for i in range(len(cnt)):
        if (yhat[i] == 0):
            cluster1_x.append(x_n[cnt[i]])
            cluster1_y.append(y_n[cnt[i]])
        else:
            cluster2_x.append(x_n[cnt[i]])
            cluster2_y.append(y_n[cnt[i]])
            
    if (want_plots==1):
        plt.hist2d(cluster1_x,cluster1_y, 40)
        ax = plt.gca()
        ax.set_ylim([0, max(y_n)])
        plt.title("Cluster 1")
        plt.xlabel(variables[1])
        plt.ylabel(variables[0])
        plt.show()
    
        plt.hist2d(cluster2_x, cluster2_y, 40)
        ax = plt.gca()
        ax.set_ylim([0, max(y_n)])
        plt.title("Cluster 2")
        plt.xlabel(variables[1])
        plt.ylabel(variables[0])
        plt.show()
    
    file1 = r.TFile(file_name, "UPDATE")
    cluster_info = r.TNtuple("clust_"+str(n_clust), "clu", "n_clust_"+str(n_clust))
    
    for i in range(len(x_n)):
        cluster_info.Fill(yhat[i])
  
    prova = model.sample(n_samples = 30000)
    Sim = r.TH2F("Sim_"+str(k0_min)+str(k1_min)+str(k2_min)+str(k3_min), "h", 40, 0, max(x_n), 40, 0, max(y_n))
    
    for i in range(len(prova[0])):
        d = prova[0][i]
        Sim.Fill(d[0],d[1])
    Sim.SetTitle("Distribuzione campionata da GMM (k0>=" + str(k0_min)+");" + variables[1]+";" + variables[0])
    cluster_info.Write()
    Sim.Write()
    file1.Close()
    return model



#versioni vecchie

def Calcolo_Variabili_Volume(k0_min):
    
    track_file = r.TFile("b000022.2.0.0.trk.root", "READ")
    tracks = track_file.Get("tracks") 
    
    VR0, VR1, VR2, VR3 = [], [], [], []

    k0s, k1s, k2s, k3s = [], [], [], [] 
    tan = []

    for track in tracks:
        plate_0 = 31
        vr0, vr1, vr2, vr3 = 0,0,0,0
        k0, k1, k2, k3 = 0,0,0,0
        for s in track.s:
            if((s.Plate()-plate_0)%4 - 0 == 0):
                vr0 += s.Volume()
                k0+=1
            elif((s.Plate()-plate_0)%4 - 1 == 0):
                vr1 += s.Volume()
                k1+=1
            elif((s.Plate()-plate_0)%4 - 2 == 0):
                vr2 += s.Volume()
                k2+=1
            elif((s.Plate()-plate_0)%4 - 3 == 0):
                vr3+= s.Volume()
                k3+=1
        if(k0!=0):
            VR0.append(vr0/k0)
        else:
            VR0.append(0)
        if(k1!=0):
            VR1.append(vr1/k1)
        else:
            VR1.append(0)
        if(k2!=0):
            VR2.append(vr2/k2)
        else:
            VR2.append(0)
        if(k3!=0):
            VR3.append(vr3/k3)
        else:
            VR3.append(0)
        k0s.append(k0)
        k1s.append(k1)
        k2s.append(k2)
        k3s.append(k3)
        tan.append(np.sqrt(track.s[0].TX()*track.s[0].TX() + track.s[0].TY()*track.s[0].TY() ))
    
    file = r.TFile("Variabili_Volume3.root", "RECREATE")   
    t2s = r.TNtuple("Tupla", "VR_i","VR0_av:VR1_av:VR2_av:VR3_av:tan:k0:k1:k2:k3")
    t2n = r.TNtuple("Tupla_k0>=" + str(k0_min), "VR_i","VR0_av:VR1_av:VR2_av:VR3_av:tan:k0:k1:k2:k3")

    for i in range(len(VR0)):
        t2s.Fill(VR0[i], VR1[i], VR2[i], VR3[i], tan[i],  k0s[i], k1s[i], k2s[i], k3s[i] )
        if (k0s[i] >= k0_min):
            t2n.Fill(VR0[i], VR1[i], VR2[i], VR3[i], tan[i],  k0s[i], k1s[i], k2s[i], k3s[i])

    tracks.AddFriend(t2s)
    tracks.Write()
    t2s.Write()
    t2n.Write()
    file.Close()
    
    return 1


def Clustering_GM(k0_min, want_plots=0):
    
    file = uproot.open("Variabili_Volume3.root")
    vrs2 = file['Tupla_k0>='+str(k0_min)]
    
    VR0_av = vrs2.arrays('VR0_av')[:]
    tan = vrs2.arrays('tan')[:]

    t, vr0 = ak.to_numpy(tan), ak.to_numpy(VR0_av)
    t, vr0 = list(t), list(vr0)

    vr0_n, t_n = [], []
    for i in range(len(vr0)):
        vr0_n.append(vr0[i][0])
        t_n.append(t[i][0])

    X, cnt, c0_s = [], [], []
    for i in range(len(t_n)):
        X.append([t_n[i], vr0_n[i]])
        cnt.append(i)
    
    from sklearn.mixture import GaussianMixture
    model = GaussianMixture(n_components=2, covariance_type = 'full', weights_init = [0.5, 0.5], random_state = 5)
    model.fit(X)
    
    yhat = model.predict(X)
    clusters = np.unique(yhat)

    cluster1_tan, cluster1_vr0 = [], [] 
    cluster2_tan, cluster2_vr0 = [], [] 

    for i in range(len(cnt)):
        if (yhat[i] == 0):
            cluster1_tan.append(t_n[cnt[i]])
            cluster1_vr0.append(vr0_n[cnt[i]])
        else:
            cluster2_tan.append(t_n[cnt[i]])
            cluster2_vr0.append(vr0_n[cnt[i]])
            
    if (want_plots==1):
        plt.hist2d(cluster1_tan,cluster1_vr0, 40)
        ax = plt.gca()
        ax.set_ylim([0, max(vr0_n)])
        plt.title("Cluster 1")
        plt.xlabel("tan(theta)")
        plt.ylabel("VR0_av")
        plt.show()
    
        plt.hist2d(cluster2_tan, cluster2_vr0, 40)
        ax = plt.gca()
        ax.set_ylim([0, max(vr0_n)])
        plt.title("Cluster 2")
        plt.xlabel("tan(theta)")
        plt.ylabel("VR0_av")
        plt.show()
    
    file = r.TFile("Variabili_Volume3.root", "UPDATE")
    cluster_info = r.TNtuple("clust_info", "", "n_clust")
    
    for i in range(len(t_n)):
        cluster_info.Fill(yhat[i])
  
    prova = model.sample(n_samples = 30000)
    Sim = r.TH2F("Sim_"+str(k0_min), "h", 40, 0, 1, 40, 0, 30000)
    
    for i in range(len(prova[0])):
        d = prova[0][i]
        Sim.Fill(d[0],d[1])
    Sim.SetTitle("Distribuzione campionata da GMM (k0>=" + str(k0_min)+"); tan(#theta); VR0_av")
    cluster_info.Write()
    Sim.Write()
    file.Close()
    return model


def Clustering_GM2(file_name, variables, tuple_name, want_plots=0, n_clust=0):
    
    file2 = uproot.open(file_name)
    vrs2 = file2[tuple_name]
    
    y = vrs2.arrays(variables[0])[:]
    x = vrs2.arrays(variables[1])[:]

    x, y = ak.to_numpy(x), ak.to_numpy(y)
    x, y = list(x), list(y)

    x_n, y_n = [], []
    for i in range(len(x)):
        x_n.append(x[i][0])
        y_n.append(y[i][0])

    X, cnt, c0_s = [], [], []
    for i in range(len(x_n)):
        X.append([x_n[i], y_n[i]])
        cnt.append(i)
    
    from sklearn.mixture import GaussianMixture
    model = GaussianMixture(n_components=2, covariance_type = 'spherical', random_state = 10)
    model.fit(X)
    
    yhat = model.predict(X)
    clusters = np.unique(yhat)

    cluster1_x, cluster1_y = [], [] 
    cluster2_x, cluster2_y = [], [] 

    for i in range(len(cnt)):
        if (yhat[i] == 0):
            cluster1_x.append(x_n[cnt[i]])
            cluster1_y.append(y_n[cnt[i]])
        else:
            cluster2_x.append(x_n[cnt[i]])
            cluster2_y.append(y_n[cnt[i]])
            
    if (want_plots==1):
        plt.hist2d(cluster1_x,cluster1_y, 40)
        ax = plt.gca()
        ax.set_ylim([0, max(y_n)])
        plt.title("Cluster 1")
        plt.xlabel(variables[1])
        plt.ylabel(variables[0])
        plt.show()
    
        plt.hist2d(cluster2_x, cluster2_y, 40)
        ax = plt.gca()
        ax.set_ylim([0, max(y_n)])
        plt.title("Cluster 2")
        plt.xlabel(variables[1])
        plt.ylabel(variables[0])
        plt.show()
    
    file1 = r.TFile(file_name, "UPDATE")
    cluster_info = r.TNtuple("clust_"+str(n_clust), "clu", "n_clust_"+str(n_clust))
    
    for i in range(len(x_n)):
        cluster_info.Fill(yhat[i])
  
    prova = model.sample(n_samples = 30000)
    Sim = r.TH2F("Sim", "h", 40, 0, max(x_n), 40, 0, max(y_n))
    
    for i in range(len(prova[0])):
        d = prova[0][i]
        Sim.Fill(d[0],d[1])
    Sim.SetTitle("Distribuzione campionata da GMM ")
    cluster_info.Write()
    Sim.Write()
    file1.Close()
    return model




def Calcolo_Variabili_Volume3(file_name, trk_name, k0_min = 0, k1_min = 0, k2_min = 0, k3_min = 0):
    
    track_file = r.TFile(trk_name, "READ")
    tracks = track_file.Get("tracks") 
    
    VR0, VR1, VR2, VR3 = [], [], [], []

    k0s, k1s, k2s, k3s = [], [], [], [] 
    tan = []

    for track in tracks:
        plate_0 = 31
        vr0, vr1, vr2, vr3 = 0,0,0,0
        k0, k1, k2, k3 = 0,0,0,0
        for s in track.s:
            if((s.Plate()-plate_0)%4 - 0 == 0):
                vr0 += s.Volume()
                k0+=1
            elif((s.Plate()-plate_0)%4 - 1 == 0):
                vr1 += s.Volume()
                k1+=1
            elif((s.Plate()-plate_0)%4 - 2 == 0):
                vr2 += s.Volume()
                k2+=1
            elif((s.Plate()-plate_0)%4 - 3 == 0):
                vr3+= s.Volume()
                k3+=1
        if(k0!=0):
            VR0.append(vr0/k0)
        else:
            VR0.append(0)
        if(k1!=0):
            VR1.append(vr1/k1)
        else:
            VR1.append(0)
        if(k2!=0):
            VR2.append(vr2/k2)
        else:
            VR2.append(0)
        if(k3!=0):
            VR3.append(vr3/k3)
        else:
            VR3.append(0)
        k0s.append(k0)
        k1s.append(k1)
        k2s.append(k2)
        k3s.append(k3)
        tan.append(np.sqrt(track.s[0].TX()*track.s[0].TX() + track.s[0].TY()*track.s[0].TY() ))
    
    file = r.TFile(file_name, "RECREATE")   
    t2s = r.TNtuple("Tupla", "VR_i","VR0_av:VR1_av:VR2_av:VR3_av:tan:k0:k1:k2:k3")

    for i in range(len(VR0)):
        t2s.Fill(VR0[i], VR1[i], VR2[i], VR3[i], tan[i],  k0s[i], k1s[i], k2s[i], k3s[i] )
	
    tracks.AddFriend(t2s)
    new_tracks = tracks.CopyTree("k0 >=" + str(k0_min) + " && k1>= "+ str(k1_min) + " && k2>= " + str(k2_min) + " && k3>= "  + str(k3_min))
    new_tracks.Write("V_tracks")
    tn = t2s.CopyTree("k0 >=" + str(k0_min) + " && k1>= "+ str(k1_min) + " && k2>= " + str(k2_min) + " && k3>= "  + str(k3_min) ) 	

    tn.Write("Volume_info")
 
    file.Close()
    
    return 1


def Calcolo_Variabili_Volume_FIN(file_name, trk_name, file_name2, k0_min = 0, k1_min = 0, k2_min = 0, k3_min = 0):
    
    track_file = r.TFile(file_name, "READ")
    tracks = track_file.Get(trk_name) 
    
    VR0, VR1, VR2, VR3 = [], [], [], []

    k0s, k1s, k2s, k3s = [], [], [], [] 
    tan = []

    for track in tracks:
        plate_0 = 31
        vr0, vr1, vr2, vr3 = 0,0,0,0
        k0, k1, k2, k3 = 0,0,0,0
        for s in track.s:
            if((s.Plate()-plate_0)%4 - 0 == 0):
                vr0 += s.Volume()
                k0+=1
            elif((s.Plate()-plate_0)%4 - 1 == 0):
                vr1 += s.Volume()
                k1+=1
            elif((s.Plate()-plate_0)%4 - 2 == 0):
                vr2 += s.Volume()
                k2+=1
            elif((s.Plate()-plate_0)%4 - 3 == 0):
                vr3+= s.Volume()
                k3+=1
        if(k0!=0):
            VR0.append(vr0/k0)
        else:
            VR0.append(0)
        if(k1!=0):
            VR1.append(vr1/k1)
        else:
            VR1.append(0)
        if(k2!=0):
            VR2.append(vr2/k2)
        else:
            VR2.append(0)
        if(k3!=0):
            VR3.append(vr3/k3)
        else:
            VR3.append(0)
        k0s.append(k0)
        k1s.append(k1)
        k2s.append(k2)
        k3s.append(k3)
        tan.append(np.sqrt(track.s[0].TX()*track.s[0].TX() + track.s[0].TY()*track.s[0].TY() ))
    
    file2 = r.TFile(file_name2, "RECREATE")   
    t2s = r.TNtuple("Tupla", "VR_i","VR0_av:VR1_av:VR2_av:VR3_av:tan:k0:k1:k2:k3")

    for i in range(len(VR0)):
        t2s.Fill(VR0[i], VR1[i], VR2[i], VR3[i], tan[i],  k0s[i], k1s[i], k2s[i], k3s[i] )
	
    tracks.AddFriend(t2s)
    new_tracks = tracks.CopyTree("k0 >=" + str(k0_min) + " && k1>= "+ str(k1_min) + " && k2>= " + str(k2_min) + " && k3>= "  + str(k3_min))
    new_tracks.Write("V_tracks")
    tn = t2s.CopyTree("k0 >=" + str(k0_min) + " && k1>= "+ str(k1_min) + " && k2>= " + str(k2_min) + " && k3>= "  + str(k3_min) ) 	

    tn.Write("Volume_info")
 
    file2.Close()
    track_file.Close()
    
    return 1


def make_classification_123(PCA_tup, fit_func, N_gaus=4, tr): 
    
    params = fit_func.GetParameters()
    file_pca2 = r.TFile("PCA2.root", "RECREATE")
    low_value, high_value = 0., 0.
    
    pca_2 = r.TNtuple("pca_2", "", "VR123:Z_c")
    Z_c = 0.
    cn_s = []
    for i in range(N_gaus+1):
        cn_s.append([])
        
    for track in PCA_tup:
        PCA_value = track.VR123
        
        if (PCA_value>tr):
            cn_s[-1].append(PCA_value)
            Z_c = 4
            pca_2.Fill(PCA_value, Z_c)
            continue
        
        random_number = np.random.uniform(0,1)
        probs = []
        for i in range(N_gaus):
            probs.append(g_func(PCA_value, params[0+int(3*i)], params[1+int(3*i)], params[2+int(3*i)])/fit_func.Eval(PCA_value) )
        o_ps = sorted(probs)
        indexes = []
        for p in probs:
            indexes.append(o_ps.index(p))
        
        if (random_number <= o_ps[0]):
            check = True
            pos = 0
            for j1 in range(len(indexes)):
                if (indexes[j1] == 0): 
                    pos = j1
            for j2 in range(len(indexes)):
                if (j2 != pos and indexes[j2] == 0):
                    check = False
            if (check):
                cn_s[pos].append(PCA_value)
                if (pos == 2 or pos == 3):
                    Z_c = 2
                elif (pos == 1):
                    Z_c = 3
                else:
                    Z_c = 4
                
                
        for k in range(N_gaus-1):
            
            if (k == 0):
                low_value = o_ps[0]
            elif (k == 1):
                low_value = o_ps[0] + o_ps[k]
            elif (k == 2):
                low_value = high_value
            high_value = low_value + o_ps[k+1]
        
            if (k < N_gaus-2):
                if (random_number > low_value and random_number <= high_value):
                    check = True
                    pos = 0
                    for j1 in range(len(indexes)):
                        if (indexes[j1] == k+1): 
                            pos = j1
                    for j2 in range(len(indexes)):
                        if (j2 != pos and indexes[j2] == k+1):
                            check = False
                    if (check):
                        cn_s[pos].append(PCA_value)
                        if (pos == 2 or pos == 3):
                            Z_c = 2
                        elif (pos == 1):
                            Z_c = 3
                        else:
                            Z_c = 4
            elif (k == N_gaus-2):  #funziona per N_gaus = 4
                if (random_number > low_value):
                    check = True
                    pos = 0
                    for j1 in range(len(indexes)):
                        if (indexes[j1] == k+1): 
                            pos = j1
                    for j2 in range(len(indexes)):
                        if (j2 != pos and indexes[j2] == k+1):
                            check = False
                    if (check):
                        cn_s[pos].append(PCA_value)
                        if (pos == 2 or pos == 3):
                            Z_c = 2
                        elif (pos == 1):
                            Z_c = 3
                        else:
                            Z_c = 4
        
        pca_2.Fill(PCA_value, Z_c)     
    pca_2.Write("pca_2")
    file_pca2.Close()
    return cn_s


def make_classification_012(PCA_tup, fit_func, N_gaus=4, tr): 
    
    params = fit_func.GetParameters()
    file_pca2 = r.TFile("PCA2.root", "RECREATE")
    low_value, high_value = 0., 0.
    
    pca_2 = r.TNtuple("pca_2", "", "VR123:Z_c")
    Z_c = 0.
    cn_s = []
    for i in range(N_gaus+1):
        cn_s.append([])
        
    for track in PCA_tup:
        PCA_value = track.VR012
        
        if (PCA_value>tr):
            cn_s[-1].append(PCA_value)
            Z_c = 4
            pca_2.Fill(PCA_value, Z_c)
            continue
        
        random_number = np.random.uniform(0,1)
        probs = []
        for i in range(N_gaus):
            probs.append(g_func(PCA_value, params[0+int(3*i)], params[1+int(3*i)], params[2+int(3*i)])/fit_func.Eval(PCA_value) )
        o_ps = sorted(probs)
        indexes = []
        for p in probs:
            indexes.append(o_ps.index(p))
        
        if (random_number <= o_ps[0]):
            check = True
            pos = 0
            for j1 in range(len(indexes)):
                if (indexes[j1] == 0): 
                    pos = j1
            for j2 in range(len(indexes)):
                if (j2 != pos and indexes[j2] == 0):
                    check = False
            if (check):
                cn_s[pos].append(PCA_value)
                if (pos == 2 or pos == 3):
                    Z_c = 2
                elif (pos == 1):
                    Z_c = 3
                else:
                    Z_c = 4
                
                
        for k in range(N_gaus-1):
            
            if (k == 0):
                low_value = o_ps[0]
            elif (k == 1):
                low_value = o_ps[0] + o_ps[k]
            elif (k == 2):
                low_value = high_value
            high_value = low_value + o_ps[k+1]
        
            if (k < N_gaus-2):
                if (random_number > low_value and random_number <= high_value):
                    check = True
                    pos = 0
                    for j1 in range(len(indexes)):
                        if (indexes[j1] == k+1): 
                            pos = j1
                    for j2 in range(len(indexes)):
                        if (j2 != pos and indexes[j2] == k+1):
                            check = False
                    if (check):
                        cn_s[pos].append(PCA_value)
                        if (pos == 2 or pos == 3):
                            Z_c = 2
                        elif (pos == 1):
                            Z_c = 3
                        else:
                            Z_c = 4
            elif (k == N_gaus-2):  #funziona per N_gaus = 4
                if (random_number > low_value):
                    check = True
                    pos = 0
                    for j1 in range(len(indexes)):
                        if (indexes[j1] == k+1): 
                            pos = j1
                    for j2 in range(len(indexes)):
                        if (j2 != pos and indexes[j2] == k+1):
                            check = False
                    if (check):
                        cn_s[pos].append(PCA_value)
                        if (pos == 2 or pos == 3):
                            Z_c = 2
                        elif (pos == 1):
                            Z_c = 3
                        else:
                            Z_c = 4
        
        pca_2.Fill(PCA_value, Z_c)     
    pca_2.Write("pca_2")
    file_pca2.Close()
    return cn_s








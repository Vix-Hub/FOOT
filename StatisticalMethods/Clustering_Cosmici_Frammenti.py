
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

def g_func(x, p1, p2, p3):
    return p1* np.exp(- (x - p2)**2 / (2*p3**2) )

def Calcolo_Variabili_Volume_New(file_name, trk_name):
    
    track_file = r.TFile(trk_name, "READ")
    tracks = track_file.Get("tracks") 
    
    VR0, VR1, VR2, VR3 = [], [], [], []

    k0s, k1s, k2s, k3s = [], [], [], [] 
    tan, tanfs = [], []
    
    file = r.TFile(file_name, "RECREATE") 
    
    vr0n = np.zeros(1, dtype = np.float32)
    vr1n = np.zeros(1, dtype = np.float32)
    vr2n = np.zeros(1, dtype = np.float32)
    vr3n = np.zeros(1, dtype = np.float32)
    tann = np.zeros(1, dtype = np.float32)
    tanf = np.zeros(1, dtype = np.float32)
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
    output_tree.Branch("tanf", tanf, "tanf/F")
    

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
        tanfs.append(np.sqrt(track.sf[0].TX()*track.sf[0].TX() + track.sf[0].TY()*track.sf[0].TY() ))
        
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
        tanf[0] = tanfs[i]
        
        
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


def make_classification_123(PCA_tup, fit_func, N_gaus=4, tr=-1.3): 
    
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


def make_classification_012(PCA_tup, fit_func, N_gaus=4, tr=-1.5): 
    
    params = fit_func.GetParameters()
    file_pca2 = r.TFile("PCA012.root", "RECREATE")
    low_value, high_value = 0., 0.
    
    pca_2 = r.TNtuple("pca_2", "", "VR123:Z_012")
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
                if (pos == 1):
                    Z_c = 2
                elif (pos == 2):
                    Z_c = 3
                elif (pos == 3 or pos == 0):
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
                        if (pos == 1):
                            Z_c = 2
                        elif (pos == 2):
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
                        if (pos == 1):
                            Z_c = 2
                        elif (pos == 2):
                            Z_c = 3
                        else:
                            Z_c = 4
        
        pca_2.Fill(PCA_value, Z_c)     
    pca_2.Write("pca_012")
    file_pca2.Close()
    return cn_s


def make_classification_123_v2(PCA_tup, fit_func, file_name, N_gaus=4, tr=-1.3): 
    
    params = fit_func.GetParameters()
    file_pca2 = r.TFile(file_name, "RECREATE")
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
                if (pos == 1 or pos == 2):
                    Z_c = 2
                elif (pos == 0):
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
                        if (pos == 1 or pos == 2):
                            Z_c = 2
                        elif (pos == 0):
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
                        if (pos == 1 or pos == 2):
                            Z_c = 2
                        elif (pos == 0):
                            Z_c = 3
                        else:
                            Z_c = 4
        
        pca_2.Fill(PCA_value, Z_c)     
    pca_2.Write("pca_2")
    file_pca2.Close()
    return cn_s


def make_classification_123_X(PCA_tup, fit_func, file_name, tupname, N_gaus=4, tr=-1.3, Z0=1, Z1=2, Z2=3, Z3=4): 
    
    params = fit_func.GetParameters()
    file_pca2 = r.TFile(file_name, "RECREATE")
    low_value, high_value = 0., 0.
    
    pca_2 = r.TNtuple(tupname, "", "VR123:Z_c:VR123b")
    Z_c = 0.
    cn_s = []
    for i in range(N_gaus+1):
        cn_s.append([])
        
    for track in PCA_tup:
        PCA_value = track.VR123
        PCA_value2 = track.VR123b
        
        if (PCA_value>tr):
            cn_s[-1].append(PCA_value)
            Z_c = Z3
            pca_2.Fill(PCA_value, Z_c, PCA_value2)
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
                if (pos == 1):
                    Z_c = Z1
                elif (pos == 2):
                    Z_c = Z2
                elif (pos == 0):
                    Z_c = Z0
                else:
                    Z_c = Z3
                
                
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
                        if (pos == 1):
                    	    Z_c = Z1
                        elif (pos == 2):
                            Z_c = Z2
                        elif (pos == 0):
                            Z_c = Z0
                        else:
                            Z_c = Z3
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
                        if (pos == 1):
                    	    Z_c = Z1
                        elif (pos == 2):
                            Z_c = Z2
                        elif (pos == 0):
                    	    Z_c = Z0
                        else:
                    	    Z_c = Z3
        
        pca_2.Fill(PCA_value, Z_c)     
    pca_2.Write(tupname)
    file_pca2.Close()
    return cn_s


def make_classification_012_X(PCA_tup, fit_func, file_name, tupname, N_gaus=4, tr=-1.5, Z0=2, Z1=3, Z2=4, Ztr=10): 
    
    params = fit_func.GetParameters()
    file_pca2 = r.TFile(file_name, "RECREATE")
    low_value, high_value = 0., 0.
    
    pca_2 = r.TNtuple(tupname, "", "VR012:Z_012")
    Z_c = 0.
    cn_s = []
    for i in range(N_gaus+1):
        cn_s.append([])
        
    for track in PCA_tup:
        PCA_value = track.VR012
        
        if (PCA_value>tr):
            cn_s[-1].append(PCA_value)
            Z_c = Ztr
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
                if (pos == 1):
                    Z_c = Z0
                elif (pos == 2):
                    Z_c = Z1
                elif (pos == 0 or pos == 3):
                    Z_c = Z2
                
                
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
                        if (pos == 1):
                    	    Z_c = Z0
                        elif (pos == 2):
                            Z_c = Z1
                        elif (pos == 0 or pos == 3):
                            Z_c = Z2
                
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
                        if (pos == 1):
                            Z_c = Z0
                        elif (pos == 2):
                            Z_c = Z1
                        elif (pos == 0 or pos == 3):
                            Z_c = Z2
    
        
        pca_2.Fill(PCA_value, Z_c)     
    pca_2.Write(tupname)
    file_pca2.Close()
    return cn_s


def Calcolo_Variabili_Volume_MC(file_name, trk_name):
    
    track_file = r.TFile(trk_name, "READ")
    tracks = track_file.Get("tracks") 

    k0s, k1s, k2s, k3s = [], [], [], [] 
    tan, tanfs = [], []
    
    file = r.TFile(file_name, "RECREATE") 
    
    tann = np.zeros(1, dtype = np.float32)
    tanf = np.zeros(1, dtype = np.float32)
    k0n = np.zeros(1, dtype = np.intc)
    k1n = np.zeros(1, dtype = np.intc)
    k2n = np.zeros(1, dtype = np.intc)
    k3n = np.zeros(1, dtype = np.intc)
    
    output_tree = tracks.CloneTree(0)
    output_tree.Branch("k0", k0n, "k0/I")
    output_tree.Branch("k1", k1n, "k1/I")
    output_tree.Branch("k2", k2n, "k2/I")
    output_tree.Branch("k3",k3n, "k3/I")
    output_tree.Branch("tan", tann, "tan/F")
    output_tree.Branch("tanf", tanf, "tanf/F")
    

    for track in tracks:
        plate_0 = 31
        k0, k1, k2, k3 = 0,0,0,0
        for s in track.s:
            if((s.Plate()-plate_0)%4 - 0 == 0):
                k0+=1
            elif((s.Plate()-plate_0)%4 - 1 == 0):
                k1+=1
            elif((s.Plate()-plate_0)%4 - 2 == 0):
                k2+=1
            elif((s.Plate()-plate_0)%4 - 3 == 0):
                k3+=1
        k0s.append(k0)
        k1s.append(k1)
        k2s.append(k2)
        k3s.append(k3)
        tan.append(np.sqrt(track.s[0].TX()*track.s[0].TX() + track.s[0].TY()*track.s[0].TY() ))
        tanfs.append(np.sqrt(track.sf[0].TX()*track.sf[0].TX() + track.sf[0].TY()*track.sf[0].TY() ))
        
    for i in range(len(k0s)):
        tracks.GetEntry(i)

        k0n[0] = k0s[i]
        k1n[0] = k1s[i]
        k2n[0] = k2s[i]
        k3n[0] = k3s[i]
        
        tann[0] = tan[i]
        tanf[0] = tanfs[i]
        
        
        output_tree.Fill()
        
    output_tree.Write("tracks_n")
    file.Close()
    track_file.Close()
    
    return 1


def make_classification_01_X(PCA_tup, fit_func, file_name, tupname, N_gaus=2, Z0=1, Z1=2): 
    
    params = fit_func.GetParameters()
    file_pca2 = r.TFile(file_name, "RECREATE")
    low_value, high_value = 0., 0.
    
    pca_2 = r.TNtuple(tupname, "", "VP01:Z_c")
    Z_c = 0.
    cn_s = []
    for i in range(N_gaus+1):
        cn_s.append([])
        
    for track in PCA_tup:
        PCA_value = track.VP01
        
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
                if (pos == 1):
                    Z_c = Z1
                elif (pos == 0):
                    Z_c = Z0
                
                
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
                        if (pos == 1):
                    	    Z_c = Z1
                        elif (pos == 0):
                            Z_c = Z0
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
                        if (pos == 1):
                    	    Z_c = Z1
                        elif (pos == 0):
                    	    Z_c = Z0
        
        pca_2.Fill(PCA_value, Z_c)     
    pca_2.Write(tupname)
    file_pca2.Close()
    return cn_s


def make_classification_01_X2(PCA_tup, fit_func, file_name, tupname, N_gaus=2, Z0=1, Z1=2, tr=2., Ztr=2): 
    
    params = fit_func.GetParameters()
    file_pca2 = r.TFile(file_name, "RECREATE")
    low_value, high_value = 0., 0.
    
    pca_2 = r.TNtuple(tupname, "", "VP01:Z_c")
    Z_c = 0.
    cn_s = []
    for i in range(N_gaus+1):
        cn_s.append([])
        
    for track in PCA_tup:
        PCA_value = track.VP01
        if (PCA_value >= tr):
            Z_c = Ztr
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
                if (pos == 1):
                    Z_c = Z1
                elif (pos == 0):
                    Z_c = Z0
                
                
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
                        if (pos == 1):
                    	    Z_c = Z1
                        elif (pos == 0):
                            Z_c = Z0
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
                        if (pos == 1):
                    	    Z_c = Z1
                        elif (pos == 0):
                    	    Z_c = Z0
        
        pca_2.Fill(PCA_value, Z_c)     
    pca_2.Write(tupname)
    file_pca2.Close()
    return cn_s



def make_classification_013_X(PCA_tup, fit_func, file_name, tupname, N_gaus=4, tr=-1.5, Z0=2, Z1=3, Z2=4): 
    
    params = fit_func.GetParameters()
    file_pca2 = r.TFile(file_name, "RECREATE")
    low_value, high_value = 0., 0.
    
    pca_2 = r.TNtuple(tupname, "", "VR012:Z_013")
    Z_c = 0.
    cn_s = []
    for i in range(N_gaus+1):
        cn_s.append([])
        
    for track in PCA_tup:
        PCA_value = track.VR013
        
        if (PCA_value>tr):
            cn_s[-1].append(PCA_value)
            Z_c = 11
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
                if (pos == 0):
                    Z_c = Z0
                elif (pos == 1):
                    Z_c = Z1
                elif (pos == 2):
                    Z_c = Z2
                
                
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
                        if (pos == 0):
                    	    Z_c = Z0
                        elif (pos == 1):
                            Z_c = Z1
                        elif (pos == 2):
                            Z_c = Z2
                
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
                        if (pos == 0):
                            Z_c = Z0
                        elif (pos == 1):
                            Z_c = Z1
                        elif (pos == 2):
                            Z_c = Z2
    
        
        pca_2.Fill(PCA_value, Z_c)     
    pca_2.Write(tupname)
    file_pca2.Close()
    return cn_s


def make_classification_123_ALL(PCA_tup, fit_func, file_name, tupname, N_gaus=6, Z0=2, Z1=3, Z2=4, Z3=5, Z4=6): 
    
    params = fit_func.GetParameters()
    file_pca2 = r.TFile(file_name, "RECREATE")
    low_value, high_value = 0., 0.
    
    pca_2 = r.TNtuple(tupname, "", "VR123:Z_c:VR123b")
    Z_c = 0.
    cn_s = []
    for i in range(N_gaus+1):
        cn_s.append([])
        
    for track in PCA_tup:
        PCA_value = track.VR123
        PCA_value2 = track.VR123b  
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
                if (pos == 0):
                    Z_c = Z1
                elif (pos == 2):
                    Z_c = Z2
                elif (pos == 0):
                    Z_c = Z0
                elif(pos == 3):
                    Z_c = Z3
                else:
                    Z_c = Z4
                
                
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
                        if (pos == 1):
                            Z_c = Z1
                        elif (pos == 2):
                            Z_c = Z2
                        elif (pos == 0):
                            Z_c = Z0
                        elif(pos == 3):
                            Z_c = Z3
                        else:
                            Z_c = Z4
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
                        if (pos == 1):
                            Z_c = Z1
                        elif (pos == 2):
                            Z_c = Z2
                        elif (pos == 0):
                            Z_c = Z0
                        elif(pos == 3):
                            Z_c = Z3
                        else:
                            Z_c = Z4
        
        pca_2.Fill(PCA_value, Z_c)     
    pca_2.Write(tupname)
    file_pca2.Close()
    return cn_s

def make_classification_013_X(PCA_tup, fit_func, file_name, tupname, N_gaus=4, tr=-1.5, Z0=2, Z1=3, Z2=4): 
    
    params = fit_func.GetParameters()
    file_pca2 = r.TFile(file_name, "RECREATE")
    low_value, high_value = 0., 0.
    
    pca_2 = r.TNtuple(tupname, "", "VR012:Z_013")
    Z_c = 0.
    cn_s = []
    for i in range(N_gaus+1):
        cn_s.append([])
        
    for track in PCA_tup:
        PCA_value = track.VR013
        
        if (PCA_value>tr):
            cn_s[-1].append(PCA_value)
            Z_c = 11
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
                if (pos == 0):
                    Z_c = Z0
                elif (pos == 1):
                    Z_c = Z1
                elif (pos == 2):
                    Z_c = Z2
                
                
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
                        if(pos==0):
                            Z_c = Z0
                        elif(pos==1):
                            Z_c = Z1
                        elif(pos==2):
                            Z_c = Z2
                
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
                        if (pos == 0):
                            Z_c = Z0
                        elif (pos == 1):
                            Z_c = Z1
                        elif (pos == 2):
                            Z_c = Z2
    
        
        pca_2.Fill(PCA_value, Z_c)     
    pca_2.Write(tupname)
    file_pca2.Close()
    return cn_s

def do_PCA_VP01(tracks_2, VERBOSE=0, a=2400, b=0.95, brick_id="GSI3", file_name01='01.root', tup_name01='tup01'):
    file_pca = r.TFile("PCA_01.root", "RECREATE")

    frag_cut = "VR0_av >= " + str(a) + "*(1 + TMath::Exp(" + str(b) + " *s[0].Theta()*s[0].Theta()))"

    cond = " && k1>0 && k2<2 && k3<2 && VR1_av<7500"
    new_cut2 = frag_cut + cond
    campione_pca = tracks_2.CopyTree(new_cut2)
    campione_pca.Write("pca")
    file_pca.Close()

    principal = r.TPrincipal(2, "ND")
    file_pca = r.TFile("PCA_01.root", "READ")
    info_pca = file_pca.Get("pca")

    for track in info_pca:
        vr0, vr1 = track.VR0_av, track.VR1_av
        vrs = np.zeros(2)
        vrs[0] = vr0
        vrs[1] = vr1
        principal.AddRow(vrs)
        
    principal.MakePrincipals()
    principal.MakeCode()
    r.gInterpreter.ProcessLine('.L pca.C+')
    r.gSystem.Load("pca_C.so")
    vr123s = []
    for track in info_pca:
        vr2, vr3 = track.VR0_av, track.VR1_av
        vrs = np.zeros(2)
        vrs[0] = vr2
        vrs[1] = vr3
        princ = np.zeros(2)
        principal.X2P(vrs, princ)
        vr123s.append(princ)

    vr123 = []
    for i in vr123s:
        vr123.append(i[0])
    file_pca.Close()

    file_pca_2 = r.TFile("PCA3.root", "RECREATE")
    pca_1 = r.TNtuple("pca_2", "", "VP01")

    for i in range(len(vr123)):
        pca_1.Fill(vr123[i])
    pca_1.Write("pca_01")

    kn = r.TCanvas()
    pca_1.Draw("VP01>>01(80, -6.5, 5.5)")

    hi = r.gDirectory.Get("01")
    hi.SetTitle("VP_{01} [Lower Population]; VP_{01}; Entries")
    hi.Draw("")

    g2 = r.TF1("g2", "gaus(0) + gaus(3)", -6.5, 6.)
    g2.SetParameters(400, 0.47, .6, 30, -2, 1.)
    g2.SetParLimits(4, -1.9, -1.8)

    hi.Fit("g2", "L","", -6.5, 4.)
    hi.Draw("PE")

    params = g2.GetParameters()

    comp1 = r.TF1("comp1", "[0]*TMath::Gaus(x, [1], [2])", -10.5, 4.)
    comp1.SetParameters(params[0], params[1], params[2])
    comp1.SetLineColor(4)
    comp1.Draw("SAME")

    comp1_2 = r.TF1("comp1_2", "[0]*TMath::Gaus(x, [1], [2])", -10.5, 4.)
    comp1_2.SetParameters(params[3], params[4], params[5])
    comp1_2.SetLineColor(2)
    comp1_2.SetLineStyle(2)
    comp1_2.Draw("SAME")
    legend = r.TLegend(0.6,0.65,0.88,0.85)
    legend.SetTextFont(0)
    legend.SetTextSize(0.04)
    legend.AddEntry("chi2 / NDF", "#chi^{2} / NDF = " + str(round(g2.GetChisquare(), 2)) + " / " + str(g2.GetNDF()), "" )
    legend.AddEntry("prob", "Prob = " + str(round(g2.GetProb(), 5)), "")
    legend.Draw("SAME")

    t1 = r.TText(0.8, 20, brick_id)
    t1.SetTextColor(1)
    t1.SetTextSize(20)
    t1.Draw("SAME")

    if (VERBOSE==1): 
        file_out.cd()
        kn.Write("VP01")

    prova = make_classification_01_X(pca_1, g2, '01.root', '01_c')
    prova = make_classification_01_X(pca_1, g2, '01.root', '01_c', 2, 1)


def do_PCA_VP123(tracks_2, VERBOSE=0, a=2400, b=0.95, brick_id="GSI3", file_name123="123.root", tup_name123 ="tup123"):
    file_pca = r.TFile("PCA_123.root", "RECREATE")

    frag_cut = "VR0_av >= " + str(a) + "*(1 + TMath::Exp(" + str(b) + " *s[0].Theta()*s[0].Theta()))"

    cond = " && k1>0 && k2>1 && k3>1"
    new_cut2 = frag_cut + cond
    campione_pca = tracks_2.CopyTree(new_cut2)
    campione_pca.Write("pca")
    file_pca.Close()

    principal = r.TPrincipal(3, "ND")
    file_pca = r.TFile("PCA_123.root", "READ")
    info_pca = file_pca.Get("pca")

    for track in info_pca:
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
    for track in info_pca:
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
        
    file_pca.Close()
    file_pca_2 = r.TFile("PCA2.root", "RECREATE")
    pca_1 = r.TNtuple("pca_1", "", "VR123:VR123b")

    for i in range(len(vr123)):
        pca_1.Fill(vr123[i], vr123b[i])

    file_pca = r.TFile("PCA_123.root", "READ")
    info_pca = file_pca.Get("pca")
    info_pca.AddFriend(pca_1)
    c = r.TCanvas()
    info_pca.Draw("VR123>>histo123(100)", "", "colz")

    h123 = r.gDirectory.Get("histo123")
    h123.SetTitle("VP_{123} (Frag_Cut, k_{1}>0, k_{2}, k_{3}>1);VP_{123};Entries")

    r.gStyle.SetOptStat(0)
    g1 = r.TF1("g1", "gaus(0) + gaus(3) + gaus(6) + gaus(9) + gaus(12) + gaus(15)", -12.5, 2.)
    g1.SetParameters(20, -9., .5, 50, -6.3, 0.5, 40, -3., 0.5, 70, -1.5)
    g1.SetParameter(11, 0.2)
    g1.SetParameter(12, 80)
    g1.SetParameter(13, -0.35)
    g1.SetParameter(14, 0.5)
    g1.SetParameter(15, 2500)
    g1.SetParameter(16, .5)
    g1.SetParameter(17, 0.8)
    #g1.SetParameters(20, -9., .5, 50, -6.3, 0.5, 40, -3., 0.5, 70, -1.5)
    g1.SetLineWidth(2)

    h123.Fit("g1", "L","", -11.5, 2.)
    h123.Draw("PE")

    params = g1.GetParameters()

    comp1 = r.TF1("comp1", "[0]*TMath::Gaus(x, [1], [2])", -10.5, 2.) 
    comp1.SetParameters(params[0], params[1], params[2])
    comp1.SetLineColor(4)  #r.kCyan+2
    comp1.SetLineStyle(2)
    comp1.Draw("SAME")

    comp2 = r.TF1("comp2", "[0]*TMath::Gaus(x, [1], [2])", -10.5, 2.) #Z=5
    comp2.SetParameters(params[3], params[4], params[5])
    comp2.SetLineStyle(2)
    comp2.SetLineColor(95) #r.kGreen+2
    comp2.Draw("SAME")

    comp3 = r.TF1("comp3", "[0]*TMath::Gaus(x, [1], [2])", -10.5, 2.) #shoulder
    comp3.SetParameters(params[6], params[7], params[8])
    comp3.SetLineColor(8)  #r.kOrange
    comp3.SetLineStyle(2)
    comp3.Draw("SAME")

    comp4 = r.TF1("comp4", "[0]*TMath::Gaus(x, [1], [2])", -10.5, 2.) #Z=2
    comp4.SetParameters(params[9], params[10], params[11])
    comp4.SetLineColor(r.kPink+3)
    comp4.SetLineStyle(2)
    comp4.Draw("SAME")

    comp5 = r.TF1("comp5", "[0]*TMath::Gaus(x, [1], [2])", -10.5, 2.) #Z=3
    comp5.SetParameters(params[12], params[13], params[14])
    comp5.SetLineColor(r.kMagenta)
    comp5.SetLineWidth(2)
    comp5.Draw("SAME")

    comp6 = r.TF1("comp6", "[0]*TMath::Gaus(x, [1], [2])", -10.5, 2.) #Z=4
    comp6.SetParameters(params[15], params[16], params[17])
    comp6.SetLineColor(r.kCyan+2)
    comp6.SetLineStyle(2)
    comp6.Draw("SAME")

    comps = [comp1, comp2, comp3, comp4, comp5, comp6]
    for comp in comps:
        comp.SetLineStyle(9)

    legend = r.TLegend(0.4,0.65,0.68,0.85)
    legend.SetTextFont(0)
    legend.SetTextSize(0.04)
    legend.AddEntry(comp1, "Z = 2")
    legend.AddEntry(comp2, "Z = 3")
    legend.AddEntry(comp3, "Z = 4")
    legend.AddEntry(comp4, "Z = 5")
    legend.AddEntry(comp5, "Z >=6")
    legend.AddEntry(comp6, "Z >=6")
    legend.Draw("SAME")

    legend2 = r.TLegend(0.15,0.7,0.4,0.85)
    legend2.SetTextFont(0)
    legend2.SetTextSize(0.04)
    legend2.AddEntry("Entries", str("Entries = "+str(h123.GetEntries())), "" )
    legend2.AddEntry("chi2 / NDF", "#chi^{2} / NDF = " + str(round(g1.GetChisquare(), 2)) + " / " + str(g1.GetNDF()), "" )
    legend2.AddEntry("prob", "Prob = " + str(round(g1.GetProb(), 4)), "")
    legend2.Draw("SAME")
    g1.Draw("SAME")
    c.SetLogy(1)

    if (VERBOSE ==1):
        file_out.cd()	
        c.Write("VP123")
    prova = make_classification_123_ALL(info_pca, g1, file_name123, tup_name123, 6, 2, 3, 4, 5, 6)


def do_PCA_VP012(tracks_2, VERBOSE=0, a=2400, b=0.95, brick_id="GSI3", file_name012="012.root", tup_name012 ="tup012"):
    file_pca = r.TFile("PCA_012.root", "RECREATE")
    
    frag_cut = "VR0_av >= " + str(a) + "*(1 + TMath::Exp(" + str(b) + " *s[0].Theta()*s[0].Theta()))"

    cond = " && k1>0 && k2>1"
    new_cut2 = frag_cut + cond

    campione_pca = tracks_2.CopyTree(new_cut2)
    campione_pca.Write("pca")

    file_pca.Close()

    principal = r.TPrincipal(3, "ND")

    file_pca = r.TFile("PCA_012.root", "READ")
    info_pca = file_pca.Get("pca")

    for track in info_pca:
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
    vr123s = []
    for track in info_pca:
        vr2, vr3, vr0 = track.VR0_av, track.VR1_av, track.VR2_av
        vrs = np.zeros(3)
        vrs[0] = vr2
        vrs[1] = vr3
        vrs[2] = vr0
        princ = np.zeros(3)
        principal.X2P(vrs, princ)
        vr123s.append(princ)

    vr123 = []
    for i in vr123s:
        vr123.append(i[0])
        
    file_pca.Close()

    file_pca_2 = r.TFile("PCA3.root", "RECREATE")

    pca_1 = r.TNtuple("pca_2", "", "VR012")

    for i in range(len(vr123)):
        pca_1.Fill(vr123[i])
    pca_1.Write("pca_012")

    kn = r.TCanvas()
    pca_1.Draw("VR012>>012(100, -10, 5.5)")

    hi = r.gDirectory.Get("012")
    hi.SetTitle("VP_{012} (k_{1}>0 & k_{2}>1); VP_{012}; Entries")
    hi.Draw("")

    g2 = r.TF1("g2", "gaus(0) + gaus(3) + gaus(6)", -10.5, 6.)
    g2.SetParameters(120, -5., .5, 75, -3.8, 0.5, 50, -1.6, 0.3)
    hi.Fit("g2", "L","", -10., -1.3)
    hi.Draw("PE")


    t1 = r.TText(-9., 1000, brick_id)
    t1.SetTextColor(1)
    t1.SetTextSize(20)
    t1.Draw("SAME")

    legend = r.TLegend(0.15,0.7,0.4,0.85)
    legend.SetTextFont(0)
    legend.SetTextSize(0.04)
    legend.AddEntry("Entries", "Entries = " + str(hi.GetEntries()), "" )
    legend.AddEntry("chi2 / NDF", "#chi^{2} / NDF = " + str(round(g2.GetChisquare(), 2)) + " / " + str(g2.GetNDF()), "" )
    legend.AddEntry("prob", "Prob = " + str(round(g2.GetProb(), 4)), "")
    legend.Draw("SAME")

    kn.SetLogy(1)

    if (VERBOSE==1):
        file_out.cd()	
        kn.Write("VP012")

    prova = make_classification_012(pca_1, g2, 3, -1.5)
    prova2 = make_classification_012_X(pca_1, g2, file_name012, tup_name012, 3, -1.5, 2, 3, 4)

def do_PCA_VP013(tracks_2, VERBOSE=0, a=2400, b=0.95, brick_id="GSI3", file_name013="013.root", tup_name013 ="tup013"):
    file_pca = r.TFile("PCA_013.root", "RECREATE")

    frag_cut = "VR0_av >= " + str(a) + "*(1 + TMath::Exp(" + str(b) + " *s[0].Theta()*s[0].Theta()))"

    cond = " && k1>0 && k3>1"
    new_cut2 = frag_cut + cond

    campione_pca = tracks_2.CopyTree(new_cut2)
    campione_pca.Write("pca")

    file_pca.Close()

    principal = r.TPrincipal(3, "ND")

    file_pca = r.TFile("PCA_013.root", "READ")
    info_pca = file_pca.Get("pca")

    for track in info_pca:
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
    vr123s = []
    for track in info_pca:
        vr2, vr3, vr0 = track.VR0_av, track.VR1_av, track.VR3_av
        vrs = np.zeros(3)
        vrs[0] = vr2
        vrs[1] = vr3
        vrs[2] = vr0
        princ = np.zeros(3)
        principal.X2P(vrs, princ)
        vr123s.append(princ)

    vr123 = []
    for i in vr123s:
        vr123.append(i[0])
        
    file_pca.Close()

    file_pca_2 = r.TFile("PCA3.root", "RECREATE")

    pca_1 = r.TNtuple("pca_2", "", "VR013")

    for i in range(len(vr123)):
        pca_1.Fill(vr123[i])
    pca_1.Write("pca_013")

    kn = r.TCanvas()
    pca_1.Draw("VR013>>013(100, -10, 5.5)")

    hi = r.gDirectory.Get("013")
    hi.SetTitle("VP_{013} (k_{1}>0 & k_{3}>1); VP_{013}; Entries")
    hi.Draw("")

    g2 = r.TF1("g2", "gaus(0) + gaus(3) + gaus(6)", -10.5, 6.)
    g2.SetParameters(50, -6.5, .5, 50, -4.5, 0.5, 30, -2.5, 0.5)
    hi.Fit("g2", "L","", -10., -2.)
    hi.Draw("PE")


    t1 = r.TText(-9., 1000, brick_id)
    t1.SetTextColor(1)
    t1.SetTextSize(20)
    t1.Draw("SAME")


    legend = r.TLegend(0.15,0.7,0.4,0.85)
    legend.SetTextFont(0)
    legend.SetTextSize(0.04)
    legend.AddEntry("Entries", "Entries = "+str(hi.GetEntries()), "" )
    legend.AddEntry("chi2 / NDF", "#chi^{2} / NDF = " + str(round(g2.GetChisquare(), 2)) + " / " + str(g2.GetNDF()), "" )
    legend.AddEntry("prob", "Prob = " + str(round(g2.GetProb(), 4)), "")
    legend.Draw("SAME")

    kn.SetLogy(1)

    if (VERBOSE==1):
        file_out.cd()	
        kn.Write("VP013")

    prova2 = make_classification_013_X(pca_1, g2, file_name013, tup_name013, 3, -1.5, 2, 3, 4)
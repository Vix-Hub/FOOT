import ROOT as r 
import fedrarootlogon
import numpy as np 
import random

# script to estimate the combinatorial bkg already present in MC true
target = "C2H4"
inFile = r.TFile("output_Oxy200MeV_" + str(target) + "_clust.root", "READ")
alpha_tup = inFile.Get("alpha_tup")

outFile = r.TFile("output_Oxy200MeV_" + str(target) + "_clust_loc.root", "RECREATE")
dtup = r.TNtuple("dtup", "", "dtheta:dtheta_rand")

MIN_LAST_PLATE = 31 + 4*2
#MIN_LAST_PLATE = 0

# find unique events
events = []
for entry in alpha_tup:
    events.append(entry.iev)

events = list(np.unique(events))
couples = [] #keep angles for each event
for i in range(len(events)):
    couples.append([])

print(events)

for j, event in enumerate(events):
    for entry in alpha_tup:
        if (entry.iev == event):
            couples[j].append([entry.tx, entry.ty, entry.lastlayer+1])

for ievent in range(len(events)):
    current_couple = couples[ievent]
    for i in range(len(current_couple)-1):
        if (current_couple[i][2]<= MIN_LAST_PLATE):
            continue
        tx1, ty1 = current_couple[i][0], current_couple[i][1]
        for j in range(i+1, len(current_couple)):
            if (current_couple[j][2]<= MIN_LAST_PLATE):
                continue
            tx2, ty2 = current_couple[j][0], current_couple[j][1]
            dtheta = np.sqrt((tx1-tx2)**2 + (ty1-ty2)**2)

            if (dtheta<0.025):
                print(ievent, dtheta)

            random_event = i
            while random_event == i:
                random_event = random.randint(0, len(events)-1)

            itrk = random.randint(0, len(couples[random_event])-1)
            tx3, ty3 = couples[random_event][itrk][0], couples[random_event][itrk][1]
            dtheta_rand = np.sqrt((tx1-tx3)**2 + (ty1-ty3)**2)

            dtup.Fill(dtheta, dtheta_rand)

    if (ievent%1000==0):
        print(" Completed " + str(100.*ievent/len(events)) + " %")

dtup.Draw("dtheta>>h(50, 0, 0.5)", "")
h = r.gDirectory.Get("h")
h.SetLineWidth(2)
h.SetLineColor(r.kRed)

dtup.Draw("dtheta_rand>>h2(50, 0, 0.5)")
h2 = r.gDirectory.Get("h2")
h2.SetLineWidth(2)
h2.SetLineColor(r.kBlue)

hs = r.THStack("hs", "")
hs.Add(h)
hs.Add(h2)

hs.SetTitle(" Estimated Bkg in MC True; #Theta_{#alpha#alpha};Entries")

leg = r.TLegend()
leg.AddEntry(h, " True #alpha Couples ")
leg.AddEntry(h2, " Random #alpha Couples ")

outFile.cd()
dtup.Write()
hs.Write()
leg.Write("leg")
outFile.Close()
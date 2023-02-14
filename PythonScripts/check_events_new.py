import ROOT as r 
import fedrarootlogon
import numpy as np 
import time 

def CalculateImpactParameter(slx, sly, slz, slftx, slfty, s0x, s0y, s0z):
    dxp = slx - (slz-s0z)*slftx -s0x
    dyp = sly - (slz-s0z)*slfty - s0y
    b = np.sqrt((dxp)*(dxp)+(dyp)*(dyp))
    return b

def FindClosestMatch(event_couples, s0x, s0y, s0z, s0tx, s0ty, s0plate, couple, EVERBOSE):

    #print("Using Track With coordinates " + str([s0x, s0y, s0z]))

    bs, inv_bs = [1000], [1000]
    gaps = [1000]
    used_couples = [(0, 0)]
    for j, ev_couple in enumerate(event_couples):
        if (ev_couple[0:2] == couple[0:2]): #if it is the same track
            continue
        if (ev_couple[10] >= s0plate or abs(s0plate-ev_couple[10])>5):  #last plate condition 
            continue
        slx, sly, slz = ev_couple[5], ev_couple[6], ev_couple[7]
        slftx, slfty = ev_couple[8], ev_couple[9]
        slplate = ev_couple[-1]

        if (j == 0):
            bs[j] = CalculateImpactParameter(slx, sly, slz, slftx, slfty, s0x, s0y, s0z)
            inv_bs[j] = CalculateImpactParameter(s0x, s0y, s0z, s0tx, s0ty, slx, sly, slz)
            gaps[j] = abs(s0plate-slplate)
            used_couples[j] = (ev_couple[0], ev_couple[1])
        else:
            bs.append(CalculateImpactParameter(slx, sly, slz, slftx, slfty, s0x, s0y, s0z))
            inv_bs.append(CalculateImpactParameter(s0x, s0y, s0z, s0tx, s0ty, slx, sly, slz))
            used_couples.append( (ev_couple[0], ev_couple[1]) )
            gaps.append(abs(s0plate-slplate))
        if (EVERBOSE == 100):
            print(" Calculated b between tracks " + str(couple[:2]) + " and " + str(ev_couple[:2]) + " = " + str(CalculateImpactParameter(slx, sly, slz, slftx, slfty, s0x, s0y, s0z)))
            print( ", inverse b = " + str(inv_bs[-1]) + "\n")
    
    idx = bs.index(min(bs))
    correct_couple = used_couples[idx]
    gap = gaps[idx]

    if (EVERBOSE == 100):
        print(" ")
        print(" IDX " + str(idx))
        print(" ")
        
        print(" event couples ")
        for el in event_couples:
            print(el)
        print(" ")
        #print(" bs ")
        #for el in bs:
        #    print(el)
        #print(" ")
        

    if (len(bs)==len(event_couples)):
        return event_couples[idx], min(bs), gap, inv_bs[idx]
    else:
        if ( (len(bs)==1 and bs[0]==1000) or (min(bs)==1000)):
            return couple, min(bs), gap, inv_bs[idx]
        else:
            for ev in event_couples:
                if ((ev[0], ev[1]) == correct_couple):
                    return ev, min(bs), gap, inv_bs[idx]


InFile = r.TFile("events.root", "READ")
evt_tree = InFile.Get("events") #contains S2 track info and MC event info 


MC_ID = 2
DEBUG_MC_EVENT = 199
DEBUG_S0_ID, DEBUG_S0_PLATE = 206150, 31
EVERBOSE = -99

TrackFileName = "b00000" + str(MC_ID) + ".0.1.2.trk.root" #S2 track file name
VertexFileName = "vertices_improved_fast_3_new_temp_2.root" #Vertex File name

TrackFile = r.TFile(TrackFileName, "READ")
tracks = TrackFile.Get("tracks")

VertexFile = r.TFile(VertexFileName, "READ")
vrec = VertexFile.Get("EdbVertexRec") #Vertex Object
n_vertices = vrec.eVTX.GetEntries()

zero_vertices, flag1_cases, flag0_cases = 0, 0, 0
tracks.BuildIndex("s[0].Plate()", "s[0].ID()")

vertex_couples = []
MC_events, track_couples = [], [] 

outName = "check_events_new.root"
if (EVERBOSE==100):
    outName = "temp_check.root"
OutFile = r.TFile(outName, "RECREATE")
out_tup = r.TNtuple("tup", "Events Info", "ev_ID:ev_flag:last_b:max_b:s0X:s0Y:s0_plate:s0_id:s0Z:s0TX:s0TY:maxbgap:b_inv")

#save trid of vertex tracks
t0 = time.time()
for i in range(n_vertices):
    vertex = vrec.eVTX.At(i)
    for j in range(vertex.N()):
        track = vertex.GetTrack(j)
        vertex_couples.append((int(track.GetSegmentFirst().Plate()), int(track.GetSegmentFirst().ID())))

#save first seg info and mcevent
for track in tracks:
    MC_events.append(track.s[0].MCEvt())
    track_couples.append((track.s[0].Plate(), track.s[0].ID(), track.s[0].X(), track.s[0].Y(), track.s[0].Z(), track.s[track.nseg-1].X(), track.s[track.nseg-1].Y(), track.s[track.nseg-1].Z(), track.s[track.nseg-1].TX(), track.s[track.nseg-1].TY(), track.s[track.nseg-1].Plate(), track.s[0].TX(), track.s[0].TY()))
t1 = time.time()    
print("Saved info in " + str(t1-t0) + " s")

c = 0

for entry in evt_tree:
    c += 1 
    flag = 0
    closest_bs, closest_gaps, closest_inv_bs = [], [], []
    last_b, max_b, max_b_gap, inv_b = 0, 0, 0, 0
    s0plate, s0id = entry.s0plate, entry.s0id
    if (s0plate==DEBUG_S0_PLATE and s0id==DEBUG_S0_ID and EVERBOSE==100):
        print(" Debug track is in evt_tree ")
    s0x, s0y, s0z, s0tx, s0ty = entry.s0X, entry.s0Y, entry.s0Z, entry.s0TX, entry.s0TY
    start_s0x, start_s0y, start_s0z, start_s0tx, start_s0ty = entry.s0X, entry.s0Y, entry.s0Z, entry.s0TX, entry.s0TY
    start_s0plate, start_s0id = s0plate, s0id
    couple = (s0plate, s0id, s0x, s0y, s0z, -99, -99, -99, -99, -99, -99, s0tx, s0ty)
    n = entry.N
    mcevt = entry.MCEvt
    if (mcevt!=DEBUG_MC_EVENT and EVERBOSE==100):
        continue
    if (EVERBOSE==10):
        print(mcevt)
    vertices = []
    for j in range(n):
        if (entry.vIDs[j] != -99):
            vertices.append(entry.vIDs[j])
        elif(entry.vIDs[j]==-99 and n==1):
            zero_vertices += 1 #count how many times the MC evt is not linked to any vertex
        

    if (len(vertices)>0):
        event_couples = [] #first seg info of all tracks connected to same MC event
        tries = []
        indices = [i for i, x in enumerate(MC_events) if x == mcevt]
        for idx in indices:
            event_couples.append(track_couples[idx])

        if (EVERBOSE == 100):
            print(" ################# ")
            print(" Studying track " + str(couple))
            print(" Vertices connected to debug event " + str(vertices))
            print(" Found  " + str(len(event_couples)) + " tracks related to same MC event ")

        closest_couple, closest_b, closest_gap, closest_inv_b = FindClosestMatch(event_couples, s0x, s0y, s0z, s0tx, s0ty, s0plate, couple, EVERBOSE)
        tries.append(closest_couple[0:2])
        closest_bs.append(closest_b)
        closest_inv_bs.append(closest_inv_b)
        closest_gaps.append(closest_gap)

        #print(closest_couple)

        
        if (EVERBOSE==100):
            print(" First Find Close Match Results: " + str(closest_couple) + " " + str(closest_b))

        while (not (closest_couple[0:2] in vertex_couples) ):
            s0x, s0y, s0z = closest_couple[2], closest_couple[3], closest_couple[4]
            s0tx, s0ty = closest_couple[11], closest_couple[12]
            s0plate = closest_couple[0]
            closest_couple, closest_b, closest_gap, closest_inv_b = FindClosestMatch(event_couples, s0x, s0y, s0z, s0tx, s0ty, s0plate, closest_couple, EVERBOSE)
            if (EVERBOSE==100):
                print(" Find Close Match Results: " + str(closest_couple) + " " + str(closest_b))
            closest_bs.append(closest_b)
            closest_gaps.append(closest_gap)
            closest_inv_bs.append(closest_inv_b)

            if (closest_couple[0:2] in tries): #if it starts circling back to past ones it means the track is not connected to a vertex
                flag = 1
                flag1_cases += 1 
                last_b = -99
                max_b = -99
                inv_b = -99
                if (EVERBOSE==100):
                    print(" Exiting because tries ")
                    print(tries)
                break
            tries.append(closest_couple[0:2])
            if (EVERBOSE==100):
                print(" Second Find Close Match Results: " + str(closest_couple) + " " + str(closest_b))
                print(" Flag " + str(flag))
                if ((closest_couple[0:2] in vertex_couples)):
                    print(" Track with first seg info " + str(closest_couple[0:2]) + " is in a vertex -> exit ")
        
        if (flag == 0):
            last_b = closest_b
            max_b = max(closest_bs)
            max_b_gap = closest_gaps[ closest_bs.index(max_b)]
            inv_b = closest_inv_bs[closest_bs.index(max_b)]
            flag0_cases += 1

    else:
        flag = 2
        last_b = -99
        max_b = -99
        max_b_gap = -99
        inv_b = -99
    
    out_tup.Fill(mcevt, flag, last_b, max_b, start_s0x, start_s0y, start_s0plate, start_s0id, start_s0z, start_s0tx, start_s0ty, max_b_gap, inv_b)

    if (c%1000 == 0):
        print(" Completed " + str(100.*c/evt_tree.GetEntries()) + " % (Time: " + str(time.time()-t1) + " s)")


OutFile.cd()
out_tup.Write()
OutFile.Close()

print(" Total Number of Events checked: " + str(evt_tree.GetEntries()))
print(" Total Number of Events with NO vertex associated: " + str(zero_vertices))
print(" Total Number of Events where it was NOT possible to find an existing vertex from the track: " + str(flag1_cases))
print(" Total Number of Events where it was possible to find an existing vertex from the track: " + str(flag0_cases))


TrackFile.Close()
VertexFile.Close()
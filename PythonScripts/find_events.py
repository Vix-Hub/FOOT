import ROOT as r 
import numpy as np 
import fedrarootlogon 
import time

MC_ID = 2
DEBUG_TRK_ID = 13234
EVERBOSE = 1
Nmax = 100
TrackFileName = "b00000" + str(MC_ID) + ".0.1.2.trk.root" #S2 track file name
VertexFileName = "vertices_improved_fast_3_Z.root" #Vertex File name
#OutFileName = "find_events.root"

TrackFile = r.TFile(TrackFileName, "READ")
VertexFile = r.TFile(VertexFileName, "READ")

vrec = VertexFile.Get("EdbVertexRec") #Vertex Object

n_vertices = vrec.eVTX.GetEntries()
vertex_couples = []
vertex_ids, vertex_MC_evts = [], []

outfile3 = open("vertex_couples.txt", "w")
outfile3.write(" First Segment Plate and ID in S2 for vertex tracks \n")

print(" --- Saving Vertex Tracks Info --- ")

for i in range(n_vertices):
    vertex = vrec.eVTX.At(i)
    for j in range(vertex.N()):
        track = vertex.GetTrack(j)

        if (not (track.GetSegment(0).MCEvt() in vertex_MC_evts)):
            vertex_MC_evts.append(track.GetSegment(0).MCEvt()) #add the MCevt even if the track does not reach S2
            vertex_ids.append(vertex.ID()) #need to be same length as previous list

        if(track.Track()==DEBUG_TRK_ID or track.ID()==DEBUG_TRK_ID):
            print(" Here for debug track with segment last plate " + str(track.GetSegmentLast().Plate()) + " and first seg aid(1) " + str(track.GetSegmentFirst().Aid(1)) + " and last seg aid (1) " + str(track.GetSegmentLast().Aid(1)) +  "\n")
        if (track.GetSegmentLast().Plate()>31 and track.GetSegmentFirst().Aid(1)==0 and track.GetSegmentLast().Aid(1)==0): #if the track reaches S2 and is not broken and is a primary fragment
            for iseg in range(track.N()):
                if (track.GetSegment(iseg).Plate()>=31):
                    if(track.Track()==DEBUG_TRK_ID or track.ID()==DEBUG_TRK_ID):
                        print(" Appending " + str((track.GetSegment(iseg).Plate(), track.GetSegment(iseg).ID())) + " for debug track")
                    vertex_couples.append((track.GetSegment(iseg).Plate(), track.GetSegment(iseg).ID()))
                    outfile3.write( str( track.GetSegment(iseg).Plate()) + " " + str(track.GetSegment(iseg).ID()))
                    outfile3.write("\n")
                    vertex_ids.append(vertex.ID())
                    vertex_MC_evts.append(track.GetSegment(iseg).MCEvt())
                    break  #saving info of the first segments in S2 from vertex tracks

outfile3.close()

print("Found " + str(len(vertex_couples)) + " tracks connected to vertices in S2")

tracks = TrackFile.Get("tracks")

MC_events_list = []

print(" --- Saving S2 Tracks Info --- ")
counter = 0
t0 = time.time()

outfile = open("events2.txt", "w")
outfile.write(" Track ID - s0plate - s0id - s0X - s0Y - s0Z - s0TX - s0TY - MCevt - vertices \n")
outfile2 = open("events.txt", "w")

#out_tree = r.TTree("events", "Events Tree")

#s0plates = np.zeros(1, np.intc)
#s0ids = np.zeros(1000, np.intc)
#vIDs = np.zeros(1000, np.intc)
#N = np.zeros(1, np.intc)
#MCevts = np.zeros(1, np.intc)

#out_tree.Branch("s0plate", s0plates, "s0plate/I")
#out_tree.Branch("N", N, "N/I") #number of vertices conected to same MC event 
#out_tree.Branch("s0id", s0ids, "s0id/I") #first segment info
#out_tree.Branch("vIDs", vIDs, "vIDs[N]/I")
#out_tree.Branch("MCEvt", MCevts, "MCEvt/I")

for track in tracks:

    if (tracks.s[0].Plate()<30): #using S1-S2 but only concerned about S2 tracks
        continue

    if (track.s[0].Aid(1)==0 and track.s[track.nseg-1].Aid(1)==0): #if track is a primary fragment, not "broken"

        if ( (track.s[0].Plate(), track.s[0].ID()) in vertex_couples): # if the track is in a vertex (check on first segment)
            continue
        else:
            MC_events_list.append(track.s[0].MCEvt())

            #s0ids[0] = int(track.s[0].ID())
            #s0plates[0] = int(track.s[0].Plate())
            #MCevts[0] = int(track.s[0].MCEvt())
            
            if (EVERBOSE == 1):
                vtxs = []
                for k, event in enumerate(vertex_MC_evts):
                    if (event == track.s[0].MCEvt()): #find the vertices related to the same MC event
                        vtxs.append(vertex_ids[k])
            
                vtxs = np.unique(vtxs)
                #N = len(list(vtxs))
                #for i in range(N):
                    #vIDs[i] = int(vtxs[i])

                outfile2.write(" Debug Track " + str(track.t.ID()) + " Connected to Event " + str(track.s[0].MCEvt()) + " and to Vertices " + str(np.unique(vtxs)))
                outfile.write( str(track.t.ID()) + " " + str(track.s[0].Plate()) + " " + str(track.s[0].ID()) + " " + str(track.s[0].X()) + " " + str(track.s[0].Y()) + " " + str(track.s[0].Z()) + " " + str(track.s[0].TX()) + " " + str(track.s[0].TY()) + " " +  str(track.s[0].MCEvt()) + " ")
                if (len(vtxs)>0):
                    for k in range(len(vtxs)):
                        outfile.write(str(vtxs[k]) + " ")
                else:
                    outfile.write("-99")

                outfile.write("\n")

                if (len(vtxs)>0):
                    outfile2.write(" draw_vertex_trackV(" + str(vtxs[0]) + ", " + str(track.s[0].Plate()) + ", " + str(track.s[0].ID()) +") ")
                outfile2.write("\n")
    counter += 1
    if (counter%5000 == 0):
        print(" Completed " + str(100.*counter/tracks.GetEntries()) + " %, Time: " + str(time.time()-t0) + " s")

print(" Found a total of " + str(len(MC_events_list)) + " events to check! ")
print(" ")
outfile.close()
#print(MC_events_list)

import ROOT as r 
import fedrarootlogon
import numpy as np 
from debug_libraries import FillVertexTracks

#script to evaluate distance between flag=1 tracks (classified in check_missing_merged.py) with vertices to which they are linked

vertexFile = r.TFile("vertices.root", "READ")
vrec = vertexFile.Get("EdbVertexRec")
n_vertices = vrec.eVTX.GetEntries()

checkFile = r.TFile("check_merged.root", "READ") #should sample from tracks also not in final vertex tree
tup = checkFile.Get("tup")

MC_events, track_couples, track_coordinates, tracks_recovered = [], [], [], []
temp_vIDs = []

# save which tracks are recovered in following postvertexing steps
VertexFile2 = r.TFile("vertices_improved_fast_3_new_temp_2.root", "READ")
vrec2 = VertexFile2.Get("EdbVertexRec")
vtx_tracks_2 = FillVertexTracks(vrec2, 1, 1)

# save FLAG=1 track info
for trackflag in tup:
    if (trackflag.ev_flag==1):
        MC_events.append(trackflag.ev_ID)
        track_couples.append((trackflag.s0_plate, trackflag.s0_id))
        track_coordinates.append((trackflag.s0X, trackflag.s0Y, trackflag.s0Z, trackflag.s0TX, trackflag.s0TY))
        if ((trackflag.s0_plate, trackflag.s0_id) in vtx_tracks_2):
            tracks_recovered.append(1)
        else:
            tracks_recovered.append(0)
    
unique_MC_events = list(np.unique(MC_events))
MC_event_vertices = []  #list of vertices with tracks coming from the corresponding MC Event in unique_MC_Events list
MC_event_vertices_coordinates = []

for i in range(len(unique_MC_events)):
    MC_event_vertices.append([])
    MC_event_vertices_coordinates.append([])

#save first seg info of vertex tracks
for i in range(n_vertices):
    vertex = vrec.eVTX.At(i)
    appended_vids = []
    for j in range(vertex.N()):
        track = vertex.GetTrack(j)

        if (vertex.ID() not in appended_vids):
            linking_mc = track.GetSegmentFirst().MCEvt()
            if ( not linking_mc in MC_events): #not all MC events have been classified as flag=1
                continue
            position = unique_MC_events.index(linking_mc)
            MC_event_vertices[position].append(vertex.ID())
            MC_event_vertices_coordinates[position].append((vertex.VX(), vertex.VY(), vertex.VZ()))
            appended_vids.append(vertex.ID())

    if (i%1000 == 0):
        print(" Vertex Step - Completed " + str(100.*i/n_vertices) + " %")

# MC_event_vertices -> list of vertex IDs connected to given MC Event in MC_events
# MC_event_vertices_coordinates -> list of vertex (vx,vy,vz) of vertices connected to given MC Event

# now calculate distance from vertices 
outFile = r.TFile("ip_flag1.root", "RECREATE")
out_tup = r.TNtuple("tup", "IP of Flag=1 tracks with vertices", "b:DZ:is_rec")

for j, unconn_track in enumerate(track_coordinates):
    x, y, z, tx, ty = unconn_track[0], unconn_track[1], unconn_track[2], unconn_track[3], unconn_track[4]

    mcevent = MC_events[j]
    position = unique_MC_events.index(mcevent) 

    Nvtx = len(MC_event_vertices[position])
    bs, dzs = [], []
    for i in range(Nvtx):
        vx, vy, vz = MC_event_vertices_coordinates[position][i][0], MC_event_vertices_coordinates[position][i][1], MC_event_vertices_coordinates[position][i][2]
        dz = z - vz
        xp = x - dz*tx - vx
        yp = y - dz*ty - vy
        b = np.sqrt(xp*xp + yp*yp)
        bs.append(b)
        dzs.append(dz)
    out_tup.Fill(min(bs), dzs[bs.index(min(bs))], tracks_recovered[j])

outFile.cd()
out_tup.Write()
outFile.Close()

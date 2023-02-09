import ROOT as r 
import fedrarootlogon
import numpy as np 

def IsTrackInVertex(vertexFileName, DEBUG_S0_PLATE=31, DEBUG_S0_ID=100):
    
    VtxFile = r.TFile(vertexFileName, "READ")
    vrec = VtxFile.Get("EdbVertexRec")

    n_vertices = vrec.eVTX.GetEntries()

    found, vID = 0, -99

    for i in range(n_vertices):
        vertex = vrec.eVTX.At(i)
        for j in range(vertex.N()):
            track = vertex.GetTrack(j)
            for iseg in range(track.N()):
                if ((int(track.GetSegment(iseg).Plate()), int(track.GetSegment(iseg).ID()))==(DEBUG_S0_PLATE, DEBUG_S0_ID)):
                    found = 1
                    vID = vertex.ID()

    return found, vID

# ----------------------------------- #

def GetVertexTracksS0Info(vID, vertexFileName):
    VtxFile = r.TFile(vertexFileName, "READ")
    vrec = VtxFile.Get("EdbVertexRec")

    n_vertices = vrec.eVTX.GetEntries()
    couples = []
    for i in range(n_vertices):
        vertex = vrec.eVTX.At(i)
        if (vertex.ID()==vID):
            for j in range(vertex.N()):
                track = vertex.GetTrack(j)
                couples.append((track.GetSegmentFirst().Plate(), track.GetSegmentFirst().ID()))
    
    return couples

# ----------------------------------- #

def GetEventTracksS0Info(MCevt, tracksFileName):
    TrkFile = r.TFile(tracksFileName, "READ")
    tracks = TrkFile.Get("tracks")

    couples = []
    for track in tracks:
        if (track.s[0].MCEvt()==MCevt):
            couples.append((track.s[0].Plate(), track.s[0].ID()))
    
    return couples

# ------------------------------ # 

def IsVertexInTree(vID, vertexFileName):
    VtxFile = r.TFile(vertexFileName, "READ")
    vrec = VtxFile.Get("EdbVertexRec")
    n_vertices = vrec.eVTX.GetEntries()
    for i in range(n_vertices):
        vertex = vrec.eVTX.At(i)
        if (vID==vertex.ID()):
            return 1
    return 0

# -------------------------------- # 

def FindNewVtxID(vtxFileName1, vtxFileName2, vID):

    VtxFile1 = r.TFile(vtxFileName1, "READ")
    VtxFile2 = r.TFile(vtxFileName2, "READ")
    vrec1 = VtxFile1.Get("EdbVertexRec")
    vrec2 = VtxFile2.Get("EdbVertexRec")

    tracks_info = []
    n_vertices1 = vrec1.eVTX.GetEntries()
    for i in range(n_vertices1):
        vertex = vrec1.eVTX.At(i)
        if (vertex.ID() == vID):
            for j in range(vertex.N()):
                track = vertex.GetTrack(j)
                tracks_info.append((track.GetSegmentFirst().Plate(), track.GetSegmentFirst().ID()))

    n_vertices2 = vrec2.eVTX.GetEntries()
    for i in range(n_vertices2):
        vertex = vrec2.eVTX.At(i)
        for j in range(vertex.N()):
            track = vertex.GetTrack(j)
            if ((track.GetSegmentFirst().Plate(), track.GetSegmentFirst().ID()) in tracks_info):
                print(" Found Shared Track " + str(track.Track()) + " with vertex in " + vtxFileName2  + " ID " + str(vertex.ID()))
        
    return 1
# ---------------------------------------------- #

def CalculateImpactParameter(slx, sly, slz, slftx, slfty, s0x, s0y, s0z):
    dxp = slx - (slz-s0z)*slftx -s0x
    dyp = sly - (slz-s0z)*slfty - s0y
    b = np.sqrt((dxp)*(dxp)+(dyp)*(dyp))
    return b

# ---------------------------------------------- #

def CalculateImpactParameterBetweenTracks(TrkFileName, DEBUG_S0_PLATE, DEBUG_S0_ID, DEBUG_S0_PLATE2, DEBUG_S0_ID2):
    TrackFile = r.TFile(TrkFileName, "READ")
    tracks = TrackFile.Get("tracks")
    tracks.BuildIndex("s[0].Plate()", "s[0].ID()")

    tracks.GetEntryWithIndex(int(DEBUG_S0_PLATE), int(DEBUG_S0_ID))

    slx, sly, slz = tracks.s[tracks.nseg-1].X(), tracks.s[tracks.nseg-1].Y(), tracks.s[tracks.nseg-1].Z()
    sltx, slty = tracks.s[tracks.nseg-1].TX(), tracks.s[tracks.nseg-1].TY()

    tracks.GetEntryWithIndex(int(DEBUG_S0_PLATE2), int(DEBUG_S0_ID2))
    s0x, s0y, s0z = tracks.s[0].X(), tracks.s[0].Y(), tracks.s[0].Z()

    print("Calculated b " + str(CalculateImpactParameter(slx, sly, slz, sltx, slty, s0x, s0y, s0z)))
    return CalculateImpactParameter(slx, sly, slz, sltx, slty, s0x, s0y, s0z)

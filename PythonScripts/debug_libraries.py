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

def GetVertex(VtxFileName, vID):
    VtxFile = r.TFile(VtxFileName, "READ")
    vrec = VtxFile.Get("EdbVertexRec")
    n_vertices = vrec.eVTX.GetEntries()
    for i in range(n_vertices):
        vertex = vrec.eVTX.At(i)
        if (vertex.ID()==vID):
            return vertex

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

# ---------------------------------------------- #

def CheckTrackSegmentsAid(TrkFileName, DEBUG_S0_PLATE, DEBUG_S0_ID):
    TrackFile = r.TFile(TrkFileName, "READ")
    tracks = TrackFile.Get("tracks")
    tracks.BuildIndex("s[0].Plate()", "s[0].ID()")

    tracks.GetEntryWithIndex(int(DEBUG_S0_PLATE), int(DEBUG_S0_ID))

    for segment in tracks.s:
        print(" Seg ID " + str(segment.ID()) + " Seg Plate " + str(segment.Plate()) + " Seg Aid[0] " + str(segment.Aid(0)) + " Seg Aid [1] " + str(segment.Aid(1)))
    
    return 1

# ---------------------------------------------- #

def CheckTrackSegmentsAid_FromVrec(VtxFileName, DEBUG_S0_PLATE, DEBUG_S0_ID):
    VtxFile = r.TFile(VtxFileName, "READ")
    vrec = VtxFile.Get("EdbVertexRec")
    n_vertices = vrec.eVTX.GetEntries()
    for i in range(n_vertices):
        vertex = vrec.eVTX.At(i)
        for j in range(vertex.N()):
            track = vertex.GetTrack(j)
            found = 0
            for iseg in range(track.N()):
                    segment = track.GetSegment(iseg)
                    if (segment.Plate()==int(DEBUG_S0_PLATE) and segment.ID()==int(DEBUG_S0_ID)):
                        found = 1
            if (found == 1):
                for iseg in range(track.N()):
                    segment = track.GetSegment(iseg)
                    if (segment.Plate()>=31):
                        print(" Seg ID " + str(segment.ID()) + " Seg Plate " + str(segment.Plate()) + " Seg Aid[0] " + str(segment.Aid(0)) + " Seg Aid [1] " + str(segment.Aid(1)))
    
    return 1

# -------------------------------------------------- #

def CalcDistFromTrack_toVertex(VtxFileName, TrkFileName, vID, s0_plate, s0_id):

    VtxFile = r.TFile(VtxFileName, "READ")
    vrec = VtxFile.Get("EdbVertexRec")
    n_vertices = vrec.eVTX.GetEntries()
    v = 0
    for i in range(n_vertices):
        vertex = vrec.eVTX.At(i)
        if (vertex.ID()==vID):
            v = vertex

    TrackFile = r.TFile(TrkFileName, "READ")
    tracks = TrackFile.Get("tracks")
    tracks.BuildIndex("s[0].Plate()", "s[0].ID()")

    tracks.GetEntryWithIndex(int(s0_plate), int(s0_id))
    dz = tracks.s[0].Z() - v.VZ()
    x = tracks.s[0].X() - dz*tracks.s[0].TX()
    y = tracks.s[0].Y() - dz*tracks.s[0].TY()
    dx = x - vertex.VX()
    dy = y - vertex.VY()

    return np.sqrt(dx*dx + dy*dy)

# -------------------------------------------------- #

def Count_Merged_Vertex_Tracks(VtxFileName):

    VtxFile = r.TFile(VtxFileName, "READ")
    vrec = VtxFile.Get("EdbVertexRec")
    n_vertices = vrec.eVTX.GetEntries()
    merged_trks = 0
    affected_vertices = 0

    for i in range(n_vertices):
        vertex = vrec.eVTX.At(i)
        counted = 0
        for j in range(vertex.N()):
            track = vertex.GetTrack(j)
            start_MCevt = track.GetSegmentFirst().MCEvt()
            for iseg in range(track.N()):
                seg = track.GetSegment(iseg)
                current_MCevt = seg.MCEvt()
                if (current_MCevt != start_MCevt):
                    merged_trks += 1 
                    if (counted == 0):
                        affected_vertices += 1
                        counted = 1

    return merged_trks, affected_vertices


# -------------------------------------------------- #


def PrintVertexTracksInfo(VtxFileName, vID, modality):
    
    VtxFile = r.TFile(VtxFileName, "READ")
    vrec = VtxFile.Get("EdbVertexRec")
    n_vertices = vrec.eVTX.GetEntries()

    if (modality == 1):
        for i in range(n_vertices):
            vertex = vrec.eVTX.At(i)
            if (vertex.ID()!=vID):
                continue
            for j in range(vertex.N()):
                track = vertex.GetTrack(j)
                print(" Track # " + str(j) + " TrackTrack " + str(track.Track()) + " First Plate " + str(track.GetSegmentFirst().Plate()) + " First Seg ID " + str(track.GetSegmentFirst().ID()) + " " + str(track.GetSegmentLast().Plate()) + " " + str(track.GetSegmentLast().ID()) + " \n")
    else:
        for i in range(n_vertices):
            vertex = vrec.eVTX.At(i)
            if (vertex.ID() == vID):
                for j in range(vertex.N()):
                    track = vertex.GetTrack(j)
                    print(" Track # " + str(j) + " TrackTrack " + str(track.Track()) + " First Plate " + str(track.GetSegmentFirst().Plate()) + " First Seg ID " + str(track.GetSegmentFirst().ID()) + " \n")

    return 1 


# -----------------------------------------------  #

def PrintTrack(TrkFileName, s0plate, s0id, EVERBOSE=1):
    
    TrackFile = r.TFile(TrkFileName, "READ")
    tracks = TrackFile.Get("tracks")
    tracks.BuildIndex("s[0].Plate()", "s[0].ID()")

    tracks.GetEntryWithIndex(int(s0plate), int(s0id))

    for iseg, seg in enumerate(tracks.s):
        print(" Seg # " + str(iseg) + " Plate " + str(seg.Plate()) + " ID " + str(seg.ID()) + " MCEvt " + str(seg.MCEvt()) + " W-70 " + str(seg.W()-70) + "\n X, Y, Z, TX, TY ")
        print([seg.X(), seg.Y(), seg.Z(), seg.TX(), seg.TY()])
        print(" seg.firstplate " + str(seg.Vid(0)) + " TRUE lastplate " + str(seg.Vid(1)))
        print("\n")
        if(EVERBOSE==0):
            break

    return 1


# --------------------------------------------  #

def FillVertexMergedTracks(vrec):
    n_vertices = vrec.eVTX.GetEntries()
    out_list = []
    for i in range(n_vertices):
        vertex = vrec.eVTX.At(i)
        for j in range(vertex.N()):
            track = vertex.GetTrack(j)
            if (track.GetSegmentFirst().Plate()<30 and track.GetSegmentLast().Plate()>30):
                out_list.append((track.GetSegmentFirst().Plate(), track.GetSegmentFirst().ID()))
    return out_list
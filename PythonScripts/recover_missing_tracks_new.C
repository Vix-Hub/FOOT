
// script to link found tracks to vertices 
const int DEBUG_VID = 21;
const int EVERBOSE = -100;

const float MAX_B = 50;
const float MAX_DXDY = 500;
const float MAX_DZ = 1300*5;
const float MAX_DZ_SHIFT = 300;


EdbVertex* GetVertex(EdbVertexRec* vrec, int vID);
void FillNewVtxTree(TTree* new_vtxtree, EdbVertexRec* vrec_new, TTree* vtxtree, bool MC);

void recover_missing_tracks_new() {

    // Reading Info From Files
    TFile *VtxFile = TFile::Open("vertices_final_Z.root", "READ");
    EdbVertexRec* vrec = (EdbVertexRec*)VtxFile->Get("EdbVertexRec");
    TTree* oldvtx = (TTree*)VtxFile->Get("vtx");

    TFile *LinkFile = TFile::Open("recover.root","READ");
    TNtuple* tup = (TNtuple*)LinkFile->Get("tup");

    float vID=0, s0plate=0, s0id=0, s0flag=0;
    float b=0, dx=0, dy=0, dz=0;
    tup->SetBranchAddress("vID", &vID); tup->SetBranchAddress("s0plate", &s0plate); tup->SetBranchAddress("s0id",&s0id);
    tup->SetBranchAddress("Zflag", &s0flag); tup->SetBranchAddress("b", &b); tup->SetBranchAddress("DX", &dx);
    tup->SetBranchAddress("DY", &dy); tup->SetBranchAddress("DZ", &dz);


    TFile* trackFile = TFile::Open("b000333.0.1.2.trk_merged_new_2.root", "READ");
    TTree* tracksTree = (TTree*)trackFile->Get("tracks");
    tracksTree->BuildIndex("s[0].Plate()", "s[0].ID()");

    EdbSegP *trk=nullptr; //important to initialize to NULL otherwise seg fault!!!
    TClonesArray *segments = new TClonesArray("EdbSegP");
    TClonesArray *fitted_segments = new TClonesArray("EdbSegP");
    
	int nseg=0, trid=0;
    tracksTree->SetBranchAddress("trid",&trid);
    tracksTree->SetBranchAddress("nseg",&nseg);
    tracksTree->SetBranchAddress("t.",&trk);

    tracksTree->SetBranchAddress("s", &segments);
    tracksTree->SetBranchAddress("sf", &fitted_segments);
    tracksTree->SetBranchAddress("nseg", &nseg);

    //OUTFILE
    TFile* OutFile = TFile::Open("vertices_final_new.root", "RECREATE");
    EdbVertexRec *vrec_new = new EdbVertexRec();
    TTree* new_vtx = new TTree();



    TObjArray* VTX_new = new TObjArray();
    int z1counter = 0;

    // Get List of Vertices that should be modified
    std::vector<int> vIDs;

    for(int ientry=0; ientry<tup->GetEntries(); ientry++) {
        tup->GetEntry(ientry);
        if (b<MAX_B && TMath::Abs(dx)<MAX_DXDY && TMath::Abs(dy)<MAX_DXDY && dz>0 && dz<MAX_DZ && s0flag<6) vIDs.push_back(vID);
    }

    // Get List of Unique Vertices 
    auto it = std::unique(vIDs.begin(), vIDs.end());
    vIDs.erase(it, vIDs.end());

    // Get all tracks corresponding to the same vertex
    std::vector<std::vector<std::pair<int, int>>> trk_list;

    for (auto vertexID : vIDs) {
        //std::cout << vertexID << std::endl;
        std::vector<std::pair<int, int>> temptracks;

        for (int i = 0; i < tup->GetEntries(); i++) {
            tup->GetEntry(i);

            if (vID == vertexID && b < MAX_B && std::abs(dx) < MAX_DXDY && std::abs(dy) < MAX_DXDY &&
                dz > 0 && dz < MAX_DZ && s0flag < 6) {
                temptracks.push_back(std::make_pair(s0plate, s0id));
            }
        }

        trk_list.push_back(temptracks);

        if (vertexID == DEBUG_VID && EVERBOSE==100) {
            std::cout << " List of tracks for vertex " << DEBUG_VID << " is ";
            for (const auto& track : trk_list.back()) {
                std::cout << "(" << track.first << ", " << track.second << ") ";
            }
            std::cout << std::endl;
        }
    }

    // Remove duplicate tracks within each vertex
    for (auto& tracks_index : trk_list) {
        std::sort(tracks_index.begin(), tracks_index.end());
        tracks_index.erase(std::unique(tracks_index.begin(), tracks_index.end()), tracks_index.end());
    }

    // Add the tracks to the vertices 

    for (size_t j = 0; j < vIDs.size(); ++j) {
        int vID = vIDs[j];

        if (vID!=DEBUG_VID && EVERBOSE==100) continue;

        TObjArray tracks_to_add;
        std::vector<std::pair<int, int>> couples = trk_list[j];
        EdbVertex* vertex = GetVertex(vrec, vID);
        if (EVERBOSE==100) std::cout << " getting vertex " << vertex->ID() << " vID " << vID << endl;
        int tempz1counter = 0;

        for (int it = 0; it < vertex->N(); it++) {
            tracks_to_add.Add(vertex->GetTrack(it));
        }

        for (const auto& couple : couples) {

            tracksTree->GetEntryWithIndex(couple.first, couple.second);

            EdbTrackP* temptrack = new EdbTrackP();
            //temptrack->Copy(*trk);
            temptrack->SetID(trk->ID());
            temptrack->SetTrack(trk->Track());

            // Loop on segments associated with the track
            for (int jseg = 0; jseg < nseg; jseg++) {
                EdbSegP* seg = (EdbSegP*)segments->At(jseg);
                EdbSegP* segf = (EdbSegP*)fitted_segments->At(jseg);

                EdbSegP* seg_new = new EdbSegP(); //ausiliary  //ausiliary approach leads to seg fault due to missing COV matrix?
                seg_new->SetX(seg->X()); seg_new->SetY(seg->Y()); seg_new->SetZ(seg->Z());
                seg_new->SetTX(seg->TX()); seg_new->SetTY(seg->TY()); seg_new->SetW(seg->W());
                seg_new->SetFlag(seg->Flag()); seg_new->SetID(seg->ID()); seg_new->SetPlate(seg->Plate());
                seg_new->SetCOV(seg->COV());

                EdbSegP* segf_new = new EdbSegP(); //ausiliary
                segf_new->SetX(segf->X()); segf_new->SetY(segf->Y()); segf_new->SetZ(segf->Z());
                segf_new->SetTX(segf->TX()); segf_new->SetTY(segf->TY()); segf_new->SetW(segf->W());
                segf_new->SetFlag(segf->Flag()); segf_new->SetID(segf->ID()); segf_new->SetPlate(segf->Plate());
                segf_new->SetCOV(segf->COV());

                temptrack->AddSegment(seg_new);
                temptrack->AddSegmentF(segf_new);

                temptrack->SetSegmentsTrack(temptrack->ID());
                temptrack->SetCounters();
            }

            if (EVERBOSE==100) std::cout << " temptrack N " << temptrack->N() << endl;
            tracks_to_add.Add(temptrack);

            if (temptrack->GetSegmentFirst()->Flag() == 1) {
                tempz1counter++;
            }
        }

        if (j%1000==0) std::cout << " Completed " << 100.*j/vIDs.size() << " % " << endl;
        
        if (vertex->N() == 0) {
            continue;
        }

        if (EVERBOSE==100){
            cout << "Tracks to add size: " << tracks_to_add.GetEntriesFast() << endl;
            for (int i = 0; i < tracks_to_add.GetEntriesFast(); i++) {
                EdbTrackP* track_temp = (EdbTrackP*)(tracks_to_add.At(i));
                if (track_temp) {
                    cout << "Track " << i << ": ID=" << track_temp->GetSegmentFirst()->ID() << ", Plate=" << track_temp->GetSegmentFirst()->Plate() << ", Flag=" << track_temp->GetSegmentFirst()->Flag() << endl;
                } else {
                    cout << "Error: Null track at index " << i << " in tracks_to_add." << endl;
                }
            }
        }
        
        EdbVertex* vertex_new = vrec->Make1Vertex(tracks_to_add, vertex->VZ());
        vertex_new->SetID(vertex->ID());

        if (vertex_new && EVERBOSE==100) {
            for (int it=0; it < vertex_new->N(); it++) {
                EdbTrackP* track_temp = vertex_new->GetTrack(it);
                cout << " vertex new trk " << track_temp->GetSegmentFirst()->ID() << " " << track_temp->GetSegmentFirst()->Plate() << " " << track_temp->GetSegmentFirst()->Flag() << endl;
            }
        }
        
        if (std::abs(vertex->VZ() - vertex_new->VZ()) < MAX_DZ_SHIFT) {
            VTX_new->Add(vertex_new);
            z1counter += tempz1counter;
        } else {
            VTX_new->Add(vertex);
        }
        
        //tracks_to_add.Delete(); leads to seg fault

    }

    // add remaining vertices
    for (int ivtx=0; ivtx<vrec->eVTX->GetEntries(); ivtx++) {
        EdbVertex* vertex = (EdbVertex*)vrec->eVTX->At(ivtx);

        auto it = std::find(vIDs.begin(), vIDs.end(), vertex->ID());
        if (it == vIDs.end() ) VTX_new->Add(vertex);
    }


   OutFile->cd();

   cout << " VTX new entries " << VTX_new->GetEntries() << endl;
   vrec_new->eVTX = (TObjArray*)VTX_new->Clone();
   vrec_new->Write("EdbVertexRec");

   FillNewVtxTree(new_vtx, vrec_new, oldvtx, 0);
   new_vtx->Write("vtx");
   OutFile->Close();



}








EdbVertex* GetVertex(EdbVertexRec* vrec, int vID) {
    int n_vertices = vrec->eVTX->GetEntries();
    for (int i = 0; i < n_vertices; ++i) {
        EdbVertex* vertex = static_cast<EdbVertex*>(vrec->eVTX->At(i));
        if (vertex->ID() == vID) {
            return vertex;
        }
    }
    return nullptr;  // Return nullptr if the vertex with the specified ID is not found
}



void FillNewVtxTree(TTree* new_vtxtree, EdbVertexRec* vrec_new, TTree* vtxtree, bool MC) {
    const int maxdim = 5000;
    
    Int_t vID;
    Float_t vx, vy, vz, maxaperture, probability;
    Int_t n, v_flag, v_plate;
    Int_t IDTrack[maxdim], TrackTrack[maxdim];
    
    Int_t s0id[maxdim], s0plate[maxdim];
    
    Float_t X[maxdim], Y[maxdim], Z[maxdim], TX[maxdim], TY[maxdim];
    Float_t Theta[maxdim], impactparameter[maxdim];
    
    Int_t nseg_vec[maxdim], npl[maxdim], nholes[maxdim];
    Int_t nseg_S1[maxdim], npl_S1[maxdim], nholes_S1[maxdim];
    
    Int_t plate[maxdim], maxgap[maxdim], incoming[maxdim], Z_flag[maxdim];
    
    Int_t MC_Charge_first[maxdim], MC_Charge_last[maxdim], MC_Charge_S2[maxdim];
    Int_t MC_evID_first[maxdim], MC_evID_last[maxdim];
    Int_t MC_trackID_first[maxdim], MC_trackID_last[maxdim];
    Int_t MC_mother_first[maxdim], MC_mother_last[maxdim];
    Int_t MC_pdgcode_first[maxdim], MC_pdgcode_last[maxdim];
    Int_t MC_firstplate[maxdim], MC_lastplate[maxdim];
    
    Int_t n_10, n_5, nseg_10, nseg_5, n_outgoing, n_S2;

    new_vtxtree->Branch("vID", &vID, "vID/I");
    new_vtxtree->Branch("vx", &vx, "vx/F");
    new_vtxtree->Branch("vy", &vy, "vy/F");
    new_vtxtree->Branch("vz", &vz, "vz/F");
    new_vtxtree->Branch("vplate", &v_plate, "vplate/I");
    new_vtxtree->Branch("v_flag", &v_flag, "v_flag/I");
    new_vtxtree->Branch("maxaperture", &maxaperture, "maxaperture/F");
    new_vtxtree->Branch("probability", &probability, "probability/F");
    new_vtxtree->Branch("n", &n, "n/I");
    new_vtxtree->Branch("IDTrack", IDTrack, "IDTrack[n]/I");
    new_vtxtree->Branch("TrackTrack", TrackTrack, "TrackTrack[n]/I");

    new_vtxtree->Branch("s0plate", s0plate, "s0plate[n]/I");
    new_vtxtree->Branch("s0id", s0id, "s0id[n]/I");

    new_vtxtree->Branch("X", X, "X[n]/F");
    new_vtxtree->Branch("Y", Y, "Y[n]/F");
    new_vtxtree->Branch("Z", Z, "Z[n]/F");
    new_vtxtree->Branch("TX", TX, "TX[n]/F");
    new_vtxtree->Branch("TY", TY, "TY[n]/F");
    new_vtxtree->Branch("Theta", Theta, "Theta[n]/F");
    new_vtxtree->Branch("nseg", nseg_vec, "nseg[n]/I");
    new_vtxtree->Branch("npl", npl, "npl[n]/I");
    new_vtxtree->Branch("nholes", nholes, "nholes[n]/I");
    new_vtxtree->Branch("nseg_S1", nseg_S1, "nseg_S1[n]/I");
    new_vtxtree->Branch("nholes_S1", nholes_S1, "nholes_S1[n]/I");
    new_vtxtree->Branch("plate", plate, "plate[n]/I");
    new_vtxtree->Branch("maxgap", maxgap, "maxgap[n]/I");
    new_vtxtree->Branch("incoming", incoming, "incoming[n]/I");
    new_vtxtree->Branch("impactparameter", impactparameter, "impactparameter[n]/F");
    new_vtxtree->Branch("Z_flag2", Z_flag, "Z_flag2[n]/I");

    new_vtxtree->Branch("MC_Charge_first", MC_Charge_first, "MC_Charge_first[n]/I");
    new_vtxtree->Branch("MC_Charge_last", MC_Charge_last, "MC_Charge_last[n]/I");
    new_vtxtree->Branch("MC_Charge_S2", MC_Charge_S2, "MC_Charge_S2[n]/I");
    new_vtxtree->Branch("MC_evID_first", MC_evID_first, "MC_evID_first[n]/I");
    new_vtxtree->Branch("MC_evID_last", MC_evID_last, "MC_evID_last[n]/I");
    new_vtxtree->Branch("MC_trackID_first", MC_trackID_first, "MC_trackID_first[n]/I");
    new_vtxtree->Branch("MC_trackID_last", MC_trackID_last, "MC_trackID_last[n]/I");
    new_vtxtree->Branch("MC_mother_first", MC_mother_first, "MC_mother_first[n]/I");
    new_vtxtree->Branch("MC_mother_last", MC_mother_last, "MC_mother_last[n]/I");
    new_vtxtree->Branch("MC_pdgcode_first", MC_pdgcode_first, "MC_pdgcode_first[n]/I");
    new_vtxtree->Branch("MC_pdgcode_last", MC_pdgcode_last, "MC_pdgcode_last[n]/I");
    new_vtxtree->Branch("MC_firstplate", MC_firstplate, "MC_firstplate[n]/I");
    new_vtxtree->Branch("MC_lastplate", MC_lastplate, "MC_lastplate[n]/I");

    new_vtxtree->Branch("n_outgoing", &n_outgoing, "n_outgoing/I");
    new_vtxtree->Branch("n_10", &n_10, "n_10/I");
    new_vtxtree->Branch("n_5", &n_5, "n_5/I");
    new_vtxtree->Branch("nseg_10", &nseg_10, "nseg_10/I");
    new_vtxtree->Branch("nseg_5", &nseg_5, "nseg_5/I");
    new_vtxtree->Branch("n_S2", &n_S2, "n_S2/I");

    int vplate=0;
    vtxtree->SetBranchAddress("vplate", &vplate);

    for (Int_t ivtx = 0; ivtx < vrec_new->eVTX->GetEntries(); ++ivtx) {
        EdbVertex* vertex = (EdbVertex*)(vrec_new->eVTX->At(ivtx));
        vtxtree->GetEntry(ivtx);
        v_plate = vplate;

        vID = vertex->ID();
        vx = vertex->VX();
        vy = vertex->VY();
        vz = vertex->VZ();
        n = vertex->N();
        v_flag = vertex->Flag();
        
        maxaperture = vertex->MaxAperture();
        probability = vertex->V()->prob();

        n_10 = 0; n_5 = 0; nseg_10 = 0; nseg_5 = 0; n_outgoing = 0; n_S2 = 0;

        for (Int_t itrk = 0; itrk < vertex->N(); itrk++) {
            EdbTrackP* track = vertex->GetTrack(itrk);
            EdbSegP* seg = track->GetSegmentFirst();

            if (vertex->GetVTa(itrk)->Zpos()==1) n_outgoing++;

            s0plate[itrk] = seg->Plate();
            s0id[itrk] = seg->ID();

            IDTrack[itrk] = track->Track();
            TrackTrack[itrk] = track->Track();
            nseg_vec[itrk] = track->N();
            npl[itrk] = track->GetSegmentLast()->Plate() - seg->Plate();

            if (track->Npl() >= 5) {
                n_5++;
                if (track->Npl() >= 10)
                    n_10++;
            }

            if (track->N() >= 5) {
                nseg_5++;
                if (track->N() >= 10)
                    nseg_10++;
            }

            if (track->GetSegmentLast()->Plate() > 30)
                n_S2++;

            
            X[itrk] = seg->X();
            Y[itrk] = seg->Y();
            Z[itrk] = seg->Z();

            if (seg->Z() > vertex->VZ()) { incoming[itrk] = 1; vertex->GetVTa(itrk)->SetZpos(1); }
            else { incoming[itrk] = 0; vertex->GetVTa(itrk)->SetZpos(0); }

            TX[itrk] = seg->TX();
            TY[itrk] = seg->TY();
            Theta[itrk] = track->Theta();
            plate[itrk] = seg->Plate();
            nholes[itrk] = track->N0();
            maxgap[itrk] = track->CheckMaxGap();
            nseg_S1[itrk] = npl_S1[itrk] = nholes_S1[itrk] = 0;

            if (plate[itrk] < 31) {
                for (Int_t iseg = 0; iseg < track->N(); ++iseg) {
                    if (track->GetSegment(iseg)->Plate() <= 30)
                        nseg_S1[itrk]++;
                }
                npl_S1[itrk] = 30 - plate[itrk] + 1;
                nholes_S1[itrk] = npl_S1[itrk] - nseg_S1[itrk];
            }

            impactparameter[itrk] = vertex->GetVTa(itrk)->Imp();
            Z_flag[itrk] = seg->Flag();

            if (MC) {
                Z_flag[itrk] = track->GetSegmentLast()->W() - 70;
            }

            for (Int_t iseg = 0; iseg < track->N(); ++iseg) {
                if (track->GetSegment(iseg)->Plate() > 30) {
                    MC_Charge_S2[itrk] = track->GetSegment(iseg)->W() - 70;
                }

                MC_evID_first[itrk] = track->GetSegmentFirst()->MCEvt();
                MC_evID_last[itrk] = track->GetSegmentLast()->MCEvt();
                MC_trackID_first[itrk] = track->GetSegmentFirst()->MCTrack();
                MC_trackID_last[itrk] = track->GetSegmentLast()->MCTrack();
                MC_pdgcode_first[itrk] = track->GetSegmentFirst()->Aid(0);
                MC_pdgcode_last[itrk] = track->GetSegmentLast()->Aid(0);
                MC_Charge_first[itrk] = track->GetSegmentFirst()->W() - 70;
                MC_Charge_last[itrk] = track->GetSegmentLast()->W() - 70;
                MC_mother_first[itrk] = track->GetSegmentFirst()->Aid(1);
                MC_mother_last[itrk] = track->GetSegmentLast()->Aid(1);
                MC_firstplate[itrk] = track->GetSegmentFirst()->Vid(0);
                MC_lastplate[itrk] = track->GetSegmentLast()->Vid(0);
            }
        }

        new_vtxtree->Fill();
    }
}

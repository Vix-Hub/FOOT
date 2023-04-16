#include<algorithm>
#include<utility>
#include<string>


#define IDBRICK 2
#define EVERBOSE -99

const int DEBUG_S0_PLATE=31; //31
const int DEBUG_S0_ID=1395973; //91690

const float xmin = 35000;//0;
const float xmax = 90000;//125000;
const float ymin = 25000;//0;
const float ymax = 75000;//100000;

const int PLMIN=0;
const int PLMAX=66; 
const int START_PLATE_S2 = 30;

EdbCell2 gridtr_S1[START_PLATE_S2+1]; 
EdbCell2 gridtr_S2[int(PLMAX-START_PLATE_S2)+1];
EdbCell2 gridtr_ALL[1+PLMAX]; //first cell left left empty

TObjArray *arrTRK = new TObjArray();   // original tracks
int MERGED=0;
float Z_LAYER[PLMAX+1]={0};

void bubbleSort(std::vector<double>& v);
int findValuePosition(const std::vector<double>& myVector, double valueToFind);
void FillTracksCells(TObjArray &arrt);
double CalcDist(float x1, float y1, float z1, float x2, float y2, float z2, float tx1, float ty1);
int IsElementInVector(const std::vector<string>& myVector, string value);
EdbTrackP* FindClosestCandidate(const int nplates, EdbTrackP* start_trk, TClonesArray *segments_new, TClonesArray *fitted_segments_new, const double MAX_B, int &added_segs);
void bubbleSort_NEW(std::vector<double>& v, std::vector<double>& v2);

double CalcDist(float x1, float y1, float z1, float x2, float y2, float z2, float tx1, float ty1){
    //prendo il primo segmento della traccia, lo proietto alla z del vertice e calcolo la distanza
    double dz=z1-z2;
    double x=x1-dz*tx1;
    double y=y1-dz*ty1;
    
    double dx=x-x2;
    double dy=y-y2;
    
    return TMath::Sqrt(dx*dx+dy*dy);
    
}

int IsElementInVector(const std::vector<string>& myVector, string value) {
    for (int i = 0; i < myVector.size(); i++) {
    if (myVector[i] == value) {
      return i;
    }
  }

  return -1; // return -1 if value is not found
}




EdbTrackP* FindClosestCandidate(const int nplates, EdbTrackP* start_trk, TClonesArray *segments_new, TClonesArray *fitted_segments_new, const double MAX_B, int &added_segs) {
    
    EdbSegP* start_seg = (EdbSegP*)start_trk->GetSegmentFirst(); //i segmenti fittati non hanno i piatti salvati bene
    EdbSegP* start_segf = (EdbSegP*)start_trk->GetSegmentFFirst();
    float xy[2] = {0,0};
    if (EVERBOSE==100) cout << " Start Seg Coordinates " << start_seg->X() << " " << start_seg->Y() << " " << start_seg->Z() << " " << start_seg->TX() << " " << start_seg->TY() << " Plate ID " << start_seg->Plate() << " " << start_seg->ID() << endl;
    int s0_plate = start_trk->GetSegmentFirst()->Plate() - START_PLATE_S2;
    float b=0, r=500, b_back=0;
    
    std::vector<double> impact_parameters, sorted_IPs, impact_parameters_back, impact_parameters_mean, sorted_IPs_for;
    std::vector<int> merge_s0plates, merge_s0ids, merge_s0plate_cand, merge_s0id_cand;
    TObjArray candidate_tracks, tr_grid_s1;

    for (int ipl=(s0_plate)-nplates; ipl<=s0_plate; ipl++) {  //look for candidates in 4 plates

            xy[0] = start_segf->X();
            xy[1] = start_segf->Y();

            if (ipl+START_PLATE_S2-1<PLMIN) continue; //avoid plate < 0
            if (ipl+START_PLATE_S2-1>START_PLATE_S2) continue;  // avoid plate >30 for S1 grid
            double z_pos = Z_LAYER[ipl+START_PLATE_S2-1];
            //cout << " z_pos " << z_pos << endl;
            xy[0] = xy[0] - (start_segf->Z()-z_pos)*start_segf->TX();  //this can be further away if angles are large
            xy[1] = xy[1] - (start_segf->Z()-z_pos)*start_segf->TY();

            //cout << " nplate S1 " << ipl+START_PLATE_S2-1 << endl;
            if (ipl+START_PLATE_S2-1<PLMIN) continue;
            int n_trk_cell = gridtr_S1[ipl+START_PLATE_S2-1].SelectObjectsC(xy, r, tr_grid_s1); //tracce in S1
            if(EVERBOSE==101) cout << " Checking Grid corresponding to plate # " << ipl+START_PLATE_S2-1 << ", found " << n_trk_cell << " objects " << endl;
            if (EVERBOSE==101) cout << " i was looking with xy " << xy[0] << " " << xy[1] << endl;
            for (int itrk1=0; itrk1<n_trk_cell; itrk1++) {
                EdbTrackP* end_trk = (EdbTrackP*)tr_grid_s1.At(itrk1);
                EdbSegP* end_seg = (EdbSegP*)end_trk->GetSegmentLast();
                EdbSegP* end_segf = (EdbSegP*)end_trk->GetSegmentFLast();
                b = CalcDist(end_segf->X(), end_segf->Y(), end_segf->Z(), start_segf->X(), start_segf->Y(), start_segf->Z(), end_segf->TX(), end_segf->TY());
                b_back = CalcDist(start_segf->X(), start_segf->Y(), start_segf->Z(), end_segf->X(), end_segf->Y(), end_segf->Z(), start_segf->TX(), start_segf->TY());
                if (EVERBOSE==100) cout << " Calculated b with track with last seg " << end_seg->Plate() << " " << end_seg->ID() << " : " << b << endl;
                if (b<MAX_B && b_back<MAX_B) {
                    impact_parameters.push_back(b); 
                    impact_parameters_back.push_back(b_back);
                    impact_parameters_mean.push_back((double)(b+b_back)/2);
                    merge_s0plate_cand.push_back(end_trk->GetSegmentFirst()->Plate());
                    merge_s0id_cand.push_back(end_trk->GetSegmentFirst()->ID());
                    candidate_tracks.Add((EdbTrackP*)end_trk);
                }
                
            }

            tr_grid_s1.Clear();
    }

    sorted_IPs = impact_parameters_mean;
    sorted_IPs_for = impact_parameters;
    bubbleSort_NEW(sorted_IPs, sorted_IPs_for); //sort by (b+b_back)/2 
    float current_b=0, current_b_back=0, current_b_mean=0;
    int pos=0;
    //cout << " after b calc sorted IP size " << sorted_IPs.size() <<  endl;

    if (sorted_IPs.size()>0) {  //if at least one candidate was found
        for (int icand=0; icand<sorted_IPs.size(); icand++) {
            current_b = sorted_IPs_for[icand];
            pos = findValuePosition(impact_parameters, current_b);
            EdbTrackP *to_merge_trk = (EdbTrackP*)candidate_tracks.At(pos);
            if (to_merge_trk==NULL) continue;
            current_b_back = impact_parameters_back[pos];
            //std::string tr_string = Form("%i_%i", merge_s0plate_cand[pos], merge_s0id_cand[pos]);
            if (to_merge_trk->Flag()!=300 && current_b<MAX_B && current_b_back<MAX_B) { //if the candidate was not used

                MERGED += 1;
                to_merge_trk->SetFlag(300);
                added_segs = added_segs + to_merge_trk->N();
                //used_tracks[merge_s0plate_cand[pos]].push_back(tr_string);
                int current_length = segments_new->GetEntries();
                int back_position = to_merge_trk->N()-1;
                for (int k = 0; k<to_merge_trk->N(); k++) {
                    new((*segments_new)[k+current_length]) EdbSegP(*to_merge_trk->GetSegment(back_position));
                    new((*fitted_segments_new)[k+current_length]) EdbSegP(*to_merge_trk->GetSegmentF(back_position));
                    back_position = back_position - 1;
                }
                candidate_tracks.Clear();
                return to_merge_trk;

            } 
            else continue; //go to the next best candidate
        }
    }

    return NULL;


}



int connect_tracks() {

    int S0 = 1, SL = 2;
    int start_s2_position = 0;
    int MC = 0;
    if (IDBRICK < 10) MC = 1;

    TString file_name = Form("b%06i.0.%i.%i.trk_trasl.root", IDBRICK, S0, SL);
    if (MC == 1) file_name = Form("b%06i.0.%i.%i.trk.root", IDBRICK, S0, SL);

    TFile *trkfile = TFile::Open(file_name, "READ");
    TTree* tracks = (TTree*) trkfile->Get("tracks");

    // End of First Section Cuts
    int PLATE_MIN=1, PLATE_MAX=31, NSEG_MIN=3;

    // Start of Second Section 
    int PLATE_MAX_S2 = 40, thr = 3000, START_PLATE_S2_ALT=31;

    // Merge Cuts
    float DT_MAX = 0.1, B_MAX=100;

    // Read Tracks (need XYZ s0 and sL and sfTX sfTY)

    EdbSegP *trk=0;
    TClonesArray *segments = new TClonesArray("EdbSegP");
    TClonesArray *fitted_segments = new TClonesArray("EdbSegP");
    
	int nseg=0, trid=0;
    tracks->SetBranchAddress("trid",&trid);
    tracks->SetBranchAddress("nseg",&nseg);
    tracks->SetBranchAddress("t.",&trk);

    tracks->SetBranchAddress("s", &segments);
    tracks->SetBranchAddress("sf", &fitted_segments);
    tracks->SetBranchAddress("nseg", &nseg);

    int N = tracks->GetEntries();

    cout << " --- Filling arrTRK with tracks --- " << endl;
    int n1=0, n2=0;

    for (int itrk=0; itrk<N; itrk++) {

        tracks->GetEntry(itrk);
        EdbSegP* seg0 = (EdbSegP*)segments->At(0);
        if (seg0->Plate()==START_PLATE_S2_ALT) start_s2_position = itrk;
        if (nseg<NSEG_MIN) continue;

        EdbSegP* segL = (EdbSegP*)segments->At(nseg-1);

        EdbTrackP *temptrack = new EdbTrackP(); 
		//fill temptrack with segments
        
		for (int k = 0; k<nseg; k++) {

		  temptrack->SetID(trk->ID());
		  EdbSegP* seg = (EdbSegP*)segments->At(k);
		  EdbSegP *segf = (EdbSegP*)fitted_segments->At(k);
		  //seg->SetDZ(300);
		  //segf->SetDZ(300);
		  temptrack->AddSegment(new EdbSegP(*((EdbSegP*)(seg))));
		  temptrack->AddSegmentF(new EdbSegP(*((EdbSegP*)(segf))));
		  temptrack->SetSegmentsTrack(temptrack->ID()); //track segments association
		  temptrack->SetCounters();	

        }

        if (seg0->Plate()<PLATE_MAX_S2 && seg0->Plate()>=START_PLATE_S2_ALT) { arrTRK->Add(temptrack); n2 = n2 +1; }
        else if (segL->Plate()>PLMIN && segL->Plate()<=PLATE_MAX) { arrTRK->Add(temptrack); n1 = n1 + 1;}

    }

    cout << " --- Filling Z LAYER --- " << endl;
    double z_value=0;
    // Fill Z Layer (to be improved)
    for (int i=2; i<=START_PLATE_S2; i++) {
        z_value = 0;
        for (int itrk=0; itrk<start_s2_position; itrk++){
                tracks->GetEntry(itrk);
		        //fill temptrack with segments
		        for (int k = 0; k<nseg; k++) {
		            EdbSegP* seg = (EdbSegP*)segments->At(k);
                    if (seg->Plate()==i) { z_value = seg->Z(); break;}
                }
                if (z_value != 0) break;
        }
        Z_LAYER[i] = z_value;
        //cout << " i " << i << endl;
    }
    for (int i=START_PLATE_S2+1; i<=PLMAX; i++) {
        z_value = 0;
        for (int itrk=start_s2_position; itrk<N; itrk++){
                tracks->GetEntry(itrk);
		        //fill temptrack with segments
		        for (int k = 0; k<nseg; k++) {
		            EdbSegP* seg = (EdbSegP*)segments->At(k);
                    if (seg->Plate()==i) { z_value = seg->Z(); break;}
                }
                if (z_value != 0) break;
        }
        Z_LAYER[i] = z_value;
        //cout << " i " << i << endl;
    }
    for (int i=1; i<=PLMAX;i++) {
        cout << " ZLAYER " << i << " " << Z_LAYER[i] << endl;
    }
    cout << endl;
    

    cout << " arrTRK entries " << arrTRK->GetEntries() << endl; 
    cout << " S2 tracks: " << n2 << ", S1 tracks: " << n1 << endl;

    cout << " --- Filling Cells --- " << endl;
    FillTracksCells(*arrTRK);
    cout << " --- Cells Ready --- " << endl;

    TObjArray tr_grid_s1, tr_grid_s2;
    EdbTrackP *start_trk, *end_trk, *merged_trk;
    EdbSegP *start_seg, *end_seg;

    int N_trks_S2 = gridtr_S2[1].SelectObjects(tr_grid_s2);
    int s0_plate_s2=0, s0_id_s2=0, n_trk_cell=0;

    float xy[2] = {0,0};
    float r = 2., b=1000;
    
    cout << " Found " << N_trks_S2 << " trks " << endl;

    TString merge_file_name = Form("b%06i.0.%i.%i.trk_merged.root", IDBRICK, S0, SL);
    if (EVERBOSE==100||EVERBOSE==101) merge_file_name = Form("b%06i.0.%i.%i.trk_merged_EVERBOSE.root", IDBRICK, S0, SL);
    TFile *merge_file = TFile::Open(merge_file_name, "RECREATE");
    TTree *tracks_new = new TTree();

    EdbSegP *trk_new=0;
    TClonesArray *segments_new = new TClonesArray("EdbSegP");
    TClonesArray *fitted_segments_new = new TClonesArray("EdbSegP");

    TClonesArray *segments_new_ordered = new TClonesArray("EdbSegP");
    TClonesArray *fitted_segments_new_ordered = new TClonesArray("EdbSegP");

    int nseg_new=0, trid_new=0, npl_new=0;
    tracks_new->Branch("trid", &trid_new, "trid/I"); 
    tracks_new->Branch("nseg", &nseg_new, "nseg/I"); 
    tracks_new->Branch("npl", &npl_new, "npl/I");
    tracks_new->Branch("t.", "EdbSegP", &trk_new); 

    tracks_new->Branch("s", &segments_new_ordered); 
    tracks_new->Branch("sf", &fitted_segments_new_ordered); 

    TStopwatch t;
    t.Start();

    //cout << " here " << endl;

    for (int iplS2=1; iplS2<PLATE_MAX_S2-START_PLATE_S2; iplS2++) {
        int N_trks_S2 = gridtr_S2[iplS2].SelectObjects(tr_grid_s2);
        cout << " Starting plate " << iplS2 << " with " << N_trks_S2 << " tracks " <<  endl;

        for (int itrk2=0; itrk2<N_trks_S2; itrk2++) {

            EdbTrackP *start_trk = (EdbTrackP*)tr_grid_s2.At(itrk2);
            //if (itrk2>100) break;
            if (EVERBOSE==100 && (start_trk->GetSegmentFirst()->ID()!=DEBUG_S0_ID || start_trk->GetSegmentFirst()->Plate()!=DEBUG_S0_PLATE)) continue;
            int backpos = start_trk->N()-1;
            // Add segments of S2 track to arrays to save later
            for (int iseg=0; iseg<start_trk->N(); iseg++) {
                new((*segments_new)[iseg]) EdbSegP(*start_trk->GetSegment(backpos));
                new((*fitted_segments_new)[iseg]) EdbSegP(*start_trk->GetSegment(backpos));
                backpos = backpos - 1; 
            }

            nseg_new = start_trk->N();
            int added_segs=0;

            // Look for pieces to merge in S1
            EdbTrackP* to_merge_trk = FindClosestCandidate(3+iplS2, start_trk, segments_new, fitted_segments_new, B_MAX, added_segs); //+iplS2 è un modo per far sì che cerchi sempre almeno fino al piatto 26
            EdbTrackP* ausiliary=NULL;
            //cout << " FindClose 1 " << endl;
            if (EVERBOSE == 100) cout << " added segs after first search " << added_segs << endl;
            if (EVERBOSE==100 && start_trk->GetSegmentFirst()->ID()==DEBUG_S0_ID && start_trk->GetSegmentFirst()->Plate()==DEBUG_S0_PLATE) {
                if (to_merge_trk) cout << " Found Candidate " << endl;
                else cout << " No Candidate Found! " << endl;
            }
            if (to_merge_trk!=NULL) {
                while(to_merge_trk!=NULL) {
                    //cout << " entered with to merge_trk plate " << to_merge_trk->GetSegmentFirst()->Plate() <<  endl;
                    ausiliary = FindClosestCandidate(4, to_merge_trk, segments_new, fitted_segments_new, B_MAX, added_segs);
                    //cout << " FindClose 2 " << endl;
                    if (EVERBOSE == 100) cout << " added segs after second search " << added_segs << endl;
                    if (ausiliary==NULL) { break; }
                    to_merge_trk->Clear();
                    for (int i=0; i<ausiliary->N(); i++) {
                        EdbSegP *seg = ausiliary->GetSegment(i);
                        EdbSegP* seg_new = new EdbSegP();
                        seg_new->SetX(seg->X());
                        seg_new->SetY(seg->Y());
                        seg_new->SetZ(seg->Z());
                        seg_new->SetTX(seg->TX());
                        seg_new->SetTY(seg->TY());
                        seg_new->SetW(seg->W());
                        seg_new->SetFlag(seg->Flag());
                        seg_new->SetID(seg->ID());
                        seg_new->SetPlate(seg->Plate());
                        to_merge_trk->AddSegment(seg_new);

                        EdbSegP *segf = ausiliary->GetSegmentF(i);
                        EdbSegP* segf_new = new EdbSegP();
                        segf_new->SetX(segf->X());
                        segf_new->SetY(segf->Y());
                        segf_new->SetZ(segf->Z());
                        segf_new->SetTX(segf->TX());
                        segf_new->SetTY(segf->TY());
                        segf_new->SetW(segf->W());
                        segf_new->SetFlag(segf->Flag());
                        segf_new->SetID(segf->ID());
                        segf_new->SetPlate(segf->Plate());
                        to_merge_trk->AddSegmentF(segf_new);
                    }
                }       
            }


            nseg_new = nseg_new + added_segs -1; //Added S1 segments
            if (EVERBOSE==100)  cout << " nseg_new " << nseg_new << " added segs " << added_segs << endl;
            
            if(to_merge_trk==NULL) to_merge_trk = start_trk;
            trid_new = to_merge_trk->ID(); //trid is first track ID
            int final_plate = start_trk->GetSegmentLast()->Plate();
            int start_plate = to_merge_trk->GetSegmentFirst()->Plate();
            npl_new = final_plate - start_plate;

            //EdbTrackP *new_trk = new EdbTrackP();
            int counter = 0;
            for (int i=segments_new->GetEntries()-1; i>0; i--) { //do it backwards to start from S1
                /*new_trk->SetID(trid_new);
                new_trk->AddSegment((EdbSegP*)segments_new->At(i));
                new_trk->GetSegment(new_trk->N()-1)->SetDZ(300);
                new_trk->SetSegmentsTrack(new_trk->ID()); //track segments association
                new_trk->SetCounters();	*/
                EdbSegP *s = (EdbSegP*)segments_new->At(i);
                EdbSegP *sf = (EdbSegP*)fitted_segments_new->At(i);
                new((*segments_new_ordered)[counter]) EdbSegP(*s);
                new((*fitted_segments_new_ordered)[counter]) EdbSegP(*sf);
                counter = counter + 1;
            }
            //new_trk->SetCounters();
            
            trk_new = (EdbSegP*)segments_new->At(segments_new->GetEntries()-1);
            //cout << " hello it's me, new_trk " << new_trk->Theta() << " " << new_trk->X() << " npl " << new_trk->Npl() << " " << new_trk->MCEvt() <<endl;
            //cout << " hello it's me, trk_new " << trk_new->Theta() << " " << trk_new->X() << endl;
            //trk_new->ForceCOV(new_trk->COV());    //leads to segmentation violation maybe recalculate COV?
            trk_new->SetVid( 0, tracks_new->GetEntries() ); 

            tracks_new->Fill();

            segments_new->Clear();
            fitted_segments_new->Clear();
            segments_new_ordered->Clear();
            fitted_segments_new_ordered->Clear();
            //candidate_tracks.Clear();
            //new_trk->Clear();

            if (itrk2%1000==0) { cout << " Completed " << 100.*itrk2/N_trks_S2 << " %, Iteration Time: " <<  t.RealTime() << " s " << endl; t.Reset(); t.Start();}  
        }

        cout << " Completed " << iplS2 << " plates in S2 " << endl;
        tr_grid_s2.Clear();

    }
    
    

    cout << " MERGED " << MERGED << endl;

    //add not used S1 tracks (with flag != 300) 
    for (int ipl=1; ipl<=START_PLATE_S2; ipl++) {
        int N_trk_S1 = gridtr_S1[ipl].SelectObjects(tr_grid_s1);
        for (int itrk1=0; itrk1<N_trk_S1; itrk1++) {
            EdbTrackP* s1_trk = (EdbTrackP*)tr_grid_s1.At(itrk1);
            if (s1_trk->Flag()==300) continue;
            for (int k = 0; k<s1_trk->N(); k++) {
                    new((*segments_new_ordered)[k]) EdbSegP(*s1_trk->GetSegment(k));
                    new((*fitted_segments_new_ordered)[k]) EdbSegP(*s1_trk->GetSegmentF(k));
            }

            nseg_new = s1_trk->N();
            npl_new = s1_trk->Npl();
            trid_new = s1_trk->ID(); //trid is first track ID

            //trk_new->Copy(*s1_trk);
            trk_new = s1_trk->GetSegment(0);
            trk_new->SetVid(0, tracks_new->GetEntries() ); 
            tracks_new->Fill();

            segments_new_ordered->Clear();
            fitted_segments_new_ordered->Clear();
        }
        tr_grid_s1.Clear();
    }


    merge_file->cd();
    tracks_new->Write("tracks");
    merge_file->Close();
    
    
    return 1;
}


void FillTracksCells(TObjArray &arrt){

    const int cellsize=1000;
    const int ncellsX = (int)(xmax-xmin)/cellsize;//35;
    const int ncellsY = (int)(ymax-ymin)/cellsize;//50;
    const int maxpercell = 100000; 
    for(int i=1; i<=START_PLATE_S2; i++) gridtr_S1[i].InitCell(ncellsX, xmin, xmax, ncellsY, ymin, ymax, maxpercell );
    for(int i=1; i<=PLMAX-START_PLATE_S2; i++) { gridtr_S2[i].InitCell(ncellsX, xmin, xmax, ncellsY, ymin, ymax, maxpercell ); }

    
    
    int ntr = arrt.GetEntries();
    for(int itr=0; itr<ntr; itr++){
        EdbTrackP *t = (EdbTrackP*)(arrt.At(itr));
        
        float x_tr = t->GetSegmentFirst()->X();
        float y_tr = t->GetSegmentFirst()->Y();

        float x_tr_end = t->GetSegmentLast()->X();
        float y_tr_end = t->GetSegmentLast()->Y();

        int plate_tr = t->GetSegmentFirst()->Plate();
        int plate_tr_end = t->GetSegmentLast()->Plate();
        
        float theta_tr = t->Theta();
        
        if (plate_tr > START_PLATE_S2) { gridtr_S2[plate_tr-START_PLATE_S2].AddObject( x_tr, y_tr, (TObject*)t ); } //organizzo le tracce in S2 in base al piatto di inizio
        else gridtr_S1[plate_tr_end].AddObject(x_tr_end, y_tr_end, (TObject*)t);  // tracce in S1 in base al piatto in cui finiscono
       
    }

}

void bubbleSort(std::vector<double>& v) {
    int n = v.size();
    for (int i = 0; i < n-1; i++) {
        for (int j = 0; j < n-i-1; j++) {
            if (v[j] > v[j+1]) {
                std::swap(v[j], v[j+1]);
            }
        }
    }
}


int findValuePosition(const std::vector<double>& myVector, double valueToFind) {
  for (int i = 0; i < myVector.size(); i++) {
    if (myVector[i] == valueToFind) {
      return i;
    }
  }

  return -1; // return -1 if value is not found
}


void bubbleSort_NEW(std::vector<double>& v, std::vector<double>& v2) {
    int n = v.size();
    for (int i = 0; i < n-1; i++) {
        for (int j = 0; j < n-i-1; j++) {
            if (v[j] > v[j+1]) {
                std::swap(v[j], v[j+1]);
                std::swap(v2[j], v2[j+1]);
            }
        }
    }
}
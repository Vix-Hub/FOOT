#define IDBRICK 2
#define EVERBOSE -100
#define NSTACKS 7
#define DELETE_TEMP_FILE 1
#define STOP_AT_FIRST_MERGE 0
#define BRAGGPLATE 26
#define IS_SECOND_STEP 0
#define LASTPLATEMIN -999 //27 for S1-S2 pieces
#define FIRSTPLATEMAX 30 //39 for S1-S2 pieces

// Sections to Merge, from 1 to 7
const int S0 = 1;
const int SL = 2;
const int LASTLAYER[NSTACKS+1]={1,30,66,76,83,90,110,120}; //esposizione Oxy@200MeV/n 2019
//int LASTLAYER[N_STACKS+1]={1,30,66,76,83,90,120,140}; //esposizione Oxy@400MeV/n 2019

// Cuts to Apply to Tracks
const int NSEG_MIN = 2;

// Need Multiple Trees to save all tracks
const int N_TREES = 20;

// Merge Cuts
float B_MAX = 100;
float DT_MAX = 0.07;

int NPLATES_S1 = 4; // number of plates in which to look for candidates in S1
int NPLATES_S2 = 8; // number of plates in which to look for candidates in S2
int NPLATES_S3 = 4; // number of plates in which to look for candidates in S3

// Debugging
const int DEBUG_S0_PLATE=31; //31
const int DEBUG_S0_ID=1767; //91690
const int DEBUG_S0_PLATE_S1=1;
const int DEBUG_S0_ID_S1 = 26102;
const int DEBUG_SL_PLATE_S1 = 30;

// X-Y Cut
const float xmin = 1000;//0;
const float xmax = 101000;//125000;
const float ymin = 1000;//0;
const float ymax = 101000;//100000;

const int PLMIN=0;  // First Brick Plate
const int PLMAX=LASTLAYER[SL];  // Last Brick Plate

// Grid, Array Definitions
EdbCell2 gridtr_ALL[1+PLMAX]; //first cell left left empty
EdbPVRec *ali = new EdbPVRec();
TObjArray *arrTRK = new TObjArray();   // original tracks
int MERGED=0, MERGED_FIRST_STEP=0;
float Z_LAYER[PLMAX+1]={0};
TTree *tracks_new = new TTree();
TTree *tracks_new2 = new TTree();
TNtuple *merge_offsets = 0;

// Function Declarations
void bubbleSort(std::vector<double>& v);
void checkpatterns();
int findValuePosition(const std::vector<double>& myVector, double valueToFind);
void FillTracksCells(TObjArray &arrt);
double CalcDist(float x1, float y1, float z1, float x2, float y2, float z2, float tx1, float ty1);
int IsElementInVector(const std::vector<string>& myVector, string value);
void bubbleSort_NEW(std::vector<double>& v, std::vector<double>& v2);
void FillZ_LAYER_SET();
double CalcDist(float x1, float y1, float z1, float x2, float y2, float z2, float tx1, float ty1);
int IsElementInVector(const std::vector<string>& myVector, string value);
EdbTrackP* FindClosestCandidate(int nplates, EdbTrackP* start_trk, TClonesArray *segments_new, TClonesArray *fitted_segments_new, const double MAX_B, int &added_segs, const double MAX_DT);
void MergeTrees(TString tempfilename, bool delete_tempfile);
void merge_trees(const char* input_file_name, const char* output_file_name, int num_trees_per_step);
int getIndex(int num);

int connect_tracks_new() {

    int MC = 0;
    if (IDBRICK < 10) MC = 1;

    TString file_name = Form("b%06i.0.%i.%i.trk.root", IDBRICK, S0, SL);
    if (MC == 1) file_name = Form("b%06i.0.%i.%i.trk.root", IDBRICK, S0, SL);

    if (IS_SECOND_STEP) {
        B_MAX = 150; DT_MAX = 0.1;
        file_name = Form("b%06i.0.%i.%i.trk_long.root", IDBRICK, S0, SL);  
    }


    TFile *trkfile = TFile::Open(file_name, "READ");
    TTree* tracks = (TTree*) trkfile->Get("tracks");
    int N_TOT = tracks->GetEntries();
    
    // Read Tracks 
    EdbSegP *trk=0;
    TClonesArray *segments = new TClonesArray("EdbSegP");
    TClonesArray *fitted_segments = new TClonesArray("EdbSegP");
    
	int nseg=0, trid=0;
    tracks->SetBranchAddress("trid",&trid);
    tracks->SetBranchAddress("nseg",&nseg);
    tracks->SetBranchAddress("t.",&trk);
    tracks->SetBranchAddress("s", &segments);
    tracks->SetBranchAddress("sf", &fitted_segments);

    int N = tracks->GetEntries();
    cout << " --- Filling arrTRK with tracks --- " << endl;
    int n1=0, n2=0;

    for (int itrk=0; itrk<N; itrk++) {

        tracks->GetEntry(itrk);
        EdbSegP* seg0 = (EdbSegP*)segments->At(0);
        if (nseg<NSEG_MIN || seg0->Plate()>PLMAX) continue;
        
        EdbSegP* segL = (EdbSegP*)segments->At(nseg-1);
        EdbTrackP *temptrack = new EdbTrackP(); //fill temptrack with segments
        
		for (int k = 0; k<nseg; k++) {
		  temptrack->SetID(trk->ID());
		  EdbSegP* seg = (EdbSegP*)segments->At(k);
		  EdbSegP *segf = (EdbSegP*)fitted_segments->At(k);
		  temptrack->AddSegment(new EdbSegP(*((EdbSegP*)(seg))));
		  temptrack->AddSegmentF(new EdbSegP(*((EdbSegP*)(segf))));
		  temptrack->SetSegmentsTrack(temptrack->ID()); //track segments association
		  temptrack->SetCounters();	

        }
        arrTRK->Add(temptrack);
    }

    cout << " arrTRK Entries: " << arrTRK->GetEntries() << endl;

    cout << " --- Filling Z LAYER --- " << endl;
    checkpatterns();
    FillZ_LAYER_SET();
    cout << endl;

    cout << " --- Filling Cells --- " << endl;
    FillTracksCells(*arrTRK);

    cout << " --- Start Merging --- " << endl;

    TObjArray tr_grid_starting;

    float xy[2] = {0,0};
    float r = 2000., b=1000;
    int counter = 0;

    TString merge_file_name_temp = Form("b%06i.0.%i.%i.trk_long_temp.root", IDBRICK, S0, SL);
    if (IS_SECOND_STEP) merge_file_name_temp = Form("b%06i.0.%i.%i.trk_long_temp_2.root", IDBRICK, S0, SL);
    if (EVERBOSE==100||EVERBOSE==101) merge_file_name_temp = Form("b%06i.0.%i.%i.trk_long_EVERBOSE.root", IDBRICK, S0, SL);
    
    //accessing section offsets
    float Xoff=0, Yoff=0, TXoff=0, TYoff=0;

    TFile *offsets_file = TFile::Open(Form("%i_S%i_S%i_offsets_all.root", IDBRICK, S0, SL), "READ");
    merge_offsets = (TNtuple*)offsets_file->Get("merge_offsets");
    merge_offsets->SetBranchAddress("Xoff", &Xoff);
    merge_offsets->SetBranchAddress("Yoff", &Yoff);
    merge_offsets->SetBranchAddress("TXoff", &TXoff);
    merge_offsets->SetBranchAddress("TYoff", &TYoff);
    
    TFile *merge_file = TFile::Open(merge_file_name_temp, "RECREATE");

    EdbSegP *trk_new=0;
    TClonesArray *segments_new = new TClonesArray("EdbSegP");
    TClonesArray *fitted_segments_new = new TClonesArray("EdbSegP");
    int nseg_new=0, trid_new=0, npl_new=0;
    // Create Tree
    tracks_new->Branch("trid", &trid_new, "trid/I"); 
    tracks_new->Branch("nseg", &nseg_new, "nseg/I"); 
    tracks_new->Branch("npl", &npl_new, "npl/I");
    tracks_new->Branch("t.", "EdbSegP", &trk_new); 

    tracks_new->Branch("s", &segments_new); 
    tracks_new->Branch("sf", &fitted_segments_new); 


    TStopwatch t, total_time;
    t.Start();
    total_time.Start();
    int ntree = 0;

    for (int iipl=1; iipl<=PLMAX; iipl++) {

        if (EVERBOSE==50 && iipl>5) break; 
        int N_trks_starting = gridtr_ALL[iipl].SelectObjects(tr_grid_starting);
        cout << " Starting plate " << iipl << " with " << N_trks_starting << " tracks " <<  endl;


        for (int itrk=0; itrk<N_trks_starting; itrk++) {

            EdbTrackP *start_trk = (EdbTrackP*)tr_grid_starting.At(itrk);
            if (EVERBOSE==100 && ( (start_trk->GetSegmentFirst()->ID()!=DEBUG_S0_ID || start_trk->GetSegmentFirst()->Plate()!=DEBUG_S0_PLATE) && (start_trk->GetSegmentFirst()->Plate()!=DEBUG_S0_PLATE_S1 || start_trk->GetSegmentFirst()->Plate()!=DEBUG_S0_ID_S1))) continue;
            if (start_trk->Flag()==300) continue; //avoid merging same track more than once
            int start_plate = start_trk->GetSegmentFirst()->Plate();
            int last_plate = start_trk->GetSegmentLast()->Plate();
            int search = 1; //used to limit search in some cases (for example if I only want S1-S2)
            //add current segments to new track to save
            for (int iseg=0; iseg<start_trk->N(); iseg++) {
                new((*segments_new)[iseg]) EdbSegP(*start_trk->GetSegment(iseg));
                new((*fitted_segments_new)[iseg]) EdbSegP(*start_trk->GetSegment(iseg));
            }
            nseg_new = start_trk->N();
            int added_segs = 0;
	        int start_nplates = NPLATES_S1;

            EdbTrackP *to_merge_trk = NULL;
            if (start_plate>FIRSTPLATEMAX || last_plate<LASTPLATEMIN) search = 0;
	        if (last_plate>BRAGGPLATE && last_plate<=LASTLAYER[1]) start_nplates += NPLATES_S2;
            if (search) EdbTrackP* to_merge_trk = FindClosestCandidate(start_nplates, start_trk, segments_new, fitted_segments_new, B_MAX, added_segs, DT_MAX); //+iplS2 è un modo per far sì che cerchi sempre almeno fino al piatto 26

            EdbTrackP* ausiliary=NULL;
            if (EVERBOSE == 100 ) cout << " added segs after first search " << added_segs << endl;
            if (EVERBOSE==100 && start_trk->GetSegmentFirst()->ID()==DEBUG_S0_ID && start_trk->GetSegmentFirst()->Plate()==DEBUG_S0_PLATE) {
                if (to_merge_trk) cout << " Found Candidate " << endl;
                else cout << " No Candidate Found! " << endl;
            }
            if (to_merge_trk!=NULL) {
                MERGED_FIRST_STEP += 1;
                while(to_merge_trk!=NULL) {
                    if (STOP_AT_FIRST_MERGE) break;
                    //cout << " entered with to merge_trk plate " << to_merge_trk->GetSegmentFirst()->Plate() <<  endl;
                    int nplates = NPLATES_S1;
                    int last_plate = to_merge_trk->GetSegmentLast()->Plate();
                    if (last_plate>BRAGGPLATE && last_plate<=LASTLAYER[1]) nplates += NPLATES_S2;
                    else if (last_plate>LASTLAYER[1] && last_plate<=LASTLAYER[2]) nplates = NPLATES_S1 + NPLATES_S2;
                    else if (last_plate>LASTLAYER[2]) nplates = NPLATES_S3; 
                    
                    ausiliary = FindClosestCandidate(nplates, to_merge_trk, segments_new, fitted_segments_new, B_MAX, added_segs, DT_MAX);
                    //cout << " FindClose 2 " << endl;
                    if (EVERBOSE == 100) cout << " added segs after second search " << added_segs << endl;
                    if (ausiliary==NULL) { break; }
                    to_merge_trk->Clear();
                    for (int i=0; i<ausiliary->N(); i++) {
                        EdbSegP *seg = ausiliary->GetSegment(i);
                        EdbSegP* seg_new = new EdbSegP();
                        seg_new->SetX(seg->X()); seg_new->SetY(seg->Y()); seg_new->SetZ(seg->Z());
                        seg_new->SetTX(seg->TX()); seg_new->SetTY(seg->TY()); seg_new->SetW(seg->W());
                        seg_new->SetFlag(seg->Flag()); seg_new->SetID(seg->ID()); seg_new->SetPlate(seg->Plate());
                        seg_new->SetVid(seg->Vid(0), seg->Vid(1));
                        seg_new->SetAid(seg->Aid(0), seg->Aid(1));
                        to_merge_trk->AddSegment(seg_new);

                        EdbSegP *segf = ausiliary->GetSegmentF(i);
                        EdbSegP* segf_new = new EdbSegP();
                        segf_new->SetX(segf->X()); segf_new->SetY(segf->Y()); segf_new->SetZ(segf->Z());
                        segf_new->SetTX(segf->TX()); segf_new->SetTY(segf->TY()); segf_new->SetW(segf->W());
                        segf_new->SetFlag(segf->Flag()); segf_new->SetID(segf->ID()); segf_new->SetPlate(segf->Plate());
                        segf_new->SetVid(segf->Vid(0), segf->Vid(1));
                        segf_new->SetAid(segf->Aid(0), segf->Aid(1));
                        to_merge_trk->AddSegmentF(segf_new);
                    }
                }       
            }


            nseg_new = nseg_new + added_segs ; //Added S1 segments
            if (EVERBOSE==100)  cout << " nseg_new " << nseg_new << " added segs " << added_segs << endl;
            
            if(to_merge_trk==NULL) to_merge_trk = start_trk;
            trid_new = to_merge_trk->ID(); //trid is first track ID
            EdbSegP *final_seg = (EdbSegP*)segments_new->At(segments_new->GetEntries()-1);
            int final_plate = final_seg->Plate();
            npl_new = final_plate - start_plate;
            
            trk_new = (EdbSegP*)segments_new->At(0);
            //cout << " hello it's me, new_trk " << new_trk->Theta() << " " << new_trk->X() << " npl " << new_trk->Npl() << " " << new_trk->MCEvt() <<endl;
            //cout << " hello it's me, trk_new " << trk_new->Theta() << " " << trk_new->X() << endl;
            //trk_new->ForceCOV(new_trk->COV());    //leads to segmentation violation maybe recalculate COV?
            trk_new->SetVid( 0, tracks_new->GetEntries() ); 
            tracks_new->Fill();
            

            segments_new->Clear();
            fitted_segments_new->Clear();

            if (itrk%1000==0) { cout << " Completed " << 100.*itrk/N_trks_starting << " %, Iteration Time: " <<  t.RealTime() << " s " << endl; t.Reset(); t.Start();}  
        }

        counter += N_trks_starting;

        if (counter>N_TOT/N_TREES) {
            counter = 0;
            ntree += 1;

            merge_file->cd();
            TString name = Form("tracks%d", ntree);
            tracks_new->Write(name);

            tracks_new->Delete();
            TTree *tracks_new2 = new TTree();

            tracks_new2->Branch("trid", &trid_new, "trid/I"); 
            tracks_new2->Branch("nseg", &nseg_new, "nseg/I"); 
            tracks_new2->Branch("npl", &npl_new, "npl/I");
            tracks_new2->Branch("t.", "EdbSegP", &trk_new); 

            tracks_new2->Branch("s", &segments_new); 
            tracks_new2->Branch("sf", &fitted_segments_new);
            tracks_new = tracks_new2;
            
        }

        cout << " Completed Plate " << iipl << endl;
        tr_grid_starting.Clear();

    }

    merge_file->cd();
    tracks_new->Write("tracks");
    merge_file->Close();

    cout << " MERGED " << MERGED << endl;
    cout << " Total Time " << total_time.RealTime() << " s " << endl; total_time.Reset(); total_time.Start();

    cout << " --- Merging the tracks trees --- " << endl;
    TString merge_file_name = Form("b%06i.0.%i.%i.trk_long.root", IDBRICK, S0, SL);
    if (IS_SECOND_STEP) merge_file_name = Form("b%06i.0.%i.%i.trk_long_2.root", IDBRICK, S0, SL);
    merge_trees(merge_file_name_temp, merge_file_name, 2);

    cout << " Total Time " << total_time.RealTime() << " s " << endl;
    if (DELETE_TEMP_FILE) std::remove(merge_file_name_temp);

    return 1;
}



void FillZ_LAYER_SET(){

    TFile *setfile = TFile::Open(Form("b%06d.0.0.0.set.root",IDBRICK));
    EdbScanSet *set = (EdbScanSet*) setfile->Get("set");
    
    int j=0;
    for(int i=0; i<PLMAX; i++) {
        Z_LAYER[i] = set->eB.GetPlate(i)->Z();
        cout << " Layer : " << i << " " << Z_LAYER[i] << endl;
    } 
    
}

void checkpatterns(){
    //check for patterns, if there is one missing add it
    int np = ali->Npatterns();
    
    TFile *setfile = TFile::Open(Form("b%06d.0.0.0.set.root",IDBRICK));
    cout << "I'm checking patterns with set " << "b" << IDBRICK << ".0.0.0.set.root" << "\t np: " << np << "\tRunning on " << PLMAX << " plates" << endl;
    EdbScanSet *set = (EdbScanSet*) setfile->Get("set");
    
    int j=0;
    for(int i=0; i<PLMAX; i++) {
        float zmissingPID = set->eB.GetPlate(i)->Z(); //note, set->eB->GetPlate(i) uses the PID, set->GetPlate(i) uses the number of plate, be careful!
        EdbPattern *pat = new EdbPattern( 0., 0.,zmissingPID);
        cout <<"missing pattern " << i << " now adding it. Z= " << zmissingPID << endl;
        pat->SetID(i);
        pat->SetScanID(0);
        ali->AddPatternAt(pat,i);
        // } //end if
        
    } //end for loop
    
    cout << "New n patters: " << ali->Npatterns() << endl;
    
}


void FillTracksCells(TObjArray &arrt){

    const int cellsize=1000;
    const int ncellsX = (int)(xmax-xmin)/cellsize;//35;
    const int ncellsY = (int)(ymax-ymin)/cellsize;//50;
    const long int maxpercell = 50000; 
    for(int i=1; i<=PLMAX; i++) gridtr_ALL[i].InitCell(ncellsX, xmin, xmax, ncellsY, ymin, ymax, maxpercell );
    
    int ntr = arrt.GetEntries();
    for(int itr=0; itr<ntr; itr++){
        EdbTrackP *t = (EdbTrackP*)(arrt.At(itr));
        
        float x_tr = t->GetSegmentFirst()->X();
        float y_tr = t->GetSegmentFirst()->Y();
        int plate_tr = t->GetSegmentFirst()->Plate();
        
        gridtr_ALL[plate_tr].AddObject(x_tr, y_tr, (TObject*)t);
        if (EVERBOSE==100 && plate_tr==DEBUG_S0_PLATE && t->GetSegmentFirst()->ID()==DEBUG_S0_ID) cout << " Adding S2 Candidate to grid with plate_tr " << plate_tr << endl;
        
    }

}


EdbTrackP* FindClosestCandidate(int nplates, EdbTrackP* start_trk, TClonesArray *segments_new, TClonesArray *fitted_segments_new, const double MAX_B, int &added_segs, const double MAX_DT) {
    
    EdbSegP* start_seg = (EdbSegP*)start_trk->GetSegmentLast(); //use segments (not fitted segments) for plate info
    EdbSegP* start_segf = (EdbSegP*)start_trk->GetSegmentFLast();
    float xy[2] = {0,0};
    if (EVERBOSE==100) cout << " Start Seg Coordinates " << start_seg->X() << " " << start_seg->Y() << " " << start_seg->Z() << " " << start_seg->TX() << " " << start_seg->TY() << " Plate ID " << start_seg->Plate() << " " << start_seg->ID() << endl;
    int s0_plate = start_trk->GetSegmentLast()->Plate();
    int start_Section = getIndex(s0_plate); // Establish to which section the current last segment belongs to
    float b=0, r=2000., b_back=0, dtx=0, dty=0;
    int next_Section = 0;
    
    
    std::vector<double> impact_parameters, sorted_IPs, impact_parameters_back, impact_parameters_mean, sorted_IPs_for;
    std::vector<int> merge_s0plates, merge_s0ids, merge_s0plate_cand, merge_s0id_cand;
    TObjArray candidate_tracks, tr_grid_candidates;

    for (int ipl=s0_plate+1; ipl<=s0_plate+nplates; ipl++) {  //look for candidates in 4 plates

            xy[0] = start_segf->X(); //last segment coordinates
            xy[1] = start_segf->Y();

            if (ipl<PLMIN) continue; //avoid plate < 0
            if (ipl>PLMAX) continue;  // avoid plate >PLMAX
            next_Section = getIndex(ipl); 

            if (next_Section != start_Section) {
                merge_offsets->GetEntry(next_Section-1); //get corresponding offsets
            } else {
                Xoff = 0; Yoff = 0; TXoff = 0; TYoff = 0;
            }

            double z_pos = Z_LAYER[ipl-1];
            //cout << " z_pos " << z_pos << endl;
            xy[0] = xy[0] - (start_segf->Z()-z_pos)*start_segf->TX();  //this can be further away if angles are large
            xy[1] = xy[1] - (start_segf->Z()-z_pos)*start_segf->TY();

            int n_trk_cell = gridtr_ALL[ipl].SelectObjectsC(xy, r, tr_grid_candidates); //tracce in S1
            if(EVERBOSE==100) cout << " Checking Grid corresponding to plate # " << ipl << ", found " << n_trk_cell << " objects " << endl;
            if (EVERBOSE==100) cout << " i was looking with xy " << xy[0] << " " << xy[1] << endl;
            for (int itrk1=0; itrk1<n_trk_cell; itrk1++) {
                EdbTrackP* end_trk = (EdbTrackP*)tr_grid_candidates.At(itrk1);
                EdbSegP* end_seg = (EdbSegP*)end_trk->GetSegmentFirst();
                EdbSegP* end_segf = (EdbSegP*)end_trk->GetSegmentFFirst();

                b_back = CalcDist(end_segf->X()+Xoff, end_segf->Y()+Yoff, end_segf->Z(), start_segf->X(), start_segf->Y(), start_segf->Z(), end_segf->TX()+TXoff, end_segf->TY()+TYoff);
                b = CalcDist(start_segf->X(), start_segf->Y(), start_segf->Z(), end_segf->X()+Xoff, end_segf->Y()+Yoff, end_segf->Z(), start_segf->TX()+TXoff, start_segf->TY()+TYoff);
                
                dtx = - (end_segf->TX()+TXoff) + start_segf->TX();
                dty = - (end_segf->TY()+TYoff) + start_segf->TY();
                if (EVERBOSE==100) cout << " Calculated b with track with first seg " << end_seg->Plate() << " " << end_seg->ID() << " : " << b << ", b_back: " << b_back << endl;
                if (b<MAX_B && TMath::Abs(dtx)<MAX_DT && TMath::Abs(dty)<MAX_DT) {
                    impact_parameters.push_back(b);
                    impact_parameters_back.push_back(b_back);
                    impact_parameters_mean.push_back((double)(b+b_back)/2);
                    merge_s0plate_cand.push_back(end_trk->GetSegmentFirst()->Plate());
                    merge_s0id_cand.push_back(end_trk->GetSegmentFirst()->ID());
                    candidate_tracks.Add((EdbTrackP*)end_trk);
                }
                
            }

            tr_grid_candidates.Clear();
    }

    sorted_IPs = impact_parameters_mean;
    sorted_IPs_for = impact_parameters;
    bubbleSort(sorted_IPs_for);
    //bubbleSort_NEW(sorted_IPs, sorted_IPs_for); //another option is to sort by (b+b_back)/2 
    float current_b=0, current_b_back=0, current_b_mean=0;
    int pos=0;
    //cout << " after b calc sorted IP size " << sorted_IPs.size() <<  endl;

    if (sorted_IPs_for.size()>0) {  //if at least one candidate was found
        for (int icand=0; icand<sorted_IPs_for.size(); icand++) {
            current_b = sorted_IPs_for[icand];
            pos = findValuePosition(impact_parameters, current_b);
            EdbTrackP *to_merge_trk = (EdbTrackP*)candidate_tracks.At(pos);
            if (to_merge_trk==NULL) continue;
            //current_b_back = impact_parameters_back[pos];
            if (to_merge_trk->Flag()!=300 && current_b<MAX_B) { //if the candidate was not used, the cut was added before

                MERGED += 1;
                to_merge_trk->SetFlag(300);
                added_segs = added_segs + to_merge_trk->N();
                int current_length = segments_new->GetEntries();
                for (int k = 0; k<to_merge_trk->N(); k++) {
                    new((*segments_new)[k+current_length]) EdbSegP(*to_merge_trk->GetSegment(k));
                    new((*fitted_segments_new)[k+current_length]) EdbSegP(*to_merge_trk->GetSegmentF(k));
                }
                candidate_tracks.Clear();
                return to_merge_trk;

            } 
            else continue; //go to the next best candidate
        }
    }

    return NULL;


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



void merge_trees(const char* input_file_name, const char* output_file_name, int num_trees_per_step) {
    TFile* input_file = new TFile(input_file_name, "READ");
    if (!input_file || input_file->IsZombie()) {
        cerr << "Error: failed to open input file " << input_file_name << endl;
        return;
    }

    TFile* output_file = new TFile(output_file_name, "RECREATE");
    if (!output_file || output_file->IsZombie()) {
        cerr << "Error: failed to create output file " << output_file_name << endl;
        return;
    }

    // Get the list of trees in the input file
    TList* input_tree_list = input_file->GetListOfKeys();
    input_tree_list->Sort();

    // Merge the trees in steps
    int num_trees = input_tree_list->GetSize();
    int step = 1;

    TList* tree_list = new TList();
    for (int i = 0; i < num_trees; i++) {
        TObject* obj = input_file->Get(input_tree_list->At(i)->GetName());
        if (TClass::GetClass(obj->ClassName())->InheritsFrom(TTree::Class())) {
            if (obj) {
                tree_list->Add(obj);
            } else {
                cerr << "Error: failed to get tree " << obj->GetName() << endl;
            }
        } else {
            cerr << "Warning: object " << obj->GetName() << " is not a TTree" << endl;
        }
    }

    if (tree_list->GetSize() == 0) {
        cerr << "Error: no trees to merge" << endl;
        delete tree_list;
        return;
    }

    TTree* output_tree = TTree::MergeTrees(tree_list);
    if (!output_tree) {
        cerr << "Error: failed to merge trees" << endl;
        delete tree_list;
        return;
    }

    output_tree->SetDirectory(output_file);
    output_file->cd();
    output_tree->Write("tracks");
    delete output_tree;
    delete tree_list;

    output_file->Close();
    input_file->Close();

    delete output_file;
    delete input_file;
}


int getIndex(int num) {
    for (int i = 0; i < NSTACKS; i++) {
        if (LASTLAYER[i] < num && num <= LASTLAYER[i+1]) {
            return i;
        }
    }
    return -1; // num is outside the range of the array
}
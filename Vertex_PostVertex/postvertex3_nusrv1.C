/*
 scp /Users/Giuliana/Desktop/Uni/FOOT/GSI_data/postvertex3.C scanner@nusrv9.na.infn.it:/home/scanner/foot/RECO_MC_TEST_giul/b000003
 
 scp /Users/Giuliana/Desktop/Uni/FOOT/GSI_data/postvertex3.C scanner@nusrv9.na.infn.it:/home/scanner/foot/2019_GSI/GSI2/b000222/
 
 scp /Users/Giuliana/Desktop/Uni/FOOT/GSI_data/postvertex3.C scanner@nusrv9.na.infn.it:/home/scanner/foot/2019_GSI/GSI1/b000111/
 
 scp /Users/Giuliana/Desktop/Uni/FOOT/GSI_data/postvertex3.C scanner@nusrv9.na.infn.it:/home/scanner/foot/2019_GSI/GSI3/b000333/
 
 Usage:
 
 root -l postvertex3.C
 
 scp /Users/Giuliana/Desktop/foot/crosssection_eval.c scanner@nusrv9.na.infn.it:/home/scanner/foot/RECO_MC_TEST_giul/
 
 */
//

//Quando cambio energia devo cambiare:
//BRICKID
//BRAGGPLATE
//int LASTLAYER[N_STACKS+1]

#define BRICKID 2
#define BRAGGPLATE 26 //26 per oxy 200 (GSI 1 e 2), 30 per oxy 400 (GSI 3 e 4)
#define FAST 0.5 //1 prepara solo AnaFake //2 preparetrks // 3 parte da AnaFake e fa il resto // 4 fa un solo vertice (quindi se non si è fatto FAST 1 e FAST 2 non si può fare FAST 3)
#define EVERBOSE 101 //1 per Print EdbCell2 // 2 found beam //3 found dau // 4 Np // 5 merge 2p vtx // 6 remove beam // 7 unione tracce // 8 merge tracks // 9 merge vtx //10 print // 11 PrepareTrk // 12 Cosmici //13 findclosetracks //100 stampa un evento particolare secondo gli ID specificati dopo //101 timing
//100 per un evento particolare // 14 show % // 15 merge checks // 16 grid occupancy checks
#define DEBUG_MCEVT -99 //
#define DEBUG_VTXID  6660 //447 //10464 //6557  //88
#define DEBUG_TRKID -99 //
#define NITROGEN_SEARCH 1//1 // se 1 ricerca O->N+p se 2 ricerca SOLO O->N+p
#define PULIZIA_EXTRA 0 // non funziona
#define DIRECTION 1 // 9 = indietro, 1c = avanti
#define N_STACKS 7
#define DEBUG_S0_PLATE_END 10 //16
#define DEBUG_S0_ID_END 1003510 //563496
#define DEBUG_S0_PLATE_ST 31
#define DEBUG_S0_ID_ST 7388
#define ALLOW_LONGER_GAPS 1 
#define UNISCI_TRACCE_LARGER_CUT 1
int LASTLAYER[N_STACKS+1]={1,30,66,76,83,90,110,120}; //esposizione Oxy@200MeV/n 2019
//int LASTLAYER[N_STACKS+1]={1,30,66,76,83,90,120,140}; //esposizione Oxy@400MeV/n 2019


const int NPLATE = LASTLAYER[N_STACKS];
#define NPLATE_1 121
const int PLMIN=0;//LASTLAYER[0];
const int PLMAX=66; // per tutto brick LASTLAYER[N_STACKS] //conta da 0
const float BEAMTHETA=0.03;
float BEAMTHETA_SMALL=0;
const float COSMICSCUT=0.40; //era 0.45 11 dec
float Z_LAYER[NPLATE_1]={0}; // z of each plate   <- 0 = 0

//const int NMISSINGPL=0; //numero di piatti mancanti // 1 per GSI2
//const int missing_PID[NMISSINGPL+1] = {84}; //piatto mancante GSI2
//const float missing_PID_z[NMISSINGPL+1] = {103658.45}; //piatto mancante GSI2

EdbCell2 gridseg_OXY[BRAGGPLATE+2]; // contiene tutte le tracce identificati come fascio incidente dall'angolo che passano in una cella
EdbCell2 gridtr_OXY[BRAGGPLATE+2]; // contiene tutte le tracce identificati come fascio incidente dall'angolo che iniziano in una cella
EdbCell2 gridtr_DAU[PLMAX+1];   //contiene tutte le tracce identificati come dau INIZIO
EdbCell2 gridtr_DAU_end[PLMAX+1];   //contiene tutte le tracce identificati come dau FINE
EdbCell2 gridvtx[BRAGGPLATE+2];   //contiene tutti i vtx
int mostfrequentevent[500000]={0};
int OXY_FOUND[500000]={0};

int MC=0;

int N_BEAM_FOUND=0, N_BEAM_TRUE=0, N_BEAM_SPLITTED=0, N_BEAM_SPLITTED_OK=0;
int N_DAU_FOUND=0, N_DAU_TRUE=0;
int N_NP_FOUND=0, N_NP_FOUND_OK=0;
int N_NPP_FOUND=0, N_NPP_FOUND_OK=0;
int N_BEAM_REMOVED =0;
int MERGEVTX=0;
int REMOVECOSMIC=0, REMOVECOSMIC_OK=0;
int REM_BEAMDAU=0, REM_BEAMDAU_OK=0;

int COUNT_ATTACHED_VTX=0, COUNT_SPLIT=0, COUNT_SPLIT_OK=0, COUNT_ATTACHED_OK=0, COUNT_ATTACHED_OK_NOBKG=0, COUNT_STACCA_OK=0, THREEP_REMOVE_BKG_OK=0, THREEP_ATTACHED_OK=0, COUNT_3P_ATTACHED_VTX=0, NV_2P=0;
int EXTEND_TRACK=0, EXTEND_TRACK_OK=0, TRACKS_REMOVED=0;

float maximp_beam=0, maximp_dau=0, maximp_Np=0,maximp_unione=0, max_deltatheta_unione=0, Max_acceptable_IP=0, Cut_Theta_cosmic=0, maxbeam_dx=0, maxbeam_dy=0;
int minseg_Np=0;

//min e max dell'area che prendo in considerazione per edbcell2
const float xmin = 35000;//0;
const float xmax = 90000;//125000;
const float ymin = 25000;//0;
const float ymax = 75000;//100000;

//passive material thickness in S1
float THICKNESS = 2000;

//verbose option
int CHECK_OXY = 0;

TObjArray vtxPat[BRAGGPLATE+2];   //vertices grouped by patterns
EdbPVRec *ali = new EdbPVRec();
EdbVertexRec *mygEVR = new EdbVertexRec();
TObjArray *arrTRK=0;   // original tracks
TObjArray *arrVTX=0;   // original vertices
TObjArray *merged_arrVTX=0;   // vertices after 2prong merge
TObjArray *new_varr=0; //final vertices array
TObjArray *final_varr=0; //vertices array after merging
TObjArray *end_tracksS1=0; //vtx tracks at the end of S1 with no candidates


TH1F *h_IP_merge = new TH1F("h_IPmerge","Impact Parameter; Impact Parameter (#mum); events",500,0,500);
TH2F *h_DXY_merge = new TH2F("h_DXY_merge","#DeltaX vs #DeltaY;#DeltaX (#mum); #DeltaY (#mum)",50,-500,500,50,-500,500);
TH2F *h_DTXTY_merge = new TH2F("h_DTXTY_merge","#Delta#theta_{X} vs #Delta#theta_{Y};#Delta#theta_{X} (rad); #Delta#theta_{Y} (rad)",100,-0.3,0.3,100,-0.3,0.3);
TH2F *h_IP_DTh_merge = new TH2F("h_IP_DTh_merge","IP vs #Delta#theta;IP; #Delta#theta",500,0,500,300,0,0.3);

// function declarations
void checkpatterns();
void FillZ_LAYER();
TObjArray* AnalyseFakeVtxs(TObjArray &arrv, EdbVertexRec *vrec);
void CreateTree(TTree *new_vtxtree, TObjArray *varr);
void FillTracksCells(TObjArray &arrt);
void FillVtxPlate(TObjArray &arrv);
void PrepareTRK();
void ReadTreeTracksVTA(TObjArray &arrt, TObjArray &arrv);
TObjArray *FindCloseTracks(TObjArray *varr, EdbVertexRec *vrec);
void FindNitrogen(TObjArray *varr, EdbVertexRec *vrec, float maximp_Np, float maximp_dau);
TObjArray *UnisciVertici(TObjArray *varr, EdbVertexRec *vrec);
TObjArray* UnisciTracce(TObjArray *varr, float maximp, float max_deltatheta);
void UnisciTracceLongerGaps(TObjArray *end_tracksS1, float maximp, float max_deltatheta);
int Get_vtx_plate(float vz);
void MergeTrack(EdbTrackP* firsttrack, EdbTrackP* secondtrack);
EdbVertex *AddTrackToVertex2(EdbVertexRec *vrec, EdbVertex *eVertex, EdbTrackP *eTr, int zpos);
EdbTrackP *AnalyseCloseBeamTracks(EdbVertex *vertex, int ipl, float maximp, float max_dx, float max_dy);
EdbVertex * AddBeamToVertex(EdbVertex *vertex, EdbTrackP *found_cand, EdbVertexRec* vrec);
EdbVertex * SelectOneBeam(EdbVertex *vertex, EdbVertexRec *vrec);
int DaughtersSearch( EdbVertex *vertex, TObjArray &tracks, EdbVertexRec* vrec, float maximp);
int FindPlate(EdbTrackP *cand, int plate);
double CalcDist3D(EdbVertex *vertex, EdbTrackP *track);
double CalcDist3D(EdbVertex *vertex, EdbSegP *seg);
double CalcDist2DX(EdbVertex *vertex, EdbTrackP *track);
double CalcDist2DX(EdbVertex *vertex, EdbSegP *seg);
double CalcDist2DY(EdbVertex *vertex, EdbSegP *seg);
double CalcDist2DY(EdbVertex *vertex, EdbTrackP *track);
int ExtendTrack(EdbTrackP *trkend, TObjArray &tracks, float maximp, float max_deltatheta);
double CalcDist(EdbVertex *vertex, EdbTrackP *track);
double CalcDist(EdbSegP *seg1, EdbSegP *seg2);
double CalcIP(EdbVertex *vertex, int itrk);
double CalcIP(EdbVertex *vertex, EdbTrackP *track);
double CalcIP(EdbVertex *vertex, EdbTrackP *track, int zpos);
double CalcIP(EdbSegP *s, double x, double y, double z);
int AddDauToVertex(EdbVertex *vertex, EdbTrackP *found_cand, EdbVertexRec *vrec);
bool SplitTrackAtZ( EdbTrackP &t, EdbTrackP &t1, EdbTrackP &t2, float vz );
Float_t GetThetaKink(EdbSegP *s1, EdbSegP *s2);
EdbTrackP *FindProton(EdbTrackP *beam, int iseg, float maximp, int minseg);
EdbVertex * MakeNpVertex(EdbTrackP *beam, EdbTrackP *proton, EdbVertexRec *vrec);
void FillVertexCells(TObjArray &arrVTX);
void LoopVertexCells_ToMerge(TObjArray &varr, EdbVertexRec *vrec);
EdbVertex* MergeVertices(TObjArray *varr,EdbVertexRec *vrec);



//---------------------------------------------------------------------


template<typename T> void safe_delete(T*& a) {
    delete a;
    a = NULL;
}



int postvertex3_new_temp()
{
    //GSI 1 (mc) 11 (dati) -> Oxy@200 MeV/n su C target 1mm
    //GSI 2 (mc) 22 (dati) -> Oxy@200 MeV/n su C2H4 target 2mm
    
    if(BRICKID==111){
        MC=0;
        maximp_beam= 60.0; //era 65  11 dec
        maximp_dau = 60.0; //era 95  11 dec
        maximp_Np = 20.0;
        minseg_Np = 10;
        maximp_unione = 70.0; //era 100  11 dec
        max_deltatheta_unione=0.1;
        Max_acceptable_IP=100;//era 150  11 dec
        Cut_Theta_cosmic=0.04; //era 0.03
        BEAMTHETA_SMALL=0.03; //era 0.03
		THICKNESS = 1000;
    }
    else if(BRICKID==1){
        MC=11;
        maximp_beam= 50.0;//60.0;
        maximp_dau = 100.0;//75.0;
        maximp_Np = 10.0;//55.0;
        minseg_Np = 10;
        maximp_unione = 100.0; //55.0;
        max_deltatheta_unione=0.1;
        Max_acceptable_IP=80;
        Cut_Theta_cosmic=0.02;
        BEAMTHETA_SMALL=0.01;
		THICKNESS = 1000;
    }
    else if(BRICKID==222 || BRICKID==444){
        MC=0;
        maximp_beam= 65.0;////era 60  11 dec
        maximp_dau = 100.0;////era 95  11 dec
        maximp_Np = 20.0; ////era 40  11 dec
        minseg_Np = 10;
        maximp_unione =100.0; //65, Era 100 11 dec
        max_deltatheta_unione=0.1;
        Max_acceptable_IP=100;//era 150  11 dec
        Cut_Theta_cosmic=0.02;
        BEAMTHETA_SMALL=0.03;
		THICKNESS = 2000;
    }
    else if(BRICKID==2 || BRICKID==4){
        MC=11;
        maximp_beam= 25.0;//era 60 // era 75 //31 gen 50
        maximp_dau = 100.0; //era 90 //31 gen 65
		maxbeam_dx = 20.;
		maxbeam_dy = 20.;
        maximp_Np = 10.0; //era 60
        minseg_Np = 10;
        maximp_unione = 100.0; //28 mar
        max_deltatheta_unione=0.1;
        Max_acceptable_IP=130;
        Cut_Theta_cosmic=0.01;
        BEAMTHETA_SMALL=0.01;
		THICKNESS = 2000;
    }
    else if(BRICKID==333){
        MC=0;
        maximp_beam= 70.0;
        maximp_dau = 100.0;
        maximp_Np = 35.0;
        minseg_Np = 6;
        maximp_unione = 100.0;
        max_deltatheta_unione=0.1;
        Max_acceptable_IP=150;
        Cut_Theta_cosmic=0.03; //
        BEAMTHETA_SMALL=0.03; //
		THICKNESS = 1000;
    }
    else if(BRICKID==3){
        MC=11;
        maximp_beam= 55.0;//60.0;
        maximp_dau = 60.0;//75.0;
        maximp_Np = 15.0;//55.0;
        minseg_Np = 6;
        maximp_unione = 55.0;
        max_deltatheta_unione=0.1;
        Max_acceptable_IP=80; //135
        Cut_Theta_cosmic=0.02;
        BEAMTHETA_SMALL=0.01;
		THICKNESS = 1000;
    }
    
    if(PLMAX<=BRAGGPLATE) cout << "Attenzione, braggplate > plmax porta a errori! Setta braggplate == plmax" << endl;
    
    arrVTX = new TObjArray();
    merged_arrVTX = new TObjArray();
    new_varr = new TObjArray();
    final_varr = new TObjArray();
    
    TFile *inputfile_vtx;
    inputfile_vtx = TFile::Open("vertices.root","READ"); //vertices_AnaFake
    if (FAST==1 || FAST==100) inputfile_vtx = TFile::Open("vertices_merged.root", "READ");
    if(FAST>1 && FAST!=100) inputfile_vtx = TFile::Open("vertices_AnaFake.root","READ");
    else if(NITROGEN_SEARCH==2) inputfile_vtx = TFile::Open("vertices_improved.root","READ");
    if (inputfile_vtx == NULL) cout<<"ERROR: inputfile_vtx not found"<<endl;
    
    EdbVertexRec *vrec = (EdbVertexRec*) inputfile_vtx->Get("EdbVertexRec");
    vrec->SetPVRec(ali);
    arrVTX = vrec->eVTX;
    arrTRK = new TObjArray();
    
    EdbScanCond *myscancond = new EdbScanCond();
    ali->SetScanCond(myscancond);
    EdbDataProc *dproc = new EdbDataProc();
    TCut tracksel("");//"t.X()>40000&&t.X()<85000 && t.Y()>36000&&t.Y()<62000"); //AREA SEGNALE
    //TCut tracksel("t.X()>26000&&t.X()<97000 && t.Y()>16000&&t.Y()<86000"); //AREA SCANSIONE GSI2 //"nseg>4&&t.X()>40000&&t.X()<85000 && t.Y()>30000&&t.Y()<70000"); //nseg>4 npl>=5 // il taglio nseg>4 torna anche dopo per le tracce di S1. Cerca "TGL4" //TAGLIO CUT FAST
    //if(FAST>0) tracksel.SetTitle("t.X()>57000&&t.X()<68000 && t.Y()>40000&&t.Y()<50000");
    dproc->InitVolume(100, tracksel);
    ali = dproc->PVR();
    checkpatterns(); //new functions to check for missing patterns
    
    ali->FillCell(30,30,0.009,0.009);
    arrTRK = ali->eTracks;
    
    TString output_vtxname;
    output_vtxname=Form("vertices_improved.root"); //vertices_improved
    if(NITROGEN_SEARCH==2) output_vtxname=Form("vertices_improved_Np.root");
    if (FAST==0.5) output_vtxname=Form("vertices_merged.root");
    if(FAST==1) output_vtxname=Form("vertices_AnaFake.root");
    if(FAST>=3) output_vtxname=Form("vertices_improved_fast_%d_new_temp.root", FAST);
    if(FAST==100) output_vtxname=Form("vertices_AnaFake_%d.root", FAST);
    
    cout << "Creating new vertices file: " << output_vtxname << endl;
    TFile *fvtx;
    if(FAST!=2) fvtx = new TFile(output_vtxname,"RECREATE"); //OUTPUT FILE NAME
    cout << "Starting from " << arrVTX->GetEntries() << " vertices" << endl;
    TTree *new_vtxtree = new TTree("vtx", "vtx");
    
    TStopwatch t_tot;
    t_tot.Start();
    
    FillZ_LAYER();

    if (FAST==0.5) {  //extend tracks connected to vertices before checking for fake ones

        FillTracksCells(*arrTRK);
        FillVtxPlate(*arrVTX);
        end_tracksS1 = UnisciTracce(arrVTX, maximp_unione, max_deltatheta_unione);
        merged_arrVTX = arrVTX;
    }
    
    if((NITROGEN_SEARCH<2&&FAST==0)||FAST==1||FAST==4||FAST==100){
        cout << "AnalyseFakeVtxs" << endl;
        merged_arrVTX = AnalyseFakeVtxs(*arrVTX, vrec);
    }
    else if(FAST>1) merged_arrVTX=arrVTX; //c'era FAST>=1.. perché?
    
    cout << "merged_arrVTX " << merged_arrVTX->GetEntries() << endl;

    if(FAST==1||FAST==100||FAST==0.5){
        CreateTree(new_vtxtree, merged_arrVTX);
        mygEVR->eVTX = merged_arrVTX;
        cout << "At the end I have " << merged_arrVTX->GetEntries() << "\t" << new_vtxtree->GetEntries() << endl; //<< "\tTime: " << t_tot.RealTime() << " s\t" << t_tot.RealTime()/60 << " min " << endl;
    }
    if(FAST!=1 && FAST!=100 && FAST!=0.5){
        cout << "FillTracksCells" << endl; //"\tTime: " << t_tot.RealTime() << " s\t" << t_tot.RealTime()/60 << " min " << endl;; //prima lo facevo all'inizio di tutto
        FillTracksCells(*arrTRK);
        
        cout << "FillVtxPlate" << endl; //"\tTime: " << t_tot.RealTime() << " s\t" << t_tot.RealTime()/60 << " min " << endl;
        FillVtxPlate(*merged_arrVTX);
    }
    if(FAST==0||FAST==2||FAST==4){
        cout << "PrepareTRK" << endl; //"\tTime: " << t_tot.RealTime() << " s\t" << t_tot.RealTime()/60 << " min " << endl;
        PrepareTRK();
        if(FAST!=4) return 0;
    }
    
    if(FAST==0||FAST==3){
        cout << "ReadTreeTracksVTA" << endl; // "\tTime: " << t_tot.RealTime() << " s\t" << t_tot.RealTime()/60 << " min " << endl;
        ReadTreeTracksVTA(*arrTRK, *merged_arrVTX);
    }
    cout << endl << "merged_arrVTX->GetEntries(): " << merged_arrVTX->GetEntries() << endl;


    if(FAST==0||(FAST>=3 && FAST<100)){
        if(NITROGEN_SEARCH<2){

            if (EVERBOSE==100 && FAST==3 && CHECK_OXY==1) {
                for(int i=0; i<merged_arrVTX->GetEntries(); i++) 
                {
                    cout << " OXY_FOUND for vertex # " << i << " is " << OXY_FOUND[i] << endl;
                }
            }

            cout << "FindCloseTracks" << endl;
            new_varr=FindCloseTracks(merged_arrVTX, vrec);
        }
        cout << "new_varr->GetEntries(): " << new_varr->GetEntries() << endl;


        if(NITROGEN_SEARCH>=1){
            cout << "FindNitrogen" << endl; //"\tTime: " << t_tot.RealTime() << " s\t" << t_tot.RealTime()/60 << " min " << endl;
            if(NITROGEN_SEARCH==2) new_varr=arrVTX;
            FindNitrogen(new_varr, vrec, maximp_Np, maximp_dau);

            if (EVERBOSE==100 && FAST==4) {
            EdbVertex* v = (EdbVertex*)(new_varr->At(0));
            for (int itrk=0; itrk<v->N(); itrk++) 
            {
                cout << " VTA_ZPos " << v->GetVTa(itrk)->Zpos() << " for track " << (int) (v->GetTrack(itrk)->Track()) << endl;
            }
        }
        }
        cout << "new_varr->GetEntries(): " << new_varr->GetEntries() << endl;
        if(NITROGEN_SEARCH<2){
            cout << "UnisciVertici" << endl;
            final_varr = UnisciVertici(new_varr, vrec); // Oltre a unirli elimino anche alcuni fake
        }
        else final_varr=new_varr;
        cout << "final_varr->GetEntries(): " << final_varr->GetEntries() << endl;

        if (EVERBOSE==100 && FAST==4 && final_varr->GetEntries()>0) {
            EdbVertex* v = (EdbVertex*)(final_varr->At(0));
            for (int itrk=0; itrk<v->N(); itrk++) 
            {
                cout << " VTA_ZPos " << v->GetVTa(itrk)->Zpos() << " for track " << (int) (v->GetTrack(itrk)->Track()) << endl;
            }
        }
                
        cout << "UnisciTracce" << endl; //"\tTime: " << t_tot.RealTime() << " s\t" << t_tot.RealTime()/60 << " min " << endl;
        end_tracksS1 = UnisciTracce(final_varr, maximp_unione, max_deltatheta_unione);
        if (UNISCI_TRACCE_LARGER_CUT) end_tracksS1 = UnisciTracce(final_varr, 2.*maximp_unione, 2.*max_deltatheta_unione);
        
        if (ALLOW_LONGER_GAPS == 1) UnisciTracceLongerGaps(end_tracksS1, maximp_unione, max_deltatheta_unione);    
                
        cout << "CreateTree" << endl; //"\tTime: " << t_tot.RealTime() << " s\t" << t_tot.RealTime()/60 << " min " << endl;
        CreateTree(new_vtxtree, final_varr);

        if (EVERBOSE==100 && FAST==4 && final_varr->GetEntries()>0) {
            EdbVertex* v = (EdbVertex*)(final_varr->At(0));
            for (int itrk=0; itrk<v->N(); itrk++) 
            {
                cout << " VTA_ZPos " << v->GetVTa(itrk)->Zpos() << " for track " << (int) (v->GetTrack(itrk)->Track()) << endl;
            }
        }
               
        mygEVR->eVTX = final_varr;
        
    }
    
    if(FAST!=2){
        fvtx->cd();
        mygEVR->Write();
        new_vtxtree->Write();
        fvtx->Close(); //close the file where vertices are saved
    }
    if(NITROGEN_SEARCH<2) {
        if(FAST<=1){
            cout << "2prong vertices deleted: " << COUNT_ATTACHED_VTX  << endl;
            if(COUNT_ATTACHED_VTX>0 && MC>0) cout << "\t tracks correctly matched: " << COUNT_ATTACHED_OK << "\t" << (float)COUNT_ATTACHED_OK/(float)COUNT_ATTACHED_VTX*100 << "\%" << endl;
            if(COUNT_ATTACHED_VTX>0 && MC>0) cout << "\t tracks correctly matched excluding bkg: " << COUNT_ATTACHED_OK_NOBKG << "\t" << (float)COUNT_ATTACHED_OK_NOBKG/(float)COUNT_ATTACHED_VTX*100 << "\%"<< endl;
            
            cout << "3prong vertices deleted: " << COUNT_3P_ATTACHED_VTX  << endl;
            if(COUNT_3P_ATTACHED_VTX>0 && MC>0) cout << "\t THREEP_REMOVE_BKG_OK: " << THREEP_REMOVE_BKG_OK  << "\t" << (float)THREEP_REMOVE_BKG_OK/(float)COUNT_3P_ATTACHED_VTX*100 << "\%"<< endl;
            if(COUNT_3P_ATTACHED_VTX>0 && MC>0)cout << "\t THREEP_ATTACHED_OK: " << THREEP_ATTACHED_OK << "\t" << (float)THREEP_ATTACHED_OK/(float)COUNT_3P_ATTACHED_VTX*100 << "\%"<< endl;
            
            cout << "vertices splitted: " << COUNT_SPLIT;
            if(COUNT_SPLIT>0 && MC>0) cout << "\tCOUNT_SPLIT_OK: " << COUNT_SPLIT_OK << "\t" << (float)COUNT_SPLIT_OK/(float)COUNT_SPLIT*100 << "\%";
            cout << endl;
            
            cout << "TRACKS REMOVED BECAUSE IP> " << Max_acceptable_IP << ": " << TRACKS_REMOVED << endl;
            
        }
        cout << "COSMICS REMOVED: " << REMOVECOSMIC;
        if(REMOVECOSMIC>0 && MC>0) cout << "\t OK: " << REMOVECOSMIC_OK << "\t" << (float)REMOVECOSMIC_OK/REMOVECOSMIC*100 << "\%";
        cout << endl;
        
        cout << "BEAMDAU VERTICES REMOVED: " << REM_BEAMDAU;
        if(MC>0) cout << "\tBEAMDAU OK " << REM_BEAMDAU_OK;
        cout << endl;
        
        
        if(FAST>=3 || FAST==0){
            cout << "N_BEAM_FOUND\t" << N_BEAM_FOUND << endl;
            if(N_BEAM_FOUND>0 && MC>0) cout << "\tN_BEAM_TRUE\t" << N_BEAM_TRUE << "\t(" << (float)N_BEAM_TRUE/(float)N_BEAM_FOUND*100 << "\%)";
            cout << endl;
            cout << "N_BEAM_SPLITTED: " << N_BEAM_SPLITTED << endl;
            if(N_BEAM_SPLITTED>0 && MC>0) cout << "\tOK: " << N_BEAM_SPLITTED_OK << "\t(" << (float)N_BEAM_SPLITTED_OK/(float)N_BEAM_SPLITTED*100 << "\%)";
            cout << endl;
            
            cout << "N_DAU_FOUND\t" << N_DAU_FOUND << endl;
            if(N_DAU_FOUND>0 && MC>0) cout << "\tN_DAU_TRUE\t" << N_DAU_TRUE << "\t(" << (float)N_DAU_TRUE/(float)N_DAU_FOUND*100 << "\%)" << endl;
            cout << "N_BEAM_REMOVED\t" << N_BEAM_REMOVED << endl;
        }
    }
    if((NITROGEN_SEARCH>=1&&FAST>=3) || (NITROGEN_SEARCH==2)) {
        cout << "N_NP_FOUND\t" << N_NP_FOUND;
        if(MC>0) cout << "\tN_NP_FOUND OK\t" << N_NP_FOUND_OK;
        cout << endl;
        cout << "N_NPP_FOUND\t" << N_NPP_FOUND;
        if(MC>0) cout << "\tN_NPP_FOUND OK\t" << N_NPP_FOUND_OK;
        cout << endl;
    }
    if(FAST>=3 && FAST<100) {
        cout << "EXTEND_TRACK: " << EXTEND_TRACK << endl;
        if(EXTEND_TRACK>0 && MC>0) cout << "\tEXTEND_TRACK_OK: "  << EXTEND_TRACK_OK <<  " (" << (float)EXTEND_TRACK_OK/EXTEND_TRACK*100 << "\%)" << endl;
        cout << "MERGED VTX: " << MERGEVTX << endl;
    }
    if(FAST!=1) cout << "At the end I have " << final_varr->GetEntries() << " vertices" << "\tTime: " << t_tot.RealTime() << " s\t" << t_tot.RealTime()/60 << " min " << endl;

    
    return 0;
    
}

//---------------------------------------------------------//

void checkpatterns_antonio(){
    //check for patterns, if there is one missing add it
    int np = ali->Npatterns();
    
    TFile *setfile = TFile::Open(Form("b%06d.0.0.0.set.root",BRICKID));
    EdbScanSet *set = (EdbScanSet*) setfile->Get("set");
    
    for(int i=0; i<np; i++) {
        EdbPattern *p = ali->GetPattern(i);
        if(!p) {
            cout<<"missing pattern "<<i<<" now adding it "<<endl;
            float zmissingPID = set->eB.GetPlate(i)->Z(); //note, set->eB->GetPlate(i) uses the PID, set->GetPlate(i) uses the number of plate, be careful!
            EdbPattern *pat = new EdbPattern( 0., 0.,zmissingPID);
            pat->SetID(i);
            pat->SetScanID(0);
            ali->AddPatternAt(pat,i);
        }//end if
    } //end for loop
}

//---------------------------------------------------------//

void checkpatterns(){
    //check for patterns, if there is one missing add it
    int np = ali->Npatterns();
    
    TFile *setfile = TFile::Open(Form("b%06d.0.0.0.set.root",BRICKID));
    cout << "I'm checking patterns with set " << "b" << BRICKID << ".0.0.0.set.root" << "\t np: " << np << "\tRunning on " << PLMAX << " plates" << endl;
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

//---------------------------------------------------------//

void FillZ_LAYER(){
    int np = ali->Npatterns();
    cout << "N patterns: " << np << endl;
    int layer=np;
    for(int i=0; i<PLMAX; i++) {
        EdbPattern *p = ali->GetPattern(i);
        if(DIRECTION==9){
            Z_LAYER[layer]=p->Z();
            cout << "layer: " << layer  << "\t" << Z_LAYER[layer] << endl;
            layer--;
        }
        else if(DIRECTION==1) {
            Z_LAYER[i]=p->Z();
            cout << "layer: " << i << "\t" << Z_LAYER[i] << endl;
        }
    }
}


//---------------------------------------------------------//


void PrepareTRK(){
    
    //storing indices of vertices start and end tracks is associated to, if present, in a TTree. Otherwise store -1
    //creating file and TTree
    TString vtatracks_filename;
    if(FAST==4) vtatracks_filename=Form("vtainfotracks_%d.root", FAST);
    else vtatracks_filename=Form("vtainfotracks.root");
    TFile *vtatracksfile = new TFile(vtatracks_filename,"RECREATE");
    TTree *vtatracks = new TTree("vtatracks","vtatracks");
    
    int trid=0, VertexS=0, VertexE=0;
    float temp_theta=-99;
    
    vtatracks->Branch("trid",&trid,"trid/I"); //first segment ID
    vtatracks->Branch("Theta",&temp_theta,"Theta/F"); //theta
    vtatracks->Branch("VertexS",&VertexS,"VertexS/I");
    vtatracks->Branch("VertexE",&VertexE,"VertexE/I");

    TString oxy_found_filename;
    if (FAST==4) oxy_found_filename = Form("oxyfound_%d.root", FAST);
    else oxy_found_filename = Form("oxyfound.root");
    TFile* oxyfoundfile = new TFile(oxy_found_filename, "RECREATE");
    TTree* oxyinfo = new TTree("oxyinfo", "oxyinfo");

    int vertexID=0, OXYFOUND=0;
    oxyinfo->Branch("vID", &vertexID, "vID/I");
    oxyinfo->Branch("OXY_FOUND", &OXYFOUND, "OXY_FOUND/I");
    
    map<int,int> frequencyEvent;
    map<int,int>::iterator it;
    int ntracks_event = 0;
    
    int r[2] = {1,1}; //raggio xy in cui voglio cercare: in questo caso una cella sola è sufficiente
    TObjArray tr_grid_dau, tr_grid_oxy;
    
    for(int ipl = PLMIN; ipl<=BRAGGPLATE; ipl++){
        TObjArray varr = vtxPat[ipl];
        int nv = varr.GetEntries();
        for(int ivtx=0; ivtx<nv; ivtx++){
            EdbVertex *v = (EdbVertex*)(varr.At(ivtx));
            v->ResetTracks();
            int idvtx = v->ID();
            vertexID = idvtx;
            if(EVERBOSE==100 && idvtx==DEBUG_VTXID) cout << "prepare trk belonging to vtx " << idvtx << "\t" << Get_vtx_plate(v->VZ()) << endl;
            OXY_FOUND[idvtx]=0;
            mostfrequentevent[idvtx]=0;
            
            for (int itrk = 0; itrk < v->N(); itrk++){
                EdbTrackP *track = v->GetTrack(itrk);
                int trkplate = track->GetSegmentFirst()->Plate();
                VertexS=-1;
                VertexE=-1;
                temp_theta=-99;
                trid=-99;
                if(EVERBOSE==100 && ((track->Track()==DEBUG_TRKID) || track->MCEvt()==DEBUG_MCEVT||idvtx==DEBUG_VTXID)) cout << endl << "Searching match with track " << track->ID() << "\t" << track->Track() << endl;
                temp_theta=track->Theta();
                //controllo se c'è già oxy entrante
                if(temp_theta<BEAMTHETA && v->GetVTa(itrk)->Zpos()==0 && (track->N()>2||trkplate<4)){ OXY_FOUND[idvtx]++; //segno in un vettore (in ivtx-esima posizione) se il vertice ivtx ha già il fascio incidente
                }

                if (EVERBOSE==100 && idvtx==DEBUG_VTXID) cout << " OXY_FOUND FOR VTX " << idvtx << " : " << OXY_FOUND[idvtx] << endl;
                if (EVERBOSE==100 && idvtx==DEBUG_VTXID) cout << " temp_theta " << temp_theta << ", VTA_Zpos " << v->GetVTa(itrk)->Zpos() << " trackN " << track->N() << " trkplate " << trkplate << ", trackTrack " << track->Track() << endl;

                //resetto VTA
                trid = track->GetSegmentFirst()->ID();
                float xy[2] = {track->X(),track->Y()}; //coordinate della cella in cui voglio cercare
                
                EdbTrackP *tr2 = new EdbTrackP();
                int tr2_found=0;
                
                //cerco tra tutte le tracce della cella qual è quella che corrisponde a "track".
                //cerco la traccia tra i dau
                if(EVERBOSE==100 && ((track->Track()==DEBUG_TRKID) || track->MCEvt()==DEBUG_MCEVT||idvtx==DEBUG_VTXID)) cout << "check trkplate " << trkplate << "\t" << xy[0] << " " << xy[1] << endl;
                int ntrks_grid_dau = gridtr_DAU[trkplate].SelectObjectsC( xy, r, tr_grid_dau);
                if(EVERBOSE==100 && ((track->Track()==DEBUG_TRKID) || track->MCEvt()==DEBUG_MCEVT||idvtx==DEBUG_VTXID)) cout << "check ntrks_grid_dau " << ntrks_grid_dau << endl;
                for (int itrk2 = 0; itrk2 < ntrks_grid_dau; itrk2++){
                    tr2 = (EdbTrackP*) tr_grid_dau.At(itrk2);
                    if(tr2->GetSegmentFirst()->ID()==trid) {
                        if(EVERBOSE==11 || ( EVERBOSE==100 && (idvtx==DEBUG_VTXID || (( track->Track()==DEBUG_TRKID || tr2->Track()==DEBUG_TRKID)) || track->MCEvt()==DEBUG_MCEVT||tr2->MCEvt()==DEBUG_MCEVT))){
                            cout << endl << "trk okk " << trkplate << " " << tr2->GetSegmentFirst()->Plate() << " | " << track->MCEvt() << " " << tr2->MCEvt() << " | " << track->ID() << " " << tr2->ID() << " | " << track->Track() << " " << tr2->Track() << " | " << trid << " " << tr2->GetSegmentFirst()->ID() << " | " << track->Theta() << " " << tr2->Theta() << " | " << track->X() << " " << tr2->X() << " | " << track->Y() << " " << tr2->Y() << "\t";
                        }
                        
                        EdbVTA *vta = new EdbVTA();
                        vta = v->GetVTa(itrk);
                        if(vta) {
                            tr2->AddVTA(vta);
                            if(tr2->VertexS()) VertexS = tr2->VertexS()->ID();
                            if(tr2->VertexE()) VertexE = tr2->VertexE()->ID();
                            
                            if(EVERBOSE==11 || ( EVERBOSE==100 && (idvtx==DEBUG_VTXID || (( track->Track()==DEBUG_TRKID || tr2->Track()==DEBUG_TRKID)) || track->MCEvt()==DEBUG_MCEVT||tr2->MCEvt()==DEBUG_MCEVT))){
                                if(tr2->VertexS()) cout << "vtxS: " << VertexS << "\t";
                                if(tr2->VertexE()) cout << "vtxE: " << VertexE << "\t";
                            }
                        }
                        tr2_found=1;
                        break;
                    }
                }
                tr_grid_dau.Clear();
                
                //}
                if(temp_theta<=BEAMTHETA && trkplate<=BRAGGPLATE+1 && track->N()>4){//Se la traccia ha l'angolo del fascio faccio il loop anche tra gli oxy
                    int ntrks_grid_oxy = gridtr_OXY[trkplate].SelectObjectsC( xy, r, tr_grid_oxy);
                    for (int itrk2 = 0; itrk2 < ntrks_grid_oxy; itrk2++){
                        tr2 = (EdbTrackP*) tr_grid_oxy.At(itrk2);
                        if(tr2->GetSegmentFirst()->ID()==trid) { //check
                            //se ho trovato la traccia resetto VTA e vado avanti
                            EdbVTA *vta = new EdbVTA();
                            vta = v->GetVTa(itrk);
                            if(vta) {
                                tr2->AddVTA(vta);
                                
                                if(tr2->VertexS()) VertexS = tr2->VertexS()->ID();
                                if(tr2->VertexE()) VertexE = tr2->VertexE()->ID();
                                
                                //                                if(EVERBOSE==100 && (idvtx==DEBUG_VTXID || (track->Track()==DEBUG_TRKID) ||track->MCEvt()==DEBUG_MCEVT)){
                                //                                    cout << tr2->MCEvt() << "\t" << tr2->ID() << "\t" << tr2->Flag() << "\t";
                                //                                    cout << "vtxS: " << VertexS << "\t";
                                //                                    cout << "vtxE: " << VertexE << "\t";
                                //                                }
                            }
                            tr2_found=1;
                            break;
                        }
                    }
                    tr_grid_oxy.Clear();
                }
                if(track==0 && track->Flag()==0) track->SetFlag(9); //se non ho trovato la traccia da nessuna parte assegno flag 9
                
                if(MC==1){ //per il MC conto le frequenze dell'eventID (MC) del primo segmento di ogni traccia
                    int eventtrack = track->GetSegmentFirst()->MCEvt();
                    frequencyEvent[eventtrack]++;
                }
                
                if(EVERBOSE==100 && ((track->Track()==DEBUG_TRKID) || track->MCEvt()==DEBUG_MCEVT||idvtx==DEBUG_VTXID)) cout << "check fill " << trid << "\t" << temp_theta << "\t" << VertexS << "\t" << VertexE << endl;
                
                
                vtatracks->Fill();
                
            }
            //esco dal loop delle tracce e per gli ev MC mi salvo l'evento più frequente nel vertice ivtx in ivtx-esima posizione del vettore mostfrequentevent
            
            OXYFOUND = OXY_FOUND[idvtx];
            oxyinfo->Fill();

            if(MC==1){
                for (it = frequencyEvent.begin(); it!=frequencyEvent.end();it++){
                    if(it->second > ntracks_event){
                        ntracks_event = it->second;
                        mostfrequentevent[idvtx] = it->first;
                    }
                }
            }
        }
    }
    
    //writing obtained info
    vtatracksfile->cd();
    vtatracks->Write();
    vtatracksfile->Close();

    oxyfoundfile->cd();
    oxyinfo->Write();
    oxyfoundfile->Close();
}


//---------------------------------------------------------------------
void ReadTreeTracksVTA(TObjArray &arrt, TObjArray &arrv){
    //need to reassociate the vtas...
    TString vtatracks_filename;
    if(FAST==4) vtatracks_filename=Form("vtainfotracks_%d.root", FAST);
    else vtatracks_filename=Form("vtainfotracks.root");
    TFile *vtatracksfile = TFile::Open(vtatracks_filename,"READ");
    TTree *vtatracks = (TTree*) vtatracksfile->Get("vtatracks");
    
    //setting branches
    int trid, VertexS, VertexE;
    float theta=0;
    vtatracks->SetBranchAddress("trid",&trid);
    vtatracks->SetBranchAddress("Theta",&theta); //theta
    vtatracks->SetBranchAddress("VertexS",&VertexS);
    vtatracks->SetBranchAddress("VertexE",&VertexE);
    vtatracks->BuildIndex("trid");
    
    TStopwatch t;
    t.Start();

    if(FAST==4 && EVERBOSE==100) cout << " arrt.GetEntries " << arrt.GetEntries() << " arrv.GetEntries " << arrv.GetEntries() << endl;
    
    for (int itrk; itrk < arrt.GetEntries(); itrk++){
        //getting track and vertices
        EdbTrackP *mytrack = (EdbTrackP*) arrt.At(itrk);
        int trid2=mytrack->GetSegmentFirst()->ID();
        vtatracks->GetEntryWithIndex(trid2);
        
        if(trid2==mytrack->GetSegmentFirst()->ID()&& theta==mytrack->Theta()){
            if(EVERBOSE==100 && ((mytrack->Track()==DEBUG_TRKID)||mytrack->MCEvt()==DEBUG_MCEVT)) cout << trid2 << "\t" << mytrack->MCEvt() << "\t" << mytrack->Track() << "\t" << mytrack->GetSegmentFirst()->ID() << "\t" << mytrack->Theta() << "\t" << theta << "\t";
            
            if(FAST==4) {
                if (VertexS == DEBUG_VTXID || VertexS == DEBUG_MCEVT) VertexS=0;
                else if (VertexE == DEBUG_VTXID || VertexE == DEBUG_MCEVT) VertexE=0;
                else continue;
            }
            
            if (VertexS >= 0){
                EdbVertex *vtxS = (EdbVertex*) arrv.At(VertexS);
                if (EVERBOSE==100) cout << " In ReadTreeTracks Vertex S" << VertexS << " Vertex ID " << vtxS->ID() << endl; 
                EdbVTA *VTAS = new EdbVTA(mytrack, vtxS);
                VTAS->SetZpos(1);
                //adding them to tracks and vertices
                vtxS->AddVTA(VTAS);
                mytrack->AddVTA(VTAS);
            }//end vertexS condition
            if (VertexE >= 0){
                EdbVertex *vtxE = (EdbVertex*) arrv.At(VertexE);
                if (EVERBOSE==100) cout << " In ReadTreeTracks Vertex E" << VertexE << " Vertex ID " << vtxE->ID() << endl;
                //building VTAs
                EdbVTA *VTAE = new EdbVTA(mytrack, vtxE);
                VTAE->SetZpos(0);
                vtxE->AddVTA(VTAE);
                mytrack->AddVTA(VTAE);
            } //end vertexE condition
            if(EVERBOSE==100 && ((mytrack->Track()==DEBUG_TRKID) || mytrack->MCEvt()==DEBUG_MCEVT)) cout << VertexS << "\t" << VertexE << endl;
        }
        
        if (itrk%20000 == 0 && (EVERBOSE==14||EVERBOSE==16||EVERBOSE==101)) {
            cout << " ReadTreeTracks: " << 100.*itrk/arrt.GetEntries() << " %, Iteration Time: " << t.RealTime() << " s\t" << t.RealTime()/60 << " min " << endl;
            t.Reset();
            t.Start();
        }
        
    }//end loop over tracks
    
    
}

//---------------------------------------------------------//

void FillTracksCells(TObjArray &arrt){
    const int cellsize=1000;
    const int ncellsX = (int)(xmax-xmin)/cellsize;//35;
    const int ncellsY = (int)(ymax-ymin)/cellsize;//50;
    cout << "CELLS: " << ncellsX << "\t" << ncellsY << "\t(cell size " << cellsize << ")" << endl;
    const int maxpercell = 100; // è il numero massimo di oggetti per ogni cella. Non restituisce errore se ce ne sono di più, per questo bisogna controllare con la stampa sotto.
    for(int i=1; i<=BRAGGPLATE+1; i++) gridseg_OXY[i].InitCell(ncellsX, xmin, xmax, ncellsY, ymin, ymax, maxpercell );
    for(int i=1; i<=BRAGGPLATE+1; i++) gridtr_OXY[i].InitCell(ncellsX, xmin, xmax, ncellsY, ymin, ymax, maxpercell );
    for(int i=1; i<=PLMAX; i++) gridtr_DAU[i].InitCell(ncellsX, xmin, xmax, ncellsY, ymin, ymax, maxpercell );
    for(int i=1; i<=PLMAX; i++) gridtr_DAU_end[i].InitCell(ncellsX, xmin, xmax, ncellsY, ymin, ymax, maxpercell );
    
    
    int ntr = arrt.GetEntries();
    for(int itr=0; itr<ntr; itr++){
        EdbTrackP *t = (EdbTrackP*)(arrt.At(itr));
        
        //    for(int ipl=PLMIN; ipl<=BRAGGPLATE; ipl++) //we cannot access track segments from a given z, trying to loop over base tracks
        //if(t->Flag()!=0){ //VORREI RIMUOVERE COSMICI
        float x_tr = t->X();
        float y_tr = t->Y();
        //if(x_tr>=xmin&&x_tr<=xmax&&y_tr>=ymin&&y_tr>=ymax){
        int plate_tr = t->GetSegmentFirst()->Plate();
        int plate_tr_end = t->GetSegmentLast()->Plate();
        
        float theta_tr = t->Theta();
        if(EVERBOSE==100 && ((t->Track()==DEBUG_TRKID) || (t->MCEvt()==DEBUG_MCEVT && MC==1))) cout << "track: " << t->Track() << "\t" << t->MCEvt() << "\t" << t->MCTrack() << "\t" << t->N() << "\t" << x_tr << "\t" << y_tr << "\t" << plate_tr << endl;
        
        if(theta_tr<=BEAMTHETA && plate_tr<=BRAGGPLATE+1 && (x_tr>40000 && x_tr<80000 && y_tr>35000 && y_tr<65000) &&  (t->N()>3||(t->GetSegmentLast()->Plate()<=3 && t->N()==t->GetSegmentLast()->Plate()))){ //seleziono gli ossigeni dall'angolo e Seleziono solo le tracce che partono in S1 prima del picco di bragg + almeno 5 segmenti o che inizi prima del piatto 5 e non abbia buchi
            EdbVertex *vtempS = t->VertexS();
            EdbVertex *vtempE = t->VertexE();
            if((!vtempS) && (!vtempE)){ //controllo che la traccia non abbia vertici entranti o uscenti
                
                gridtr_OXY[plate_tr].AddObject( x_tr, y_tr, (TObject*)t ); //qui ci sono tutte le tracce che iniziano in un piatto
                if(EVERBOSE==100 && ((t->Track()==DEBUG_TRKID) || (t->MCEvt()==DEBUG_MCEVT && MC==1))) cout << "\tfill OXY" << endl;
                
                for (int iseg=0; iseg < t->N(); iseg++){ //salvo la traccia in corrispondenza di tutti i segmenti perché l'idea è di spezzare tracce passanti attraverso il vertice
                    //float z = Z_LAYER[ipl];
                    EdbSegP *seg = t->GetSegment(iseg);
                    float x = seg->X();
                    float y = seg->Y();
                    int ipl = seg->Plate();
                    if(ipl<=BRAGGPLATE+1){
                        gridseg_OXY[ipl].AddObject( x, y, (TObject*)t );
                        if(EVERBOSE==100 && ((t->Track()==DEBUG_TRKID) || (t->MCEvt()==DEBUG_MCEVT && MC==1))) cout << "\tfill segOXY: seg" << ipl << "\t" << seg->ID() << "\t" << seg->MCEvt() << "\t" << seg->MCTrack() << "\t" << seg->Plate() << endl;
                    }
                } //end loop plates
            }// if((!vtempS) && (!vtempE)){
        } // end BEAMTHETA
        //else{
        
        gridtr_DAU[plate_tr].AddObject( x_tr, y_tr, (TObject*)t ); //nel caso dei dau voglio solo le tracce che effettivamente iniziano lì
        gridtr_DAU_end[plate_tr_end].AddObject(x_tr, y_tr, (TObject*)t); //nel caso dei dau voglio solo le tracce che effettivamente finiscono lì
        if(t->GetSegmentFirst()->Plate()==DEBUG_S0_PLATE_ST && t->GetSegmentFirst()->ID()==DEBUG_S0_ID_ST && EVERBOSE==15) {
            cout << " Candidate for track is added to gridtr_DAU with x_tr " << x_tr << " and y_tr " << y_tr << endl;
        }
        
        if(EVERBOSE==100 && ((t->Track()==DEBUG_TRKID) || (t->MCEvt()==DEBUG_MCEVT && MC==1))) cout << "\tfill DAU" << endl;
        //} // end else
        //}//end flag
    }//end loop tracks
    
    if(EVERBOSE==1){
        for(int i=1; i<BRAGGPLATE; i++) {
            cout << endl << endl << i << endl;
            gridseg_OXY[i].PrintStat();
            gridtr_OXY[i].PrintStat();
            gridtr_DAU[i].PrintStat();
            gridtr_DAU_end[i].PrintStat();
        }
    }
    if(EVERBOSE==16) {
        int check = 0;
        for(int i=1; i<PLMAX; i++) {
            TH2F* h2 = gridtr_DAU[i].DrawH2();
            int Nmax = h2->GetMaximum();
            if (Nmax>=maxpercell) {cout << " ERROR! Too many objects inside cell # " << i << endl; check=1;}
        }
        if (check == 0)  cout << " Grid Occupancy OK " << endl;
        
        //TFile* f2 = new TFile("grid200.root", "RECREATE");
        //h2->Write("h_gr");
        //f2->Close();
    }
    
    cout << endl;
    
}


//---------------------------------------------------------//

TObjArray* AnalyseFakeVtxs(TObjArray &arrv, EdbVertexRec *vrec){
    EdbTrackP* firsttrack;
    EdbTrackP* secondtrack;
    int COUNT_STRANGE=0;
    int nv = arrv.GetEntries();
    TObjArray *merged_varr = new TObjArray(); // Metterò qui tutti i vertici senza quelli eliminati
    int temp_trkidmc=0;
    int vID=0;
    
    if(EVERBOSE==5.5) cout << "ARRAY INIZIO " << nv << endl;
    
    for(int iv=0; iv<nv; iv++){
        EdbVertex *vertex = (EdbVertex*)(arrv.At(iv));
        vertex->ResetTracks();
    }
    
    for(int iv=0; iv<nv; iv++){
        EdbVertex *vertex = (EdbVertex*)(arrv.At(iv));
        vID=vertex->ID();
		int mcevt = (int) (vertex->GetTrack(0)->MCEvt());
        
        if((FAST==4||FAST==100) && (DEBUG_VTXID!=-99 || mcevt!= -99)){ //TAGLIO CUT FAST 2
            if(vID!=DEBUG_VTXID && mcevt!= DEBUG_MCEVT) {
                vertex->SetFlag(-99);
                continue;
            } 
            else cout << "vertex id " << vID << " found with flag " << vertex->Flag() <<  endl;
        }
        int First_MCev=0, Second_MCev=0, First_MCtr=0, Second_MCtr=0, Removed_MCEv=0;
        
        if(vertex->N()==2){//only 2 prong vertices
            firsttrack = (EdbTrackP*) vertex->GetTrack(0); //in two prongs with flag 1, first track is always the most upstream
            secondtrack = (EdbTrackP*) vertex->GetTrack(1);
            if(EVERBOSE==5.5) cout << "N=2" << endl;
            
            int First_ID = firsttrack->Track();
            int Second_ID = secondtrack->Track();
            int First_N = firsttrack->N();
            int Second_N = secondtrack->N();
            float First_Theta = firsttrack->Theta();
            float Second_Theta = secondtrack->Theta();
            int First_Plate1 = firsttrack->GetSegmentFirst()->Plate();
            int Second_Plate1 = secondtrack->GetSegmentFirst()->Plate();
            
            //riattacco tracce erroneamente viste come 2 prong vtx
            if (vertex->Flag()==1 || vertex->Flag()==4){
                if(First_N>3){
                    if(vertex->GetVTa(0)->Zpos()!=vertex->GetVTa(1)->Zpos() && TMath::Abs(First_Theta-Second_Theta)<0.1){
                        
                        //if((firsttrack->Flag()!=0 && secondtrack->Flag()!=0)|| MC==11){
                        if(1==1){
                            MergeTrack(firsttrack,secondtrack);
                            vertex->SetFlag(-99);
                            //if(EVERBOSE==5.5 || (EVERBOSE==100 && (vID==DEBUG_VTXID || (( First_ID==DEBUG_TRKID||Second_ID==DEBUG_TRKID) || firsttrack->MCEvt()==DEBUG_MCEVT || secondtrack->MCEvt()==DEBUG_MCEVT)))) cout << "vertex " << vID << " deleted. Tracks " << First_ID << " " << Second_ID << " merged" << endl;
                            
                            COUNT_ATTACHED_VTX++;
                            
                            if(MC==11){
                                First_MCev = firsttrack->GetSegmentLast()->MCEvt();
                                First_MCtr = firsttrack->GetSegmentLast()->MCTrack();
                                Second_MCev = secondtrack->GetSegmentFirst()->MCEvt();
                                Second_MCtr = secondtrack->GetSegmentFirst()->MCTrack();
                                if(EVERBOSE==5.5 || (EVERBOSE==100 && ((firsttrack->MCEvt()==DEBUG_MCEVT || secondtrack->MCEvt()==DEBUG_MCEVT)))) cout<< vID<< "\t"<< First_N << " - " << First_MCev << " - " << First_MCtr << " " << First_Theta << "\t" << Second_N << " - " << Second_MCev << " - " << Second_MCtr << endl;
                                
                                if(First_MCev == Second_MCev && First_MCtr == Second_MCtr) {
                                    if(EVERBOSE==5.5 || (EVERBOSE==100 && ((firsttrack->MCEvt()==DEBUG_MCEVT || secondtrack->MCEvt()==DEBUG_MCEVT)))) cout << First_MCev << "\t" << First_MCtr << "\t" << vertex->VZ() << endl;
                                    COUNT_ATTACHED_OK++;
                                    if(First_MCev!=-999) COUNT_ATTACHED_OK_NOBKG++;
                                }
                            }
                            
                            EdbVertex *v1 = secondtrack->VertexE();
                            if(v1){
                                v1=AddTrackToVertex2(vrec, v1, firsttrack,vertex->GetVTa(0)->Zpos());
                                if(v1){
                                    //if(EVERBOSE==5.5 || (EVERBOSE==100 && ((firsttrack->MCEvt()==DEBUG_MCEVT || secondtrack->MCEvt()==DEBUG_MCEVT)||((First_ID==DEBUG_TRKID || Second_ID==DEBUG_TRKID)))))  cout << "new vtx " << v1->ID() << " N: " << v1->N() << endl;
                                    
                                    int its=0;
                                    for(its=0; its<v1->N(); its++){
                                        EdbTrackP *temptrack = v1->GetTrack(its);
                                        //cout << temptrack->Track() << " (" << temptrack->N() << ")\t" << secondtrack->Track() << " (" << Second_N << ")" << endl;
                                        if(temptrack->Track()==secondtrack->Track())
                                            break;
                                    }
                                    if(v1->N()>2) v1=vrec->RemoveTrackFromVertex(v1, its);
                                    else v1->SetFlag(-99);
                                }
                            }
                        }
                    }
                }
            }
            
            int split=0;
            //Elimino vertici fatti da due ossigeni o ossigeno + fondo o fondo+fondo
            
            if(vertex->Flag()==0 || vertex->Flag()==3){
                
                if(First_Theta<=BEAMTHETA && Second_Theta<=BEAMTHETA){ //due ossigeni
                    vertex->SetFlag(-99);
                    if(EVERBOSE==100 && (vID==DEBUG_VTXID||vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT)) cout << "vertex deleted rigo 457" << endl;
                    
                    split=1;
                }
                else if((First_Theta<=BEAMTHETA&&Second_N<=3) || (Second_Theta<=BEAMTHETA && First_N<=3)){ //ossigeno + bkg
                    vertex->SetFlag(-99);
                    if(EVERBOSE==100 && (vID==DEBUG_VTXID||vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT)) cout << "vertex deleted rigo 463" << endl;
                    
                    split=1;
                }
            }
            if((First_N<=3 && Second_N<=3)||(First_N==2||Second_N==2)){ //bkg+bkg
                vertex->SetFlag(-99);
                if(EVERBOSE==100 && (vID==DEBUG_VTXID||vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT)) cout << "vertex deleted rigo 470" << endl;
                
                split=1;
            }
            
            if(split==1){
                COUNT_SPLIT++;
                
                if(MC==11){
                    
                    First_MCev = firsttrack->GetSegmentFirst()->MCEvt();
                    First_MCtr = firsttrack->GetSegmentFirst()->MCTrack();
                    Second_MCev = secondtrack->GetSegmentFirst()->MCEvt();
                    Second_MCtr = secondtrack->GetSegmentFirst()->MCTrack();
                    
                    if(First_MCev!=Second_MCev||First_MCev==-999||Second_MCev==-999) COUNT_SPLIT_OK++;
                    
                    //if(EVERBOSE==5.5 || (EVERBOSE==100 && ((firsttrack->MCEvt()==DEBUG_MCEVT || secondtrack->MCEvt()==DEBUG_MCEVT)||((firsttrack->Track()==DEBUG_TRKID || secondtrack->Track()==DEBUG_TRKID))))) cout << vID << "\t" << vertex->Flag() << "\t1st: " << First_MCev << " " << First_MCtr << " " <<  First_Plate1 << "\t\t2nd: " << Second_MCev << " " << Second_MCtr << " " <<  Second_Plate1 << endl;
                    
                }
            }
        }
        else if(vertex->N()==3){
            int NPl[3]={0,0,0};
            float Th[3]={0,0,0};
            int zpos[3]={0,0,0};
            int warning=0, bad_itrk=-99, beam_itrk=-99, attach_itrk=-99;
            if(EVERBOSE==5.5) cout << "N=3" << endl;
            for (int itrk = 0; itrk < vertex->N(); itrk++){
                EdbTrackP *track = vertex->GetTrack(itrk);
                NPl[itrk]=track->N();
                Th[itrk]=track->Theta();
                zpos[itrk]=vertex->GetVTa(itrk)->Zpos();
                if(zpos[itrk]==0 && Th[itrk]<=BEAMTHETA) {
                    beam_itrk=itrk;
                    //cout << "beam trk: " << zpos[itrk] << "\t" << Th[itrk] << "\t" << NPl[itrk] << endl;
                }
                else if(zpos[itrk]==0 && Th[itrk]>BEAMTHETA && NPl[itrk]<=2) {
                    warning++;
                    bad_itrk = itrk;
                    //cout << "bad trk: " << zpos[itrk] << "\t" << Th[itrk] << "\t" << NPl[itrk] << endl;
                }
                else if(zpos[itrk]==1 && Th[itrk]<=BEAMTHETA){
                    warning++;
                    attach_itrk=itrk;
                    //cout << "attach trk: " << zpos[itrk] << "\t" << Th[itrk] << "\t" << NPl[itrk] << endl;
                }
            }
            if(warning==2){
                
                vertex->SetFlag(-99);
                if(EVERBOSE==100 && (vID==DEBUG_VTXID||vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT) ) cout << "vertex deleted rigo 520" << endl;
                
                COUNT_3P_ATTACHED_VTX++;
                
                //EdbTrackP *firsttrack, *secondtrack, *badtrack;
                EdbTrackP *badtrack;
                
                if(beam_itrk!=-99) firsttrack = vertex->GetTrack(beam_itrk);
                if(attach_itrk!=-99) secondtrack = vertex->GetTrack(attach_itrk);
                if(bad_itrk!=-99) badtrack = vertex->GetTrack(bad_itrk);
                
                if(MC==11){
                    if(beam_itrk!=-99) First_MCev = firsttrack->GetSegmentLast()->MCEvt();
                    if(beam_itrk!=-99) First_MCtr = firsttrack->GetSegmentLast()->MCTrack();
                    if(attach_itrk!=-99) Second_MCev = secondtrack->GetSegmentFirst()->MCEvt();
                    if(attach_itrk!=-99) Second_MCtr = secondtrack->GetSegmentFirst()->MCTrack();
                    if(bad_itrk!=-99) Removed_MCEv = badtrack->GetSegmentFirst()->MCEvt();
                    
                    if(Removed_MCEv==-999||First_MCev==-999||Second_MCev==-999)
                        THREEP_REMOVE_BKG_OK++;
                    //cout << "A " << First_MCev << "\t" << iv << "\t" << vID << endl;
                }
                
                //if(EVERBOSE==5 || (EVERBOSE==100 && ((firsttrack->MCEvt()==DEBUG_MCEVT || secondtrack->MCEvt()==DEBUG_MCEVT)||((firsttrack->Track()==DEBUG_TRKID || secondtrack->Track()==DEBUG_TRKID)))))  cout << "3p removed " << Removed_MCEv << "\t"<< badtrack->N() << "\t" << badtrack->Theta() << "\t1st: " << First_MCev << " " << First_MCtr << " " <<  First_Plate1 << "\t\t2nd: " << Second_MCev << " " << Second_MCtr << " " <<  Second_Plate1 << endl;
                
                
                if(firsttrack && secondtrack){
                    if(secondtrack->GetSegmentFirst()->Plate() > firsttrack->GetSegmentLast()->Plate()){
                        MergeTrack(firsttrack, secondtrack);
                        
                        if(First_MCev==Second_MCev && First_MCtr==Second_MCtr) {
                            THREEP_ATTACHED_OK++;
                            
                            //if(EVERBOSE==5 || (EVERBOSE==100 && ((firsttrack->MCEvt()==DEBUG_MCEVT || secondtrack->MCEvt()==DEBUG_MCEVT)||((firsttrack->Track()==DEBUG_TRKID || secondtrack->Track()==DEBUG_TRKID)))))  cout << "tracks " << First_MCev << " " << First_MCtr << " and " << Second_MCev << " " << Second_MCtr << "merged " << endl;
                            
                        }
                    }
                    
                }
            }
        }
        
        //possibile crash!
        if(vertex->Flag()!=-99){
            for (int itrk = 0; itrk < vertex->N(); itrk++){
                EdbTrackP *track = vertex->GetTrack(itrk);
                if(track->GetSegmentFirst()->Z()<vertex->VZ() && track->GetSegmentLast()->Z()>vertex->VZ() && vertex->GetVTa(itrk)->Z_pos()==0 && track->N()<8){
                    COUNT_STRANGE++;
                    vertex->SetFlag(-99);
					if(EVERBOSE==100 && (vID==DEBUG_VTXID||vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT)) { cout << "vertex deleted rigo 935" << endl; cout << " Track was " << itrk << " Seg " << track->GetSegmentFirst()->ID() << " " << track->GetSegmentFirst()->Plate() << " " << track->GetSegmentFirst()->Z() << " Last seg " << track->GetSegmentLast()->Plate() << " " << track->GetSegmentLast()->Z() << " vertex z " << vertex->VZ() << endl;}
                    //if(EVERBOSE==5.5 || (EVERBOSE==100 && (vID==DEBUG_VTXID || (( First_ID==DEBUG_TRKID||Second_ID==DEBUG_TRKID) || firsttrack->MCEvt()==DEBUG_MCEVT || secondtrack->MCEvt()==DEBUG_MCEVT))))             cout << vertex->ID() << " strange" << endl;
                    
                    if((vertex->GetVTa(itrk)->Zpos()==1 && vertex->VZ()-track->GetSegmentLast()->Z()>200)|| (vertex->GetVTa(itrk)->Zpos()==0 && track->GetSegmentFirst()->Z()-vertex->VZ()>200)){
                        cout << iv << "\t" << vertex->ID() << "\t" << vertex->VZ() << "\t" << track->GetSegmentFirst()->Z() << "\t" << track->GetSegmentLast()->Z()<< "\t" << track->N() << "\t" << track->Track() << "\t" << track->Track() << endl;
                        //procedura che crasha se c'è più di una traccia rimossa WARNING!
                        //EdbVTA *vta1 = (EdbVTA*)vertex->GetVTa(itrk);
                        //vertex = vrec->RemoveVTAFromVertex(*vertex, *vta1);
                        break;
                    }
                }}}
        
        
        if(vertex->VZ()<Z_LAYER[1]-THICKNESS || vertex->VZ()>Z_LAYER[BRAGGPLATE]){
			if(EVERBOSE==100 && (vID==DEBUG_VTXID||vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT)) cout << "vertex deleted rigo 949 because vertex VZ = " << vertex->VZ() << " Z_LAYER[0] " << Z_LAYER[0] << " Z_LAYER[1] " << Z_LAYER[1] << endl;
            vertex->SetFlag(-99);

            //if(EVERBOSE==5.5 || (EVERBOSE==100 && (vID==DEBUG_VTXID || (( First_ID==DEBUG_TRKID||Second_ID==DEBUG_TRKID) || firsttrack->MCEvt()==DEBUG_MCEVT || secondtrack->MCEvt()==DEBUG_MCEVT))))             cout << "vz " << vertex->VZ() << "\t" << Z_LAYER[BRAGGPLATE] << endl;
        }
              
        
    }//chiudo ciclo sui vtx
    
    
    //salvo il nuovo vertice di array che esclude quelli eliminati
    if(EVERBOSE==5.5) cout << "ARRAY FINE " << arrv.GetEntries() << endl;
    int FLAG99=0;
    for(int iv=0; iv<arrv.GetEntries(); iv++){
        EdbVertex *vertex = (EdbVertex*)(arrv.At(iv));
        if(vertex->Flag()!=-99){
            merged_varr->Add(vertex);
            if(EVERBOSE==100 && ( (vertex->GetTrack(0)->Track()==DEBUG_TRKID) || vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT ||  vID==DEBUG_VTXID)){                cout << "Final vtx " << iv << endl;
                for(int itrk=0; itrk<vertex->N(); itrk++){
                    cout << "\t\t" << vID << "\t" << vertex->GetTrack(itrk)->MCEvt() << "\t" << vertex->GetTrack(itrk)->Track() << "\t" << vertex->GetVTa(itrk)->Imp() << "\t" << vertex->Flag() << endl;
                }
            }
        }
        else if(EVERBOSE==100 && vID==DEBUG_VTXID){
            cout << "delete vtx" << endl;
            for(int itrk=0; itrk<vertex->N(); itrk++){
                cout << "\t\t" << vID << "\t" <<  vertex->GetTrack(itrk)->MCEvt() << "\t" << vertex->GetVTa(itrk)->Imp() << "\t" << vertex->Flag() << endl;
            }
        }
        else FLAG99++;
    }
    
    cout << "FLAG99 " << FLAG99 << endl;
    
    cout << "COUNT_STRANGE " << COUNT_STRANGE << endl;
    
    return merged_varr;
}

//---------------------------------------------------------//

void MergeTrack(EdbTrackP* firsttrack, EdbTrackP* secondtrack){
    //first track diventerà firsttrack+second track.
    
    int Second_N = secondtrack->N();
    int Second_NF = secondtrack->NF();
    
    int First_Plate1 = firsttrack->GetSegmentFirst()->Plate();
    int First_Plate_last = firsttrack->GetSegmentLast()->Plate();
    int Second_Plate1 = secondtrack->GetSegmentFirst()->Plate();
    int Second_Plate_last = secondtrack->GetSegmentLast()->Plate();
    
    if(EVERBOSE==8 || (EVERBOSE==100 && ((firsttrack->MCEvt()==DEBUG_MCEVT || secondtrack->MCEvt()==DEBUG_MCEVT)||((firsttrack->Track()==DEBUG_TRKID || secondtrack->Track()==DEBUG_TRKID))))) cout << endl << "\tMerge tr: id " << firsttrack->Track()  << " mc " << firsttrack->MCEvt() << " " << firsttrack->MCTrack() << "\t" << secondtrack->MCEvt() << " " << secondtrack->MCTrack() << "\tPl " << First_Plate1 << "\t" << First_Plate_last << "\t" << Second_Plate1 << "\t" << Second_Plate_last << "\tFlag " << firsttrack->Flag() << "\t" << secondtrack->Flag() << "\tN " << firsttrack->N() << "\t" << Second_N << "\t"<< endl;
    
    //associo carica
    if(First_Plate_last>LASTLAYER[1] && First_Plate_last<=LASTLAYER[2] && Second_Plate_last>LASTLAYER[1] && Second_Plate_last<=LASTLAYER[2]){
        //if(secondtrack->N()>firsttrack->N() || firsttrack->Flag()>=9 || firsttrack->Flag()==0)
        if(secondtrack->N()>firsttrack->N() || firsttrack->Flag()>=9 || firsttrack->Flag()<=0) firsttrack->SetFlag(secondtrack->Flag());
        else firsttrack->SetFlag(firsttrack->Flag());
    }
    else if(Second_Plate_last>LASTLAYER[1] && Second_Plate_last<=LASTLAYER[2]){
        firsttrack->SetFlag(secondtrack->Flag());
    }
    
    for(int i=0; i<Second_N; i++){
        firsttrack->AddSegment(new EdbSegP(*secondtrack->GetSegment(i)));
        firsttrack->GetSegment(i)->SetDZ(300);
        // if(EVERBOSE==8 || (EVERBOSE==100 && ((MC==0&&firsttrack->Track()==DEBUG_TRKID) || firsttrack->MCEvt()==DEBUG_MCEVT || (MC==0&&secondtrack->Track()==DEBUG_TRKID) || secondtrack->MCEvt()==DEBUG_MCEVT))) cout << "\t\t\tadd seg: " << secondtrack->GetSegment(i)->Plate() << "\t" << secondtrack->GetSegment(i)->X() << "\t" << secondtrack->GetSegment(i)->Y() << "\t" << secondtrack->GetSegment(i)->TX() << "\t" << secondtrack->GetSegment(i)->TY() << endl;
    }
    for(int i=0; i<Second_NF; i++){
        firsttrack->AddSegmentF(new EdbSegP(*(secondtrack->GetSegmentF(i))));
        firsttrack->GetSegmentF(i)->SetDZ(300);
    }
    
    firsttrack->SetCounters();
    
    secondtrack->GetSegmentFirst()->SetProb(9);
    
    for(int i=0; i<firsttrack->N(); i++){
        
        if(EVERBOSE==8 || (EVERBOSE==100 && ((MC==0&&firsttrack->Track()==DEBUG_TRKID) || firsttrack->MCEvt()==DEBUG_MCEVT))) cout << "\t\tseg: " << firsttrack->GetSegment(i)->Plate() << "\t" << firsttrack->GetSegment(i)->X() << "\t" << firsttrack->GetSegment(i)->Y() << "\t" << firsttrack->GetSegment(i)->TX() << "\t" << firsttrack->GetSegment(i)->TY() << "\t" << firsttrack->GetSegment(i)->Theta() << endl;
        
    }
    
    if(EVERBOSE==8 || (EVERBOSE==100 && ((MC==0&&firsttrack->Track()==DEBUG_TRKID) || firsttrack->MCEvt()==DEBUG_MCEVT))) cout << "\tMerge trX: " << firsttrack->GetSegmentFirst()->Plate() << "\t" << firsttrack->GetSegmentLast()->Plate() << "\t" << firsttrack->N() << "\t" << firsttrack->Npl() << "\t" << firsttrack->Flag() << endl;
    
    
    //if(EVERBOSE==8) cout << "\tMerge tr: " << First_Plate1 << "\t" << Second_Plate1 << "\t" << firsttrack->N() << "\t" << Second_N << endl;
    
}

//---------------------------------------------------------//

void FillVtxPlate(TObjArray &arrv){
    for(int i=PLMIN; i<=BRAGGPLATE; i++) { TObjArray* newobj = new TObjArray(); vtxPat[i] = *newobj; }
    
    int nv = arrv.GetEntries();
    cout << "nv fill: " << nv << endl;
    for(int iv=0; iv<nv; iv++){
        EdbVertex *vertex = (EdbVertex*)(arrv.At(iv));
        int vID=vertex->ID();
        int mcevt = vertex->GetTrack(0)->MCEvt();
        if(FAST==4 && (vID!=DEBUG_VTXID&&mcevt!=DEBUG_MCEVT)) continue; //FAST 4
        if(vertex->Flag()==-99) cout << "\t FillVtxPlate 99 " << vID << endl;
        //cout << iv << "\t" << vID << endl;
        float vz = vertex->VZ();
        //FillClosestPat(v) -> ; //need to fill array of vertices with closest plate
        
        //Mod 10dec
//        int iclosestplate = -1;
//        float dzmin = 1000000000;
//        for(int i=PLMIN; i<=BRAGGPLATE; i++) {
//            float dz = TMath::Abs(vz - Z_LAYER[i]);
//            if (dz < dzmin){
//                dzmin = dz;
//                iclosestplate = i;
//            }
//        }//end of for loop over plates
        int iclosestplate = Get_vtx_plate(vertex->VZ()); //Mod 10 dec
        if(iclosestplate>= PLMIN && iclosestplate<=BRAGGPLATE) vtxPat[iclosestplate].Add(vertex);
        else vertex->SetFlag(-99);
        
        
        if(EVERBOSE==100 && ((MC==0&&vertex->GetTrack(0)->Track()==DEBUG_TRKID)||vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT || vID==DEBUG_VTXID)){
            cout << "\tvtx " << vID << " in fillvtxplate flag: " << vertex->Flag() << "\t" << vertex->VZ() << "\t" << Get_vtx_plate(vertex->VZ()) << endl;
            for(int itrk=0; itrk<vertex->N(); itrk++){
                cout << "\t\t" << vertex->GetTrack(itrk)->MCEvt() << "\t" << vertex->GetTrack(itrk)->Track() << "\t" << vertex->GetTrack(itrk)->Theta() << "\t" << vertex->GetVTa(itrk)->Imp() << "\t" << vertex->GetTrack(itrk)->Flag() << endl;
            }
        }
        
    }
}

//---------------------------------------------------------//


TObjArray *FindCloseTracks(TObjArray *varr, EdbVertexRec *vrec){
    int count=0;
    float maximp_dau_plate=0;

    TString oxy_found_filename;
    if (FAST==4) oxy_found_filename = Form("oxyfound_%d.root", FAST);
    else oxy_found_filename = Form("oxyfound.root");
    TFile* oxyfoundfile = new TFile(oxy_found_filename, "READ");
    TTree* oxyinfo = (TTree*)oxyfoundfile->Get("oxyinfo");
    if (oxyinfo->GetEntries()>0) oxyinfo->BuildIndex("vID");

    int OXYFOUND=0;
    if (oxyinfo->GetEntries()>0) oxyinfo->SetBranchAddress("OXY_FOUND", &OXYFOUND);


    TStopwatch t;
    t.Start();

    EdbVertex *newvertex;
    
    for(int ipl = PLMIN; ipl<=BRAGGPLATE; ipl++){ //BRAGGPLATE
        TObjArray varrpl = vtxPat[ipl];
        int nv = varrpl.GetEntries();
        if(EVERBOSE==10||EVERBOSE==100) cout << "\t >>> ipl " << ipl << "\tnv " << nv << endl;
        for(int iv=0; iv<nv; iv++){
            EdbVertex *vertex  = (EdbVertex*)(varrpl.At(iv));
            if (EVERBOSE==100 && FAST==4) cout << " Vertex Flag (Start FindCloseTracks) " << vertex->Flag() << " ID " << vertex->ID() << " N " << vertex->N() << endl;
            if(vertex->ID()==-99||vertex->Flag()==-99) continue;
            float vx = vertex->VX();
            float vy = vertex->VY();
            float vz = vertex->VZ();
            int vplate = Get_vtx_plate(vz);
            float xy[2] = {vx,vy};
            float r = 2000.;
            int idvtx= vertex->ID();
            if (oxyinfo->GetEntries()>0) {
                if (oxyinfo->GetEntryWithIndex(idvtx) ) OXY_FOUND[idvtx] = OXYFOUND;
                else OXY_FOUND[idvtx]=0;
            }
            if(vplate<1 && vplate>BRAGGPLATE) continue;
            
            newvertex= new EdbVertex();
            newvertex=vertex;
            
            //if(EVERBOSE==13 || (EVERBOSE==100 && newvertex->ID()==DEBUG_VTXID)) cout << "\t FindCT-A " << newvertex->ID() << "\t" << newvertex->Flag() << "\t" << ipl << "\t" << newvertex->N() << "\t" << newvertex->VZ() << endl;
            
            
            if(OXY_FOUND[idvtx]==0){
                EdbTrackP *found_cand = AnalyseCloseBeamTracks(vertex,ipl, maximp_beam, maxbeam_dx, maxbeam_dy);
                
                if(found_cand){
                    
                    newvertex = AddBeamToVertex(vertex, found_cand, vrec);
                    
                    if(!(newvertex)) continue;
                    if(EVERBOSE==100 && (newvertex->ID()==DEBUG_VTXID || vertex->ID()==DEBUG_VTXID)) cout << "\t FindCT-B " << newvertex->ID() << "\t" << newvertex->Flag() << "\t" << ipl << "\t" << newvertex->N() << "\t" << found_cand->MCTrack() << "\t" << newvertex->VZ() << endl;
                }
            }
            else if(OXY_FOUND[idvtx]>1){
                newvertex = SelectOneBeam(vertex, vrec);
                if(!(newvertex)) continue;
                if(EVERBOSE==100 && (newvertex->ID()==DEBUG_VTXID || vertex->ID()==DEBUG_VTXID)) cout << "\tFindCT remove oxy " << newvertex->ID() << "\t" << newvertex->N() << "\t" << newvertex->Flag() << endl;
                
                N_BEAM_REMOVED++;
                if(newvertex->N()==2){
                    EdbTrackP *firsttrack= newvertex->GetTrack(0);
                    EdbTrackP *secondtrack= newvertex->GetTrack(1);
                    
                    if(EVERBOSE==13 || (EVERBOSE==100 && (newvertex->ID()==DEBUG_VTXID || vertex->ID()==DEBUG_VTXID))) cout << "\tFindCT after rem " << firsttrack->MCTrack() << "\t" << firsttrack->N() <<  "\t" << secondtrack->MCTrack() << "\t" << secondtrack->N() << endl;
                    MergeTrack(firsttrack, secondtrack);
                    newvertex->SetFlag(-99);
                    if(EVERBOSE==13 || (EVERBOSE==100 && newvertex->ID()==DEBUG_VTXID)) cout << "vertex deleted rigo 695" << endl;
                    
                    if(EVERBOSE==13 || (EVERBOSE==100 && (newvertex->ID()==DEBUG_VTXID || vertex->ID()==DEBUG_VTXID))) cout << "\tFindCT-C " << newvertex->ID() << "\t" << newvertex->Flag() << "\t" << ipl << "\t" << newvertex->N() << "\t" << firsttrack->MCTrack() << "\t" << firsttrack->N() << "\t" << newvertex->VZ() << endl;
                    continue;
                }
            }
            
            if(EVERBOSE==13 || (EVERBOSE==100 && (newvertex->ID()==DEBUG_VTXID || vertex->ID()==DEBUG_VTXID))) cout << "\tFindCT-D " << newvertex->ID() << "\t" << ipl << "\t" << newvertex->N() << "\t" << newvertex->VZ() << " \t " << newvertex->Flag() << endl;
            
            TObjArray tr_grid_dau;
            int modifiedvtx=0;
            int tempid=newvertex->ID();
            int tempflag=newvertex->Flag();
            for(int iipl=ipl-3; iipl<ipl+4; iipl++){
                if(iipl>=BRAGGPLATE+1) break;
                if(iipl<0) continue;
                int ndau_grid = gridtr_DAU[iipl].SelectObjectsC( xy, r, tr_grid_dau);
                
                if(EVERBOSE==13 || (EVERBOSE==100 && (newvertex->ID()==DEBUG_VTXID || vertex->ID()==DEBUG_VTXID))) cout << "\tiipl: " << iipl << "\t" <<  ndau_grid << endl;
                if(iipl==ipl+1||iipl==ipl-1) maximp_dau_plate=maximp_dau/2;
                else maximp_dau_plate=maximp_dau;
                newvertex->SetFlag(-99);
                modifiedvtx += DaughtersSearch(newvertex, tr_grid_dau, vrec, maximp_dau_plate);
				if (EVERBOSE == 100 && modifiedvtx) cout << " modified vertex " << modifiedvtx<< endl;
                if(!(newvertex)) continue;
                tr_grid_dau.Clear();
            }
            newvertex->SetID(tempid);
            if(modifiedvtx==0) { newvertex->SetFlag(tempflag); if (EVERBOSE==100) cout << " newvertex flag = " << newvertex->Flag() << endl; }
            else if(modifiedvtx>0) { //newvertex->Flag()!=-99 &&
                newvertex->SetFlag(100);
				if (EVERBOSE == 100 && FAST==4) cout << " set flag 100 to newvertex " << endl;
                //vtxPat[ipl].Add(newvertex);
            }

			if (EVERBOSE == 100 && FAST==4) cout << " newvertex flag " << newvertex->Flag() << "\t" << " vertex ID " <<  newvertex->ID()  << " vertex N " <<  newvertex->N() << "\t" << endl;
            
            //vertex->SetFlag(newvertex->Flag()); 
            //vtxPat[ipl][iv] = newvertex;
            if (EVERBOSE==100 && FAST==4) { TObjArray varrplNEW = vtxPat[8];  EdbVertex *vertexNEW  = (EdbVertex*)(varrplNEW.At(iv)); cout << " Vertex NEW Flag " << vertexNEW->Flag() << endl;}
            //{if(newvertex->Flag()!=-99 && modifiedvtx==1) new_varrpl->Add(newvertex);
            
            if(EVERBOSE==13 || (EVERBOSE==100 && (newvertex->ID()==DEBUG_VTXID || vertex->ID()==DEBUG_VTXID))) {
                cout << "\tFindCT-E " << newvertex->ID() << "\t" << newvertex->Flag() << "\t" << ipl << "\t" << newvertex->N() << "\t" << newvertex->VZ() << endl;
                for (int i=0; i<newvertex->N(); i++){
                    cout << "\t\t\t" << newvertex->GetTrack(i)->Track() << endl;
                }
            }
            
        } //end for vertices
        
        if ((BRAGGPLATE-ipl)%1 == 0 && (EVERBOSE==14||EVERBOSE==16||EVERBOSE==101)) {
            cout << " FindCloseTracks: " << ipl << "\t" << 100.*(ipl)/(BRAGGPLATE-PLMIN+1) << " %, Iteration Time: " << t.RealTime() << " s\t" << t.RealTime()/60 << " min " <<  endl;
            t.Reset();
            t.Start();
        }
        
    } //end for plates

    //if (vertex && newvertex) if(EVERBOSE==13 || (EVERBOSE==100 && (newvertex->ID()==DEBUG_VTXID || vertex->ID()==DEBUG_VTXID))) cout << "final array" << endl;
    
    TObjArray *new_varr_def = new TObjArray(); // Elimino i vertici con flag -99
    
//    for (int ivtx = 0; ivtx < varr->GetEntries(); ivtx++){
//        EdbVertex *vertex = (EdbVertex*)(varr->At(ivtx));
    for(int ipl = PLMIN; ipl<=BRAGGPLATE; ipl++){
        TObjArray varrpl = vtxPat[ipl];
        int nv = varrpl.GetEntries();
        int count_newvtx=0;
        //if(EVERBOSE==10)
            cout << "\t >>> ipl " << ipl << "\tnv " << nv << "\t";
        for(int iv=0; iv<nv; iv++){
            EdbVertex *vertex  = (EdbVertex*)(varrpl.At(iv));
			if (EVERBOSE == 100) cout << " vertex flag " << vertex->Flag() << "\t" << " vertex ID " <<  vertex->ID()  << " vertex N " <<  vertex->N() << "\t";
            if(vertex->Flag()!=-99) {
                new_varr_def->Add(vertex);
                count_newvtx++;
                if(EVERBOSE==13 || (EVERBOSE==100 && (newvertex->ID()==DEBUG_VTXID || vertex->ID()==DEBUG_VTXID))){
                    cout << vertex->ID() << "\t" << vertex->Flag() << endl;
                    for (int iii=0; iii<vertex->N(); iii++) cout << "\t\t" << vertex->GetTrack(iii)->ID() << "\t" << vertex->GetTrack(iii)->Track() << "\t" << vertex->GetTrack(iii)->TX() << "\t" << vertex->GetTrack(iii)->N() << endl;
                }
            }
        }
        //if(EVERBOSE==10)
            cout << count_newvtx << endl;
    }
    
    oxyfoundfile->Close();

    return new_varr_def;
}


//---------------------------------------------------------//

EdbTrackP *AnalyseCloseBeamTracks(EdbVertex *vertex, int ipl, float maximp, float max_dx, float max_dy){
    int save_itrk=0;
    float minimp = maximp;
	float min_dx = max_dx;
	float min_dy = max_dy;
    
    
    int vID=vertex->ID();
    float vx = vertex->VX();
    float vy = vertex->VY();
    float vz = vertex->VZ();
    
    float xy[2] = {vx,vy};
    float r = 2000.;
    TObjArray tracks;
    if(EVERBOSE==1 || (EVERBOSE==100 && vID==DEBUG_VTXID)) cout <<"coord vtx: " << xy[0] << "\t" << xy[1] << "\tPl: " << ipl << endl;
    
    int iplmin=ipl-3;
    if(BRICKID==333) iplmin=ipl-4;
    for(int iipl=ipl; iipl>iplmin; iipl--){
        if(iipl<=0) break;
        int nbeam_grid = gridseg_OXY[ipl].SelectObjectsC( xy, r, tracks);
        if(EVERBOSE==1 || (EVERBOSE==100 && vID==DEBUG_VTXID)) cout << "\t nbeam_grid: " << nbeam_grid << endl;
        
        for (int itrk = 0; itrk < tracks.GetEntries(); itrk++){
            EdbTrackP *cand = (EdbTrackP*) tracks.At(itrk);
            EdbVertex *vtempS = cand->VertexS();
            EdbVertex *vtempE = cand->VertexE();
            
            if(!(vtempS)&&(!(vtempE))){//condition A: la traccia non deve avere altri vertici entranti o uscenti
                //  if((!(vtempS)||(vtempS && vtempS->ID()!=vID))&&(!(vtempE)) && (cand->Flag()!=8&&cand->Flag()!=9)){//condition A: la traccia non deve avere altri vertici entranti o uscenti o, se ce li ha, devono avere id diverso dal vertice in esame. Inoltre non deve essere già stata associata ad altri vertici nuovi.
                if ((cand->GetSegmentFirst()->Z() < vz) && (cand->N()>4||(cand->GetSegmentLast()->Plate()<=4 && cand->N()==cand->GetSegmentLast()->Plate()))){ //beam candidate should start before vz and have at least 5 segments or if it starts before plate 5 should have all segments
                    //transverse IP to vertex
                    int iseg=FindPlate(cand, iipl);
                    if(iseg<NPLATE){
                        float dist = CalcDist3D(vertex, cand->GetSegment(iseg)); //calcolo il parametro d'impatto tra il vertice e il segmento della traccia in esame
						float distX = CalcDist2DX(vertex, cand->GetSegment(iseg));
						float distY = CalcDist2DY(vertex, cand->GetSegment(iseg));
                        if(EVERBOSE==100 && (vID==DEBUG_VTXID||(cand->ID() == DEBUG_TRKID) || cand->MCEvt()==DEBUG_MCEVT)) {
                            cout << "Search Beam: " << cand->Track() << "\t" << cand->N() << "\t" << cand->GetSegmentFirst()->Z() << "\t" << vz << endl;
                            cout << "\t" << iipl << "\t" <<  cand->GetSegment(iseg)->Plate() << "\t " << cand->GetSegment(iseg)->Z() << "\t " << cand->GetSegment(iseg)->Theta() << "\tflag " << cand->Flag() << "\tdist: " << dist << " min " << minimp << " distX "  << distX << " distY " << distY << endl;
                        }
                        if(dist < minimp  && distX<min_dx && distY<min_dy){ //found candidate, storing track and new min IP
                            minimp = dist;
                            save_itrk = itrk;
                            //cout << "\t\tdist: " << dist << "\tcoord vtx: " << xy[0] << "\t" << xy[1] << "\t" << vz << "\tPl: " << ipl << "\t" << cand->TY() << endl;
                            //save_iplate = iplate;
                        } //imp condition
                    } //if (tr->GetSegmentFirst()->Z() < vertex->VZ()){
                } //condition A
            }
        }
        
        if (minimp < maximp){
            EdbTrackP *found_candidate = (EdbTrackP*) tracks.At(save_itrk);
            if(found_candidate->Flag()==0) found_candidate->SetFlag(8); //con questa flag evito che la stessa traccia sia associata ad altri vertici
            if(EVERBOSE==100 && (vertex->ID()==DEBUG_VTXID)) cout << "found cand: " << found_candidate->ID() << "\t" << found_candidate->GetSegmentFirst()->ID() << "\t" << found_candidate->MCEvt() << "\t" << found_candidate->MCTrack() << "\t" << found_candidate->N() << "\timp " << minimp << endl;
            // if(EVERBOSE==2) cout <<  "\tcandidate found " << found_candidate->ID() << "\t" << found_candidate->MCEvt() << "\t" << found_candidate->MCTrack() << "\tIP " << minimp << "\t" << found_candidate->Theta() <<  "\t" << found_candidate->N() << endl;
            
            return found_candidate;
        }
    }
    return NULL;
    
}
//---------------------------------------------------------//

TObjArray* UnisciTracce(TObjArray *varr, float maximp, float max_deltatheta){
    
    TStopwatch t;
    t.Start();
    TObjArray *end_tracksS1 = new TObjArray();
    EdbTrackP *track;
    for (int ivtx = 0; ivtx < varr->GetEntries(); ivtx++){
        EdbVertex *vertex = (EdbVertex*)(varr->At(ivtx));
        if((EVERBOSE==100||EVERBOSE==15) && vertex->ID()==DEBUG_VTXID) cout << vertex->ID() << "\tbefore unisci: " << vertex->N() << endl;
        float vx = vertex->VX();
        float vy = vertex->VY();
        float vz = vertex->VZ();
        int vplate=Get_vtx_plate(vz);
        if(vplate<0) vplate=0;
        float xy[2] = {0,0};
        float r = 2000.;
        
        TObjArray tr_grid;
        
        for (int itrk = 0; itrk < vertex->N(); itrk++){
            track = vertex->GetTrack(itrk);
            EdbVertex *vtempE = track->VertexE();

            if (EVERBOSE==15) cout << "Track # " << itrk << " with s0 ID " << track->GetSegmentFirst()->ID() << endl;
            if (EVERBOSE == 15)
            {
                if (vtempE) cout << " Track VertexE is not NULL (ID= " << vtempE->ID() << ")" << endl;
                EdbVertex *vtempS = track->VertexS();
                if (vtempS) cout << " Track Vertex S is not NULL (ID= " << vtempS->ID() << ")" << endl;
            }
            // if(EVERBOSE==7 || (EVERBOSE==100 && ((trkend->Track()==DEBUG_TRKID)|| trkend->MCEvt()==DEBUG_MCEVT)))  cout << "\t" << cand->MCEvt() << "\t" << cand->MCTrack() << "\tnseg " << cand->N() << "\t" << cand->GetSegmentFirst()->Plate() << "\t" << cand->GetSegmentLast()->Plate() << "\t" << endl;
            //if(EVERBOSE==7) cout << "\t" << itrk << "\t" << track->Track() <<"\t" << track->MCEvt() << "\t" << track->MCTrack() << endl;
            if(EVERBOSE==7 || (EVERBOSE==100 && (vertex->ID()==DEBUG_VTXID || track->Track()==DEBUG_TRKID || track->MCEvt()==DEBUG_MCEVT))) cout << "E check\t" << itrk << "\t" << track->Track() << "\t" << track->Theta() << "\t" << vertex->GetVTa(itrk)->Zpos() << " (vtx id: " << vertex->ID() << ")" << "\t" << track->MCEvt() << "\t" << track->MCTrack() << endl;
            if(!(vtempE)){//condition A: la seconda traccia non deve entrare in un vertice
                if(vertex->GetVTa(itrk)->Zpos()==1){ //lo faccio solo per le tracce che escono dal vertice!
                    if(EVERBOSE==7 || (EVERBOSE==100 && (vertex->ID()==DEBUG_VTXID || track->Track()==DEBUG_TRKID || track->MCEvt()==DEBUG_MCEVT))) cout  << "Extend Track id " << track->ID() << "\ttr " << track->Track() << "\ttheta " << track->Theta()<< "\tincoming " << vertex->GetVTa(itrk)->Zpos() << "\tmc " << track->MCEvt() << "\t" << itrk+1 << "/" << vertex->N() <<"\tnseg " << track->N() << "\tnpl " << track->Npl() << "\tPl " << track->GetSegmentFirst()->Plate() << " " << track->GetSegmentLast()->Plate() << endl;
                    if (EVERBOSE==15) cout << "Extending Track # " << itrk << endl;
                    int lastplate=track->GetSegmentLast()->Plate();
                    int lastplate0=lastplate;
                    int end=-99, found=1, ngapmax=-99, justonce=0;
                    while(lastplate<PLMAX && found==1){
                        EdbSegP *lastseg=track->GetSegmentLast();
                        if(lastseg->Plate()>LASTLAYER[1]-4&&lastseg->Plate()<=LASTLAYER[2]) ngapmax=9;
                        else ngapmax=4;
                        if (EVERBOSE==15 && track->GetSegmentFirst()->Plate()==DEBUG_S0_PLATE_END && track->GetSegmentFirst()->ID()==DEBUG_S0_ID_END && lastseg->Plate()!=lastplate0 ) cout << " last seg plate  " << lastseg->Plate() << endl;
                        for(int iipl=lastplate+1; iipl<=lastplate+ngapmax; iipl++){
                            if(iipl>PLMAX) break;
                            xy[0] = lastseg->X()-lastseg->TX()*(lastseg->Z()-Z_LAYER[iipl]);
                            xy[1] = lastseg->Y()-lastseg->TY()*(lastseg->Z()-Z_LAYER[iipl]);
                            int n_grid = gridtr_DAU[iipl].SelectObjectsC(xy, r, tr_grid);
                            if(EVERBOSE==100 && (vertex->ID()==DEBUG_VTXID || track->Track()==DEBUG_TRKID || track->MCEvt()==DEBUG_MCEVT)) cout << "\t searching in plate " << iipl << " cell " << "\t" << xy[0] << " " << xy[1] << "\t between " << n_grid << " tracks" << endl;
                            
                            if (EVERBOSE==15 && track->GetSegmentFirst()->Plate()==DEBUG_S0_PLATE_END && track->GetSegmentFirst()->ID()==DEBUG_S0_ID_END ) {
                                cout << " Checking Grid Plate " << iipl << endl;
                                cout << " Ngapmax " << ngapmax << endl;
                                if (iipl == DEBUG_S0_PLATE_ST) {
                                    cout << " Selecting Objects from gridtr_DAU with xy[0] " << xy[0] << " xy[1] " << xy[1] << endl;
                                    cout << " Selected a total of " << tr_grid.GetEntries() << " candidates starting in plate 31 " << endl;
                                    cout << " Ngrid " << n_grid << endl;
                                    for (int itrk = 0; itrk < tr_grid.GetEntries(); itrk++){
                                        EdbTrackP *cand = (EdbTrackP*) tr_grid.At(itrk);
                                        if (cand->GetSegmentFirst()->ID()==DEBUG_S0_ID_ST) cout << " Got Candidate from Tr_Grid " << endl;
                                        //else {cout << " X " << cand->GetSegmentFirst()->X() << " Y " << cand->GetSegmentFirst()->Y() << endl;
                                        //     cout << " S0_plate " << cand->GetSegmentFirst()->Plate() << " S0_ID " << cand->GetSegmentFirst()->ID() << endl;}
                                    }
                                }
                            }
                            if(n_grid>0) found = ExtendTrack(track, tr_grid, maximp, max_deltatheta+0.01*(iipl-lastplate-1));
                            else {found=0; if (ALLOW_LONGER_GAPS == 1 && ngapmax==4) end_tracksS1->Add(track); }// allow ngapmax to be slightly larger if no candidate was found
                            if(found == 1) {
                                lastplate=track->GetSegmentLast()->Plate();
                                if(EVERBOSE==100 && (vertex->ID()==DEBUG_VTXID || track->Track()==DEBUG_TRKID || track->MCEvt()==DEBUG_MCEVT)){ cout << "\t\t track extended up to plate " << lastplate << "\tnseg new " << track->N() << endl;}
                                break;
                            }
                            tr_grid.Clear();
                        }
                    }
                }
            }
        }
        if(EVERBOSE==7 || (EVERBOSE==100 && (vertex->ID()==DEBUG_VTXID || track->Track()==DEBUG_TRKID || track->MCEvt()==DEBUG_MCEVT))) cout << "\tafter unisci: " << vertex->N() << endl;
        
        if (ivtx%1000 == 0 && (EVERBOSE==14||EVERBOSE==16||EVERBOSE==101)) {
            cout << " Unisci Tracce: " << 100.*(ivtx)/varr->GetEntries() << " %, Iteration Time: " << t.RealTime() << " s\t" << t.RealTime()/60 << " min " << endl;
            t.Reset();
            t.Start();
        }
    }
    
    return end_tracksS1;
}

//---------------------------------------------------------//


int ExtendTrack(EdbTrackP *trkend, TObjArray &tracks, float maximp, float max_deltatheta){
    
    float dist=maximp;
    float deltatheta=max_deltatheta, deltatheta_x=max_deltatheta, deltatheta_y=max_deltatheta;
    float temp_dist=1000, temp_deltatheta=1000, temp_deltatheta_x=1000, temp_deltatheta_y=1000, temp_x=1000, temp_y=1000, temp_dz=1000, temp_delta_x=1000, temp_delta_y=1000;
    int iitrk=-99, c=0;
    
    for (int itrk = 0; itrk < tracks.GetEntries(); itrk++){
        EdbTrackP *cand = (EdbTrackP*) tracks.At(itrk);
        EdbVertex *vtempS = cand->VertexS();
        
        if(cand->GetSegmentFirst()->Prob()==9) continue;
        
        if (EVERBOSE==15 && cand->GetSegmentFirst()->ID()==DEBUG_S0_ID_ST && cand->GetSegmentFirst()->Plate()==DEBUG_S0_PLATE_ST && trkend->GetSegmentFirst()->ID()==DEBUG_S0_ID_END && trkend->GetSegmentFirst()->Plate()==DEBUG_S0_PLATE_END) {
         bool cond = vtempS;
         cout << " cand vtempS " << cond << endl;
         cout << " Z new " << cand->GetSegmentFirst()->Z() << " Z old " << trkend->GetSegmentLast()->Z() << endl;
         }
        
        // if(EVERBOSE==7 || (EVERBOSE==100 && ((trkend->Track()==DEBUG_TRKID)|| trkend->MCEvt()==DEBUG_MCEVT)))  cout << "\t" << cand->MCEvt() << "\t" << cand->MCTrack() << "\tnseg " << cand->N() << "\t" << cand->GetSegmentFirst()->Plate() << "\t" << cand->GetSegmentLast()->Plate() << "\t" << endl;
        if(!(vtempS)){//condition A: la seconda traccia non deve uscire da un vertici
            temp_dist=1000; temp_deltatheta=1000; temp_deltatheta_x=1000; temp_deltatheta_y=1000; temp_x=1000; temp_y=1000; temp_dz=1000; temp_delta_x=1000; temp_delta_y=1000;
            
            if (trkend->GetSegmentLast()->Z() < cand->GetSegmentFirst()->Z()){ //next segment should start after end of first one
                
                //transverse IP to vertex
                temp_dist = CalcDist(trkend->GetSegmentLast(),cand->GetSegmentFirst());
                temp_deltatheta = trkend->GetSegmentLast()->Theta()-cand->GetSegmentFirst()->Theta();
                temp_deltatheta_x = trkend->GetSegmentLast()->TX()-cand->GetSegmentFirst()->TX();
                temp_deltatheta_y = trkend->GetSegmentLast()->TY()-cand->GetSegmentFirst()->TY();
                temp_dz=trkend->GetSegmentLast()->Z()-cand->GetSegmentFirst()->Z();
                temp_x=trkend->GetSegmentLast()->X()-temp_dz*trkend->GetSegmentLast()->TX();
                temp_y=trkend->GetSegmentLast()->Y()-temp_dz*trkend->GetSegmentLast()->TY();
                temp_delta_x=temp_x-cand->GetSegmentFirst()->X();
                temp_delta_y=temp_y-cand->GetSegmentFirst()->Y();
                //if(TMath::Abs(temp_delta_x)<50 && TMath::Abs(temp_delta_y)<50) cout << "\t\t" << temp_dist << "\t" << temp_deltatheta << "\t" << dist << "\t" << deltatheta << "\t" << temp_x << "\t" << temp_delta_x << "\t" << temp_y << "\t" << temp_delta_y << endl;
                
                if(cand->GetSegmentFirst()->Plate()>LASTLAYER[1]){
                    h_IP_merge->Fill(temp_dist);
                    h_DXY_merge->Fill(temp_delta_x, temp_delta_y);
                    h_DTXTY_merge->Fill(temp_deltatheta_x,temp_deltatheta_y);
                    h_IP_DTh_merge->Fill(temp_dist, sqrt(temp_deltatheta_x*temp_deltatheta_x+temp_deltatheta_y*temp_deltatheta_y));
                }
                
                if (EVERBOSE==15 && cand->GetSegmentFirst()->ID()==DEBUG_S0_ID_ST && cand->GetSegmentFirst()->Plate()==DEBUG_S0_PLATE_ST && cand->GetSegmentFirst()->ID()==DEBUG_S0_ID_ST && cand->GetSegmentFirst()->Plate()==DEBUG_S0_PLATE_ST && trkend->GetSegmentFirst()->ID()==DEBUG_S0_ID_END && trkend->GetSegmentFirst()->Plate()==DEBUG_S0_PLATE_END) {
                 
                 cout << " Extend Track Calculated Parameters " << endl;
                 cout << " b " << temp_dist << " DTX_temp " << temp_deltatheta_x << " DTY_temp " << temp_deltatheta_y  << " DX temp " << temp_delta_x << " DY temp " << temp_delta_y << endl;
                 
                 }
                
                
                if(TMath::Abs(temp_deltatheta_x)<max_deltatheta && TMath::Abs(temp_deltatheta_y)<max_deltatheta && TMath::Abs(temp_delta_x)<maximp && TMath::Abs(temp_delta_y)<maximp && temp_dist<dist && temp_deltatheta<max_deltatheta){ //temp_dist<dist &&
                    deltatheta=temp_deltatheta;
                    deltatheta_x=TMath::Abs(temp_deltatheta_x);
                    deltatheta_y=TMath::Abs(temp_deltatheta_y);
                    
                    dist=temp_dist;
                    iitrk=itrk;
                    if (EVERBOSE==15 && cand->GetSegmentFirst()->ID()==DEBUG_S0_ID_ST && cand->GetSegmentFirst()->Plate()==DEBUG_S0_PLATE_ST && trkend->GetSegmentFirst()->ID()==DEBUG_S0_ID_END && trkend->GetSegmentFirst()->Plate()==DEBUG_S0_PLATE_END) {
                     
                     cout << " First Check OK -> b = " << dist << endl;
                     }
                    
                    if(EVERBOSE==7) cout << "\t\tpar " << temp_dist << "\t" << temp_deltatheta << "\t" << deltatheta_x << "\t" << deltatheta_y <<  "\t" << temp_delta_x << "\t" << temp_delta_y << endl;
                    
                    if((EVERBOSE==7 && MC==11) || (EVERBOSE==100 && ((cand->Track() == DEBUG_TRKID||trkend->Track()==DEBUG_TRKID) || cand->MCEvt()==DEBUG_MCEVT|| trkend->MCEvt()==DEBUG_MCEVT))) cout << "\t\t\t\t" << cand->MCEvt() << "\t" << cand->MCTrack() << " selected " << dist << "\t" << deltatheta << endl;
                    
                } else if (EVERBOSE==15 && cand->GetSegmentFirst()->ID()==DEBUG_S0_ID_ST && cand->GetSegmentFirst()->Plate()==DEBUG_S0_PLATE_ST && trkend->GetSegmentFirst()->ID()==DEBUG_S0_ID_END && trkend->GetSegmentFirst()->Plate()==DEBUG_S0_PLATE_END){
                  
                  cout << " First Check NOT OK " << endl;
                  cout << " Used These Parameters: deltatheta_x " << deltatheta_x << " deltatheta_y " << deltatheta_y << " dist " << dist << endl;
                  
                  }
                
            } //if (tr->GetSegmentFirst()->Z() < vertex->VZ()){
        } //condition A
        
    } //for
    
    if(dist<maximp && deltatheta<max_deltatheta){
        EdbTrackP *cand = (EdbTrackP*) tracks.At(iitrk);
        
        if(EVERBOSE==8 || (EVERBOSE==100 && ((trkend->Track()==DEBUG_TRKID)|| trkend->MCEvt()==DEBUG_MCEVT))) {
            cout << endl << endl << "A: " << trkend->Track() << " " << trkend->GetSegmentFirst()->ID() << " Pl " << trkend->GetSegmentLast()->Plate() << " MC " << trkend->GetSegmentLast()->MCEvt() << " " << trkend->GetSegmentLast()->MCTrack() << "\t\tB: " << cand->Track() << " " << cand->GetSegmentFirst()->ID() << " Pl " << cand->GetSegmentFirst()->Plate() << " MC " << cand->GetSegmentFirst()->MCEvt() << " " << cand->GetSegmentLast()->MCTrack() << "\t\t" << dist << "\t" << deltatheta << endl;
        }
        if(MC==11) {
            if(trkend->GetSegmentLast()->MCEvt()==cand->GetSegmentFirst()->MCEvt() && trkend->GetSegmentLast()->MCTrack()==cand->GetSegmentFirst()->MCTrack()){
                EXTEND_TRACK_OK++;
                if(EVERBOSE==7) cout << "\tok\t";
            }
        }
        if(EVERBOSE==7 || (EVERBOSE==100 && ((trkend->Track()==DEBUG_TRKID)|| trkend->MCEvt()==DEBUG_MCEVT))) cout << "before: nseg " << trkend->N() << " first Pl " << trkend->GetSegmentFirst()->Plate() << " last Pl " << trkend->GetSegmentLast()->Plate() << " Npl " << trkend->Npl() << endl;
        
        MergeTrack(trkend, cand);
        
        if(cand->GetSegmentFirst()->Plate()==DEBUG_S0_PLATE_ST && cand->GetSegmentFirst()->ID()==DEBUG_S0_ID_ST && EVERBOSE==15) {
            cout << " Candidate in S2 has been merged to track with s0_plate " << trkend->GetSegmentFirst()->Plate() << " s0_id " << trkend->GetSegmentFirst()->ID() << endl;
        }
        
        if(EVERBOSE==7 || (EVERBOSE==100 && ((trkend->Track()==DEBUG_TRKID)|| trkend->MCEvt()==DEBUG_MCEVT))) {
            cout << endl << "new: nseg " << trkend->N() << " first Pl " << trkend->GetSegmentFirst()->Plate() << " last Pl " << trkend->GetSegmentLast()->Plate() << " Npl " << trkend->Npl() << "\t";
        }
        EXTEND_TRACK++;
        return 1;
        
    }
    else {
        return 0;
    }
    
}

//---------------------------------------------------------//

void UnisciTracceLongerGaps(TObjArray *end_tracksS1, float maximp, float max_deltatheta) {
    
    if (EVERBOSE==14||EVERBOSE==16||EVERBOSE==15) cout << " Found " << end_tracksS1->GetEntries() << " tracks possibly having distant candidates " << endl;
    TStopwatch t2;
    t2.Start();
    
    float xy[2] = {0,0};
    float r = 2000.;
    TObjArray tr_grid;
    
    for (int itrk2=0; itrk2<end_tracksS1->GetEntries(); itrk2++) {
        
        EdbTrackP* track = (EdbTrackP*)end_tracksS1->At(itrk2);
        
        int lastplate=track->GetSegmentLast()->Plate();
        int found=1, ngapmax=6;
        //cout << " Looking to extend track # " << itrk2 << ", with lastplate " << lastplate << endl;
        //cout << " track s[0] ID " << track->GetSegmentFirst()->ID() << " s[0] Plate " << track->GetSegmentFirst()->Plate() << endl;
        while(lastplate<PLMAX && found==1){
            //cout << " while " << c<< endl;
            EdbSegP *lastseg = track->GetSegmentLast();
            //cout << " lastplate " << lastplate << endl;
            for(int iipl=lastplate+5; iipl<=lastplate+ngapmax; iipl++){
                if(iipl>PLMAX) {found=0; break;}
                xy[0] = lastseg->X()-lastseg->TX()*(lastseg->Z()-Z_LAYER[iipl]);
                xy[1] = lastseg->Y()-lastseg->TY()*(lastseg->Z()-Z_LAYER[iipl]);
                int n_grid = gridtr_DAU[iipl].SelectObjectsC(xy, r, tr_grid);
                //cout << " plate " << iipl << " found " << found << endl;
                if(n_grid>0) found = ExtendTrack(track, tr_grid, 0.5*maximp, 0.5*max_deltatheta+0.01*(iipl-lastplate-1));
                else found=0;
                //cout << " Found Ngrid = " << n_grid << " candidates in plate " << iipl  << " found " << found << endl;
                if(found == 1) {
                    lastplate=track->GetSegmentLast()->Plate();
                    break;
                }
                tr_grid.Clear();
            }
        }
        if (itrk2%2000==0 && (EVERBOSE==14||EVERBOSE==16)) {
            cout << " UnisciTracceLongGaps: " << float(100.*itrk2/end_tracksS1->GetEntries()) << " %,  Iteration Time: " << t2.RealTime() << " s" << endl;
            t2.Reset();
            t2.Start();
        }
    }
    
    
    
}

//---------------------------------------------------------//

int DaughtersSearch( EdbVertex *vertex, TObjArray &tracks, EdbVertexRec* vrec, float maximp){
    
    int save_itrk=0;
    int vnew=0;
    float minimp = maximp;
    //cout <<"start dau\t" <<  vertex.N() << endl;
    int vID=vertex->ID();
    int vtemp_ID=-99, vtemp_PL=-99;
    
    float vx = vertex->VX();
    float vy = vertex->VY();
    float vz = vertex->VZ();
    vertex->ResetTracks();
    
    if(EVERBOSE==100 && vID==DEBUG_VTXID) cout << "Searching dau for vertex " << vID << "\tVZ " << vz << endl;
    
    for (int itrk = 0; itrk < tracks.GetEntries(); itrk++){
        EdbTrackP *cand = (EdbTrackP*) tracks.At(itrk);
        EdbVertex *vtemp=NULL;
        EdbSegP* cand_seg;
        vtemp_ID=-999;
        vtemp_PL=-999;
        int readd_vtx=0;
        
        if (cand->GetSegmentFirst()->Z() > vz) {
            cand_seg=cand->GetSegmentFirst(); //dau candidate should start after vz
            if(cand->VertexS()) {
                vtemp_ID=cand->VertexS()->ID();
                vtemp_PL = Get_vtx_plate(cand->VertexS()->VZ());
            }
        }
        else {
            cand_seg=cand->GetSegmentLast();
            if(cand_seg->Z() > vz) continue;
            if(cand->VertexE()) {
                vtemp_ID = cand->VertexE()->ID();
                vtemp_PL = Get_vtx_plate(cand->VertexE()->VZ());
            }
        }
        
        if(vtemp_ID>=0 && vtemp_ID!=vID && vtemp_PL>0 && vtemp_PL<BRAGGPLATE){
            //            int iipl=0;
            //            int i1=-1;
            //            int ij1=0;
            //            for(int ij=0; ij<=10; ij++){
            //                i1=i1*(-1);
            //                if(ij%2==1 && ij>=1) ij1++;
            //                iipl=cand_seg->Plate()+i1*ij1;
            //                if(iipl<0||iipl>BRAGGPLATE) continue;
            //  TObjArray varrpl2 = vtxPat[iipl];
            TObjArray varrpl2 = vtxPat[vtemp_PL];
            for(int ivtx1=0; ivtx1<varrpl2.GetEntries(); ivtx1++){
                EdbVertex *vtempsearch = (EdbVertex*)(varrpl2.At(ivtx1));
                if(vtemp_ID == vtempsearch->ID() && vtempsearch->Flag()!=-99){
                    vtemp=vtempsearch;
                    if(EVERBOSE==100 && (cand->Track() == DEBUG_TRKID || cand->MCEvt()==DEBUG_MCEVT || vID==DEBUG_VTXID || vtemp_ID==DEBUG_VTXID)){ cout << "cand " << cand->Track() << "\tflag " << cand->Flag() << "\tvtemp found " << vtemp->ID() << endl;}
                    break;
                }
            }
        }
        //}
        
        //        if(vtemp){
        //            if(EVERBOSE==100 && ((cand->Track() == DEBUG_TRKID) || cand->MCEvt()==DEBUG_MCEVT)) cout << "\t vtemp found " << vtemp->Track() << endl;
        //            if(vtemp==vID) continue;
        //        }
        
        /*  if(EVERBOSE==100 && ((cand->Track() == DEBUG_TRKID) || cand->MCEvt()==DEBUG_MCEVT || vID==DEBUG_VTXID)) {
         cout << "Dau search -> cand: " << cand->MCEvt() << " " << cand->MCTrack() << "\t" << cand->Track() << "\t" << cand->Flag() << "\tPl " << cand->GetSegmentFirst()->Plate() << "\tnseg " << cand->N() << "\t";
         if(vtemp) cout << "vtx: " << vtemp->Track() << "\t" << vtemp->N() << "\t";
         cout << endl;
         }*/
        
        //transverse IP to vertex
        float imp =  CalcIP(cand_seg, vx,vy,vz); //calcolo il parametro d'impatto tra il segmento del dau e il vertice
        if(EVERBOSE==100 && vID==DEBUG_VTXID)  cout << "\t" << cand->MCEvt() << "\t" << cand->MCTrack() << "\t" << cand->Track() << "\t" << cand->GetSegmentFirst()->Plate() << "-" << cand->GetSegmentLast()->Plate() << "\t" << cand->N() << "\t" << imp << "\t" << CalcDist(vertex,cand) << "\t" << CalcDist3D(vertex,cand) << endl;
        if(EVERBOSE==100 && ((cand->Track() == DEBUG_TRKID) || cand->MCEvt()==DEBUG_MCEVT)){ cout << "\t\tSearch Dau: imp " << imp << "\t";
            if(vtemp) cout << "vtx: " << vtemp->ID() << "\t" << vtemp->N();
            cout << endl;
        }
        if(imp < maximp_dau && (cand->N()>=2&&(cand->N0()<=1||cand->N()>3))){ //found candidate, storing track and new min IP
            int alreadybusy=-99, trackfound=-99;
            
            if(vtemp) {
                if(vtemp->ID()!=vID){
                    for(int its=0; its<vtemp->N(); its++){
                        EdbTrackP *temptrack= vtemp->GetTrack(its);
                        if(cand->GetSegmentFirst()->ID()==temptrack->GetSegmentFirst()->ID()){
                            trackfound=its;
                            if(EVERBOSE==100 && ((cand->Track() == DEBUG_TRKID) || cand->MCEvt()==DEBUG_MCEVT||vID==DEBUG_VTXID||vtemp->ID() ==DEBUG_VTXID)){ cout << "track already belonging to vertex " <<  vtemp->ID() << "\timp " << vtemp->GetVTa(its)->Imp() << "\timp new vtx " << imp << endl;}
                            if(vtemp->GetVTa(its)->Imp()<=imp) {
                                alreadybusy=its;
                                break;
                            }
                        }
                    }
                }
            }
            if(EVERBOSE==100 && (cand->Track() == DEBUG_TRKID || cand->MCEvt()==DEBUG_MCEVT||vID==DEBUG_VTXID)) {
                if(vtemp) {
                    cout << "\t\t" << vtemp->ID() << " becomes N " << vtemp->N() << "\tFlag " << vtemp->Flag() << endl;
                    for(int kkk=0; kkk<vtemp->N(); kkk++) {
                        cout << vtemp->GetTrack(kkk)->Track() << " ";
                        if (vtemp->GetTrack(kkk)->VertexS()->ID()) cout << vtemp->GetTrack(kkk)->VertexS()->ID();
                        cout << endl;
                    }
                }
            }
            
            if(alreadybusy==-99){ //Se la traccia è associata a un vertice peggiore
                //rimuovo la traccia da quel vertice
                if(vtemp){
                    if(trackfound!=-99){
                        if(EVERBOSE==100 && ((cand->Track() == DEBUG_TRKID) || cand->MCEvt()==DEBUG_MCEVT||vID==DEBUG_VTXID||vtemp->ID()==DEBUG_VTXID)) cout << "\t Track removed from vID " << vtemp->ID() << "\tN " << vtemp->N() << "\tFlag " << vtemp->Flag() << endl;
                        if(vtemp->N()>2) {
                            //vtempS = vrec->RemoveTrackFromVertex(vtempS, trackfound);
                            int vertexid=vtemp->ID();
                            vtemp->SetFlag(-99);
                            EdbVTA *vta = vtemp->GetVTa(trackfound);
                            vtemp = vrec->RemoveVTAFromVertex(*vtemp, *vta);
                            vtemp->SetFlag(15);
                            vtemp->SetID(vertexid);
                            readd_vtx=1;
                        }
                        else vtemp->SetFlag(-99);
                        
                        if(EVERBOSE==100 && (cand->Track() == DEBUG_TRKID || cand->MCEvt()==DEBUG_MCEVT||vID==DEBUG_VTXID)) {
                            if(vtemp) {
                                cout << "\t\t" << vtemp->ID() << " becomes N " << vtemp->N() << "\tFlag " << vtemp->Flag() << endl;
                                
                            }
                        }
                    }
                }
                
                //aggiungo la traccia al nuovo vertice
                vnew = AddDauToVertex(vertex, cand, vrec);
                
                if(EVERBOSE==100 && ((cand->Track() == DEBUG_TRKID) || cand->MCEvt()==DEBUG_MCEVT||vID==DEBUG_VTXID)) {cout << "\t track " << cand->MCTrack() << "\t" << cand->Track() << " " << cand->Track() << "\t" << cand->MCEvt() << "\t" << cand->N() << "\t" << cand->X() << "\t" << cand->Y() << " added (if not duplicated) " << "\t";
                    if(cand->VertexS()) cout << "-> VS " << cand->VertexS()->ID();
                    cout << endl;
                }
                if(EVERBOSE==100 && vID==DEBUG_VTXID) cout << "add dau: N prong new " << vertex->N();
                if(EVERBOSE==100 && vID==DEBUG_VTXID) cout << endl;
            }
        } //if(imp < maximp_dau && cand->N()>=2)
        if(readd_vtx==1){
            int closeplate=Get_vtx_plate(vtemp->VZ());
            if(closeplate>= PLMIN && closeplate<=BRAGGPLATE){
                vtxPat[closeplate].Add(vtemp);
                if(EVERBOSE==100 && vID==DEBUG_VTXID) {cout << "Add to varr " << vtemp->ID() << "\t" << vtemp->N() << "\t" << vtemp->Flag() << "\t" << vtxPat[closeplate].GetEntries() << endl;}
            }
            else vtemp->SetFlag(-99);
        }
    } //if (tr->GetSegmentFirst()->Z() < vertex->VZ()){
    
    if(vnew==0) return 0;
    else return 1;
    
}

//---------------------------------------------------------//

EdbVertex * AddBeamToVertex(EdbVertex *vertex, EdbTrackP *found_cand, EdbVertexRec* vrec){
    int foundtrueBeam=0;
    float imp_beam=0, imp_dau=0, imp=0;
    
    float vz = vertex->VZ();
    int idvtx = vertex->ID();
    
    //tappullo per evitare vengano riaggiunte tracce già appartenenti al vertice
    for(int itrk=0; itrk<vertex->N(); itrk++){
        EdbTrackP *checkduplicates=vertex->GetTrack(itrk);
        if(found_cand->Theta()==checkduplicates->Theta()){
            if(EVERBOSE==3 || ( EVERBOSE==100 && idvtx==DEBUG_VTXID)) cout << "track\t" << found_cand->Track() << " duplicated" << endl;
            return vertex;
        }
    }
    
    if(EVERBOSE==2 || (EVERBOSE==100 && vertex->ID()==DEBUG_VTXID)) cout << "AddBeamToVertex " << idvtx << endl << "found cand Z " << found_cand->GetSegmentLast()->Z() << " (pl "<< found_cand->GetSegmentLast()->Plate() << ")\tVZ " << vz << " (pl " << Get_vtx_plate(vz) << ")" << endl;
    //cout << " qui " << found_cand->GetSegmentLast()->Z() << endl;
    if (found_cand->GetSegmentLast()->Z() > vz){ //track ending after vertex, splitting into two
        EdbTrackP * splitted_dau = new EdbTrackP();
        EdbTrackP * splitted_beam = new EdbTrackP();
        if(EVERBOSE==100 && vertex->ID()==DEBUG_VTXID) cout << "before split " << found_cand->N() << endl;
        SplitTrackAtZ(*found_cand, *splitted_beam, *splitted_dau, vz); //splitto found_cand in splitted_beam + splitted_dau nel punto vz
        imp_beam = CalcDist3D(vertex, splitted_beam); //CalcIP(vertex, splitted_beam, 0); //
        imp_dau = CalcDist3D(vertex, splitted_dau); //CalcIP(vertex, splitted_beam, 1); //
        if(EVERBOSE==100 && vertex->ID()==DEBUG_VTXID) cout << imp_beam << "\t" << imp_dau << endl;
        if(imp_beam > maximp_beam && imp_dau > maximp_dau) return vertex;
        else{
            if(found_cand->Flag()==0) found_cand->SetFlag(9); //evito venga aggiunta nuovamente
            splitted_beam->SetPDG(8); //al fascio assegno PDG 8
            if(EVERBOSE==2 || ( EVERBOSE==100 && (((splitted_beam->Track()==DEBUG_TRKID || splitted_dau->Track()==DEBUG_TRKID))|| splitted_beam->MCEvt()==DEBUG_MCEVT || splitted_dau->MCEvt()==DEBUG_MCEVT))) {
                cout << "\t\tbeam splitted confirmed: " << splitted_beam->Track() << "\t" << splitted_beam->MCEvt() << "\t" << splitted_beam->MCTrack() <<  "\t" << splitted_beam->N() << "\t" << splitted_beam->GetSegmentFirst()->Plate() << " " << splitted_beam->GetSegmentLast()->Plate() << "\t" << splitted_beam->N() << endl;
                cout << "\t\tdau confirmed: " << splitted_dau->Track() << "\t" << splitted_dau->MCEvt() << "\t" << splitted_dau->MCTrack() << "\t" << splitted_dau->N() << "\t" << splitted_dau->GetSegmentFirst()->Plate() << " " << splitted_dau->GetSegmentLast()->Plate() << "\t" << splitted_dau->Theta() << endl;
            }
            
            if(EVERBOSE==100 && vertex->ID()==DEBUG_VTXID) cout << "creo vta " << endl;
            
            //creo i vta
            EdbVTA *vta = new EdbVTA(splitted_beam,vertex);
            vta->SetFlag(2);
            //vta->SetZpos(0);
            (splitted_beam->Z() >= vertex->VZ())? vta->SetZpos(1) : vta->SetZpos(0);
            vertex->AddVTA(vta);
            splitted_beam->AddVTA(vta);
            found_cand->AddVTA(vta);
            
            if(EVERBOSE==100 && vertex->ID()==DEBUG_VTXID) cout << "creo vta 1" << endl;
            
            EdbVTA *vta1 = new EdbVTA(splitted_dau,vertex);
            vta1->SetFlag(2);
            //vta1->SetZpos(1);
            (splitted_dau->Z() >= vertex->VZ())? vta1->SetZpos(1) : vta1->SetZpos(0);
            vertex->AddVTA(vta1);
            splitted_dau->AddVTA(vta1);
            
            if(EVERBOSE==100 && vertex->ID()==DEBUG_VTXID) cout << "fine vta " << endl;
            
            N_BEAM_FOUND++;
            N_BEAM_SPLITTED++;
            //cout << endl << ivertex << "\tvID: " << vID << "\t" << vertex->N() << "\t" << mostfrequentevent[ivertex];
            if (MC==1 && splitted_beam->MCEvt() == mostfrequentevent[idvtx] && splitted_beam->GetSegmentFirst()->MCTrack()==1) foundtrueBeam = 1;  //verifico se l'id MC della traccia è lo stesso del vtx
            else if(MC==11){
                if(splitted_beam->GetSegmentLast()->MCEvt()==splitted_dau->GetSegmentFirst()->MCEvt() && splitted_beam->GetSegmentLast()->MCTrack()!=splitted_dau->GetSegmentFirst()->MCTrack()) N_BEAM_SPLITTED_OK++;
                
                for(int itrk=0; itrk<vertex->N(); itrk++){
                    EdbTrackP *tr=vertex->GetTrack(itrk);
                    if(((found_cand->GetSegmentLast()->MCEvt() == tr->MCEvt()||found_cand->GetSegmentFirst()->MCEvt() == tr->MCEvt()) && (found_cand->Track() != tr->Track()))||(splitted_dau->GetSegmentFirst()->MCEvt() == tr->MCEvt() && splitted_dau->Track() != tr->Track())) {
                        foundtrueBeam=1;//verifico se l'id MC della traccia è lo stesso del vtx
                        break;
                    }
                }
            }
        }
    }
    else{ //l'ossigeno si ferma prima del vertice
        EdbVTA *vta = new EdbVTA(found_cand,vertex);
        vta->SetFlag(2);
        //vta->SetZpos(0);
        (found_cand->Z() >= vertex->VZ())? vta->SetZpos(1) : vta->SetZpos(0);
        
        vertex->AddVTA(vta);
        imp = CalcDist3D(vertex, found_cand); //CalcIP(vertex, found_cand); //
        found_cand->AddVTA(vta);
        N_BEAM_FOUND++;
        if(EVERBOSE==2 || (EVERBOSE==100 && (found_cand->Track() == DEBUG_TRKID) || found_cand->MCEvt()==DEBUG_MCEVT)) cout << "\t\tbeam confirmed: " << found_cand->Track() << "\t" << found_cand->MCEvt() << "\t" << found_cand->MCTrack() <<  "\t" << found_cand->N() << "\t" << found_cand->GetSegmentFirst()->Plate() << " " << found_cand->GetSegmentLast()->Plate() << endl;
        
        //cout << endl << ivertex << "\tvID: " << vID << "\t" << vertex->N() << "\t" << mostfrequentevent[idvtx];
        if (MC==1 && found_cand->MCEvt() == mostfrequentevent[idvtx] && found_cand->GetSegmentFirst()->MCTrack()==1) foundtrueBeam = 1;//verifico se l'id MC della traccia è lo stesso del vtx
        else if(MC==11){
            for(int itrk=0; itrk<vertex->N(); itrk++){
                EdbTrackP *tr=vertex->GetTrack(itrk);
                if((found_cand->GetSegmentLast()->MCEvt() == tr->MCEvt()) && (found_cand->Track() != tr->Track())) {
                    foundtrueBeam=1;//verifico se l'id MC della traccia è lo stesso di almeno una traccia del vtx
                    break;
                }
            }
        }
    }
    if ((MC==1 || MC==11) && foundtrueBeam==1) N_BEAM_TRUE++;
    
    vertex->EstimateVertexFlag();
    
    float final_imp=0;
    
    for(int i=0; i<vertex->N(); i++){
        final_imp=CalcIP(vertex, i);
        
        if(final_imp>Max_acceptable_IP && (PULIZIA_EXTRA==1)){
            int tempid = vertex->ID();
            if(vertex->N()>2) vertex=vrec->RemoveTrackFromVertex(vertex,i);
            else vertex->SetFlag(-99);
            vertex->SetID(tempid);
        }
        else{
            vertex->GetVTa(i)->SetImp(final_imp);
        }
    }
    
    return vertex;
}

//---------------------------------------------------------//

int AddDauToVertex(EdbVertex *vertex, EdbTrackP *found_cand, EdbVertexRec *vrec){
    
    float vz = vertex->VZ();
    int vID = vertex->ID();
    int idvtx = vertex->ID(); //check if this is correct, this variable was nowhere in the function!
    if(EVERBOSE==3 || ( EVERBOSE==100 && vID==DEBUG_VTXID)) cout << "Trying to add a dau to vertex id " << vID << endl;
    
    //tappullo per evitare vengano riaggiunte tracce già appartenenti al vertice
    for(int itrk=0; itrk<vertex->N(); itrk++){
        EdbTrackP *checkduplicates=vertex->GetTrack(itrk);
        if(found_cand->Theta()==checkduplicates->Theta()){
            if(EVERBOSE==3 || ( EVERBOSE==100 && vID==DEBUG_VTXID)) cout << "track\t" << found_cand->Track() << " duplicated" << endl;
            return 0;
        }
    }
    
    if(EVERBOSE==3 || ( EVERBOSE==100 && vID==DEBUG_VTXID)) cout << "\t\tdau confirmed: " << found_cand->Track() << "\t" << found_cand->MCEvt() << "\t" << found_cand->MCTrack() << "\t" << found_cand->N() << "\t" << found_cand->GetSegmentFirst()->Plate() << " " << found_cand->GetSegmentLast()->Plate() << endl;
    
    //aggiungo traccia al vertice
    EdbVTA *vta = new EdbVTA(found_cand,vertex);
    vta->SetFlag(2);
    vertex->AddVTA(vta);
    if(found_cand){
        (found_cand->Z() >= vertex->VZ())? vta->SetZpos(1) : vta->SetZpos(0);
        if( EVERBOSE==100 && vID==DEBUG_VTXID) cout << "Z check " << found_cand->Z() << "\t" << vertex->VZ() << "\t" << vta->Zpos() << endl;
    }
    //    if (found_cand->GetSegmentFirst()->Z() > vz) vta->SetZpos(1); //track starts after vertex
    //    else vta->SetZpos(0);
    
    float imp = CalcIP(vertex, found_cand);//CalcDist3D(vertex, found_cand);
    if(imp>Max_acceptable_IP){
        vertex->SetFlag(-99);
        vertex=vrec->RemoveVTAFromVertex(*vertex,*vta);
        vertex->SetID(vID);
    }
    else{
        vta->SetImp(imp);
        found_cand->AddVTA(vta);
        if(found_cand->Flag()==0) found_cand->SetFlag(9); //evito venga aggiunta nuovamente
        
        N_DAU_FOUND++;
        //cout << endl << ivertex << "\tvID: " << vID << "\t" << vertex->N() << "\t" << mostfrequentevent[ivertex];
        
        if(MC==1 && found_cand->MCEvt() == mostfrequentevent[idvtx]) N_DAU_TRUE++;//verifico se l'id MC della traccia è lo stesso del vtx
        if(MC==11){
            for(int itrk=0; itrk<vertex->N(); itrk++){
                EdbTrackP *tr=vertex->GetTrack(itrk);
                if(found_cand->MCEvt() == tr->MCEvt() && found_cand->Track() != tr->Track()) {
                    N_DAU_TRUE++;//verifico se l'id MC della traccia è lo stesso del vtx
                    break;
                }
            }
        }
    }
    return 1;
    
}

//---------------------------------------------------------//

EdbVertex * SelectOneBeam(EdbVertex *vertex, EdbVertexRec *vrec){
    float maximp=10000;
    int beam_itrk[10]={-99,-99,-99,-99,-99,-99,-99,-99,-99,-99}; //immagino e spero non possano esserci + di 10 beam collegati!
    int save_itrk=-99;
    int nbeam=0;
    if(EVERBOSE==6) cout << "vtx n " <<  vertex->N() << " vtx ID "<<vertex->ID()<<endl;
    for (int itrk = 0; itrk < vertex->N(); itrk++){
        EdbTrackP *track = vertex->GetTrack(itrk);
        if(EVERBOSE==6 || (EVERBOSE==100 && (( track->Track()==DEBUG_TRKID) || track->MCEvt()==DEBUG_MCEVT))) cout << endl << track->MCEvt() << "\t" << vertex->GetVTa(itrk)->Zpos() << "\t" << track->GetSegmentFirst()->W()-70 << "\t" << track->GetSegmentLast()->W()-70 << "\t";
        if(track->Theta()<BEAMTHETA && vertex->GetVTa(itrk)->Zpos()==0){
            beam_itrk[nbeam]=itrk;
            nbeam++;
            float imp =  CalcIP(track->GetSegmentLast(), vertex->VX(),vertex->VY(),vertex->VZ()); //calcolo il parametro d'impatto tra l'ultimo seg del beam e il vertice
            if(EVERBOSE==6 || (EVERBOSE==100 && (( track->Track()==DEBUG_TRKID) || track->MCEvt()==DEBUG_MCEVT))) cout << "Beam " << "\t" << track->Theta() << "\t" << track->GetSegmentFirst()->Plate() << "\t" << track->N();
            
            if(imp < maximp){ //found candidate, storing track and new min IP
                maximp=imp;
                save_itrk=itrk;
            }
        }
    }
    
    if(EVERBOSE==6){
        EdbTrackP *track_ok = vertex->GetTrack(save_itrk);
        cout << "Chosen beam " << save_itrk << "\t" << track_ok->MCEvt() << "\t" << track_ok->W()-70 << "\t" << track_ok->Theta() << "\t" << track_ok->GetSegmentFirst()->Plate() << "\t" << track_ok->N() <<  endl;
    }
    //we can't remove and access VTA at the same time! Proposed solution: filling an array of vta to remove
    TObjArray vta_toremove;
    for(int i=0; i<nbeam; i++){
        if(beam_itrk[i]!=-99 && beam_itrk[i]!=save_itrk){
            EdbVTA *vta = vertex->GetVTa(beam_itrk[i]);
            EdbTrackP *track = vertex->GetTrack(beam_itrk[i]);
            if(EVERBOSE==6|| (EVERBOSE==100 && (( track->Track()==DEBUG_TRKID) || track->MCEvt()==DEBUG_MCEVT))){
                cout << "Adding track to array of to be removed " << track->MCEvt() << "\t" << track->W()-70 << "\t" << track->Theta() << "\t" << track->GetSegmentFirst()->Plate() << "\t" << track->N() <<  endl;
            }
            //only adding to container now, we will remove them later
            vta_toremove.Add(vta);
        }
    }
    
    
    //actually removing now
    int vID = vertex->ID();
    for (int ivta = 0; ivta<vta_toremove.GetEntries(); ivta++){
        EdbVTA *vta = (EdbVTA*) vta_toremove.At(ivta);
        
        for(int itrk=0; itrk<vertex->N(); itrk++){
            EdbVTA *vta1 = (EdbVTA*) vertex->GetVTa(itrk);
            if(vta->GetTrack()->Track()==vta1->GetTrack()->Track()){
                vertex->SetFlag(-99);//add 8dec
                vertex = vrec->RemoveVTAFromVertex(*vertex, *vta1);
                break;
            }}
    } //end loop of vta to remove
    vertex->SetID(vID);
    
    
    if(EVERBOSE==6||(EVERBOSE==100 && vertex->ID()==DEBUG_VTXID)){
        cout << "SelectOneBeam end " << vertex->ID() << endl;
        for(int i=0; i<vertex->N(); i++){
            cout << i << "\t" << vertex->GetVTa(i)->Imp() << " " << vertex->GetTrack(i)->Track() <<  endl;
        }
    }
    return vertex;
}


//---------------------------------------------------------//

void FindNitrogen(TObjArray *varr, EdbVertexRec *vrec, float maximp_Np, float maximp_dau){
    
    float kink=0, temp_kink=0, mean_kink=0;
    TObjArray tr_grid_oxy;
    
    TStopwatch t;
    t.Start();
    
    for(int ipl = PLMIN; ipl<=BRAGGPLATE; ipl++){
        int nbeam_grid = gridtr_OXY[ipl].SelectObjects(tr_grid_oxy);
        for (int itrk = 0; itrk < nbeam_grid; itrk++){
            EdbTrackP *beam = (EdbTrackP*) tr_grid_oxy.At(itrk);
            if(FAST==4 && ((beam->Track()!=DEBUG_TRKID)||beam->MCEvt()!=DEBUG_MCEVT)) continue;
            int nseg = beam->N();
            if(beam->Npl()>15 && beam->GetSegmentLast()->Plate()>BRAGGPLATE){ //la traccia deve avere più di 15 piatti, angolo del fascio e terminare oltre il picco di bragg. Le flag evitano siano tracce già associate a qualcosa nel loop precedente
                EdbVertex *vtempS = beam->VertexS();
                EdbVertex *vtempE = beam->VertexE();
                if(!(vtempS) && !(vtempE)){ //controllo che la traccia non abbia vertici entranti o uscenti
                    if(EVERBOSE==4 || (EVERBOSE==100 && ((beam->Track()==DEBUG_TRKID) || beam->MCEvt()==DEBUG_MCEVT))) cout << "\t Entering N search" << endl;
                    kink=0;
                    mean_kink=0;
                    int isplit=0;
                    int seg_endsearch=0;
                    if(BRAGGPLATE<=LASTLAYER[1]) seg_endsearch=nseg-1;
                    else seg_endsearch=LASTLAYER[1];
                    for (int iseg = 2; iseg < seg_endsearch; iseg++){ //loop on segments
                        if(beam->GetSegment(iseg+1)->Plate()<BRAGGPLATE){
                            temp_kink=GetThetaKink(beam->GetSegment(iseg), beam->GetSegment(iseg+1));
                            mean_kink+=temp_kink;
                            if(temp_kink>=kink) {
                                kink=temp_kink;
                                isplit=iseg;
                            }
                        }
                    }
                    if(isplit!=0){
                        if(EVERBOSE==4 || (EVERBOSE==100 && ((beam->Track()==DEBUG_TRKID) || beam->MCEvt()==DEBUG_MCEVT))) cout << "NP Debug kink: " << kink << "\t" << (double)mean_kink/(nseg-3) << "\t" <<  3*(double)mean_kink/(nseg-3) << "\t" << beam->GetSegment(isplit)->W()-70 << "\t" << beam->GetSegment(isplit+1)->W()-70 << "\t" << beam->GetSegment(isplit)->MCEvt() << "\t" << beam->GetSegment(isplit+1)->MCEvt() << "\t" << beam->GetSegmentLast()->Plate() << "\t" << beam->GetSegment(isplit)->Plate() << endl;
                        if(kink>2*(double)mean_kink/(nseg-3)){
                            EdbTrackP *cand_proton = FindProton(beam, isplit, maximp_Np, minseg_Np);
                            if(cand_proton){
                                EdbVertex *newvertex= new EdbVertex();
                                newvertex = MakeNpVertex(beam, cand_proton, vrec);
                                
                                //se ho trovato un vertice non è detto che sia O->N+p, potrebbe anche essere O->X+nP (n>1). Cerco allora altri figli nei paraggi del vertice appena creato
                                if(newvertex){
                                    int vplate=Get_vtx_plate(newvertex->VZ());
                                    if(vplate==-99) {
                                        newvertex->SetFlag(-99);
                                        continue;
                                    }
                                    float xy[2] = {newvertex->X(),newvertex->Y()};
                                    float r = 2000.;
                                    for(int iipl=vplate-1; iipl<vplate+3; iipl++){
                                        TObjArray tr_grid_dau;
                                        int ndau_grid = gridtr_DAU[iipl].SelectObjectsC( xy, r, tr_grid_dau);
                                        
                                        if(EVERBOSE==4||(EVERBOSE==100 && ((beam->Track()==DEBUG_TRKID) || beam->MCEvt()==DEBUG_MCEVT))) cout << "\tSearching other dau in plate: " << iipl << " between " <<  ndau_grid << " tracks" << endl;
                                        
                                        DaughtersSearch(newvertex, tr_grid_dau, vrec, maximp_dau);
                                    }//for(int iipl=vplate-1; iipl<vplate+3; iipl++)
                                    if(vplate>= PLMIN && vplate<=BRAGGPLATE){ vtxPat[vplate].Add(newvertex);
                                        varr->Add(newvertex);
                                    }
                                    else newvertex->SetFlag(-99);
                                    if(MC==11){
                                        if(newvertex->N()>3){
                                            N_NPP_FOUND++;
                                            int temp_ok=0;
                                            if(EVERBOSE==4||(EVERBOSE==100 && ((beam->Track()==DEBUG_TRKID) || beam->MCEvt()==DEBUG_MCEVT))) cout << "\tNEW NPP: " << newvertex->GetTrack(0)->MCEvt() << "\tFlag " << newvertex->Flag() << endl;
                                            for(int it=0; it<newvertex->N(); it++){
                                                EdbTrackP *newt=newvertex->GetTrack(it);
                                                if(EVERBOSE==4||(EVERBOSE==100 && ((beam->Track()==DEBUG_TRKID) || beam->MCEvt()==DEBUG_MCEVT))) cout << "\t\t" << it << "\t" << newt->MCEvt() << "\t" << newt->Track() << "\t" << newt->W()-70 << "\t" << newvertex->GetVTa(it)->Imp() <<endl;
                                                if(it>2 && (newt->MCEvt()== newvertex->GetTrack(1)->MCEvt() || newt->MCEvt()== newvertex->GetTrack(2)->MCEvt())) temp_ok++;
                                            }
                                            if(temp_ok>=1) N_NPP_FOUND_OK++;
                                        }
                                    }//if(MC==11){
                                }//if(newvertex){
                            }//if(cand_proton)
                        } //if(kink>2*(double)mean_kink/(nseg-3))
                    }//if(isplit!=0)
                }
            }//if(beam->Npl()>20 && beam->GetSegmentLast()->Plate()>BRAGGPLATE && (beam->Flag()!=8 && beam->Flag()!=9)
        }//for (int itrk = 0; itrk < nbeam_grid; itrk++)
        tr_grid_oxy.Clear();
        
        if ((BRAGGPLATE-ipl)%1 == 0 && (EVERBOSE==14||EVERBOSE==16||EVERBOSE==101)) {
            cout << " Find Nitrogen: " << 100.*(ipl)/(BRAGGPLATE-PLMIN+1) << " %, Iteration Time: " << t.RealTime() << " s\t" << t.RealTime()/60 << " min " << endl;
            t.Reset();
            t.Start();
        }
    }
    
}


//---------------------------------------------------------//

EdbTrackP *FindProton(EdbTrackP *beam, int iseg, float maximp, int minseg){
    float minimp=maximp;
    float imp=0;
    int splitatseg=0;
    int r[2] = {2,2};
    
    EdbTrackP *proton_candidate=new EdbTrackP();
    EdbSegP* seg;
    for(int is=iseg; is<=iseg+2; is++){
        if(is<beam->N()) seg=beam->GetSegment(is);
        else break;
        
        int segPlate=seg->Plate();
        if(EVERBOSE==100 && ((beam->Track()==DEBUG_TRKID) || beam->MCEvt()==DEBUG_MCEVT)) cout << segPlate << endl;
        if(segPlate>=3&&segPlate<BRAGGPLATE){ //controllo che il segmento non sia nei primissimi piatti e sia prima del picco di bragg
            
            float xy[2] = {seg->X(),seg->Y()}; //coordinate della cella
            float imp=0;
            TObjArray tr_cand_protons;
            for(int isegplate=segPlate-1; isegplate<segPlate+2; isegplate++){
                TObjArray tr_grid_protons;
                int ntrks_grid_protons = gridtr_DAU[isegplate].SelectObjectsC( xy, r, tr_grid_protons); //la griglia selezionata corrisponde al piatto del segmento in esame
                for (int itrk2 = 0; itrk2 < ntrks_grid_protons; itrk2++){
                    EdbTrackP *cand_proton = (EdbTrackP*) tr_grid_protons.At(itrk2);
                    
                    if ((cand_proton->N()>=minseg_Np||cand_proton->GetSegmentFirst()->Plate()>=BRAGGPLATE-minseg_Np) && cand_proton->Theta()>BEAMTHETA){ //MIN N PROTON
                        EdbVertex *vtempS2 = cand_proton->VertexS();
                        EdbVertex *vtempE2 = cand_proton->VertexE();
                        if(EVERBOSE==4||(EVERBOSE==100 && ((beam->Track()==DEBUG_TRKID) || beam->MCEvt()==DEBUG_MCEVT))){
                            cout << endl << "\t Np deb proton: " << cand_proton->GetSegmentFirst()->MCEvt() << "\t" << cand_proton->Track() << "\t" << cand_proton->GetSegmentFirst()->Plate() << "\t" << cand_proton->Theta() << "\t" << cand_proton->N() << "\t" << cand_proton->Flag() << "\t";
                            if(vtempS2) cout << "vtempS2: " << vtempS2->ID() << " " << vtempS2->N() << "\t";
                            if(vtempE2) cout << "vtempE2: " << vtempE2->ID() << " " << vtempE2->N() << "\t";
                        }
                        if(!(vtempS2) && !(vtempE2)){ //non deve essere collegata a nessun vertice
                            float pv[3]={0};
                            bool parallel;
                            imp = mygEVR->CheckImpactN(seg, cand_proton->GetSegmentFirst(), pv, parallel, 2500);
                            if(EVERBOSE==4||(EVERBOSE==100 && ((beam->Track()==DEBUG_TRKID) || beam->MCEvt()==DEBUG_MCEVT))) cout << "\t imp " << imp << "\t" << minimp << endl;
                            //if(EVERBOSE==100 && seg->MCEvt()==DEBUG_MCEVT) cout << "\t imp " << imp << endl;
                            if(imp < minimp){ //found candidate, storing track and new min IP
                                proton_candidate = cand_proton;
                                minimp = imp;
                                splitatseg=is;
                                if(minimp<maximp) tr_cand_protons.Add(proton_candidate);
                            } //imp condition
                        } //if((!vtempS2) && (!vtempE2))
                    } //if (cand_proton...
                } // for itrk2
            } //for(int isegplate=segPlate-1; isegplate<segPlate+1; isegplate++){
        } //if(segPlate>3&&segPlate<BRAGGPLATE)
    }//for(int is=iseg; is<iseg+1; iseg++){
    
    // if(tr_cand_protons.GetEntries()>1) MakeNpvecVertex(beam, tr_cand_protons, vrec);
    
    if(minimp<maximp) return proton_candidate;
    else return NULL;
}


//---------------------------------------------------------//


EdbVertex * MakeNpVertex(EdbTrackP *beam, EdbTrackP *proton, EdbVertexRec *vrec){
    
    EdbVertex *NEWvertex = new EdbVertex();
    
    EdbTrackP *splitted_beam = new EdbTrackP();
    EdbTrackP *splitted_dau = new EdbTrackP();
    
    SplitTrackAtZ(*beam, *splitted_beam, *splitted_dau, proton->GetSegmentFirst()->Z()); //divido la traccia del fascio in 2 pezzi in corrispondenza di dove ho trovato il protone
    
    splitted_beam->SetPDG(8); //setto a 8 il pdg per il fascio
    splitted_dau->SetPDG(7); //setto a 7 il pdg per l'azoto N
    proton->SetPDG(1); //setto a 1 il pdg per il protone
    splitted_beam->SetFlag(8); //setto a 8 il pdg per il fascio
    splitted_dau->SetFlag(7); //setto a 7 il pdg per l'azoto N
    proton->SetFlag(1); //setto a 1 il pdg per il protone
    //    if(splitted_beam->Flag()==0) splitted_beam->SetFlag(9); //setto flag 9 per evitare siano associate di nuovo
    //    if(splitted_dau->Flag()==0) splitted_dau->SetFlag(9);//setto flag 9 per evitare siano associate di nuovo
    //    if(proton->Flag()==0) proton->SetFlag(9);//setto flag 9 per evitare siano associate di nuovo
    //
    TObjArray tr_cand_protons;
    
    //creo il nuovo vertice
    TObjArray tarr;
    tarr.Add(splitted_beam);
    tarr.Add(splitted_dau);
    tarr.Add(proton);
    NEWvertex->SetXYZ( 0,0, splitted_beam->GetSegmentLast()->Z() );
    int ntr = tarr.GetEntries();
    int BADTRK=0;
    for(int i=0; i<ntr; i++) {
        EdbTrackP *t = (EdbTrackP*)tarr.At(i);
        EdbVTA *vta = new EdbVTA(t,NEWvertex);
        vta->SetFlag(2);
        NEWvertex->AddVTA(vta);
        (t->Z() >= splitted_beam->GetSegmentLast()->Z())? vta->SetZpos(1) : vta->SetZpos(0);
        t->AddVTA(vta);
    }
    /*
     if( vrec->MakeV(*NEWvertex) )  {
     vrec->AddVertex(NEWvertex);
     }
     */
    vrec->MakeV(*NEWvertex);
    
    
    double imp[3]={0};
    for(int i=0; i<ntr; i++) {
        EdbVTA *vta = NEWvertex->GetVTa(i);
        EdbTrackP *t = NEWvertex->GetTrack(i);
        imp[i] = CalcIP(NEWvertex, i);
        if(EVERBOSE==4|| (EVERBOSE==100 && (((splitted_beam->Track()==DEBUG_TRKID || splitted_dau->Track()==DEBUG_TRKID))||(splitted_beam->GetSegmentLast()->MCEvt()==DEBUG_MCEVT || splitted_beam->GetSegmentFirst()->MCEvt()==DEBUG_MCEVT || splitted_dau->GetSegmentFirst()->MCEvt()==DEBUG_MCEVT)))) {
            
            cout <<"NP debug IP " << t->MCEvt() << "\t" << t->Track() << "\t"<< t->GetSegmentFirst()->Plate() << " - " << t->GetSegmentLast()->Plate() << ": " << imp[i] << "\t" << CalcDist3D(NEWvertex, t) << endl;
        }
        if(imp[i]>100) return NULL;
        //if(imp>maximp_Np) BADTRK++;
        vta->SetImp(imp[i]);
    }
    
    //if(BADTRK>1) return NULL;
    if(imp[1]>maximp_Np || imp[2]>maximp_Np) return NULL;
    else{
        N_NP_FOUND++;
        NEWvertex->SetID(900000+N_NP_FOUND);
        if(EVERBOSE==4 || (EVERBOSE==100 && ((beam->Track()==DEBUG_TRKID)||beam->MCEvt()==DEBUG_MCEVT))) {
            cout << "Newvertex\tIP";
            for(int i=0; i<ntr; i++) {
                cout << "\t" << NEWvertex->GetVTa(i)->Imp();
            }
            cout << endl;
        }
        if(EVERBOSE==4 || (EVERBOSE==100 && ((beam->Track()==DEBUG_TRKID)||beam->MCEvt()==DEBUG_MCEVT))) cout << "New NP-> temporary vID: " << NEWvertex->ID() << "\t" << NEWvertex->GetTrack(0)->MCEvt() <<  "\t" << NEWvertex->GetTrack(0)->Track() << "\t" << NEWvertex->GetTrack(0)->Theta() << endl;
        if (MC==11 && proton->GetSegmentFirst()->MCEvt() == splitted_dau->GetSegmentFirst()->MCEvt()) N_NP_FOUND_OK++; //se protone e N hanno stesso id l'evento è ok
        
        // else if(EVERBOSE==4 || (EVERBOSE==100 && (splitted_beam->GetSegmentLast()->MCEvt()==DEBUG_MCEVT || splitted_beam->GetSegmentFirst()->MCEvt()==DEBUG_MCEVT || splitted_dau->GetSegmentFirst()->MCEvt()==DEBUG_MCEVT)){
        //cout << "WRONG NP VTX\t" << splitted_beam->GetSegmentLast()->MCEvt() << " (" << splitted_beam->GetSegmentLast()->W()-70 << ")\t" << splitted_dau->GetSegmentFirst()->MCEvt() << " (" << splitted_dau->GetSegmentLast()->W()-70 << ")\t" << proton->GetSegmentFirst()->MCEvt() << " (" << proton->GetSegmentLast()->W()-70 << ")" << endl;}
        
        if(EVERBOSE==4|| (EVERBOSE==100 && (((splitted_beam->Track()==DEBUG_TRKID || splitted_dau->Track()==DEBUG_TRKID))||(splitted_beam->GetSegmentLast()->MCEvt()==DEBUG_MCEVT || splitted_beam->GetSegmentFirst()->MCEvt()==DEBUG_MCEVT || splitted_dau->GetSegmentFirst()->MCEvt()==DEBUG_MCEVT)))) {
            cout << "VERTEX CREATED!\t Nprong: " << NEWvertex->N() << "\t" << NEWvertex->ID() << endl;
            cout << "\t\tbeam: " << splitted_beam->MCEvt()  <<  "\t" << splitted_beam->Track()  <<  "\t" << splitted_beam->N() << "\t" << splitted_beam->GetSegmentFirst()->Plate() << " " << splitted_beam->GetSegmentLast()->Plate() << "\t" << splitted_beam->GetSegmentFirst()->W()-70 << " " << splitted_beam->GetSegmentLast()->W()-70 << endl;
            
            cout << "\t\tdau: " << splitted_dau->MCEvt() << "\t" << splitted_dau->Track()  <<  "\t" << splitted_dau->N() << "\t" << splitted_dau->GetSegmentFirst()->Plate() << " " << splitted_dau->GetSegmentLast()->Plate() << "\t" << splitted_dau->GetSegmentFirst()->W()-70 << " " << splitted_dau->GetSegmentLast()->W()-70 << endl;
            
            cout <<  "\t\tproton: " << proton->MCEvt() << "\t" << proton->Track()  <<  "\t" <<  proton->N() << "\t" << proton->GetSegmentFirst()->Plate() << " " << proton->GetSegmentLast()->Plate() << "\t" << proton->GetSegmentFirst()->W()-70 << " " << proton->GetSegmentLast()->W()-70 << "\t" << proton->Theta() << endl;
        }
        
        return NEWvertex;
    }
}

//---------------------------------------------------------//


bool SplitTrackAtZ( EdbTrackP &t, EdbTrackP &t1, EdbTrackP &t2, float vz ){
    // Code copied from EdbTrackFitter::SplitTrack to split a track in 2.
    // Only here I set the split point at z of vertex, instead of in a certain segment number
    //cout << "in SplitTrackAtZ " << vz << endl;
    
    int j=0, jj=0;
    for(int i=t.N()-1; i>=0; i--) {
        EdbSegP *seg = t.GetSegment(i);
        if ((DIRECTION==9 && seg->Z() > vz)||(DIRECTION==1 && seg->Z() < vz)){
            t1.AddSegment(   t.GetSegment(i) );
            t1.AddSegmentF(   t.GetSegmentF(i) );
            t1.GetSegment(j)->SetTrack(t.ID());
            t1.GetSegmentF(j)->SetTrack(t.ID());
            t1.GetSegment(j)->SetDZ(300);
            t1.GetSegmentF(j)->SetDZ(300);
            
            
            j++;
        }
        else if ((DIRECTION==9 && seg->Z() <= vz)||(DIRECTION==1 && seg->Z() >= vz)){
            t2.AddSegment(   t.GetSegment(i) );
            t2.AddSegmentF(   t.GetSegmentF(i) );
            t2.GetSegment(jj)->SetTrack(t.ID()+9000000);
            t2.GetSegmentF(jj)->SetTrack(t.ID()+9000000);
            t2.GetSegment(j)->SetDZ(300);
            t2.GetSegmentF(j)->SetDZ(300);
            
            jj++;
        }
    }
    t.SetCounters();
    //t1.FitTrackKFS(true, X0, 0);
    t1.SetCounters();
    t2.SetCounters();
    
    t1.SetX(t.X());
    t1.SetY(t.Y());
    t1.SetZ(t.Z());
    t1.SetTX(t.TX());
    t1.SetTY(t.TY());
    
    t1.SetM(t.M());
    t1.SetP(t.P());
    t1.SetSegmentsTrack();
    t1.SetID(t.ID());
    t1.SetTrack(t.ID());
    t1.SetMC(t.MCEvt(),t.MCTrack());
    t1.SetAid(t.Aid(0), t.Aid(1)); //setto pdgcode e motherID
    //t1.FitTrack();
    
    
    
    t2.SetSegmentsTrack();
    t2.SetM(t.M());
    t2.SetP(t.P());
    EdbSegP *seg0 = t2.GetSegment(0);
    //t1.FitTrackKFS(true, X0, 0);
    t2.SetAid(seg0->Aid(0), t1.ID()); //setto pdgcode e motherID
    t2.SetX(seg0->X());
    t2.SetY(seg0->Y());
    t2.SetZ(seg0->Z());
    t2.SetTX(seg0->TX());
    t2.SetTY(seg0->TY());
    t2.SetID(t.ID()+90000000);
    t2.SetTrack(t.ID()+90000000);
    t2.SetMC(seg0->MCEvt(),seg0->MCTrack());
    //t2.FitTrack();
    
    // if(t.ID()==1671) cout << "counters " << t1.N() << "\t" << t2.N() << endl;
    
    
    return true;
}

//---------------------------------------------------------------------

void CreateTree(TTree *new_vtxtree, TObjArray *varr){
    
    Float_t vx, vy, vz;
    Float_t maxaperture;
    Float_t probability;
    Int_t n;
    Int_t v_flag;
    Int_t vplate;
    const Int_t maxdim = 5000; //big arrays for containers of track variables
    //
    Int_t vID=0;
    Int_t IDTrack[maxdim]={0};
    Int_t TrackTrack[maxdim]={0};
    Int_t maxgap[maxdim]={0};
    Int_t nseg_vec[maxdim]={0};
    Int_t nseg_S1[maxdim]={0};
    Int_t npl[maxdim]={0};
    Int_t npl_S1[maxdim]={0};
    Int_t nholes[maxdim]={0};
    Int_t nholes_S1[maxdim]={0};
    Int_t plate[maxdim]={0};
    Float_t X[maxdim]={0};
    Float_t Y[maxdim]={0};
    Float_t Z[maxdim]={0};
    Float_t TX[maxdim]={0};
    Float_t TY[maxdim]={0};
    Float_t Theta[maxdim]={0};
    Float_t impactparameter[maxdim]={0};
    Int_t incoming[maxdim]={0}; //coming from the vertex
    Int_t Z_flag[maxdim]={0};
    Int_t MC_Charge_first[maxdim]={0};
    Int_t MC_Charge_last[maxdim]={0};
    Int_t MC_Charge_S2[maxdim]={0};
    Int_t MC_evID_last[maxdim]={0};
    Int_t MC_evID_first[maxdim]={0};
    Int_t MC_trackID_first[maxdim]={0};
    Int_t MC_trackID_last[maxdim]={0};
    Int_t MC_IDPart[maxdim]={0};
    Int_t MC_nseg[maxdim]={0};
    Int_t MC_mother_first[maxdim]={0};
    Int_t MC_mother_last[maxdim]={0};
    Int_t MC_pdgcode_first[maxdim]={0};
    Int_t MC_pdgcode_last[maxdim]={0};
    Int_t MC_firstplate[maxdim]={0};
    Int_t MC_lastplate[maxdim]={0};
    //
    new_vtxtree->Branch("vID",&vID,"vID/I");
    new_vtxtree->Branch("vx",&vx,"vx/F");
    new_vtxtree->Branch("vy",&vy,"vy/F");
    new_vtxtree->Branch("vz",&vz,"vz/F");
    new_vtxtree->Branch("vplate",&vplate,"vplate/I");
    new_vtxtree->Branch("v_flag",&v_flag,"v_flag/I");
    new_vtxtree->Branch("maxaperture",&maxaperture,"maxaperture/F");
    new_vtxtree->Branch("probability",&probability,"probability/F");
    new_vtxtree->Branch("n",&n,"n/I");
    new_vtxtree->Branch("IDTrack",&IDTrack,"IDTrack[n]/I");
    new_vtxtree->Branch("TrackTrack",&TrackTrack,"TrackTrack[n]/I");
    new_vtxtree->Branch("X",&X,"X[n]/F");
    new_vtxtree->Branch("Y",&Y,"Y[n]/F");
    new_vtxtree->Branch("Z",&Z,"Z[n]/F");
    new_vtxtree->Branch("TX",&TX,"TX[n]/F");
    new_vtxtree->Branch("TY",&TY,"TY[n]/F");
    new_vtxtree->Branch("Theta",&Theta,"Theta[n]/F");
    new_vtxtree->Branch("nseg",&nseg_vec,"nseg_vec[n]/I");
    new_vtxtree->Branch("npl",&npl,"npl[n]/I");
    new_vtxtree->Branch("nholes",&nholes,"nholes[n]/I");
    new_vtxtree->Branch("nseg_S1",&nseg_S1,"nseg_S1[n]/I");
    new_vtxtree->Branch("nholes_S1",&nholes_S1,"nholes_S1[n]/I");
    new_vtxtree->Branch("plate",&plate," plate[n]/I");
    new_vtxtree->Branch("maxgap",&maxgap,"maxgap[n]/I");
    new_vtxtree->Branch("incoming",&incoming,"incoming[n]/I");
    new_vtxtree->Branch("impactparameter",&impactparameter,"impactparameter[n]/F");
    new_vtxtree->Branch("Z_flag",&Z_flag,"Z_flag[n]/I");
    new_vtxtree->Branch("MC_Charge_first",&MC_Charge_first,"MC_Charge_first[n]/I");
    new_vtxtree->Branch("MC_Charge_last",&MC_Charge_last,"MC_Charge_last[n]/I");
    new_vtxtree->Branch("MC_Charge_S2",&MC_Charge_S2,"MC_Charge_S2[n]/I");
    new_vtxtree->Branch("MC_evID_first",&MC_evID_first,"MC_evID_first[n]/I");
    new_vtxtree->Branch("MC_evID_last",&MC_evID_last,"MC_evID_last[n]/I");
    new_vtxtree->Branch("MC_trackID_first",&MC_trackID_first,"MC_trackID_first[n]/I");
    new_vtxtree->Branch("MC_trackID_last",&MC_trackID_last,"MC_trackID_last[n]/I");
    new_vtxtree->Branch("MC_mother_first",&MC_mother_first,"MC_mother_first[n]/I");
    new_vtxtree->Branch("MC_mother_last",&MC_mother_last,"MC_mother_last[n]/I");
    new_vtxtree->Branch("MC_pdgcode_first",&MC_pdgcode_first,"MC_pdgcode_first[n]/I");
    new_vtxtree->Branch("MC_pdgcode_last",&MC_pdgcode_last,"MC_pdgcode_last[n]/I");
    new_vtxtree->Branch("MC_firstplate",&MC_firstplate,"MC_firstplate[n]/I");
    new_vtxtree->Branch("MC_lastplate",&MC_lastplate,"MC_lastplate[n]/I");
    
    
    cout << "Writing " << varr->GetEntries() << " vertices" << endl;
    for (int ivtx = 0; ivtx < varr->GetEntries(); ivtx++){
        EdbVertex *vertex = (EdbVertex*)(varr->At(ivtx));
        //RIEMPIO NEW VTX FILE:
        
        if(EVERBOSE==10 || (EVERBOSE==100 && (vertex->ID()==DEBUG_VTXID|| (vertex->GetTrack(0)->Track()==DEBUG_TRKID)||vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT))) cout << endl << ivtx << "\tWriting vertex " << vertex->ID() <<"\t" << vertex->Flag() << "\t" << vertex->N() << "\t" << vertex->VX() << "\t" << vertex->VY() << "\t" << vertex->VZ() << "\t" << vplate << endl;
        
        
        if(vertex->ID()==0 && FAST!=1) vID=ivtx+900000;
        else vID=ivtx;
        vertex->SetID(vID);
        vx=vertex->VX();
        vy=vertex->VY();
        vz=vertex->VZ();
        if (EVERBOSE==100 && FAST==4) cout << " vertex VZ " << vertex->VZ() << ", vertex Z " << vertex->Z() << endl;
        n=vertex->N();
        vplate=Get_vtx_plate(vz);
        v_flag=vertex->Flag();
        
        // if(vplate<=30){
        
        
        maxaperture = vertex->MaxAperture();
        probability = vertex->V()->prob();
        for (int itrk = 0; itrk < vertex->N(); itrk++){
            EdbTrackP *track = vertex->GetTrack(itrk);
            //tarr.Add(new EdbTrackP(*track));
            
            if(EVERBOSE==10 || (EVERBOSE==100 && ((track->Track()==DEBUG_TRKID) || track->MCEvt()==DEBUG_MCEVT))) cout << "\t" << track->MCEvt() << " (" << track->MCTrack() << ")\t" << track->Track() << "\t" << track->GetSegmentFirst()->W()-70 << endl;
            
            
            IDTrack[itrk] = track->Track(); //track->GetSegment(0)->Track();//track->Track();
            //cout << track->GetSegment(0)->Track() << "\t" << track->Track() << "\t" << track->Theta() << endl;
            TrackTrack[itrk] = track->Track();
            nseg_vec[itrk] = track->N();
            npl[itrk] = track->Npl();
            Int_t zpos = vertex->GetVTa(itrk)->Zpos();
            incoming[itrk] = zpos;
            X[itrk] = track->GetSegment(0)->X();
            Y[itrk] = track->GetSegment(0)->Y();
            Z[itrk] = track->GetSegment(0)->Z();
            TX[itrk] = track->TX();
            TY[itrk] = track->TY();
            Theta[itrk] = track->Theta();
            plate[itrk] = track->GetSegment(0)->Plate();
            nholes[itrk] = track->N0();
            maxgap[itrk] = track->CheckMaxGap();
            
            nseg_S1[itrk]=0;
            npl_S1[itrk]=0;
            nholes_S1[itrk]=0;
            
            
            if(EVERBOSE==10 || (EVERBOSE==100 && ((track->Track()==DEBUG_TRKID) || track->MCEvt()==DEBUG_MCEVT))){
                cout << "\t " << track->Track() << "\t"  << track->N() << "\tzpos " << zpos << "\t" << track->MCTrack() << "\t" << endl;
                //             for (int iseg = 0; iseg < track->N(); iseg++){ //loop on segments
                //             cout << "\t " << track->GetSegment(iseg)->Plate() << "\t"  << track->GetSegment(iseg)->Z() << "\t" << track->GetSegment(iseg)->X() << "\t" << track->GetSegment(iseg)->Y() << "\t" << track->GetSegment(iseg)->TX() << "\t" << track->GetSegment(iseg)->TY() << endl;
                //             }
            }
            
            if(plate[itrk]<31){
                for (int iseg = 0; iseg < track->N(); iseg++){ //loop on segments
                    if(track->GetSegment(iseg)->Plate()<=LASTLAYER[1]){
                        nseg_S1[itrk]++;
                    }
                }
                npl_S1[itrk]=LASTLAYER[1]-plate[itrk]+1;
                nholes_S1[itrk]=npl_S1[itrk]-nseg_S1[itrk];
            }
            impactparameter[itrk] = vertex->GetVTa(itrk)->Imp();
            //            if(impactparameter[itrk]==0) impactparameter[itrk] = CalcIP(track, TVector3(TVector3(vx,vy,vz)));
            //            else cout <<  itrk << "\t" << incoming[itrk] << "\t" << track->N() << "\t" << vertex->VX() << " " << vertex->VY() << " " << vertex->VZ() << "\t" << impactparameter[itrk] << "\t" << CalcIP(vertex, itrk)  << endl; //<< "\t" << distance(vertex, track)
            Z_flag[itrk] = track->Flag();
            for(int is=0; is<track->N(); is++){
                if(track->GetSegment(is)->Plate()>LASTLAYER[1]+1){
                    MC_Charge_S2[itrk] = (int)track->GetSegment(is)->W()-70;
                    continue;
                }
            }
            
            MC_evID_first[itrk] = track->GetSegmentFirst()->MCEvt();
            MC_evID_last[itrk] = track->GetSegmentLast()->MCEvt();
            MC_trackID_first[itrk] = track->GetSegmentFirst()->MCTrack();
            MC_trackID_last[itrk] = track->GetSegmentLast()->MCTrack();
            MC_pdgcode_first[itrk] = track->GetSegmentFirst()->Aid(0);
            MC_pdgcode_last[itrk] = track->GetSegmentLast()->Aid(0);
            MC_Charge_first[itrk] = (int)track->GetSegmentFirst()->W()-70;
            MC_Charge_last[itrk] = (int)track->GetSegmentLast()->W()-70;
            MC_mother_first[itrk] = track->GetSegmentFirst()->Aid(1);
            MC_mother_last[itrk] = track->GetSegmentLast()->Aid(1);
            MC_firstplate[itrk] = track->GetSegmentFirst()->Vid(0);
            MC_lastplate[itrk] = track->GetSegmentLast()->Vid(0);
            
            if(EVERBOSE==10) cout << "\t\t" << MC_evID_first[itrk] << "\t" << MC_evID_last[itrk] << "\t" << incoming[itrk] << "\t" << track->N() << "\t" << track->Theta() << "\t" << track->GetSegmentFirst()->Plate() << "\t" << track->GetSegmentLast()->Plate() << "\t" << impactparameter[itrk] << "\t" << CalcIP(vertex, itrk) << "\t" << CalcDist(vertex, track) << "\t" << CalcDist3D(vertex, track) << endl;
        }
        new_vtxtree->Fill();
        //}if(vplate<=30)
    }
    
}

//______________________________________________________________________________
EdbVertex *AddTrackToVertex2(EdbVertexRec *vrec, EdbVertex *eVertex, EdbTrackP *eTr, int zpos){
    EdbVTA *vta = 0;
    EdbVertex *old = 0;
    EdbVertex *eWorking = 0;
    if (!eVertex)
    {
        printf("No working vertex selected!\n");
        fflush(stdout);
        return 0;
    }
    if (!eTr)
    {
        printf("No working track selected!\n");
        fflush(stdout);
        return 0;
    }
    
    if (eWorking == 0)
    {
        eWorking = new EdbVertex();
        int ntr = eVertex->N();
        int i = 0, n = 0;
        for(i=0; i<ntr; i++)
        {
            if ((vta = vrec->AddTrack(*(eWorking), (eVertex)->GetTrack(i), (eVertex)->Zpos(i))))
            {
                (eVertex)->GetTrack(i)->AddVTA(vta);
                n++;
            }
        }
        if (n < 2)
        {
            safe_delete(eWorking);
            (eVertex)->ResetTracks();
            printf("Can't create working copy of the vertex! if (n < 2)\n");
            fflush(stdout);
            return 0;
        }
        
        if (!vrec->MakeV(*eWorking))
        {
            safe_delete(eWorking);
            eVertex->ResetTracks();
            printf("Can't create working copy of the vertex! if (!MakeV(*eWorking))\n");
            fflush(stdout);
            return 0;
        }
    }
    if ((vta = vrec->AddTrack(*eWorking, eTr, zpos)))
    {
        eTr->AddVTA(vta);
        eWorking->SetID(eVertex->ID());
    }
    else
    {
        //    printf("Track not added! May be Prob < ProbMin.\n");
        //    fflush(stdout);
        
        safe_delete(eWorking);
        eVertex->ResetTracks();
        return 0;
    }
    return eWorking;
}


//---------------------------------------------------------------------

TObjArray *UnisciVertici(TObjArray *varr, EdbVertexRec *vrec){
    
    TObjArray tr_grid_end;
    int vplate=-99;
    
    cout << "UnisciVertici starting from " << varr->GetEntries() << endl;
    TObjArray *merged_varr = new TObjArray();
    
    FillVertexCells(*varr);
    LoopVertexCells_ToMerge(*varr, vrec);
    
    int arrentries_temp = varr->GetEntries();
    cout << "UnisciVertici After LoopVertexCells_ToMerge " << arrentries_temp << endl;

    if (EVERBOSE==100 && FAST==4) {
        EdbVertex* v = (EdbVertex*)(varr->At(0));
        for (int itrk=0; itrk<v->N(); itrk++) 
        {
            cout << " VTA_ZPos " << v->GetVTa(itrk)->Zpos() << " for track " << (int) (v->GetTrack(itrk)->Track()) << endl;
        }
    }
        
    for (int ivtx = 0; ivtx < arrentries_temp; ivtx++){
        
        EdbVertex *vertex = (EdbVertex*) varr->At(ivtx);
        int vID = vertex->ID();
        
        if(EVERBOSE==12 || (EVERBOSE==100 && (vID==DEBUG_VTXID  || vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT))) {
            cout << "UVa " << ivtx << ": " << vID << "\t" << vertex->Flag() << "\t" << vertex->N() << "\t" << vertex->VX() << "\t" << vertex->VY() << "\t" << vertex->VZ() << endl;
        }
        
        if(vertex->Flag() == -99 ) continue;
        
        if(EVERBOSE==12 || (EVERBOSE==100 && (vID==DEBUG_VTXID || vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT))) {
            cout << "UVb " << ivtx << ": " << vID << "\t" << vertex->Flag() << "\t" << vertex->N() << "\t" << vertex->VX() << "\t" << vertex->VY() << "\t" << vertex->VZ() << endl;
            for(int itrk=0; itrk<vertex->N(); itrk++){
                EdbTrackP *track = vertex->GetTrack(itrk);
                cout << "\t\t" << track->MCEvt() << "\t" << track->Track() << "\t" << track->N() << "\t" << track->Theta() << " VTA_Zpos " << vertex->GetVTa(itrk)->Zpos() << endl;
            }
        }
        
        
        //RIMUOVO VERTICI CHE HANNO TROPPI DAU AD ANGOLO BEAM RISPETTO AL TOTALE
        int beamdau=0;
        int countoxy=0;
        for(int itrk=0; itrk<vertex->N(); itrk++){
            EdbTrackP *track = vertex->GetTrack(itrk);
            if(track->Theta()<BEAMTHETA_SMALL && track->VertexS()) beamdau++;
            if(MC==11 && track->GetSegmentFirst()->W()-70==8) countoxy++;
        }
        if(beamdau>=2 && vertex->N()-beamdau<=2) {
            vertex->SetFlag(-99);
            if(EVERBOSE==12 || (EVERBOSE==100 && (vID==DEBUG_VTXID || vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT))) cout << " Here I Set Flag -99 (I) " << endl;
            REM_BEAMDAU++;
            //            cout << endl << vID << endl;
            //            for(int itrk=0; itrk<vertex->N(); itrk++){
            //                EdbTrackP *track = vertex->GetTrack(itrk);
            //                cout << track->MCEvt() << "\t" << track->MCTrack() << "\t" << track->GetSegmentFirst()->W()-70 << "\t" << track->N() << "\t" << track->Theta() << endl;
            //            }
            if(MC==11 && countoxy>0) REM_BEAMDAU_OK++;
        }
        
        //RIMUOVO I COSMICI PASSANTI
        int cosmicremove=0;
        TObjArray vta_toremove;
        for(int itrk=0; itrk<vertex->N(); itrk++){ // qui elimino 2 tracce che costituiscono un CR entrambe collegate a un vtx
            EdbTrackP *track = vertex->GetTrack(itrk);
            EdbVTA *vta1 = vertex->GetVTa(itrk);
            if(vta1->Zpos()==1 && track->Theta()>COSMICSCUT){
                for(int itrk2=0; itrk2<vertex->N(); itrk2++){
                    EdbTrackP *track2 = vertex->GetTrack(itrk2);
                    EdbVTA *vta2 = vertex->GetVTa(itrk2);
                    if(vta2->Zpos()==0 && track2->Theta()>COSMICSCUT){
                        if(TMath::Abs(track->Theta()-track2->Theta())<Cut_Theta_cosmic){
                            vertex->SetFlag(-99);
                            
                            if(EVERBOSE==12 || (EVERBOSE==100 && (vID==DEBUG_VTXID || vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT))) cout << " Here I Set Flag -99 (II) " << endl;
                            if(vertex->N()>3) {
                                
                                vta_toremove.Add(vta1);
                                vta_toremove.Add(vta2);
                                cosmicremove++;
                                
                                if(EVERBOSE==12 || ( EVERBOSE==100 && vID==DEBUG_VTXID)){
                                    cout << "Removing tracks " << endl << "\t" << track->MCEvt() << "\t" << track->Track() << "\t" << track->W()-70 << "\t" << track->N() << "\t" << track->Theta() << endl << "\t" << track2->MCEvt() << "\t" << track2->Track() << "\t" << track2->W()-70 << "\t" << track2->N() << "\t" << track2->Theta() << endl;
                                }
                                if(MC==11 && (track->MCEvt()==-999 || track2->MCEvt()==-999 || (track->GetSegmentLast()->MCEvt()==track2->GetSegmentFirst()->MCEvt() && track->GetSegmentLast()->MCTrack()==track2->GetSegmentFirst()->MCTrack()))) REMOVECOSMIC_OK++;
                                //break; //Mod 12 dec
                            }
                        }
                    }
                }
                if(cosmicremove==0){ //qui vado a cercare se c'è una traccia che sembra CR se ha un corrispettivo tra le tracce esistenti (non collegata a quel vtx)
                    int plate = track->GetSegmentFirst()->Plate();
                    float xy[2] = {0,0};
                    int r[2] = {1,1}; //raggio xy in cui voglio cercare: in questo caso una cella sola è sufficiente
                    int nh=4;
                    for(int iipl=plate-3; iipl<plate; iipl++){
                        nh--;
                        if(iipl<1 || iipl>PLMAX) continue;
                        xy[0] = track->GetSegmentFirst()->X()-track->GetSegmentFirst()->TX()*(track->GetSegmentFirst()->Z()-Z_LAYER[iipl]);
                        xy[1] = track->GetSegmentFirst()->Y()-track->GetSegmentFirst()->TY()*(track->GetSegmentFirst()->Z()-Z_LAYER[iipl]);
                        int ndau_grid = gridtr_DAU_end[iipl].SelectObjectsC( xy, r, tr_grid_end);
                        for(int itrk2=0; itrk2<ndau_grid; itrk2++){
                            EdbTrackP *track2 = (EdbTrackP*) tr_grid_end.At(itrk2);
                            if((TMath::Abs(track->Theta()-track2->Theta())<Cut_Theta_cosmic) && TMath::Abs(xy[0]-track2->GetSegmentLast()->X()<nh*60) && TMath::Abs(xy[1]-track2->GetSegmentLast()->Y()<nh*60)){
                                if(vertex->N()>2) {
                                    vta_toremove.Add(vta1);
                                    cosmicremove++;
                                }
                                vertex->SetFlag(-99);
                                if(EVERBOSE==12 || (EVERBOSE==100 && (vID==DEBUG_VTXID  || vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT))) cout << " Here I Set Flag -99 (III) " << endl;
                            }
                        }
                        tr_grid_end.Clear();
                    }
                }
            }
        }
        if(cosmicremove>0){
            //actually removing now
            REMOVECOSMIC++;
            for (int ivta = 0; ivta<vta_toremove.GetEntries(); ivta++){
                EdbVTA *vta = (EdbVTA*)vta_toremove.At(ivta);
                if(vertex->N()>2){
                    for(int itrk=0; itrk<vertex->N(); itrk++){
                        EdbVTA *vta1 = (EdbVTA*)vertex->GetVTa(itrk);
                        if(vta->GetTrack()->Track()==vta1->GetTrack()->Track()){
                            vertex->SetFlag(-99);//add 8dec

                            if ( (FAST==4||(FAST==3 && vertex->ID()==DEBUG_VTXID)) && EVERBOSE==100) {
                                cout << " Before RemoveVTAFromVertex " << endl;
                                for (int itrk=0; itrk<vertex->N(); itrk++) 
                                    {
                                        cout << " VTA_ZPos " << vertex->GetVTa(itrk)->Zpos() << " for track " << (int) (vertex->GetTrack(itrk)->Track()) << endl;
                                        EdbTrackP *track_check = vertex->GetTrack(itrk);
                                        EdbVertex* vtempE_check = track_check->VertexE();
                                        EdbVertex* vtempS_check = track_check->VertexS();
                                        if (EVERBOSE==100 && FAST==4) 
                                        {
                                            if (vtempE_check) cout << " Track VertexE is not NULL (ID= " << vtempE_check->ID() << ")" << endl;
                                            if (vtempS_check) cout << " Track VertexS is not NULL (ID= " << vtempS_check->ID() << ")" << endl;
                                            cout << " \n " << endl;
                                        } 

                                    }
                                cout << " \n " << endl;
                            }

                            if (FAST==4 && EVERBOSE==100) cout << " vertex ID Before RemoveVTA " << vertex->ID() << " VZ " << vertex->VZ() << " Z " << vertex->Z() <<  endl;
                            vertex = vrec->RemoveVTAFromVertex(*vertex, *vta1);
                            if (FAST==4 && EVERBOSE==100) cout << " vertex ID After RemoveVTA " << vertex->ID() << " VZ " << vertex->VZ() << " Z " << vertex->Z() <<  endl;
                            //FixVertexVTA_Zpos(vertex); not needed 

                            vertex->SetFlag(16);

                            if ((FAST==4|| (FAST==3 && vertex->ID()==DEBUG_VTXID)) && EVERBOSE==100) {
                                cout << " After RemoveVTAFromVertex " << endl;
                                for (int itrk=0; itrk<vertex->N(); itrk++) 
                                    {
                                        cout << " VTA_ZPos " << vertex->GetVTa(itrk)->Zpos() << " for track " << (int) (vertex->GetTrack(itrk)->Track()) << endl;
                                        EdbTrackP *track_check = vertex->GetTrack(itrk);
                                        EdbVertex* vtempE_check = track_check->VertexE();
                                        EdbVertex* vtempS_check = track_check->VertexS();
                                        if (EVERBOSE==100 && FAST==4) 
                                        {
                                            if (vtempE_check) cout << " Track VertexE is not NULL (ID= " << vtempE_check->ID() << ")" << endl;
                                            if (vtempS_check) cout << " Track VertexS is not NULL (ID= " << vtempS_check->ID() << ")" << endl;
                                        } 
                                    }
                                cout << " \n " << endl;
                            }

                            


                            if(EVERBOSE==12 || ( EVERBOSE==100 && vID==DEBUG_VTXID)){
                                cout << "\tTrack " << vta->GetTrack()->Track() << " removed" << endl;
                            }
                            break; //add 12 dec
                        }
                    }
                } //if(vertex->N()>2){
                else {
                    vertex->SetFlag(-99);
                    if(EVERBOSE==12 || (EVERBOSE==100 && (vID==DEBUG_VTXID  || vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT))) cout << " Here I Set Flag -99 (IV) " << endl;
                }
            } //for (int ivta = 0;
            vertex->SetID(vID);
            vplate=Get_vtx_plate(vertex->VZ());

            if(vplate>= PLMIN && vplate<=BRAGGPLATE && vertex->Flag()!=-99){
                vtxPat[vplate].Add(vertex);
                if (FAST==4 && EVERBOSE==100) {
                    cout << " Just Before Adding to varr  " << endl;
                    for (int itrk=0; itrk<vertex->N(); itrk++) 
                        {
                            cout << " VTA_ZPos " << vertex->GetVTa(itrk)->Zpos() << " for track " << (int)(vertex->GetTrack(itrk)->Track()) << endl;
                        }
                }
                varr->Add(vertex);
            }
            vta_toremove.Clear();
            if(EVERBOSE==12 || ( EVERBOSE==100 && vID==DEBUG_VTXID)){
                vplate=Get_vtx_plate(vertex->VZ());
                cout << "\tafter cosmic tracks removed\tvID " << vertex->ID() << "\t" << vertex->Flag() << "\tvplate " << vplate << endl;
                for(int itrk=0; itrk<vertex->N(); itrk++){
                    cout << "\t\t"  << vertex->GetTrack(itrk)->MCEvt() << "\t" << vertex->GetTrack(itrk)->Track() << "\t" << vertex->GetVTa(itrk)->Imp() << "\t" << vertex->Flag() << " VTA_Zpos " << vertex->GetVTa(itrk)->Zpos() << endl;
                }
            }
        } //if(cosmicremove>0)
        
        if(EVERBOSE==12 || (EVERBOSE==100 && (vID==DEBUG_VTXID|| vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT))) {
            cout << "UVc " << ivtx << ": " << vID << "\t" << vertex->Flag() << "\t" << vertex->N() << "\t" << vertex->VX() << "\t" << vertex->VY() << "\t" << vertex->VZ() << endl;
        }
        
        int vtemp_ID=-99;
        //SE UNA TRACCIA È ASSOCIATA A DUE VERTICI (Incoming e ougoing) SCELGO IL MIGLIORE
        for(int itrk=0; itrk<vertex->N(); itrk++){
            EdbTrackP *track = vertex->GetTrack(itrk);
            Int_t incoming = vertex->GetVTa(itrk)->Zpos();
            if(incoming==0 && track->VertexS()) {
                vtemp_ID = track->VertexS()->ID();
            }
            else if(incoming==1 && track->VertexE()){
                vtemp_ID = track->VertexE()->ID();
            }
            else continue;
            if(EVERBOSE==12 || (EVERBOSE==100 && (vID==DEBUG_VTXID|| vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT))){
                cout << "traccia condivisa: " << track->ID() << " " << track->Track() <<  " " << incoming << " " << vtemp_ID << endl;}
            EdbVertex *vtemp=NULL;
            if (vtemp_ID == vertex->ID()) continue; //mod vincenzo 8 feb 
            for (int ivtx1 = 0; ivtx1 < varr->GetEntries(); ivtx1++){ //cerco l'altro vertice in varr
                
                EdbVertex *vtempsearch = (EdbVertex*) varr->At(ivtx1);
                if(vtemp_ID == vtempsearch->ID() && vtempsearch->Flag()!=-99){
                    vtemp=vtempsearch;
                    for(int itrk1=0; itrk1<vtemp->N(); itrk1++){
                        if(vtemp->GetTrack(itrk1)->Track()==track->Track()){
                            if(EVERBOSE==12 || (EVERBOSE==100 && (vID==DEBUG_VTXID|| vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT)))
                            { cout << "trk found " << vertex->GetVTa(itrk)->Imp() << " " << vtemp->GetVTa(itrk1)->Imp() << endl;}
                            
                            if(vertex->GetVTa(itrk)->Imp()>vtemp->GetVTa(itrk1)->Imp() ||( (vertex->GetVTa(itrk)->Imp()==vtemp->GetVTa(itrk1)->Imp()) && incoming==1)) { // ip con vertex è peggiore-> la rimuovo da vertex. Se la devo rimuovere da vtemp aspetto che arrivi il suo turno. A parità di IP rimuovo quella entrante nel vertice
                                if(vertex->N()>2){
                                    int vertexid=vertex->ID();
                                    vertex->SetFlag(-99);
                                    EdbVTA *vta = vertex->GetVTa(itrk);
                                    vertex = vrec->RemoveVTAFromVertex(*vertex, *vta);
                                    //FixVertexVTA_Zpos(vertex);
                                    vertex->SetFlag(20);
                                    vertex->SetID(vertexid);
                                    vplate=Get_vtx_plate(vertex->VZ());
                                    if(vplate>= PLMIN && vplate<=BRAGGPLATE){ vtxPat[vplate].Add(vertex);
                                    varr->Add(vertex);
                                    }
                                }
                                else vertex->SetFlag(-99);
                            }
                            break;
                        }
                    }
                }
            }
        }
        
        
        if(EVERBOSE==12 || (EVERBOSE==100 && (vID==DEBUG_VTXID || vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT))) {
            cout << "UVd " << ivtx << ": " << vID << "\t" << vertex->Flag() << "\t" << vertex->N() << "\t" << vertex->VX() << "\t" << vertex->VY() << "\t" << vertex->VZ() << endl;
        }
    }
    cout << "new check here varr->GetEntries(): " << varr->GetEntries() << endl;
    
    for (int ivtx = 0; ivtx < varr->GetEntries(); ivtx++){
        
        EdbVertex *vertex = (EdbVertex*) varr->At(ivtx);
        int vID = vertex->ID();
        
        vplate = Get_vtx_plate(vertex->VZ());
        
        if((EVERBOSE==100 && (vID==DEBUG_VTXID || vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT))) {
            cout << "UVe " << ivtx << ": " << vID << "\t" << vertex->Flag() << "\t" << vertex->N() << "\t" << vertex->VX() << "\t" << vertex->VY() << "\t" << vertex->VZ() << endl;
        }
        if(vertex->Flag() == -99 || vplate<1 || vplate>BRAGGPLATE) {
            vertex->SetFlag(-99);
            continue;
        }
        else merged_varr->Add(vertex);

        if((EVERBOSE==100 && (vID==DEBUG_VTXID || vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT))) {
            cout << "UVf " << ivtx << ": " << vID << "\t" << vertex->Flag() << "\t" << vertex->N() << "\t" << vertex->VX() << "\t" << vertex->VY() << "\t" << vertex->VZ() << endl;
        }
        if(EVERBOSE==12 ||(EVERBOSE==100 && (vID==DEBUG_VTXID || vertex->GetTrack(0)->MCEvt()==DEBUG_MCEVT))){
            
            cout << "\tInserisco vtx " << vID << "\t" << vertex->N() << "\t" << vertex->VX() << "\t" << vertex->VY() << "\t" << vertex->VZ() << "\t" << vertex->Flag() << endl;
            for(int itrk=0; itrk<vertex->N(); itrk++){
                EdbTrackP *track = vertex->GetTrack(itrk);
                cout << "\t\t" << track->MCEvt() << "\t" << track->Track() << "\t" << track->N() << "\t" << track->Theta() << " VTA_Zpos " << vertex->GetVTa(itrk)->Zpos() << endl;
            }
        }
        
    }
    
    cout << "UnisciVertici end " << merged_varr->GetEntries() << endl;
    return merged_varr;
    
}

//---------------------------------------------------------------------

void FillVertexCells(TObjArray &arrVTX){
    //filling cells with vertices
    //    const int ncellsX = 500; //500
    //    const int ncellsY = 600; //600
    int ncellsX,ncellsY;
    if(MC>0){
        ncellsX = 450; //500;
        ncellsY = 550; //600;
    }
    else{
        ncellsX = 400; //era 450
        ncellsY = 500; //era 550
    }
    
    const int maxpercell = 30; // è il numero massimo di oggetti per ogni cella. Non restituisce errore se ce ne sono di più, per questo bisogna controllare con la stampa sotto.
    //preparing cells
    for(int i=PLMIN; i<=BRAGGPLATE; i++) gridvtx[i].InitCell(ncellsX, xmin, xmax, ncellsY, ymin, ymax, maxpercell );
    
    const int nvertices = arrVTX.GetEntries();
    //starting loop on vertices
    for (int ivtx = 0; ivtx < nvertices; ivtx++){
        
        EdbVertex *vtx = (EdbVertex*) arrVTX.At(ivtx);
        if(vtx->Flag()==-99) continue;
        float vx = vtx->VX();
        float vy = vtx->VY();
        float vz = vtx->VZ();
        
        int vtxplate = Get_vtx_plate(vz);
        if(vtxplate>=0 && vtxplate<=BRAGGPLATE){
            gridvtx[vtxplate].AddObject( vx, vy, (TObject*)vtx);
        }
        else vtx->SetFlag(-99);
    } //end of loop, checking stats
    //for(int i=PLMIN; i<=BRAGGPLATE; i++) gridvtx[i].PrintStat();
}

//---------------------------------------------------------------------

void LoopVertexCells_ToMerge(TObjArray &varr, EdbVertexRec *vrec){//EdbVertexRec *vrec){
    TObjArray *closevertices = new TObjArray();
    //start loop in plates
    for(int ipl = PLMIN; ipl<=BRAGGPLATE; ipl++){
        //start actual loop in cells of this plate
        for(int jcell=0; jcell<gridvtx[ipl].Ncell(); jcell++) {
            int nvertices = gridvtx[ipl].Bin(jcell);
            if(nvertices<2)      continue;
            // select cells with at least two vertices
            closevertices->Clear();
            //adding close vertices to the list
            for(int k=0; k<nvertices; k++)   {
                closevertices->Add(gridvtx[ipl].GetObject(jcell,k));
                
            } //end loop over vertices in cell
            
            EdbVertex *mergedvertex = MergeVertices(closevertices, vrec);
            if(EVERBOSE==9 ||( EVERBOSE==100 && mergedvertex->ID()==DEBUG_VTXID )) {cout << "LoopVertexCells_ToMerge " << mergedvertex->ID() << "\t" << mergedvertex->VX() << "\t" << mergedvertex->VY() << "\t" << mergedvertex->VZ() << "\t" << Get_vtx_plate(mergedvertex->VZ()) << "\t" << mergedvertex->Flag() << "\t" << mergedvertex->N() << endl;
            for(int i=0; i<mergedvertex->N(); i++){
                cout << "\t\t" << mergedvertex->GetTrack(i)->Track() << endl;
            }}
            
            if (mergedvertex) {
                int closeplate=Get_vtx_plate(mergedvertex->VZ());
                if(closeplate>=PLMIN&&closeplate<=BRAGGPLATE){
                varr.Add(mergedvertex);//vrec->eVTX->Add(mergedvertex);
                }
            }
        }//end loop over cells
        
    }//end loop over plates
    
}

//---------------------------------------------------------------------

EdbVertex* MergeVertices(TObjArray *varr,EdbVertexRec *vrec){
    //merge vertices in varr into one (inspired from TestVTAGroup in FEDRA)
    EdbVertex *mergedvtx = new EdbVertex();
    int nvertices = varr->GetEntries();
    //loop over tracks within vertices
    //cout<<"Trying to merge "<<nvertices<<" into one "<<endl;
    // cout << endl;
    for (int ivtx = 0; ivtx < nvertices; ivtx++){
        EdbVertex *vtx = (EdbVertex*) varr->At(ivtx);
        int ntracks = vtx->N();
        if(EVERBOSE==9 ||( EVERBOSE==100 && vtx->ID()==DEBUG_VTXID )) cout << "MergeVertices " << vtx->ID() << "\t" << vtx->VX() << "\t" << vtx->VY() << "\t" << vtx->VZ() << "\t" << Get_vtx_plate(vtx->VZ()) << "\t" << vtx->Flag() << "\t" << ntracks << endl;
        
        for (int itrack = 0; itrack < ntracks; itrack++){
            
            EdbVTA *vta = vtx->GetVTa(itrack); //get vertex to track association
            
            EdbVTA *newvta =  new EdbVTA(vta->GetTrack(),mergedvtx);
            newvta->SetZpos(vta->Zpos());
            newvta->SetFlag(2);
            mergedvtx->AddVTA(newvta);
            if(EVERBOSE==9) cout << "\t" << vtx->GetTrack(itrack)->MCEvt() << endl;
        }
    }
    int newid=0;
    //start vertex fit
    if(vrec->MakeV(*mergedvtx))
        if( mergedvtx->V() )
            if( mergedvtx->V()->valid() ) {
                //printf("prob = %f\n", newvtx->V()->prob() );
                if( mergedvtx->V()->prob() >= -1 )   // accept new N-tracks vertex
                {
                    //vertex accepted, original vertices get flag -99
                    for(int ivtx=0; ivtx < nvertices; ivtx++)    {
                        EdbVertex *vtx = (EdbVertex*) varr->At(ivtx);
                        if(ivtx==0) newid=vtx->ID();
                        vtx->SetFlag(-99);
                        if(EVERBOSE==9 ||( EVERBOSE==100 && vtx->ID()==DEBUG_VTXID )) cout << "vertex " << vtx->ID() << " deleted rigo 1684" << endl;
                    }
                    mergedvtx->SetFlag(30);
                    mergedvtx->SetID(newid+100000);
                    
                    if(EVERBOSE==9 ||( EVERBOSE==100 && mergedvtx->ID()==DEBUG_VTXID )) cout << "merged vtx: " << mergedvtx->ID() << "\t" << mergedvtx->VX() << "\t" << mergedvtx->VY() << "\t" << mergedvtx->VZ() << "\t" << mergedvtx->Flag() << "\t" << mergedvtx->N() << endl;
                    //cout<<"\t vertex accepted :)"<<endl;
                    
                    MERGEVTX++;
                    return mergedvtx;
                }
            }
    
    if(EVERBOSE==9) cout<<"but it refused. Why the vertex was refused?"<<endl;
    if(EVERBOSE==9) cout<<mergedvtx->V()<<endl;
    if( mergedvtx->V() && EVERBOSE==9) cout<<mergedvtx->V()->valid()<<endl;
    if( mergedvtx->V() ){
        if( mergedvtx->V()->valid() && EVERBOSE==9) cout<<mergedvtx->V()->prob()<<endl;
    }
    delete mergedvtx;   // discard new vertex
    return 0;
    
}


//---------------------------------------------------------------------

int Get_vtx_plate(float vz){
    
    for(int layer=0;layer<=PLMAX;layer++){
        if(vz>Z_LAYER[layer]&&vz<Z_LAYER[layer+1]){
            //cout << "vz: " << vz << "\tlayer " << layer << "\tmin " << Z_LAYER[layer] << "\t" << Z_LAYER[layer+1] << endl;
            return layer+1;
        }
    }
    cout << "plate -99: vz: " << vz << endl;
    return -99;
}

//---------------------------------------------------------------------

Float_t GetThetaKink(EdbSegP *s1, EdbSegP *s2){
    
    float tx_p  = s1->TX();
    float ty_p  = s1->TY();
    float th_p  = acos(1/sqrt(tx_p*tx_p+ty_p*ty_p+1));
    float phi_p = atan2(ty_p,tx_p);
    
    float tx_d  = s2->TX();
    float ty_d  = s2->TY();
    float th_d  = acos(1/sqrt(tx_d*tx_d+ty_d*ty_d+1));
    float phi_d = atan2(ty_d,tx_d);
    
    float th_kink = acos( cos(th_p)*cos(th_d)+sin(th_p)*sin(th_d)*cos(phi_d-phi_p));
    
    return th_kink;
}

//---------------------------------------------------------//

double CalcIP(EdbSegP *s, double x, double y, double z){
    // Calculate IP between a given segment and a given x,y,z.
    // return the IP value.
    
    double ax = s->TX();
    double ay = s->TY();
    double bx = s->X()-ax*s->Z();
    double by = s->Y()-ay*s->Z();
    
    double a;
    double r;
    double xx,yy,zz;
    
    a = (ax*(x-bx)+ay*(y-by)+1.*(z-0.))/(ax*ax+ay*ay+1.);
    xx = bx +ax*a;
    yy = by +ay*a;
    zz = 0. +1.*a;
    r = sqrt((xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z));
    
    return r;
}

//---------------------------------------------------------//

double CalcDist(EdbVertex *vertex, EdbTrackP *track){
    //prendo il primo segmento della traccia, lo proietto alla z del vertice e calcolo la distanza
    double dz=track->Z()-vertex->VZ();
    double x=track->X()-dz*track->TX();
    double y=track->Y()-dz*track->TY();
    
    double dx=x-vertex->VX();
    double dy=y-vertex->VY();
    
    return TMath::Sqrt(dx*dx+dy*dy);
    
}

//---------------------------------------------------------//

double CalcDist(EdbSegP *seg1, EdbSegP *seg2){
    //prendo il primo segmento della traccia, lo proietto alla z del vertice e calcolo la distanza
    double dz=seg1->Z()-seg2->Z();
    double x=seg1->X()-dz*seg1->TX();
    double y=seg1->Y()-dz*seg1->TY();
    
    double dx=x-seg2->X();
    double dy=y-seg2->Y();
    
    return TMath::Sqrt(dx*dx+dy*dy);
    
}

//---------------------------------------------------------//

double CalcDist3D(EdbVertex *vertex, EdbTrackP *track) {
    
    double A  = track->TX();
    double B  = track->TY();
    double C  = 1.;
    double xt = track->X();
    double yt = track->Y();
    double zt = track->Z();
    double xv = vertex->VX();
    double yv = vertex->VY();
    double zv = vertex->VZ();
    
    
    double denom = A*A + B*B + C*C;
    double nom   = A*(xt-xv) + B*(yt-yv) + C*(zt-zv);
    double rho   = nom / denom;
    
    // point of track closest to Vertex
    double x = xt - A*rho;
    double y = yt - B*rho;
    double z = zt - C*rho;
    
    double dx = x - xv;
    double dy = y - yv;
    double dz = z - zv;
    
    
    return sqrt(dx*dx + dy*dy + dz*dz);
}

//---------------------------------------------------------//

double CalcDist3D(EdbVertex *vertex, EdbSegP *seg) {
    
    double A  = seg->TX();
    double B  = seg->TY();
    double C  = 1.;
    double xt = seg->X();
    double yt = seg->Y();
    double zt = seg->Z();
    double xv = vertex->VX();
    double yv = vertex->VY();
    double zv = vertex->VZ();
    
    
    double denom = A*A + B*B + C*C;
    double nom   = A*(xt-xv) + B*(yt-yv) + C*(zt-zv);
    double rho   = nom / denom;
    
    // point of track closest to Vertex
    double x = xt - A*rho;
    double y = yt - B*rho;
    double z = zt - C*rho;
    
    double dx = x - xv;
    double dy = y - yv;
    double dz = z - zv;
    
    
    return sqrt(dx*dx + dy*dy + dz*dz);
}

//---------------------------------------------------------//

double CalcDist2DX(EdbVertex *vertex, EdbTrackP *track) {
    
    double A  = track->TX();
    double C  = 1.;
    double xt = track->X();
    double zt = track->Z();
    double xv = vertex->VX();
    double zv = vertex->VZ();
    
    
    double denom = A*A + C*C;
    double nom   = A*(xt-xv) + C*(zt-zv);
    double rho   = nom / denom;
    
    // point of track closest to Vertex
    double x = xt - A*rho;
    double z = zt - C*rho;
    
    double dx = x - xv;
    double dz = z - zv;
    
    
    return sqrt(dx*dx + dz*dz);
}

//---------------------------------------------------------//

double CalcDist2DX(EdbVertex *vertex, EdbSegP *seg) {
    
    double A  = seg->TX();
    double C  = 1.;
    double xt = seg->X();
    double zt = seg->Z();
    double xv = vertex->VX();
    double zv = vertex->VZ();
    
    
    double denom = A*A + C*C;
    double nom   = A*(xt-xv) + C*(zt-zv);
    double rho   = nom / denom;
    
    // point of track closest to Vertex
    double x = xt - A*rho;
    double z = zt - C*rho;
    
    double dx = x - xv;
    double dz = z - zv;
    
    
    return sqrt(dx*dx + dz*dz);
}

double CalcDist2DY(EdbVertex *vertex, EdbSegP *seg) {
    
    double B  = seg->TY();
    double C  = 1.;
    double yt = seg->Y();
    double zt = seg->Z();
    double yv = vertex->VY();
    double zv = vertex->VZ();
    
    
    double denom = B*B + C*C;
    double nom   = B*(yt-yv) + C*(zt-zv);
    double rho   = nom / denom;
    
    // point of track closest to Vertex
    double y = yt - B*rho;
    double z = zt - C*rho;
    
    double dy = y - yv;
    double dz = z - zv;
    
    
    return sqrt(dy*dy + dz*dz);
}

//---------------------------------------------------------//

//---------------------------------------------------------//

double CalcDist2DY(EdbVertex *vertex, EdbTrackP *track) {
    
    double B  = track->TY();
    double C  = 1.;
    double yt = track->Y();
    double zt = track->Z();
    double yv = vertex->VY();
    double zv = vertex->VZ();
    
    
    double denom = B*B + C*C;
    double nom   = B*(yt-yv) + C*(zt-zv);
    double rho   = nom / denom;
    
    // point of track closest to Vertex
    double y = yt - B*rho;
    double z = zt - C*rho;
    
    double dy = y - yv;
    double dz = z - zv;
    
    
    return sqrt(dy*dy + dz*dz);
}

//---------------------------------------------------------//

double CalcIP(EdbVertex *vertex, int itrk){
    // calculate IP between a given track and a given vertex.
    // return the IP value.
    
    if(vertex==NULL) return -1.;
    
    EdbSegP *seg;
    EdbTrackP *track = vertex->GetTrack(itrk);
    if(vertex->GetVTa(itrk)->Zpos()==1){ //se è dau prendo il primo segmento
        seg = track->GetSegmentFFirst();
    }
    else seg = track->GetSegmentFLast(); //se è beam prendo l'ultimo segmento
    return CalcIP(seg, vertex->VX(), vertex->VY(), vertex->VZ());
}

//---------------------------------------------------------//

double CalcIP(EdbVertex *vertex, EdbTrackP *track){
    // calculate IP between a given track and a given vertex.
    // return the IP value.
    
    if(vertex==NULL) return -1.;
    int itrk=-99;
    for(int i=0; i<vertex->N(); i++){
        EdbTrackP* trk = vertex->GetTrack(i);
        if(trk->GetSegmentFirst()->ID()==track->GetSegmentFirst()->ID()){
            itrk = i;
            break;
        }
    }
    if(itrk==-99) return -99;
    EdbSegP *seg;
    if(vertex->GetVTa(itrk)->Zpos()==1){ //se è dau prendo il primo segmento
        seg = track->GetSegmentFFirst();
    }
    else seg = track->GetSegmentFLast(); //se è beam prendo l'ultimo segmento
    return CalcIP(seg, vertex->VX(), vertex->VY(), vertex->VZ());
}

//---------------------------------------------------------------------


double CalcIP(EdbVertex *vertex, EdbTrackP *track, int zpos){
    // calculate IP between a given track and a given vertex.
    // return the IP value.
    
    if(vertex==NULL) return -1.;
    
    EdbSegP *seg;
    if(zpos==1){ //se è dau prendo il primo segmento
        seg = track->GetSegmentFFirst();
    }
    else seg = track->GetSegmentFLast(); //se è beam prendo l'ultimo segmento
    return CalcIP(seg, vertex->VX(), vertex->VY(), vertex->VZ());
}
//---------------------------------------------------------------------

int FindPlate(EdbTrackP *cand, int plate){
    
    for(int i=0; i<cand->N(); i++){
        if(cand->GetSegment(i)->Plate()==plate) return i;
    }
    return 1000;
}


//---------------------------------------------------------------------

void FixVertexVTA_Zpos(EdbVertex* vertex) {

    for (int itrk=0; itrk<vertex->N(); itrk++)
    {
        EdbTrackP* track = vertex->GetTrack(itrk);
        EdbVTA* vta = vertex->GetVTa(itrk);
        (track->Z() >= vertex->VZ())? vta->SetZpos(1):vta->SetZpos(0);
    }

}



// -----------------------------------------------------------------------

//float CalcIP(EdbTrackP *tr, TVector3 V){
//    //transverse distance (IP) from track to vertex (n.d.r. tranvserse with respect to beam z direction),
//    //taken from nearest upstream segment
//    float imp = -99;
//    for (int iseg = 0; iseg < tr->N(); iseg++){
//        EdbSegP *seg = tr->GetSegmentF(iseg);
//        if (seg->Z() <= V(2)){
//            float dz = V(2) - seg->Z();
//            float ipx = seg->TX() * dz + seg->X() - V(0);
//            float ipy = seg->TY() * dz + seg->Y() - V(1);
//            imp = TMath::Sqrt(pow(ipx,2)+pow(ipy,2));
//        }
//    }
//    return imp;
//
//}
////---------------------------------------------------------//
//
//double distance(EdbVertex *v, EdbTrackP *t) {
//    // distance: spatial distance between Track and Vertex
//
//    double dx = v->VX() - t->X();
//    double dy = v->VY() - t->Y();
//    double dz = v->VZ() - t->Z();
//    double tx = t->TX();
//    double ty = t->TY();
//
//    double nom = sqrt(dx*ty - dy*tx) + sqrt(dy - dz*ty) + sqrt(dz*tx - dx);
//    double denom = sqrt(tx) + sqrt(ty) + 1.;
//    return sqrt(nom/denom);
//}

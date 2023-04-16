
#define IDBRICK 2
#define EVERBOSE 100
#define DEBUG_S0_PLATE 31
#define DEBUG_S0_ID 411417
#define DEBUG_S0_PLATE_S1 25
#define DEBUG_S0_ID_S1 411412

double CalcDistMiddle(float x1, float y1, float z1, float tx1, float ty1, float x2, float y2, float z2, float tx2, float ty2);
double CalcDist(float x1, float y1, float z1, float x2, float y2, float z2, float tx1, float ty1);



int merge_offsets() {

    int S0 = 1, SL = 2;

    TFile *trkfile = TFile::Open(Form("b%06i.0.%i.%i.trk.root", IDBRICK, S0, SL));
    TTree* tracks = (TTree*) trkfile->Get("tracks");


    // End of First Section Cuts
    int PLATE_MIN=26, PLATE_MAX=31, NSEG_MIN=3;
    float X_min=35000, X_max=90000, Y_min = 25000, Y_max=75000;
    // Start of Second Section 
    int PLATE_MAX_S2 = 40, thr = 3000, START_PLATE_S2=31;

    // Merge Cuts
    float DT_MAX = 0.5, B_MAX=500;
    float MAX_B = 1000;

    float X_off_manual = 0;
    float Y_off_manual = 0;
    // Read Tracks (need XYZ s0 and sL and sfTX sfTY)
    vector<double> SL_s0_X;
    vector<double> SL_s0_Y;
    vector<double> SL_s0_Z;
    vector<double> SL_sf0_TX;
    vector<double> SL_sf0_TY;

    vector<int> SL_s0_ID;
    vector<int> SL_s0_PLATE;
    vector<int> SL_s0_flag;

    vector<double> S0_sl_X;
    vector<double> S0_sl_Y;
    vector<double> S0_sl_Z;
    vector<double> S0_sfl_TX;
    vector<double> S0_sfl_TY;

    vector<int> S0_s0_ID;
    vector<int> S0_s0_PLATE;
    vector<int> S0_sl_PLATE;
    // MC Info
    vector<int> SL_s0_MCEvt;
    vector<int> S0_s0_MCEvt;
    vector<int> SL_s0_MCtrk;
    vector<int> S0_s0_MCtrk;

    TClonesArray *segments = new TClonesArray("EdbSegP",100);
    TClonesArray *fitted_segments = new TClonesArray("EdbSegP",100);
    int nseg=0;

    tracks->SetBranchAddress("s", &segments);
    tracks->SetBranchAddress("sf", &fitted_segments);
    tracks->SetBranchAddress("nseg", &nseg);

    int N = tracks->GetEntries();

    cout << " --- Retrieving First and Last Segment Info --- " << endl;

    for (int itrk=0; itrk<N; itrk++) {

        tracks->GetEntry(itrk);
        
        if (nseg<NSEG_MIN) continue;

        EdbSegP* seg0 = (EdbSegP*)segments->At(0);
        EdbSegP* segL = (EdbSegP*)segments->At(nseg-1);

        //cout << " seg0 plate " << seg0->Plate() << endl; 

        if (seg0->Plate()<PLATE_MAX_S2 && seg0->Plate()>=START_PLATE_S2) {

            EdbSegP* segf0 = (EdbSegP*)fitted_segments->At(0);

            //cout << " here " << endl;

            SL_s0_X.push_back(segf0->X()+X_off_manual);
            SL_s0_Y.push_back(segf0->Y()+Y_off_manual);
            SL_s0_Z.push_back(segf0->Z());

            SL_sf0_TX.push_back(segf0->TX());
            SL_sf0_TY.push_back(segf0->TY());

            SL_s0_ID.push_back(seg0->ID());
	        SL_s0_PLATE.push_back(seg0->Plate());
            SL_s0_flag.push_back(seg0->Flag());
            SL_s0_MCEvt.push_back(seg0->MCEvt());
            SL_s0_MCtrk.push_back(seg0->MCTrack());
            cout << " seg flag " << (double)seg0->Flag() << endl;

        }

        if (segL->Plate()>=PLATE_MIN && segL->Plate()<=PLATE_MAX) {

            EdbSegP* segfL = (EdbSegP*)fitted_segments->At(nseg-1);

            S0_sl_X.push_back(segfL->X());
            S0_sl_Y.push_back(segfL->Y());
            S0_sl_Z.push_back(segfL->Z());

            S0_sfl_TX.push_back(segfL->TX());
            S0_sfl_TY.push_back(segfL->TY());

            S0_s0_ID.push_back(seg0->ID());
	        S0_s0_PLATE.push_back(seg0->Plate());
            S0_sl_PLATE.push_back(segL->Plate());
            S0_s0_MCEvt.push_back(seg0->MCEvt());
            S0_s0_MCtrk.push_back(seg0->MCTrack());
            //cout << " filling s0 s0 ID with " << seg0->ID() << endl;

        }

    }

    cout << " Using " << SL_s0_X.size() << " S2 tracks and " << S0_sl_X.size() << " S1 tracks " << endl;

    cout << " --- Calculating Impact Parameter --- " << endl;

    double b_back=0, b_for=0;
    float x0=0, x1=0, y0=0, y1=0, z0=0, z1=0, tx0=0, ty0=0, tx1=0, ty1=0;
    float dtx=0, dty=0, dxp=0, dyp=0, dx=0, dy=0, dxp2=0, dyp2=0, dx_middle=0, dy_middle=0, dz=0;
    int mcevt_S2=0, mcevt_S1=0, mc_trkS2=0, mc_trkS1=0;

    TStopwatch t;
    t.Start();

    TFile *outfile = TFile::Open(Form("%i_S%i_S%i_offsets.root", IDBRICK, S0, SL), "RECREATE");
    TNtuple *merge_all = new TNtuple("t_merge", "Impact Parameter between ALL selected tracks ", "b:DTX:DTY:DXp:DYp:b_back:DXp_back:DYp_back:DXp_middle:DYp_middle");
    TNtuple *merge_best = new TNtuple("merge_best", "Only Best Candidates (b+b_back / 2)", "b:s0_id:s0_plate:s0_idS1:s0_plateS1:DXp:DYp:DX_middle:DY_middle:DXp2:DYp2:s0TX:s0TY:b_back:s0_flag");
    TNtuple *merge_best_info = new TNtuple("merge_best_info", "Only Best Candidates (b+b_back / 2), other info", "s0_MCEvt:s0_MCEvtS1:s0_MCtrk:s0_MCtrkS1:ngap");
    TNtuple *merge_comp = new TNtuple("merge_comp", "Comparing b forward, b back, b middle for best b middle candidates ", "b:b_back:b_middle");

    cout << " size " << SL_s0_X.size() << endl;
    float b_min = 0;
    int plate=0, id=0, id2=0, plate2=0, save_plate=0, save_id=0, flag=0, save_mcevt=0, save_mctrk=0, ngap=0;
    float save_dxp=0, save_dyp=0, save_dx=0, save_dy=0, save_dyp2=0, save_dxp2=0, save_bback=0, save_bf=0;
    float b_middle=0, save_bmid=0;
    int lastplate=0, save_ngap=0;

    for (int i=0; i<SL_s0_X.size(); i++) {

        x0 = SL_s0_X[i];
        y0 = SL_s0_Y[i];
        z0 = SL_s0_Z[i];

        tx0 = SL_sf0_TX[i];
        ty0 = SL_sf0_TY[i];

        plate = SL_s0_PLATE[i];
        id = SL_s0_ID[i];
        flag = SL_s0_flag[i];
        mcevt_S2 = SL_s0_MCEvt[i];
        mc_trkS2 = SL_s0_MCtrk[i];

        b_min = 10000;
        save_plate=0; 
        save_id = 0;

        for (int j=0; j<S0_sl_X.size(); j++) {

            x1 = S0_sl_X[j];
            y1 = S0_sl_Y[j];
            z1 = S0_sl_Z[j];

            tx1 = S0_sfl_TX[j];
            ty1 = S0_sfl_TY[j];

            plate2 = S0_s0_PLATE[j];
            lastplate = S0_sl_PLATE[j];
            ngap = TMath::Abs(plate-lastplate);
            id2 = S0_s0_ID[j];
            mcevt_S1 = S0_s0_MCEvt[j];
            mc_trkS1 = S0_s0_MCtrk[j];

            b_back = CalcDist(x0, y0, z0, x1, y1, z1, tx0, ty0);
            b_for = CalcDist(x1, y1, z1, x0, y0, z0, tx1, ty1);
            b_middle = CalcDistMiddle(x1, y1, z1, tx1, ty1, x0, y0, z0, tx0, ty0);

            if (EVERBOSE==100 && plate2==DEBUG_S0_PLATE_S1 && id2==DEBUG_S0_ID_S1 && plate==DEBUG_S0_PLATE && id == DEBUG_S0_ID) {
                cout << " Calculated b_back: " << b_back << endl;
            }

            dtx = tx1 - tx0;
            dty = ty1 - ty0;
            dxp = x1 - (z1-z0)*tx1 - x0;
            dyp = y1 - (z1-z0)*ty1 - y0;
            dxp2 = x0 - (z0-z1)*tx0 - x1;
            dyp2 = y0 - (z0-z1)*ty0 - y1;
            dx = x1 - x0;
            dy = y1 - y0;
	        dz = z1 - z0;

            double dz_middle = TMath::Abs(dz/2);
            double x1_middle = x1+dz_middle*tx1;
            double y1_middle = y1+dz_middle*ty1;
    
            double x2_middle = x0 - dz_middle*tx0;
            double y2_middle = y0 - dz_middle*ty0;

            dx_middle = x1_middle-x2_middle;
            dy_middle = y1_middle-y2_middle;

            if (b_for<MAX_B && b_back<MAX_B) merge_all->Fill(b_for, dtx, dty, dxp, dyp, b_back, dxp2, dyp2, dx_middle, dy_middle);
            if ( (b_for+b_back)/2<b_min) { b_min = (b_for+b_back)/2; save_id = id2; save_plate = plate2; save_dxp=dxp; save_dyp=dyp; save_dx=dx_middle; save_dy=dy_middle; save_dyp2=dyp2; save_dxp2=dxp2;
                                save_bback = b_back; save_bmid=b_middle; save_bf=b_for; save_mcevt=mcevt_S1; save_mctrk=mc_trkS1; save_ngap=ngap;}

        }

        if (i%1000==0) { cout << " Completed " << 100.*i/SL_s0_X.size() << " %, Iteration time: " << t.RealTime() << " s" << endl; t.Reset(); t.Start(); }
        merge_best->Fill(b_min, id, plate, save_id, save_plate, save_dxp, save_dyp, save_dx, save_dy, save_dxp2, save_dyp2, tx0, ty0, save_bback, flag);
        merge_comp->Fill(save_bf, save_bback, b_min);
        merge_best_info->Fill(mcevt_S2, save_mcevt, mc_trkS2, save_mctrk, ngap);
        
    }

    outfile->cd();
    merge_all->Write();
    merge_best->Write();
    merge_comp->Write();
    merge_best_info->Write();
    outfile->Close();

    return 1;


}




double CalcDistMiddle(float x1, float y1, float z1, float tx1, float ty1, float x2, float y2, float z2, float tx2, float ty2){
    //prendo il primo segmento della traccia, lo proietto alla z del vertice e calcolo la distanza
    double dz=z1-z2;
    double dz_middle = TMath::Abs(dz/2);
    double x1_middle=x1+dz_middle*tx1;
    double y1_middle=y1+dz_middle*ty1;
    
    double x2_middle = x2 - dz_middle*tx2;
    double y2_middle = y2 - dz_middle*ty2;

    double dx=x1_middle-x2_middle;
    double dy=y1_middle-y2_middle;
    
    return TMath::Sqrt(dx*dx+dy*dy);
    
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

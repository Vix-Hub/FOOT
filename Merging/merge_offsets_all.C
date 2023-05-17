
#define EVERBOSE -99
#define NSTACK 7

double CalcDistMiddle(float x1, float y1, float z1, float tx1, float ty1, float x2, float y2, float z2, float tx2, float ty2);
double CalcDist(float x1, float y1, float z1, float x2, float y2, float z2, float tx1, float ty1);
float* merge_offsets(int IDBRICK=222, int S0=1, int SL=2, int NSEG_MIN=3, int n_plates_minus=5, int n_plates_plus=5);
float* merge_offsets_couple(TTree* tracks, int IDBRICK=222, int S0=1, int SL=2, int NSEG_MIN=3, int n_plates_minus=4, int n_plates_plus=9, int angle_correction=1, float *angular_offsets=0, int check=0) ;
float fit_histogram(TH1F* hist, float* offset, float &mean0, float sigma0, float sigma, float maximum, float prob_min, int N);


// merge_offsets_all calculates all offsets between each couple of section from S0 to S1
// merge_offsets is used to calculate the offsets for each couple

void merge_offsets_all(int IDBRICK=222, int S0=1, int SL=2, int check=0) {

    TString filename = Form("%i_S%i_S%i_offsets_all.root", IDBRICK, S0, SL);
    if (check) filename = Form("%i_S%i_S%i_offsets_all_check.root", IDBRICK, S0, SL);

    TFile *outfile = TFile::Open(filename, "RECREATE");
    TNtuple *off_tup = new TNtuple("merge_offsets", "Each Entry corresponds to a couple of sections (ie S1-S2...)", "Xoff:Yoff:TXoff:TYoff");

    int NSEG_MIN = 3, n_plates_minus=4, n_plates_plus=9;

    TString trackname = Form("b%06i.0.%i.%i.trk.root", IDBRICK, S0, SL);
    if (check) trackname = Form("b%06i.0.%i.%i.trk_trasl.root", IDBRICK, S0, SL);

    TFile *trkfile = TFile::Open(trackname, "READ");
    cout << " using track file " << trackname << endl;
    TTree* tracks = (TTree*) trkfile->Get("tracks");

    float *offsets = new float[4]; // DX, DY, DTX, DTY
    float *previous_offsets = new float[4];

    for (int i = 0; i < 4; i++) { offsets[i] = 0.0; previous_offsets[i] = 0.0;}

    for (int iS=S0; iS<SL; iS++) {

        if (iS==2) n_plates_minus = 9;  //ad hoc S2 choice due to lack of passive material
        else n_plates_minus = 4;

        if (iS+1 == 2) n_plates_plus = 9;
        else n_plates_plus = 4;

        n_plates_plus = 0;
        n_plates_minus = 0; //modifica temporanea

        cout << " Calculating Offsets between S" << S0 << " and S" << SL << endl;
        float *angular_offsets = merge_offsets_couple(tracks, IDBRICK, iS, iS+1, NSEG_MIN, n_plates_minus, n_plates_plus, 1, nullptr, check); //current angular offsets
        for (int i=2; i<4; i++ ) { offsets[i] += angular_offsets[i]; previous_offsets[i] = angular_offsets[i]; } //save angular offsets
        for (int i=0; i<2; i++) { offsets[i] += angular_offsets[i]; } //temp

        /*if (!check) //mod temp
        {
            float *offsets_current = merge_offsets_couple(tracks, IDBRICK, iS, iS+1, NSEG_MIN, n_plates_minus, n_plates_plus, 1, nullptr, check); //previous offsets not yet changed
            for (int i=2; i<4; i++ ) { offsets[i] = angular_offsets[i]; previous_offsets[i] = angular_offsets[i]; } //save angular offsets
            for (int i=0; i<2; i++) { offsets[i] = offsets_current[i]; }
        }*/
        
        off_tup->Fill(offsets);

    }

    outfile->cd();
    off_tup->Write();
    outfile->Close();
}


float* merge_offsets_couple(TTree* tracks, int IDBRICK=222, int S0=1, int SL=2, int NSEG_MIN=3, int n_plates_minus=4, int n_plates_plus=9, int angle_correction=1, float* angular_offsets=0, int check=0) {

    //int LASTLAYER[NSTACK+1]={0,30,66,76,83,90,110,120};
    int LASTLAYER[NSTACK+1]={1,30,66,76,83,90,120,140};
    int START_PLATE_SL = LASTLAYER[SL-1]+1;
    int PLATE_MAX_SL = START_PLATE_SL + n_plates_plus;
    int PLATE_MAX = LASTLAYER[S0];
    int PLATE_MIN = PLATE_MAX - n_plates_minus;
    cout << " S" << S0 << " section: using tracks ending between " << PLATE_MIN << " and " << PLATE_MAX << endl;
    cout << " S" << SL << " section: using tracks starting between " << START_PLATE_SL << " and " << PLATE_MAX_SL << endl;; 
    if (angle_correction) cout << " --- Calculating Angular Offsets --- " << endl;
    else cout << " --- Calculating XY Offsets --- " << endl;

    int N = tracks->GetEntries();
    float* offsets = new float[4];

    float TX_offset=0, TY_offset=0;
    if (!angle_correction && S0!=1) { TX_offset = angular_offsets[3]; TY_offset = angular_offsets[4]; }

    // Merge Cuts
    float DT_MAX = 0.5, B_MAX=500;
    float MAX_B = 1000;

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

    for (int itrk=0; itrk<N; itrk++) {

        tracks->GetEntry(itrk);
        if (nseg<NSEG_MIN) continue;

        EdbSegP* seg0 = (EdbSegP*)segments->At(0);
        EdbSegP* segL = (EdbSegP*)segments->At(nseg-1);

        if (seg0->Plate()<=PLATE_MAX_SL && seg0->Plate()>=START_PLATE_SL) {

            EdbSegP* segf0 = (EdbSegP*)fitted_segments->At(0);

            SL_s0_X.push_back(segf0->X());
            SL_s0_Y.push_back(segf0->Y());
            SL_s0_Z.push_back(segf0->Z());

            SL_sf0_TX.push_back(segf0->TX());
            SL_sf0_TY.push_back(segf0->TY());

            SL_s0_ID.push_back(seg0->ID());
	        SL_s0_PLATE.push_back(seg0->Plate());
            SL_s0_flag.push_back(seg0->Flag());
            SL_s0_MCEvt.push_back(seg0->MCEvt());
            SL_s0_MCtrk.push_back(seg0->MCTrack());

        } 

        if (segL->Plate()>=PLATE_MIN && segL->Plate()<=PLATE_MAX) {

            EdbSegP* segfL = (EdbSegP*)fitted_segments->At(nseg-1);

            S0_sl_X.push_back(segfL->X());
            S0_sl_Y.push_back(segfL->Y());
            S0_sl_Z.push_back(segfL->Z());

            S0_sfl_TX.push_back(segfL->TX()+TX_offset);
            S0_sfl_TY.push_back(segfL->TY()+TY_offset);

            S0_s0_ID.push_back(seg0->ID());
	        S0_s0_PLATE.push_back(seg0->Plate());
            S0_sl_PLATE.push_back(segL->Plate());
            S0_s0_MCEvt.push_back(seg0->MCEvt());
            S0_s0_MCtrk.push_back(seg0->MCTrack());

        }

    }

    cout << " Using " << SL_s0_X.size() << " S" << SL << " tracks and " << S0_sl_X.size() << " S" << S0 << " tracks " << endl;

    double b_back=0, b_for=0;
    float x0=0, x1=0, y0=0, y1=0, z0=0, z1=0, tx0=0, ty0=0, tx1=0, ty1=0;
    float dtx=0, dty=0, dxp=0, dyp=0, dx=0, dy=0, dxp2=0, dyp2=0, dx_middle=0, dy_middle=0, dz=0;
    int mcevt_S2=0, mcevt_S1=0, mc_trkS2=0, mc_trkS1=0;

    TStopwatch t;
    t.Start();

    TString outname = Form("merge_offsets/%i_S%i_S%i_offsets.root", IDBRICK, S0, SL);
    if (check) outname = Form("merge_offsets/%i_S%i_S%i_offsets_check.root", IDBRICK, S0, SL);

    TFile *outfile = TFile::Open(outname, "RECREATE");
    TNtuple *merge_all = new TNtuple("t_merge", "Impact Parameter between ALL selected tracks ", "b:DTX:DTY:DXp:DYp:b_back:DXp_back:DYp_back:DXp_middle:DYp_middle:X_S1:Y_S1:X_S2:Y_S2");
    TNtuple *merge_best = new TNtuple("merge_best", "Only Best Candidates (b+b_back / 2)", "b:s0_id:s0_plate:s0_idS1:s0_plateS1:DXp:DYp:DX_middle:DY_middle:DXp2:DYp2:s0TX:s0TY:b_back:s0_flag");
    TNtuple *merge_best_info = new TNtuple("merge_best_info", "Only Best Candidates (b+b_back / 2), other info", "s0_MCEvt:s0_MCEvtS1:s0_MCtrk:s0_MCtrkS1:ngap:s0X:s0Y:DTX:DTY");
    TNtuple *merge_comp = new TNtuple("merge_comp", "Comparing b forward, b back, b middle for best b middle candidates ", "b:b_back:b_middle");

    cout << " size " << SL_s0_X.size() << endl;
    float b_min = 0;
    int plate=0, id=0, id2=0, plate2=0, save_plate=0, save_id=0, flag=0, save_mcevt=0, save_mctrk=0, ngap=0;
    float save_dxp=0, save_dyp=0, save_dx=0, save_dy=0, save_dyp2=0, save_dxp2=0, save_bback=0, save_bf=0, save_dtx=0, save_dty=0;
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

            if (b_for<MAX_B && b_back<MAX_B) {merge_all->Fill(b_for, dtx, dty, dxp, dyp, b_back, dxp2, dyp2, dx_middle, dy_middle, x1, y1, x0, y0); }
            if ( (b_for+b_back)/2<b_min) { b_min = (b_for+b_back)/2; save_id = id2; save_plate = plate2; save_dxp=dxp; save_dyp=dyp; save_dx=dx_middle; save_dy=dy_middle; save_dyp2=dyp2; save_dxp2=dxp2;
                                save_bback = b_back; save_bmid=b_middle; save_bf=b_for; save_mcevt=mcevt_S1; save_mctrk=mc_trkS1; save_ngap=ngap; save_dtx=dtx; save_dty=dty;}

        }

        if (i%5000==0) { cout << " Completed " << 100.*i/SL_s0_X.size() << " %, Iteration time: " << t.RealTime() << " s" << endl; t.Reset(); t.Start(); }
        merge_best->Fill(save_bf, id, plate, save_id, save_plate, save_dxp, save_dyp, save_dx, save_dy, save_dxp2, save_dyp2, tx0, ty0, save_bback, flag);
        merge_comp->Fill(save_bf, save_bback, b_min);
        merge_best_info->Fill(mcevt_S2, save_mcevt, mc_trkS2, save_mctrk, ngap, x0, y0, save_dtx, save_dty);
        
    }

    outfile->cd();
    merge_all->Write();
    merge_best->Write();
    merge_comp->Write();
    merge_best_info->Write();

    // access the files to determine the offsets /// avoid using different files!
    TNtuple *t_merge = merge_all;
    float Xoff=0, Yoff=0, TXoff=0, TYoff=0;

    TH1F *h_DXp = new TH1F("h_DXp", "", 300, -300, 300);
    TH1F *h_DYp = new TH1F("h_DYp", "", 300, -300, 300);
    TH1F *h_DTX = new TH1F("h_DTX", "", 500, -0.1, 0.1);
    TH1F *h_DTY = new TH1F("h_DTY", "", 500, -0.1, 0.1);

    t_merge->Draw("DXp>>h_DXp");
    t_merge->Draw("DYp>>h_DYp");
    t_merge->Draw("DTX>>h_DTX");
    t_merge->Draw("DTY>>h_DTY");

    //float fit_histogram(TH1F* hist, float* offset, float &mean0, float sigma0, float sigma, float maximum, float prob_min, int N)
    float fit_mean = 0;
    float fit_sigma = fit_histogram(h_DXp, offsets, fit_mean, 5, 20, 100, 0.001, -99); //-99 non modifica offsets
    /*h_DXp->Delete();
    TH1F *h_DXp_new = new TH1F("h_DXp_new", "", 100, fit_mean-fit_sigma, fit_mean+fit_sigma);
    t_merge->Draw("DXp>>h_DXp_new");
    fit_sigma = fit_histogram(h_DXp, offsets, fit_mean, fit_sigma, fit_sigma, 100, 0.001, 0);*/

    int bin_max = 0;
    float binc = 0, fit_prob=0, mean=0;

    /*int bin_max = h_DXp->GetMaximumBin();
    float binc = h_DXp->GetXaxis()->GetBinCenter(bin_max);
    TF1 *g1 = new TF1("g1", "gaus(0)", binc-50, binc+50);
    g1->SetParameter(0, 100);
    g1->SetParameter(1, binc);
    g1->SetParameter(2, 5);
    h_DXp->Fit("g1", "L", "", binc-10, binc+10);
    float fit_prob = g1->GetProb();
    float mean = g1->GetParameter(1);
    float fit_sigma = g1->GetParameter(2);
    if (fit_prob>0.0001)  offsets[0] = mean; 
    else { offsets[0] = binc; cout << " Low Fit Probability: " << fit_prob << endl; }*/

    // --------------- DY

    bin_max = h_DYp->GetMaximumBin();
    binc = h_DYp->GetXaxis()->GetBinCenter(bin_max);
    TF1 *g2 = new TF1("g2", "gaus(0)", binc-50, binc+50);
    g2->SetParameter(0, 100);
    g2->SetParameter(1, binc);
    g2->SetParameter(2, 5);
    h_DYp->Fit("g2", "L", "", binc-10, binc+10);
    fit_prob = g2->GetProb();
    fit_sigma = g2->GetParameter(2);
    mean = g2->GetParameter(1);
    if (fit_prob>0.0001) offsets[1] = mean;
    else { offsets[1] = binc; cout << " Low Fit Probability: " << fit_prob << endl; }


    // --------------- DTX

    bin_max = h_DTX->GetMaximumBin();
    binc = h_DTX->GetXaxis()->GetBinCenter(bin_max);
    TF1 *g3 = new TF1("g3", "gaus(0)", binc-0.02, binc+0.02);
    g3->SetParameter(0, 100);
    g3->SetParameter(1, binc);
    g3->SetParameter(2, 5);
    h_DTX->Fit("g3", "L", "", binc-0.01, binc+0.01);
    fit_prob = g3->GetProb();
    fit_sigma = g3->GetParameter(2);
    mean = g3->GetParameter(1);
    if (fit_prob>0.0001) offsets[2] = mean;
    else { offsets[2] = binc; cout << " Low Fit Probability: " << fit_prob << endl; }


    // --------------- DTY

    bin_max = h_DTY->GetMaximumBin();
    binc = h_DTY->GetXaxis()->GetBinCenter(bin_max);
    TF1 *g4 = new TF1("g4", "gaus(0)", binc-0.02, binc+0.02);
    g4->SetParameter(0, 100);
    g4->SetParameter(1, binc);
    g4->SetParameter(2, 5);
    h_DTY->Fit("g4", "L", "", binc-0.01, binc+0.01);
    fit_prob = g4->GetProb();
    fit_sigma = g4->GetParameter(2);
    mean = g4->GetParameter(1);
    if (fit_prob>0.0001) offsets[3] = mean;
    else { offsets[3] = binc; cout << " Low Fit Probability: " << fit_prob << endl; }


    TCanvas *c_DXp = new TCanvas();
    c_DXp->cd();
    h_DXp->SetTitle("DXp;DXp[#mum];Entries");
    h_DXp->Draw("colz");
    c_DXp->Write("c_DXp");
    c_DXp->Close();

    TCanvas *c_DYp = new TCanvas();
    c_DYp->cd();
    h_DYp->SetTitle("DYp;DYp[#mum];Entries");
    h_DYp->Draw("colz");
    c_DYp->Write("c_DYp");
    c_DYp->Close();

    TCanvas *c_DTX = new TCanvas();
    c_DTX->cd();
    h_DTX->SetTitle("DTX;DTX;Entries");
    h_DTX->Draw("colz");
    c_DTX->Write("c_DTX");
    c_DTX->Close();

    TCanvas *c_DTY = new TCanvas();
    c_DTY->cd();
    h_DTY->SetTitle("DTY;DTY;Entries");
    h_DTY->Draw("colz");
    c_DTY->Write("c_DTY");
    c_DTY->Close();

    // combine the canvas into one?
    TCanvas *c_merged = new TCanvas("report","offsets report",800,800);
    c_merged->Divide(2,2);
    c_merged->cd(1);
    h_DXp->DrawClone();
    c_merged->cd(2);
    h_DYp->DrawClone();
    c_merged->cd(3);
    h_DTX->DrawClone();
    c_merged->cd(4);
    h_DTY->DrawClone();
   
    c_merged->Write("report");

    outfile->Close();
    outfile->Delete();

    return offsets;


}





float* merge_offsets(int IDBRICK=222, int S0=1, int SL=2, int NSEG_MIN=3, int n_plates_minus=4, int n_plates_plus=9) {

    //int LASTLAYER[NSTACK+1]={0,30,66,76,83,90,110,120};
    int LASTLAYER[NSTACK+1]={1,30,66,76,83,90,120,140};
    int START_PLATE_SL = LASTLAYER[SL-1]+1;
    int PLATE_MAX_SL = START_PLATE_SL + n_plates_plus;
    int PLATE_MAX = LASTLAYER[S0];
    int PLATE_MIN = PLATE_MAX - n_plates_minus;
    cout << " S" << S0 << " section: using tracks ending between " << PLATE_MIN << " and " << PLATE_MAX << endl;
    cout << " S" << SL << " section: using tracks starting between " << START_PLATE_SL << " and " << PLATE_MAX_SL << endl;; 

    TFile *trkfile = TFile::Open(Form("b%06i.0.%i.%i.trk.root", IDBRICK, S0, SL));
    cout << " using track file " << Form("b%06i.0.%i.%i.trk.root", IDBRICK, S0, SL) << endl;
    TTree* tracks = (TTree*) trkfile->Get("tracks");
    int N = tracks->GetEntries();

    float* offsets = new float[4];

    // Merge Cuts
    float DT_MAX = 0.5, B_MAX=500;
    float MAX_B = 1000;

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

    for (int itrk=0; itrk<N; itrk++) {

        tracks->GetEntry(itrk);
        if (nseg<NSEG_MIN) continue;

        EdbSegP* seg0 = (EdbSegP*)segments->At(0);
        EdbSegP* segL = (EdbSegP*)segments->At(nseg-1);

        if (seg0->Plate()<=PLATE_MAX_SL && seg0->Plate()>=START_PLATE_SL) {

            EdbSegP* segf0 = (EdbSegP*)fitted_segments->At(0);

            SL_s0_X.push_back(segf0->X());
            SL_s0_Y.push_back(segf0->Y());
            SL_s0_Z.push_back(segf0->Z());

            SL_sf0_TX.push_back(segf0->TX());
            SL_sf0_TY.push_back(segf0->TY());

            SL_s0_ID.push_back(seg0->ID());
	        SL_s0_PLATE.push_back(seg0->Plate());
            SL_s0_flag.push_back(seg0->Flag());
            SL_s0_MCEvt.push_back(seg0->MCEvt());
            SL_s0_MCtrk.push_back(seg0->MCTrack());

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

        }

    }

    cout << " Using " << SL_s0_X.size() << " S" << SL << " tracks and " << S0_sl_X.size() << " S" << S0 << " tracks " << endl;
    cout << " --- Calculating Impact Parameter --- " << endl;

    double b_back=0, b_for=0;
    float x0=0, x1=0, y0=0, y1=0, z0=0, z1=0, tx0=0, ty0=0, tx1=0, ty1=0;
    float dtx=0, dty=0, dxp=0, dyp=0, dx=0, dy=0, dxp2=0, dyp2=0, dx_middle=0, dy_middle=0, dz=0;
    int mcevt_S2=0, mcevt_S1=0, mc_trkS2=0, mc_trkS1=0;

    TStopwatch t;
    t.Start();

    TFile *outfile = TFile::Open(Form("%i_S%i_S%i_offsets.root", IDBRICK, S0, SL), "RECREATE");
    TNtuple *merge_all = new TNtuple("t_merge", "Impact Parameter between ALL selected tracks ", "b:DTX:DTY:DXp:DYp:b_back:DXp_back:DYp_back:DXp_middle:DYp_middle:X_S1:Y_S1:X_S2:Y_S2");
    TNtuple *merge_best = new TNtuple("merge_best", "Only Best Candidates (b+b_back / 2)", "b:s0_id:s0_plate:s0_idS1:s0_plateS1:DXp:DYp:DX_middle:DY_middle:DXp2:DYp2:s0TX:s0TY:b_back:s0_flag");
    TNtuple *merge_best_info = new TNtuple("merge_best_info", "Only Best Candidates (b+b_back / 2), other info", "s0_MCEvt:s0_MCEvtS1:s0_MCtrk:s0_MCtrkS1:ngap:s0X:s0Y:DTX:DTY");
    TNtuple *merge_comp = new TNtuple("merge_comp", "Comparing b forward, b back, b middle for best b middle candidates ", "b:b_back:b_middle");

    cout << " size " << SL_s0_X.size() << endl;
    float b_min = 0;
    int plate=0, id=0, id2=0, plate2=0, save_plate=0, save_id=0, flag=0, save_mcevt=0, save_mctrk=0, ngap=0;
    float save_dxp=0, save_dyp=0, save_dx=0, save_dy=0, save_dyp2=0, save_dxp2=0, save_bback=0, save_bf=0, save_dtx=0, save_dty=0;
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

            if (b_for<MAX_B && b_back<MAX_B) {merge_all->Fill(b_for, dtx, dty, dxp, dyp, b_back, dxp2, dyp2, dx_middle, dy_middle, x1, y1, x0, y0); }
            if ( (b_for+b_back)/2<b_min) { b_min = (b_for+b_back)/2; save_id = id2; save_plate = plate2; save_dxp=dxp; save_dyp=dyp; save_dx=dx_middle; save_dy=dy_middle; save_dyp2=dyp2; save_dxp2=dxp2;
                                save_bback = b_back; save_bmid=b_middle; save_bf=b_for; save_mcevt=mcevt_S1; save_mctrk=mc_trkS1; save_ngap=ngap; save_dtx=dtx; save_dty=dty;}

        }

        if (i%5000==0) { cout << " Completed " << 100.*i/SL_s0_X.size() << " %, Iteration time: " << t.RealTime() << " s" << endl; t.Reset(); t.Start(); }
        merge_best->Fill(save_bf, id, plate, save_id, save_plate, save_dxp, save_dyp, save_dx, save_dy, save_dxp2, save_dyp2, tx0, ty0, save_bback, flag);
        merge_comp->Fill(save_bf, save_bback, b_min);
        merge_best_info->Fill(mcevt_S2, save_mcevt, mc_trkS2, save_mctrk, ngap, x0, y0, save_dtx, save_dty);
        
    }

    outfile->cd();
    merge_all->Write();
    merge_best->Write();
    merge_comp->Write();
    merge_best_info->Write();
    outfile->Close();

    outfile->Delete();

    // access the files to determine the offsets

    TFile *result_file = TFile::Open(Form("%i_S%i_S%i_offsets.root", IDBRICK, S0, SL), "READ");
    TNtuple *t_merge = (TNtuple*) result_file->Get("t_merge");
    float Xoff=0, Yoff=0, TXoff=0, TYoff=0;

    TFile *fit_file = TFile::Open(Form("%i_S%i_S%i_offsets_fit.root", IDBRICK, S0, SL), "RECREATE");

    TH1F *h_DXp = new TH1F("h_DXp", "", 600, -300, 300);
    TH1F *h_DYp = new TH1F("h_DYp", "", 600, -300, 300);
    TH1F *h_DTX = new TH1F("h_DTX", "", 1000, -0.1, 0.1);
    TH1F *h_DTY = new TH1F("h_DTY", "", 600, -0.1, 0.1);

    t_merge->Draw("DXp>>h_DXp");
    t_merge->Draw("DYp>>h_DYp");
    t_merge->Draw("DTX>>h_DTX");
    t_merge->Draw("DTY>>h_DTY");

    int bin_max = h_DXp->GetMaximumBin();
    float binc = h_DXp->GetXaxis()->GetBinCenter(bin_max);
    TF1 *g1 = new TF1("g1", "gaus(0)", binc-50, binc+50);
    g1->SetParameter(0, 100);
    g1->SetParameter(1, binc);
    g1->SetParameter(2, 5);
    h_DXp->Fit("g1", "L", "", binc-15, binc+15);
    float fit_prob = g1->GetProb();
    float mean = g1->GetParameter(1);
    if (fit_prob>0.0001)  offsets[0] = mean; 
    else { offsets[0] = binc; cout << " Low Fit Probability: " << fit_prob << endl; }

    bin_max = h_DYp->GetMaximumBin();
    binc = h_DYp->GetXaxis()->GetBinCenter(bin_max);
    TF1 *g2 = new TF1("g2", "gaus(0)", binc-50, binc+50);
    g2->SetParameter(0, 100);
    g2->SetParameter(1, binc);
    g2->SetParameter(2, 5);
    h_DYp->Fit("g2", "L", "", binc-15, binc+15);
    fit_prob = g2->GetProb();
    mean = g2->GetParameter(1);
    if (fit_prob>0.0001) offsets[1] = mean;
    else { offsets[1] = binc; cout << " Low Fit Probability: " << fit_prob << endl; }

    bin_max = h_DTX->GetMaximumBin();
    binc = h_DTX->GetXaxis()->GetBinCenter(bin_max);
    TF1 *g3 = new TF1("g3", "gaus(0)", binc-0.02, binc+0.02);
    g3->SetParameter(0, 100);
    g3->SetParameter(1, binc);
    g3->SetParameter(2, 5);
    h_DTX->Fit("g3", "L", "", binc-0.01, binc+0.01);
    fit_prob = g3->GetProb();
    mean = g3->GetParameter(1);
    if (fit_prob>0.0001) offsets[2] = mean;
    else { offsets[2] = binc; cout << " Low Fit Probability: " << fit_prob << endl; }

    bin_max = h_DTY->GetMaximumBin();
    binc = h_DTY->GetXaxis()->GetBinCenter(bin_max);
    TF1 *g4 = new TF1("g4", "gaus(0)", binc-0.02, binc+0.02);
    g4->SetParameter(0, 100);
    g4->SetParameter(1, binc);
    g4->SetParameter(2, 5);
    h_DTY->Fit("g4", "L", "", binc-0.01, binc+0.01);
    fit_prob = g4->GetProb();
    mean = g4->GetParameter(1);
    if (fit_prob>0.0001) offsets[3] = mean;
    else { offsets[3] = binc; cout << " Low Fit Probability: " << fit_prob << endl; }

    fit_file->cd();

    TCanvas *c_DXp = new TCanvas();
    c_DXp->cd();
    h_DXp->SetTitle("DXp;DXp[#mum];Entries");
    h_DXp->Draw("colz");
    c_DXp->Write("c_DXp");
    c_DXp->Close();

    TCanvas *c_DYp = new TCanvas();
    c_DYp->cd();
    h_DYp->SetTitle("DYp;DYp[#mum];Entries");
    h_DYp->Draw("colz");
    c_DYp->Write("c_DYp");
    c_DYp->Close();

    TCanvas *c_DTX = new TCanvas();
    c_DTX->cd();
    h_DTX->SetTitle("DTX;DTX;Entries");
    h_DTX->Draw("colz");
    c_DTX->Write("c_DTX");
    c_DTX->Close();

    TCanvas *c_DTY = new TCanvas();
    c_DTY->cd();
    h_DTY->SetTitle("DTY;DTY;Entries");
    h_DTY->Draw("colz");
    c_DTY->Write("c_DTY");
    c_DTY->Close();

    fit_file->Close();

    return offsets;


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


float fit_histogram(TH1F* hist, float* offset, float &mean0, float sigma0, float sigma, float maximum, float prob_min, int N) {
  // Find the bin with the maximum content
  int bin_max = hist->GetMaximumBin();
  float binc = hist->GetXaxis()->GetBinCenter(bin_max);

  // Fit a Gaussian function to the histogram in the vicinity of the maximum bin
  TF1* g1 = new TF1("g1", "gaus(0)", binc - sigma, binc + sigma);
  g1->SetParameter(0, maximum);
  g1->SetParameter(1, binc);
  g1->SetParameter(2, sigma0);
  hist->Fit("g1", "L", "", binc - sigma, binc + sigma);

  // Get the fit probability, mean, and sigma
  float fit_prob = g1->GetProb();
  float mean = g1->GetParameter(1);
  float fit_sigma = g1->GetParameter(2);

  // Set the offset to the fitted mean if the fit probability is above a threshold
  if (fit_prob > prob_min && N>=0) {
    offset[N] = mean;
    cout << " Good Fit Result, using " << mean << " instead of " << binc << endl;
  } else {
    offset[N] = binc;
    cout << " Low Fit Probability: " << fit_prob << endl;
  }

  // Delete the Gaussian function
  delete g1;
  
  mean0 = mean;
  return fit_sigma;
}

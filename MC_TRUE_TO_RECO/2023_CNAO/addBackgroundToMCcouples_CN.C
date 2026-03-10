/* Questo script aggiunge tracce di background dai dati al MC
 
 scp /Users/Giuliana/Desktop/Uni/FOOT/TEST_RECO_MIC/addBackgroundToMCcouples.C scanner@nusrv9.na.infn.it:/home/scanner/foot/RECO_MC_TEST_giul/b000028
 
 root -l
 .L addBackgroundToMCcouples.C
 addBackgroundToMCcouples(#FROMPLATE, #TOPLATE)
 
 e dopo far girare:
 source merge_background_MC.sh
 */

//NB Da S3 in poi modificare in modo che continui a prendere il fondo random da S1

#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include "TH1.h"
#include "TCanvas.h"
#include <fstream>
#include <sstream>
#include <vector>
#include "TFile.h"


#define CN 7 //2
#define BRICK_MC 18 //25
#define BRICK_DATA 777 //222
#define LASTLAYER_S1 33
#define LASTLAYER_S2 71
#define ADD_LONGCR 1
#define ADD_RANDOMCR 0
#define SHIFT 2 //number of plates to be shifted if there are different geometries in data and MC
#define eVERBOSE 0

const int MISSING_PLATES[3] = {27, 44, 62};

TRandom3 *grandom = new TRandom3(0);
void addBackgroundToMCcouples_plate(int plate);

void addBackgroundToMCcouples_V(int fromplate=1, int toplate=140){
    
    for(int plate=fromplate; plate<=toplate; plate++){
        addBackgroundToMCcouples_plate(plate);
    }
    
}
void smearing_theta (Float_t &TX, Float_t &TY, const float angres);

void addBackgroundToMCcouples_plate(int plate) {

    for (int i = 0; i < 3; ++i) {
        if (plate == MISSING_PLATES[i]) {
            std::cout << "Plate " << plate << " is in MISSING_PLATES. Skipping." << std::endl;
            return;
        }
    }
    
    char path[200]; // Definisci un buffer di dimensione sufficiente per contenere il percorso del file
    char tmpl_path[200];

    snprintf(tmpl_path, sizeof(tmpl_path), "/home/foot/foot_sdb/2023_CNAO/CN%d/invert_angles/b%06i/p%03i/%i.%i.0.0.cp.root", 7, 777, 3, 777, 3); //temporary fix because of different couples structure
    TFile* tmpl_file = TFile::Open(tmpl_path, "READ");
    TTree* tmpl_tree = (TTree*)tmpl_file->Get("couples");

    char newpath[100], longCRpath[100];
    TString name = "rand";
	if (ADD_LONGCR) name = "bisDZ2";    

    cout << "Making Background couples for plate " << plate << endl;
    sprintf(newpath,"/home/foot/foot_sdb/2023_CNAO/CN7_MC/b%06i/p%03i/%i%s.%i.0.0.cp.root", BRICK_MC, plate-SHIFT, BRICK_MC, name.Data(), plate-SHIFT);
    
    float cutMIN;
    float cutMAX;
    
    if(CN==7){
        cutMIN=25000;
        cutMAX=65000;
    }
    else if(CN==8){
        cutMIN=40000; //40000 
        cutMAX=80000; // 77000
        
    }

    TFile *file = NULL;
    std::stringstream ss;
    
    // Array dei parametri in cascata
    // std::vector<int> parameters = {BRICK_DATA, 111, 222, 333, 4};
    
    std::vector<int> parameters;
    parameters.push_back(BRICK_DATA);
    parameters.push_back(111);
    parameters.push_back(222);
    parameters.push_back(333);
    parameters.push_back(4);
    
    
    // Itera sui parametri
    for (Long64_t i = 0; i < parameters.size(); i++) {
        int parameter = parameters[i];
        ss.str("");
        int tempCN=0;
        if (parameter >= 100) {
            tempCN = parameter / 100;  // Ottieni la prima cifra di parameter
        } else if (parameter >= 10) {
            tempCN = parameter / 10;  // Ottieni la prima cifra di parameter
        } else {
            tempCN = parameter;  // parameter ha una sola cifra, assegna direttamente a CN
        }
        
        snprintf(path, sizeof(path), "/home/foot/foot_sdb/2023_CNAO/CN%d/invert_angles/b%06i/p%03i/%i.%i.0.0.cp.root", tempCN, parameter, plate, parameter, plate);

        std::ifstream ifile(path, std::ifstream::in);
        if (ifile) { // Verifica se il file esiste
            ifile.close(); // Chiudi il file prima di aprirlo con ROOT
            file = new TFile(path, "READ");
            break; // Esci dal ciclo se il file è stato trovato
        }
    }
    
    // Se il file non è stato trovato, stampa un messaggio di errore
    if (file == NULL) {
        std::cout << "Nessun file trovato con i parametri forniti." << std::endl;
    }
    
    
    
    
    //sprintf(path,"/home/scanner/foot/2019_CN/CN%d/b%06i/b%06i.0.1.7.trk.root", CN, BRICK);
    //sprintf(path,"/home/scanner/foot/2019_CN/CN%d/b000333/b000333.0.1.2.trk_merged_new.root", CN, BRICK);
    
    TFile *file_trk = NULL;
    
    if (file_trk == NULL) {
        ss.str("");
        //snprintf(path, sizeof(path), "/home/scanner/foot/2019_CN/CN%d/b%06i/b%06i.0.1.2.trk_merged_new_2.root", CN, BRICK_DATA, BRICK_DATA);
        //snprintf(longCRpath, sizeof(longCRpath), "/home/scanner/foot/2019_CN/CN2/b000222/prove_V/aus_tracks_sigma80.root");//, CN, BRICK_DATA, BRICK_DATA);
        if (CN == 7) snprintf(longCRpath, sizeof(longCRpath), "/home/foot/foot_sdb/2023_CNAO/CN7/invert_angles/b000777/aus_tracks_mod2.root");//, CN, BRICK_DATA, BRICK_DATA);
        if (CN == 8) snprintf(longCRpath, sizeof(longCRpath), "/home/scanner/foot/2023_CNAO/CN8/b000888/aus_tracks_mod2.root");//, CN, BRICK_DATA, BRICK_DATA);
        

        std::ifstream ifile2(longCRpath, std::ifstream::in);
        if (ifile2) { // Verifica se il file esiste
            ifile2.close(); // Chiudi il file prima di aprirlo con ROOT
            file_trk = new TFile(longCRpath, "READ");
        }
    }
    

    
    
    TTree *tree_trk = (TTree*)file_trk->Get("tracks");
    
    EdbSegP *tr=0, *seg=0, *lastseg=0, *firstseg=0;
    TClonesArray *segt  = new TClonesArray("EdbSegP");
    int nseg=0, npl=0, n0=0, in_vtx=0;
    float avg_volume_S1=0;
    tree_trk->SetBranchAddress("t.",&tr);
    tree_trk->SetBranchAddress("s",&segt); //array of measured segments
    tree_trk->SetBranchAddress("nseg",&nseg); //total number of segments
    tree_trk->SetBranchAddress("npl",&npl); //total number of plates passed by track
	tree_trk->SetBranchAddress("avg_volume_S1", &avg_volume_S1);
    
    
    
    TTree *tree = (TTree*)file->Get("couples");
    TFile* newfile = new TFile(newpath, "RECREATE");
    if (!newfile->IsOpen()) {
        std::cerr << "Errore: impossibile aprire il file " << newpath << std::endl;
        // Puoi inserire qui ulteriori istruzioni per gestire l'errore
        return; // o return false se la funzione ha un valore di ritorno
    }
    // Continua con il codice per il file aperto correttamente
    
    
    TTree *newtree = tmpl_tree->CloneTree(0);
    
    cout << "New path for bkg couples: " << newpath << endl;
    
    EdbSegP *s=0;
    tree->SetBranchAddress("s.",&s);
    int COUNT=0, COUNT_COSMIC=0;
    int nentries_trk = tree_trk->GetEntries();
    int nentries = tree->GetEntries();
    cout << "Number of entries: " << nentries << endl;
    Float_t X_coord=0, Y_coord=0;
    Float_t Xmin=0, Xmax=0;
    Float_t Ymin=0, Ymax=0;

    tree->Draw("s.eX >> htmp", "", "goff");
    TH1F *htmp = (TH1F*)gDirectory->Get("htmp");
    Xmin = htmp->GetXaxis()->GetXmin() + 2000;
    Xmax = htmp->GetXaxis()->GetXmax() + 2000;
    cout << "Xmin: " << Xmin << "\tXmax: " << Xmax << endl;

    tree->Draw("s.eY >> htmpy", "", "goff");
    TH1F *htmpy = (TH1F*)gDirectory->Get("htmpy");
    Ymin = htmpy->GetXaxis()->GetXmin() + 2000;
    Ymax = htmpy->GetXaxis()->GetXmax() + 2000;
    cout << "Ymin: " << Ymin << "\tYmax: " << Ymax << endl;
    
    
    float shiftMin = cutMIN-Xmin;
    float shiftMax = Xmax-cutMAX;
    float shift = TMath::Max(shiftMin,shiftMax);
    float gap_size = (cutMAX-shiftMax)-(cutMIN+shiftMin);
    if (gap_size < 0) gap_size = cutMAX - cutMIN;

    int num_copies = std::ceil(gap_size / shift);

    cout << "shiftMin = " << cutMIN << "-"<< Xmin << "=" << shiftMin << "\tshiftMax = "<< Xmax <<"-"<<cutMAX << "=" << shiftMax << "\tshift = " << shift << endl;
    cout << "gap_size " << gap_size << "\tnum_copies "<< num_copies << endl;
    
    int countmax=0;
    int countmin=0;
    int n=0;

	EdbSegP *newcp = new EdbSegP();
    newtree->SetBranchAddress("s.", &newcp);
	
    if(ADD_LONGCR==1){ //modV to force the same DZ as in DATA
        cout << "Long cosmic rays taken from: " << longCRpath << endl;

        float thickS1 = 2320; float thickS2 = 320;
		float theta0 = 0;
        if (CN==1) thickS1 = 1320;
        
        for (Long64_t i=0;i<nentries_trk; i++) {
            tree_trk->GetEntry(i);
            bool found_cosmic_signal = false;
            lastseg=(EdbSegP*)segt->At(nseg-1);
			firstseg=(EdbSegP*)segt->At(0);
			theta0 = firstseg->Theta();
            if(lastseg->Plate()>=LASTLAYER_S1+5&&lastseg->Plate()>=plate){ //aggiungo i cosmici flaggati in S2
                for(int iseg=0;iseg<nseg; iseg++){
                    seg=(EdbSegP*)segt->At(iseg);
                    if(seg->Plate()>LASTLAYER_S1 && seg->Plate()<=LASTLAYER_S2){
                        if(seg->Flag()==-1){//seleziono i cosmici
                            for(int iseg2=0;iseg2<nseg; iseg2++){
                                s=(EdbSegP*)segt->At(iseg2);
                                if(s->Plate()==plate){
                                    //cout << s->Plate() << "\t" << s->Flag() << "\t" << s->TX() << endl;

                                    *newcp = *s; //copy the segment
                                    newcp->SetVid(900, 900);
                                    float expected_z = thickS2*(plate - LASTLAYER_S1 - 1) + (LASTLAYER_S1)*thickS1;
                                    float dx = (newcp->Z() - expected_z)*newcp->TX();
                                    float dy = (newcp->Z() - expected_z)*newcp->TY();
									newcp->SetFlag(450);

                                    //cout << " For plate " << plate << " expected_z " << expected_z; 
                                    //cout << " old X " << newcp->X() << " dx " << dx << " old Z " << newcp->Z() << endl;

                                    //newcp->SetX(newcp->X() - dx); newcp->SetY(newcp->Y()-dy);

                                    newtree->Fill();
                                    COUNT_COSMIC++;
                                    found_cosmic_signal = true;
                                    break;
                                }
                            }
                        }
                    }
                    if (found_cosmic_signal) {
                        break;
                    }
                }
            }

			
            if(theta0>0.2 && avg_volume_S1<3000*(1 + TMath::Exp(1.4*theta0*theta0) ) && lastseg->Plate()>=plate && !found_cosmic_signal && avg_volume_S1>0){ // aggiungo i cosmici non flaggati in S2 //modV in_vtx==0 &&

							//cout << " here because avg volume is " << avg_volume_S1 << endl;
                            for(int iseg2=0;iseg2<nseg; iseg2++){
                                s=(EdbSegP*)segt->At(iseg2);
                                if(s->Plate()==plate){
                                    //cout << s->Plate() << "\t" << s->Flag() << "\t" << s->TX() << endl;

                                    *newcp = *s; //copy the segment
									newcp->SetFlag(450);
                                    newcp->SetVid(800, 800);
                                    float expected_z = (plate-1)*thickS1;
                                    if (plate>LASTLAYER_S1) expected_z = thickS2*(plate - LASTLAYER_S1 - 1) + (LASTLAYER_S1)*thickS1;
                                    
                                    float dx = (newcp->Z() - expected_z)*newcp->TX();
                                    float dy = (newcp->Z() - expected_z)*newcp->TY();

                                    if (plate>1 && eVERBOSE==2) {
                                        cout << " S1 - For plate " << plate << " expected_z " << expected_z; 
                                        cout << " old X " << newcp->X() << " dx " << dx << " old Z " << newcp->Z() << endl;
                                    }

                                    //newcp->SetX(newcp->X() - dx); newcp->SetY(newcp->Y()-dy);

                                    newtree->Fill();
                                    COUNT_COSMIC++;
                                    found_cosmic_signal = true;
                                    break;
                                }
                            }
            } 
//            if (found_cosmic_signal) {
//                continue;
//            }
        }
    }
    
    if(ADD_RANDOMCR==1){
        
        cout << "Random cosmic rays taken from: " << path << endl;
        
        int Nbkg = tree->GetEntries(Form(" ( s.eY >= %f && s.eY <= %f ) && ( (s.eX >= %f && s.eX <= %f) || (s.eX >= %f && s.eX <= %f) )", Ymin, Ymax, Xmin, cutMIN, cutMAX, Xmax));
        cout << "Number of background entries: " << Nbkg << endl;
        cout << Form(" ( s.eY >= %f && s.eY <= %f ) && ( (s.eX >= %f && s.eX <= %f) || (s.eX >= %f && s.eX <= %f) )", Ymin, Ymax, Xmin, cutMIN, cutMAX, Xmax) << endl;
        
        Float_t area = (Ymax - Ymin)* ( (cutMIN-Xmin) + (Xmax-cutMAX));
        cout << "Density of background entries: " << Nbkg/area << " um^-2 " << endl;

        int nBinsTX = 200;
        int nBinsTY = 200;

        // Find TX and TY ranges
        Float_t TXmin = -1; Float_t TXmax = 1;
        Float_t TYmin = -1; Float_t TYmax = 1;

        // 2D histogram of angular distribution
        TH2F* hTXTY = new TH2F("hTXTY","Background TX-TY distribution", nBinsTX, TXmin, TXmax, nBinsTY, TYmin, TYmax);

        // Fill histogram only for background tracks
        tree->Draw(Form("s.TX():s.TY() >> hTXTY"), Form("( s.eY >= %f && s.eY <= %f ) && ( (s.eX >= %f && s.eX <= %f) || (s.eX >= %f && s.eX <= %f) )",
                Ymin, Ymax, Xmin, cutMIN, cutMAX, Xmax), "goff");
        
        int N_sample = Nbkg/area * (Ymax - Ymin) * (Xmax - Xmin);
        cout << " sampling " << N_sample << " random background base-tracks " << endl;

        EdbScanCond* scancond = new EdbScanCond();
        scancond->SetDefault();
        scancond->SetDegrad(1); 
        scancond->SetSigma0(5, 5, 0.005, 0.005); //to check

        TMatrixD temp_cov_matrix = TMatrixD(5,5);

        *newcp = *s;

        Double_t newX = gRandom->Uniform(Xmin, Xmax);
        cout << "Entries in hTXTY: " << hTXTY->GetEntries() << endl;

        Double_t newTX, newTY;
        hTXTY->GetRandom2(newTX,newTY);

        for(int i=0; i<N_sample; i++){
            float newX = gRandom->Uniform(Xmin, Xmax);
            float newY = gRandom->Uniform(Ymin, Ymax);

            // Sample TX, TY according to the original background distribution
            Double_t newTX, newTY;
            hTXTY->GetRandom2(newTX,newTY);

            float tx = float(newTX); float ty = float(newTY);

            // Create a new segment
            newcp->SetX(newX);
            newcp->SetY(newY);
            newcp->SetTX(tx);
            newcp->SetTY(ty);

            newcp->SetVid(700,700);
            scancond->FillErrorsCov(tx, ty, temp_cov_matrix);
            newcp->SetCOV(temp_cov_matrix);
            newtree->Fill();
        }

    }

    //newtree->Print();
    newtree->AutoSave();
    newfile->Close();
    
    cout << COUNT_COSMIC << " cosmic rays basektracks saved = " << (float)COUNT_COSMIC/nentries_trk*100 << "\%" << endl;
    cout << COUNT << " basetracks saved = " << (float)COUNT/nentries*100 << "\%" << endl;
    
    
    cout << "countmax " << countmax << "\tcountmin " << countmin <<endl;

    file->Close();
    file_trk->Close();
    
    delete file;
    delete file_trk;
    delete segt;
    
    
    
}

void smearing_theta (Float_t &TX, Float_t &TY, const float angres){
    float deltaTX = grandom->Gaus(0,angres); //angular resolution, adding a gaussian offset to TX and TY
    float deltaTY = grandom->Gaus(0,angres);
    //cout<<TX<<endl;
    TX = TX + deltaTX;
    TY = TY + deltaTY;
}

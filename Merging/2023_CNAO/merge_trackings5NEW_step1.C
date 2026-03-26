/*
 scp /Users/Giuliana/Desktop/Uni/FOOT/GSI_data/merge_trackings2NEW.C scanner@nusrv9.na.infn.it:/home/scanner/foot/RECO_MC_TEST_giul/
 
 su nusrv9 lanciare dalla cartella
 
 /home/scanner/foot/RECO_MC_TEST_giul/
 
 root -l launch_mergetrackings3.C
 
 NB: bisogna prima aver fatto makescanset dell'intero brick con
 
 tracks->Scan("nseg:s[0].ID():s[0].Plate():s[0].eX:s[0].eY:s[0].eTX:s.[0]eTY:s[nseg-1].ID():s[nseg-1].Plate():s.[nseg-1].eX:s.[nseg-1].eY:s[nseg-1].eTX:s[nseg-1].eTY","s[0].Plate()<10&&s[nseg-1].Plate()>25")
 
 */

#define EVERBOSE 0
#define NSTACK 7
#define DIRECTION 1 // direction 1 = da 1 a 66. direction 9 = da 66 a 1
#define IDBRICK 777	
#define STEP 1 // 1= only correct TX,TY, 2= only correct X,Y

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include "TH1.h"
#include "TCanvas.h"
#include <stdio.h>
#include <TROOT.h>
#include "TRandom3.h"
#include "TVector.h"
#include <vector>
#include <algorithm>
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"
/*#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbCouplesTree.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbDataSet.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbScanSet.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbSegP.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbPVRec.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbPattern.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbTrackFitter.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbLayer.h"*/
#include "TEntryList.h"

using namespace std;

EdbDataProc sproc;
EdbPVRec *ali = new EdbPVRec();
TString folderpath;


void add_stack(int S, TString folderpath);
void TransformAngles(EdbSegP &segment, EdbAffine2D aff);

int S0=1;
int SL=2;

int LASTLAYER[NSTACK+1]={3,33,71,76,83,90,110,120}; //CN7
//int LASTLAYER[NSTACK+1]={0,30,66,76,83,90,110,120}; GSI 1 e GSI2

void merge_trackings5NEW_step_1(){
    
    cout << endl << endl << endl << "Brick: " << IDBRICK << endl;
    folderpath = "./";
    
    //TString newpath=Form("%sb%06d/b%06d.0.%d.%d.trk_trasl.root", folderpath.Data(), IDBRICK, IDBRICK, S0, SL);
    
    cout << " STEP " << STEP << ": will only correct TXY or XY " << endl;
    
    TString newpath=Form("%sb%06d/b%06d.0.%d.%d.trk_trasl_step_%d.root", folderpath.Data(), IDBRICK, IDBRICK, S0, SL, STEP);

    //cout << "Attenzione! Il file set deve comprendere solo gli stack che si intende unire" << endl; // per direction 9
    
    
    for(int iS=S0; iS<=SL; iS++){
        add_stack(iS, folderpath);
    }
    
    
    sproc.MakeTracksTree(ali, newpath);
    cout << "npat: " << ali->Npatterns() << endl;
    
}

void add_stack(int S, TString folderpath)
{
    
    cout << endl << endl;
    
    TString setfilename=Form("%sb%06d/b%06d.0.0.0.set.root", folderpath.Data(), IDBRICK, IDBRICK); //set di tutto il brick
    TFile *setfile = TFile::Open(setfilename);
    EdbScanSet *scanset = (EdbScanSet*) setfile->Get("set");

    TString path=Form("%sb%06d/b%06d.%d.0.0.trk.root", folderpath.Data(), IDBRICK, IDBRICK, S);
    if (S==2) path=Form("%sb%06d/b%06d.%d.0.0.trk_updateV.root", folderpath.Data(), IDBRICK, IDBRICK, S);
    TFile *file = new TFile(path, "READ");
    
    EdbBrickP brick=scanset->Brick();
    int nplate_set=scanset->Brick().Npl();
    
    EdbPlateP *p=0, *thisp=0;
    EdbAffine2D *aff, *aff_temp, *aff_alXY, *affXY, *affTXTY, *aff_alTXTY; //affine transformation from local to global;
    EdbAffine2D *invaff, *invaff2, *invaffTXTY; //affine transformation from global to local
    EdbTrackFitter *trackfitter=new EdbTrackFitter();
    
    EdbPattern *pat=0;
    
    TTree *tracks = (TTree*)file->Get("tracks");
    
    cout << "Using scanset: " << setfilename << " with " << nplate_set << " plates" << endl;
    cout << "Adding file: " << path << " with ntrks : " << tracks->GetEntries() << endl;
    
    // tracks
    
    int lastsegID = 0,nseg = 0, trid=0, npl=0, n0=0;
    float w=0;
    float Zseg=0;
    int newpid=0;
    double mean_volume=0;
    TClonesArray *segments = new TClonesArray("EdbSegP",1000);
    TClonesArray *sf = new TClonesArray("EdbSegP",1000);
    EdbSegP *tr=NULL;
    
    tracks->SetBranchAddress("trid",&trid);
    tracks->SetBranchAddress("t.",&tr);
    tracks->SetBranchAddress("s",&segments); //array of measured segments
    tracks->SetBranchAddress("sf",&sf); //array of fitted segments
    tracks->SetBranchAddress("nseg",&nseg); //total number of segments
    tracks->SetBranchAddress("npl",&npl); //total number of plates passed by track
    tracks->SetBranchAddress("n0",&n0); //number of "holes" - plates where the segments were not found. In principle n0 should always be npl-nseg
    tracks->SetBranchAddress("w",&w); //the sum of all s.eW
    
    
    EdbSegP *seg = new EdbSegP();
    EdbTrackP *newtr = new EdbTrackP();
    
    int TRACKS=0;
    int pid30=0;

	// offset corrections *NEED merge_offsets.py in brick dir

	float X_offset=0, Y_offset=0; 
	float TX_offset=0, TY_offset=0;

	float b=0, dtx=0, dty=0, dxp=0, dyp=0;
    float fit_prob=0, mean=0;
	float* row_content;

	TString offset_file_name, offset_file_name_step_1;

	offset_file_name = Form("%sb%06d/%d_S%i_S%i_offsets_all.root", folderpath.Data(), IDBRICK, IDBRICK, S0, SL); 
    offset_file_name_step_1 = Form("%sb%06d/%d_S%i_S%i_offsets_all_step_1.root", folderpath.Data(), IDBRICK, IDBRICK, S0, SL); 
	cout << " Using Offset file: " << offset_file_name << endl;
    TFile* off_file = TFile::Open(offset_file_name);
	TNtuple* tup = (TNtuple*) off_file->Get("merge_offsets"); // containing offset info

    tup->SetBranchAddress("Xoff", &X_offset);
    tup->SetBranchAddress("Yoff", &Y_offset);

    tup->SetBranchAddress("TXoff", &TX_offset);
    tup->SetBranchAddress("TYoff", &TY_offset);
    const int N = SL-S0;
    std::vector<float> Z_refs(N, 0.0f);
    std::vector<float> X_offsets(N, 0.0f);
    std::vector<float> Y_offsets(N, 0.0f);
    std::vector<float> TX_offsets(N, 0.0f);
    std::vector<float> TY_offsets(N, 0.0f);

    //std::vector<float> Z_refs, X_offsets, Y_offsets, TX_offsets, TY_offsets;

    for(int i=0; i<N; i++) {
        tup->GetEntry(i);
        X_offsets[i] = X_offset;
        Y_offsets[i] = Y_offset;
        TX_offsets[i] = TX_offset;
        TY_offsets[i] = TY_offset;
        p = scanset->Brick().GetPlate(LASTLAYER[i+1]);
        float Zref = p->Z(); //Z_refs.push_back(Zref);
        Z_refs[i] = Zref;
    }

    if (STEP == 2) { //need to apply position correction from step 1 and angular correction from step 0 
        TFile* off_file_step1 = TFile::Open(offset_file_name_step_1);
        TNtuple* tup_step_1 = (TNtuple*) off_file_step1->Get("merge_offsets"); // containing offset info

        float X_offset_step_1=0, Y_offset_step_1=0; 
        tup_step_1->SetBranchAddress("Xoff", &X_offset_step_1);
        tup_step_1->SetBranchAddress("Yoff", &Y_offset_step_1);

        for(int i=0; i<N; i++) {
            tup_step_1->GetEntry(i);
            X_offsets[i] = X_offset_step_1;
            Y_offsets[i] = Y_offset_step_1;
            cout << " here X off set " << X_offset_step_1 << endl;
        }
    }

	cout << " X offsets [0] " << X_offsets[0] << endl;
    cout << " TX offsets [0] " << TX_offsets[0] << endl;

    if(S>1){

        //X_offsets.push_back(X_offset); Y_offsets.push_back(Y_offset);
        //TX_offsets.push_back(TX_offset); TY_offsets.push_back(TY_offset); 

        p = scanset->Brick().GetPlate(LASTLAYER[S-1]-1); //prendo direttamente makescansetunico. voglio piatto 31
        for(int ipid=0; ipid<scanset->Brick().Npl(); ipid++){
            //cout << ipid << "\t" << scanset->GetID(ipid)->GetPlate() << "\t" << LASTLAYER[S-1]+1 << endl;
            if(scanset->GetID(ipid)->GetPlate()==LASTLAYER[S-1]+1) {
                pid30=ipid;
                break;
            }
        }
        aff = p->GetAffineXY();
        affTXTY = p->GetAffineTXTY();
        cout << "Interruzione al piatto " << p->ID() << "\t" << pid30 << "\t" << p->Z() << endl;
        
        cout << "Aff X Y: ";
        aff->Print();
        cout << endl << "Aff TX TY: ";
        affTXTY->Print();

		off_file->Close();
    }
   
    int count=0;
    float dz = 0;
    
    for (int itrk = 0; itrk < tracks->GetEntries(); itrk++){ //tracks->GetEntries()
        
        if(itrk%10000==0) cout << "S " << S << "\t" << (float)itrk/tracks->GetEntries()*100 << "\%" << endl;
        tracks->GetEntry(itrk);
        newtr = new EdbTrackP(nseg);
        newtr->SetNpl(npl);
        newtr->SetN0(n0);
        newtr->SetID(trid);
        newtr->SetW(w);
        newtr->SetFlag(tr->Flag());
        
        if(IDBRICK<50){
            newtr->SetMC(tr->MCEvt(), tr->MCTrack());//tr2->W()-70);
        }
        newtr->SetAid(tr->Aid(0), tr->Aid(1));
        newtr->SetVid(tr->Vid(0), tr->Vid(1));
        
        for(int iseg=0; iseg<nseg; iseg++){
            
            seg=(EdbSegP *)segments->At(iseg);
            seg->SetTrack(trid);
            seg->SetPlate(seg->Plate());
            seg->SetP(-999);
            //seg->SetDZ(300);
            // if(S!=2)
            mean_volume+=seg->Volume();
            Zseg = seg->Z();

            if(DIRECTION==1){
                
                seg->SetPID(seg->Plate()-1);

                if(S>1) {
                //thisp = scanset->Brick().GetPlate(seg->Plate()-1); //makescansetunico

                seg->Transform(aff);
                TransformAngles(*seg, *affTXTY);
                Zseg = p->Z()+seg->Z();
                float total_Xshift=0, total_Yshift=0, total_TXshift=0, total_TYshift=0;
                for (int ns=0; ns<S-S0; ns++) {
                    
                    if (STEP == 1) {
                        total_Xshift += (Zseg - Z_refs[ns])*TX_offsets[ns];
                        total_Yshift += (Zseg - Z_refs[ns])*TY_offsets[ns];
                    }
                    else if (STEP == 2) {
                        total_Xshift += X_offsets[ns] + (Zseg - Z_refs[ns])*TX_offsets[ns];
                        total_Yshift += Y_offsets[ns] + (Zseg - Z_refs[ns])*TY_offsets[ns]; 
                    }
                    total_TXshift += TX_offsets[ns];
                    total_TYshift += TY_offsets[ns];

                }
                //newpid = pid30+seg->PID();
                //seg->SetPID(newpid);
				seg->SetX(seg->X() + total_Xshift);
				seg->SetY(seg->Y() + total_Yshift);
				seg->SetTX(seg->TX() + total_TXshift);
				seg->SetTY(seg->TY() + total_TYshift);
                if(count<10) {
                    cout << seg->ID() << "\t" << seg->PID() << "\t" << seg->Plate() << "\t" << seg->Z() << "\t"  << p->Z() <<  "\t-->" << Zseg << endl;
                    count ++;
                }
                seg->SetZ(Zseg);


                //cout << "NEW: " << seg->PID() << "\t" << seg->Plate() << "\t" << seg->Z() << endl;
            }
        }
            
            if(IDBRICK>=50) seg->SetFlag(seg->Flag());
            if(IDBRICK<50){
                seg->SetFlag(10);//tr2->W()-70);
                seg->SetMC(seg->MCEvt(), seg->MCTrack());//tr2->W()-70);
            }
            
            
            pat = ali->GetPattern( seg->PID() );
            if(!pat) {
                pat = new EdbPattern( 0., 0., seg->Z() );
                pat->SetID(seg->PID());
                pat->SetScanID(IDBRICK);
                //pat->SetScanID(seg->ScanID());
                ali->AddPatternAt(new EdbPattern(*pat),seg->PID());
                cout << "EdbDataProc::ReadTracksTree: WARNING: no pattern with pid " << seg->PID() << " ( plate " << seg->Plate() << ") Z: " << seg->Z() << " -> creating new one!" << endl;
                
            }
            //cout << seg->PID() << "\t" << seg->Plate() << "\t" << p->Z() << "\t" << seg->Z() << endl;
            newtr->AddSegment(new EdbSegP(*seg));
            
            if(iseg==0) {
                newtr->SetPID(seg->PID());
                newtr->SetPlate(seg->Plate());
                //cout << seg->PID() << "\t" << newtr->PID() << "\t\t" << seg->Plate() << "\t" << newtr->Plate() << endl;
            }
            
            int temp_pid=seg->PID();
            int tempid=seg->ID();

            
            ///////segmenti fittati /////
            
            seg=(EdbSegP *)sf->At(iseg);
            seg->SetTrack(trid);
            seg->SetP(-999);
            seg->SetID(tempid);
            
            if(DIRECTION==1 && S>1) {
                seg->Transform(aff);
                TransformAngles(*seg, *affTXTY);
                seg->SetZ(Zseg);
                float total_Xshift=0, total_Yshift=0, total_TXshift=0, total_TYshift=0;
                for (int ns=0; ns<S-S0; ns++) {
                    if (STEP == 1) {
                        total_Xshift += (Zseg - Z_refs[ns])*TX_offsets[ns];
                        total_Yshift += (Zseg - Z_refs[ns])*TY_offsets[ns];
                    }
                    else if (STEP == 2) {
                        total_Xshift += X_offsets[ns] + (Zseg - Z_refs[ns])*TX_offsets[ns];
                        total_Yshift += Y_offsets[ns] + (Zseg - Z_refs[ns])*TY_offsets[ns]; 
                    }
                    total_TXshift += TX_offsets[ns];
                    total_TYshift += TY_offsets[ns];
                }

				seg->SetX(seg->X() + total_Xshift);
				seg->SetY(seg->Y() + total_Yshift);
				seg->SetTX(seg->TX() + total_TXshift);
				seg->SetTY(seg->TY() + total_TYshift);
            }
            
            seg->SetPlate(0);
            seg->SetPID(temp_pid);
            
            if(IDBRICK>=50) seg->SetFlag(10);
            if(IDBRICK<50){
                seg->SetFlag(10);//tr2->W()-70);
                seg->SetMC(seg->MCEvt(), seg->MCTrack());//tr2->W()-70);
            }
            newtr->AddSegmentF(new EdbSegP(*seg));
            
            if(iseg==0) {
                newtr->SetCOV(seg->COV());
                newtr->SetX(seg->X());
                newtr->SetY(seg->Y());
                newtr->SetZ(seg->Z());
                newtr->SetTX(seg->TX());
                newtr->SetTY(seg->TY());
                //cout << "F " << seg->PID() << "\t" << newtr->PID() << "\t\t" << seg->Plate() << "\t" << newtr->Plate() << endl;
            }

        }
        
        
        // if(S!=2) {
        mean_volume=(double)mean_volume/nseg;
        newtr->SetVolume(mean_volume);
        // }
        //        if(newtr->Theta()>=0.2 && mean_volume<= 13687*newtr->Theta()+13781 && nseg==2) {
        //            newtr->SetFlag(0); //setto flag cosmici
        //        }
        //        else newtr->SetFlag(10);
        
        
        
        //newtr->FitTrack();
        ali->AddTrack(new EdbTrackP(*newtr));
        newtr->Clear();
        TRACKS++;
        // }
    }
    
    cout << TRACKS << " tracks added from Stack " << S << endl;
    
}


void TransformAngles(EdbSegP& segment, EdbAffine2D aff) {

	float new_tx = aff.A11()*segment.TX() + aff.A12()*segment.TY() + aff.B1();
	float new_ty = aff.A21()*segment.TX() + aff.A22()*segment.TY() + aff.B2();
	segment.SetTX(new_tx);
	segment.SetTY(new_ty);


}
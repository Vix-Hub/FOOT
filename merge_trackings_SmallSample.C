/*
 scp /Users/Giuliana/Desktop/Uni/FOOT/GSI_data/merge_trackings2NEW.C scanner@nusrv9.na.infn.it:/home/scanner/foot/RECO_MC_TEST_giul/
 
 su nusrv9 lanciare dalla cartella
 
 /home/scanner/foot/RECO_MC_TEST_giul/
 
 root -l launch_mergetrackings2.C
 
 NB: bisogna prima aver fatto makescanset dell'intero brick con
 
 tracks->Scan("nseg:s[0].ID():s[0].Plate():s[0].eX:s[0].eY:s[0].eTX:s.[0]eTY:s[nseg-1].ID():s[nseg-1].Plate():s.[nseg-1].eX:s.[nseg-1].eY:s[nseg-1].eTX:s[nseg-1].eTY","s[0].Plate()<10&&s[nseg-1].Plate()>25")
 
 */

#define EVERBOSE 0
#define NSTACK 7
#define DIRECTION 1 // direction 1 = da 1 a 66. direction 9 = da 66 a 1
#define IDBRICK 3
#define PID_NOHOLE 0

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
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbCouplesTree.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbDataSet.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbScanSet.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbSegP.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbPVRec.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbPattern.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbTrackFitter.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbLayer.h"
#include "TEntryList.h"

using namespace std;

EdbDataProc sproc;
EdbPVRec *ali = new EdbPVRec();

TString folderpath;


void add_stack(int S);

int S0=1;
int SL=2;

//int LASTLAYER[NSTACK+1]={0,30,66,76,83,90,110,120}; //GSI 1 e GSI2
int LASTLAYER[NSTACK+1]={0,30,66,76,83,90,110,140}; //GSI 3 e GSI4

void merge_trackings(){
    
    cout << endl << endl << endl << "Brick: " << IDBRICK << endl;
    
    if(IDBRICK<10) folderpath = "/home/scanner/foot/RECO_MC_TEST_giul/";
    else if(IDBRICK==11||IDBRICK==111) folderpath = "/home/scanner/foot/2019_GSI/GSI1/";
    else if(IDBRICK==22||IDBRICK==222) folderpath = "/home/scanner/foot/2019_GSI/GSI2/";
    else if(IDBRICK==33||IDBRICK==333) folderpath = "/home/scanner/foot/2019_GSI/GSI3/";
    else if(IDBRICK==44||IDBRICK==444) folderpath = "/home/scanner/foot/2019_GSI/GSI4/";
    
    TString newpath=Form("%sb%06d_SmallSample/b%06d/b%06d.0.%d.%d.trk.root", folderpath.Data(), IDBRICK, IDBRICK, IDBRICK, S0, SL);
    

    //cout << "Attenzione! Il file set deve comprendere solo gli stack che si intende unire" << endl; // per direction 9
    
    
    for(int iS=S0; iS<=SL; iS++){
        add_stack(iS);
    }
    
    
    sproc.MakeTracksTree(ali, newpath);
    cout << "npat: " << ali->Npatterns() << endl;
    
}

void add_stack(int S)
{
    
    cout << endl << endl;
    
    TString setfilename=Form("%sb%06d/b%06d.0.0.0.set.root", folderpath.Data(), IDBRICK, IDBRICK); //set di tutto il brick
    TFile *setfile = TFile::Open(setfilename);
    EdbScanSet *scanset = (EdbScanSet*) setfile->Get("set");

    TString path=Form("%sb%06d/b%06d.%d.0.0.trk.root", folderpath.Data(), IDBRICK, IDBRICK, S);
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
    
    if(S>1){
    	cout << "LASTLAYER[S-1]" << endl;
        p = scanset->Brick().GetPlate(LASTLAYER[S-1]); //prendo direttamente makescansetunico. voglio piatto 31
        for(int ipid=0; ipid<scanset->Brick().Npl(); ipid++){
            cout << ipid << "\t" << scanset->GetID(ipid)->GetPlate() << "\t" << LASTLAYER[S-1]+1 << endl;
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
    }
   
    int count=0;
    
    
    for (int itrk = 0; itrk < tracks->GetEntries(); itrk++){ //tracks->GetEntries()
        
        if(itrk%10000==0) cout << "S " << S << "\t" << (float)itrk/tracks->GetEntries()*100 << "\%" << endl;
        tracks->GetEntry(itrk);
        newtr = new EdbTrackP(nseg);
        newtr->SetNpl(npl);
        newtr->SetN0(n0);
        newtr->SetID(trid);
        newtr->SetW(w);
        newtr->SetFlag(tr->Flag());
        
        if(IDBRICK<10){
            newtr->SetMC(tr->MCEvt(), tr->MCTrack());//tr2->W()-70);
        }
        newtr->SetAid(tr->Aid(0), tr->Aid(1));
        newtr->SetVid(tr->Vid(0), tr->Vid(1));
        
        for(int iseg=0; iseg<nseg; iseg++){
            
            seg=(EdbSegP *)segments->At(iseg);
            seg->SetTrack(trid);
            seg->SetPlate(seg->Plate());
            seg->SetP(-999);
			seg->SetDZ(300);
            // if(S!=2)
            mean_volume+=seg->Volume();
            Zseg = seg->Z();
            if(DIRECTION==1){
                
                if(PID_NOHOLE==0) seg->SetPID(seg->Plate()-1);

                if(S>1) {
                //thisp = scanset->Brick().GetPlate(seg->Plate()-1); //makescansetunico
                seg->Transform(aff);
                seg->Transform(affTXTY);
                Zseg = p->Z()+seg->Z();
                    if(PID_NOHOLE==1){
                        newpid = pid30+seg->PID();
                        seg->SetPID(newpid);
                    }
                if(count<10) {
                    cout << seg->ID() << "\t" << seg->PID() << "\t" << seg->Plate() << "\t" << seg->Z() << "\t"  << p->Z() <<  "\t-->" << Zseg << endl;
                    count ++;
                }
                seg->SetZ(Zseg);


                //cout << "NEW: " << seg->PID() << "\t" << seg->Plate() << "\t" << seg->Z() << endl;
            }
        }
            
            if(IDBRICK>=10) seg->SetFlag(seg->Flag());
            if(IDBRICK<10){
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
                seg->Transform(affTXTY);
                seg->SetZ(Zseg);
            }
            
            seg->SetPlate(0);
            seg->SetPID(temp_pid);
            
            if(IDBRICK>=10) seg->SetFlag(10);
            if(IDBRICK<10){
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



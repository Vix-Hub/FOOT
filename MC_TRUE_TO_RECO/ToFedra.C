//Questo script trasforma la simulazione MC in file di fedra, come se fossero stati acquisiti al mic
//tool to load hit from emulsion to FEDRA (13 aprile 2018) //updated to work with the new TTreeReader structure (6 March)
//to use it, go in a directory and create the folders b000012/p001 to b000012/p029
//then launch it from the directory mother of b000012

/*USAGE:
scp /Users/Giuliana/Desktop/Uni/FOOT/GSI_AnalisiMC/ToFedra.C scanner@nusrv9.na.infn.it:/home/scanner/foot/RECO_MC_TEST_giul
pw: neutrino.01
//once do:
g++ -shared -o EMU_Track_c.so EMU_Track_c.o GSI_AnaMCDict.o

/*run the script in the folder /home/scanner/foot/RECO_MC_TEST_giul
 
 root -l
 .except
 gROOT->ProcessLine(".L EMU_Track_c.so");
 .L ToFedra.C+
 ToFedra()
 
 oppure
 
 root -l launch_ToFedra.C
 
 */

#define NR 4 // numero di piatti con diverso refresh in S2

#include <stdio.h>
#include <TROOT.h>
#include "TRandom3.h"
#include "TVector.h"
#include <vector>

#include "GSI_AnaMC.h"
#include "/mnt/sdb/opera/ana/valeri/sw_oracle_dev/fedra/include/EdbCouplesTree.h"

//gROOT->ProcessLine(".L EMU_Track_c.so");
//gROOT->ProcessLine(".L GSI_AnaMCDict.so");


using namespace TMath;
TRandom *grandom = new TRandom3(0); //creating every time a TRandom3 is a bad idea
void smearing_theta (Float_t &TX, Float_t &TY, const float angres);
void smearing_XY (Float_t &X, Float_t &Y, const float posres);
bool efficiency(const float emuefficiency);

int ToFedra(){
    
    const float emuefficiency = 0.8;//usual 0.90;
    const float angres = 0.006; // usual 5 milliradians
    const float posres = 0.02; // usual 0.01 micron

    const float ngrains = 70; //the same number for all the couples, so they have the same weigth.
    const int minplate = 0;
    
    sprintf(fascio,"Oxy400n"); //Oxy200n //He4_700 //C200 // Oxy400n //C700n
    sprintf(materiale_S1, "Carbonio"); // Carbonio //Polietilene
    sprintf(materiale_S3, "LxWWPbS3"); // C2H4WWPbS3 //LxWWPbS3 //40371014144010
    geo = 32; //32 per albox esposizione 2019 // 7 per esposizione 2020
    
    
    //**********************OPENING INPUT FILE***************************
    char inputfilename[100];
    filename(inputfilename, fascio, materiale_S1, geo);
    
    cout << "Last Layers: ";
    for(int l=1; l<=N_STACKS; l++){
        cout << "\t" << LASTLAYER[l];
    }
    cout << endl;
    
    const int nplates = LASTLAYER[N_STACKS];
    cout << "Using file: " << inputfilename << endl; // << "Output file: " << outputfilename << endl << endl;
    TFile *f = new TFile(inputfilename);
    if(!f) { cout << "Input file not exists!"<< endl; return -99;}
    
    cout << "Creating cp.root files for brick " << idbrick << endl;
    
    EMU_Event *evt = NULL;
    EMU_BaseTrack bt;
    EMU_VolumeTrack vt;
    
    //getting tree and arrays
    TTree *EventTree = (TTree*)f->Get("EventTree");
    EventTree->SetBranchAddress("event", &evt);
    
    const Int_t nevents = EventTree->GetEntries();
    
    //    TTreeReader reader("cbmsim",inputfile);
    //    TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
    //    TTreeReaderArray<BoxPoint> emulsionhits(reader,"BoxPoint");
    
    //TTree* cbmsim = (TTree*)inputfile->Get("cbmsim");
    Float_t tx = 0, ty=0, xem= 0, yem = 0;
    int ihit = 0, ievent = 0, nseg=0;
    int nfilmhit = 0;
    int trackID = 0, motherID = 0, pdgcode = 0, eventID=0, flukacharge=0;
    int firstplate=0, lastplate=0;
    float momentum=0, Ekin=0;
    Int_t Flag = 0;

    //***********************CREATING FEDRA TREES**************************
    //gInterpreter->AddIncludePath("/afs/cern.ch/work/a/aiuliano/public/fedra/include");
    vector<EdbCouplesTree *> vec_ect;
    EdbCouplesTree *ect[nplates+1]; //= new EdbCouplesTree();
    
    cout << "nevents: " << nevents << endl;
    cout << "nplates: " << nplates+1 << endl;
    // EdbCouplesTree *ect[nplates] = NULL;
    //I still need to prepare the other ones, even if empty
    for (int i= 0; i < minplate; i++){
        ect[i] = new EdbCouplesTree();
        vec_ect.push_back(ect[i]);
    }
    for (int i = minplate; i <= nplates; i++){
        cout << "Plate: " << i+1 << endl;
        ect[i] = new EdbCouplesTree();
        ect[i]->InitCouplesTree("couples",Form("b%06d/p%03d/%d.%d.0.0.cp.root",idbrick, i+1,idbrick, i+1),"RECREATE");
        vec_ect.push_back(ect[i]);
        //delete ect;
    }
    bool savehit;
    //************************STARTING LOOP ON SIMULATION******************
    //    while (reader.Next()){
    //        for (const BoxPoint& emupoint:emulsionhits){
    
    
    for(int iev = 0; iev<14375; iev++){ //nevents //19375
        EventTree->GetEntry(iev);
        //vt = evt->GetVolumeTrack(0);
        //if(vt.GetX()>-0.1&&vt.GetX()<0.1&&vt.GetY()>-0.1&&vt.GetY()<0.1){ //producing a smaller file, commented by default
//        if(evt->GetVolumeTracks().size()==1) savehit=true;
       // if(evt->GetID()==9703) savehit=true;
       // else savehit=false;
        savehit=true;
        //if (savehit) {cout << iev << " " << evt->GetID() << " " << evt->GetVolumeTracks().size() << " " << endl;}
        
        for(unsigned int iv=0; iv<evt->GetVolumeTracks().size(); iv++){ // evt->GetVolumeTracks().size()
            // for(vector<EMU_VolumeTrack>::iterator iv = evt->GetVolumeTracks().begin();
            //  iv!=evt->GetVolumeTracks().end(); iv++){

            vt = evt->GetVolumeTrack(iv);
            
            //  cout << iv->GetX() << endl;
            //cout << vt.GetIDPart() << endl;
            nseg=vt.GetBaseTracks().size();
            //            eventID=vt.GetEvent();
            //            trackID=vt.GetIDPart();
            //            motherID = vt.GetParentIDPart();
            //            pdgcode = vt.GetFlukaID();
            firstplate=vt.GetBaseTrack(0).GetLayer()+1;
            lastplate=vt.GetBaseTrack(vt.GetBaseTracks().size()-1).GetLayer()+1;
            eventID= evt->GetID();
            flukacharge=vt.GetCharge();
            //cout << iev << " " << eventID << " " << trackID << endl;
            //cout << "\t" << iv << "\t" << eventID << " " << trackID << " " << nseg << " " << firstplate << " " << lastplate << endl;
            
            
            for(unsigned int ib=0; ib<vt.GetBaseTracks().size(); ib++){
                
                bt = vt.GetBaseTrack(ib);
                
                xem = bt.GetX()* 1E+4 + 62500;
                yem = bt.GetY()* 1E+4 + 49500;
                tx = bt.GetPX()/bt.GetPZ();
                ty = bt.GetPY()/bt.GetPZ();
                //momentum = bt.GetP();
                momentum = bt.GetEkin();
                //trackID = bt.GetIDPart();
                motherID = bt.GetParentIDPart();
                pdgcode = bt.GetFlukaID();
                trackID = vt.GetID();
                Flag= 0;//vt.GetIDPart();
                Ekin=bt.GetEkin();
                
                nfilmhit = bt.GetLayer(); //getting number of the film
                //cout << eventID << " " << trackID << " " << lastplate << " " << nfilmhit << endl;
                //if (savehit){cout << "\t\t\t" << nfilmhit << endl;}
                //
                //                //*************EXCLUDE HITS FROM BEING SAVED*******************
                //cout << nfilmhit << "\t";
                //if ((nfilmhit > LASTLAYER[N_STACKS])||(!efficiency(emuefficiency))) savehit = false;
                ////                //double charge;
                ////                //if(charge == 0.) savehit = false; //we do not track neutral particles
                //////                //saving the hits for a plate in the corresponding couples (only one layer saved, the other has ID + 10000)
                //////                //Inserting real data effects in the simulation (now COMMENTED OUT)
                //if(!efficiency(emuefficiency)) savehit = false; //inserting some holes due to emulsion inefficiency
                
                savehit=efficiency(emuefficiency);
                smearing_theta(tx,ty,angres);
                smearing_XY(xem,yem,posres);

                if(bt.GetP()<0.1) savehit=false;
                
//                //Z=0 non salvo i segmenti in R1, R2 e R3
//                if(flukacharge==0&&bt.GetStack()==2){
//                    int temp=-1;
//                    for (int t=1; t<NR; t++){
//                        temp=temp*(-1);
//                        if ((nfilmhit+1+t+temp)%NR==0){
//                            savehit=false; //non salva segmenti in R1, R2 e R3
//                            //cout << iev << "\t" << iv << "\t" << pdgcode << "\tnot saving plate " << nfilmhit << endl;
//                        }
//                    }
//                }
                
                //Z=1 non salvo i segmenti in R2 e R3
                if(flukacharge==1&&bt.GetStack()==2){
                    int temp=-1;
                    for (int t=2; t<NR; t++){
                        temp=temp*(-1);
                        if ((nfilmhit+1+t+temp)%NR==0){
                            savehit=false; //non salva segmenti in R2 e R3
                            //cout << iev << "\t" << iv << "\t" << pdgcode << "\tnot saving plate " << nfilmhit+1 << endl;
                        }
                    }
                }
                
                ///
                //                 //**************SAVING HIT IN FEDRA BASE-TRACKS****************
                if (savehit && nfilmhit>=minplate){
                    //cout << " salvo " << endl;
                    vec_ect.at(nfilmhit)->eS->Set(ihit,xem,yem,tx,ty,1,Flag);
                    vec_ect.at(nfilmhit)->eS->SetMC(eventID, trackID); //objects used to store MC true information
                    vec_ect.at(nfilmhit)->eS->SetAid(pdgcode, motherID); //forcing areaID member to store mother MC track information and MC total number of plates
                    vec_ect.at(nfilmhit)->eS->SetVid(firstplate, lastplate); //forcing ViewID member to store first and last plate track information
                    vec_ect.at(nfilmhit)->eS->SetW(70+flukacharge); //need a high weight to do tracking
                    //vec_ect.at(nfilmhit)->eS->SetP(momentum); //need a high weight to do tracking
                    vec_ect.at(nfilmhit)->eS->SetDZem(Ekin); //set Ekin

                    vec_ect.at(nfilmhit)->Fill();
                    ihit++; //hit entry, increasing as the tree is filled
                }
                //else cout << " NON salvo " << endl;

            }//end of loop on emulsion points
        }
        ievent++;
       // } // end loop su selection if
    } //end of loop on tree
    
    cout << "ihit: " << ihit << endl << "ievent: " << ievent << endl;
    
    for (int iplate = minplate; iplate <= nplates; iplate++){
        vec_ect.at(iplate)->Close();
    }
    
    return 0;
}

void smearing_theta (Float_t &TX, Float_t &TY, const float angres){
    float deltaTX = grandom->Gaus(0,angres); //angular resolution, adding a gaussian offset to TX and TY
    float deltaTY = grandom->Gaus(0,angres);
    //cout<<TX<<endl;
    TX = TX + deltaTX;
    TY = TY + deltaTY;
}

void smearing_XY (Float_t &X, Float_t &Y, const float posres){
    
    float deltaX = grandom->Gaus(0,posres); //psition resolution, adding a gaussian offset to X and Y
    float deltaY = grandom->Gaus(0,posres);
    X = X + deltaX;
    Y = Y + deltaY;
}


bool efficiency(const float emuefficiency){ //for now, just a constant, to be replaced with an efficiency map (probably with the angle)
    
    float prob = grandom->Uniform(0,1);
    //cout << prob << "\t";
    if (prob < emuefficiency) return true; //efficiency larger than probability, we take the event
    else return false;
}



//Fluka ID:
//-7 = frammenti da evaporazione nucleone
//-6 = 4He
//-5 = 3He
//-4 = 3H = Trizio
//-3 = 2H = deuterio
//-2 = heavy ion
//
//1 = protone
//3 = elettrone
//4 = positrone
//7 = fotone
//8 = neutrone
//13 = pione+
//14 = pione-

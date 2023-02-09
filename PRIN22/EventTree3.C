//   In a ROOT session, you can do:
//      root> .L EventTree.C
//      root> EventTree t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

/*Usage
 
 root -l EventTree.C
 main()
*/

#define EventTree_cxx
#define EVERBOSE 0
#include "EventTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>

int main(){
    EventTree tr;
    tr.Loop();
    return 0;
}

void EventTree::Loop()
{
    if (fChain == 0) return;
    
    sprintf(fascio,"P200");
    geo=1; //
    
    const float L_THRESH=0.00002; // 200 nanometri in cm
    const float L_THRESH_LONG=0.001; // 10 micron in cm
    const float TH_THRESH=0.01; //rad
    const int N_CHARGE=8;
    
	// DEF NUMBERS
    
    int MULT=0, MULT_visible=0, mult_AgBr=0, RECOVTX=0;
    int Z_atom=0, Z_reco=0, ch_reco=0, ZTree=0;
    float length, energy, kink, sqrt1suL;
    int p=0, D=0, Tr=0, He3=0, He4=0, HeavIo=0, others=0;
    int nbeam=0, beam_interagenti=0, beam_interagenti_AgBr=0, beam_interagenti_visibili=0, beam_interagenti_AgBr_visibili=0, nbeam_passanti=0, prodotti=0, EXIT_LASTPLATE=0;
    int justone=0, justone_v=0, ipbeam=0;
    int longtrack=0, largekinktrack=0, longandkink=0, startinemu=0, startinplastic=0;
    int TOT_vtxemu=0, TOT_vtxpl=0;
	int emutrk=0, MULT_visible_long=0;

	const int nregFirstEm1 = 47;
    const int nregLastEm1 = 48;
    const int nregFirstEm2 = 49;
    const int nregLastEm2 = 248;
    const int nregFirstPl1 = 357; //259;
    const int nregFirstPl2 = 358; //260;
    const int nregLastPl2 = 458; //360;
    const int nregLastPB = 565; //no
    const float ZLAST = 8; //1.96;

	double ax=0, ay=0, az=0, ix=0, iy=0, iz=0, fax=0, fay=0, faz=0, pax=0, pay=0,paz=0,fx=0,fy=0,fz=0;
	int ev_ID=0, pa_ID=0, A=0, event_inemu=0, event_inplast=0, interac=0, first=0;
    
    Long64_t events = fChain->GetEntriesFast();

	cout << " --- Looping over " << events << " entries " << endl;
    
    TFile *fout = new TFile("outfile.root","RECREATE");
    fout->cd();

    TTree* treevar = new TTree("treevar","treevar");
    treevar->Branch("Z", &ZTree);
    treevar->Branch("Kink", &kink);
    treevar->Branch("sqrt1suL", &sqrt1suL);
	treevar->Branch("n", &MULT, "n/I");
	treevar->Branch("n_vis", &MULT_visible, "n_vis/I");
	treevar->Branch("L", &length, "L/F");
	treevar->Branch("iax", &ax, "iax/D");
	treevar->Branch("iay", &ay, "iay/D");
	treevar->Branch("iaz", &az, "iaz/D");

	treevar->Branch("fax", &fax, "fax/D");
	treevar->Branch("fay", &fay, "fay/D");
	treevar->Branch("faz", &faz, "faz/D");

	treevar->Branch("pax", &pax, "pax/D");
	treevar->Branch("pay", &pay, "pay/D");
	treevar->Branch("paz", &paz, "paz/D");


	treevar->Branch("A", &A, "A/I");
	treevar->Branch("ev_ID", &ev_ID, "ev_ID/I");
	treevar->Branch("pa_ID", &pa_ID, "pa_ID/I");

	treevar->Branch("ix", &ix, "ix/D");
	treevar->Branch("iy", &iy, "iy/D");
	treevar->Branch("iz", &iz, "iz/D");

	treevar->Branch("fx", &fx, "fx/D");
	treevar->Branch("fy", &fy, "fy/D");
	treevar->Branch("fz", &fz, "fz/D");

	treevar->Branch("atom", &Z_atom, "Z_atom/I");;
	
	treevar->Branch("Ekin", &energy, "Ekin/F");
	treevar->Branch("emu_trk", &emutrk, "emu_trk/I");
	treevar->Branch("interac", &interac, "interac/I");
	treevar->Branch("first", &first, "first/I");


	Double_t cosdir[3];
	Double_t cosdirparent[3];

	TH1I *h_mult_all = new TH1I("h_mult_all","Multiplicity distribution; multiplicity; entries", 20, 0, 20);

	TH1I *h_mult_visible = new TH1I("h_mult_visible","Multiplicity distribution; multiplicity; entries", 20, 0, 20);
	TH1I *h_mult_visible_long = new TH1I("h_mult_visible_long","Multiplicity distribution; multiplicity; entries", 20, 0, 20);

	int save_beam = 0;
	
    
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<events;jentry++) { //events
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        
        justone=0; justone_v=0;
        MULT=0; MULT_visible=0; mult_AgBr=0, MULT_visible_long=0;
        Z_atom=0, Z_reco=0, ch_reco=0;
        longtrack=0; largekinktrack=0; longandkink=0; startinemu=0; startinplastic=0;
        ipbeam=0; emutrk=0;
		event_inemu=0; event_inplast=0;
		ev_ID = jentry + 1;
		interac = 0;
        
        nb = fChain->GetEntry(jentry);   nbytes += nb;
		 
        if(EVERBOSE==1) if(TRn>0) cout << endl << ientry << "\t" << TRn << endl;

		if(TRn>1) interac = 1;
        
        for (Long64_t ip=0; ip<TRn; ip++) {
            
            if(TRfz[ip]<0) continue;
            
            //beam
            if(TRpaid[ip]==0 && TRiz[ip]<0){
                nbeam++;
                if(TRfz[ip]>ZLAST) {
					nbeam_passanti++;
       				ipbeam=ip;
				}

				//cout << "ev id " << ev_ID << endl;
				//cout << "ipbeam " << ipbeam << endl;

				if (save_beam == 1) { 

					energy=(sqrt(TRipx[ip]*TRipx[ip]+TRipy[ip]*TRipy[ip]+TRipz[ip]*TRipz[ip]+TRmass[ip]*TRmass[ip])-TRmass[ip])*1000.;
                    kink=0;
                        

                    sqrt1suL=sqrt(1/length);
                    ZTree=TRcha[ip];
					A = TRbar[ip];

					pa_ID = TRpaid[ip];
				
					//angles
					GetDirectionCosine(TRipx[ip], TRipy[ip], TRipz[ip], cosdir);
					ax = TMath::ACos(cosdir[0]);
					ay = TMath::ACos(cosdir[1]);
					az = TMath::ACos(cosdir[2]);

					GetDirectionCosine(TRfpx[ip], TRfpy[ip], TRfpz[ip], cosdir);
					fax = TMath::ACos(cosdir[0]);
					fay = TMath::ACos(cosdir[1]);
					faz = TMath::ACos(cosdir[2]);

					ix = TRix[ip];
					iy = TRiy[ip];
					iz = TRiz[ip];

					fx = TRfx[ip];
					fy = TRfy[ip];
					fz = TRfz[ip];

					treevar->Fill();

				}

            }

            else if((TRpaid[ip]==1) && ( (TRreg[ip]>=nregFirstEm2 && TRreg[ip]<=nregLastEm2) || (TRreg[ip]>=nregFirstPl2 && TRreg[ip]<=nregLastPl2) ) ){ 
//daughters //|| TRpaid[ip]==0

				//if (ev_ID == 26) cout << "Events " << events << ", ev_ID " << ev_ID << ", TRreg[ip] " <<  TRreg[ip] << " nregFirstEM2 " << nregFirstEm2 << ", nregLastPl2 " << nregLastPl2 << ", TRcha " << TRcha[ip] << endl;


                if(justone==0) {beam_interagenti++; justone++;}
                MULT++;

				if (ip==1) first = 1;
				else first = 0;
                
                
                if(TRcha[ip]!=0) { // cariche
                    Z_atom+=TRcha[ip];
                    length=TRtrlen[ip]; // in cm

                    if(TRfz[ip]>ZLAST) { //trk che escono dal rivelatore
                        length=ZLAST+3.;
                    }


                    if(length>L_THRESH){
                        if(justone_v==0) {beam_interagenti_visibili++; justone_v++;}
                        MULT_visible++;
						if (length>0.0950) MULT_visible_long++;
                        prodotti++;
                        
                        energy=(sqrt(TRipx[ip]*TRipx[ip]+TRipy[ip]*TRipy[ip]+TRipz[ip]*TRipz[ip]+TRmass[ip]*TRmass[ip])-TRmass[ip])*1000.;
                        kink=Eval_Kink(TRfpx[ipbeam],TRfpy[ipbeam],TRfpz[ipbeam],TRipx[ip],TRipy[ip],TRipz[ip]);

                        sqrt1suL=sqrt(1/length);
                        ZTree=TRcha[ip];
						A = TRbar[ip];

						pa_ID = TRpaid[ip];
				
						//angles
						GetDirectionCosine(TRipx[ip], TRipy[ip], TRipz[ip], cosdir);
						ax = TMath::ACos(cosdir[0]);
						ay = TMath::ACos(cosdir[1]);
						az = TMath::ACos(cosdir[2]);

						GetDirectionCosine(TRfpx[ip], TRfpy[ip], TRfpz[ip], cosdir);
						fax = TMath::ACos(cosdir[0]);
						fay = TMath::ACos(cosdir[1]);
						faz = TMath::ACos(cosdir[2]);

						//parent angles
						GetDirectionCosine(TRfpx[ipbeam],TRfpy[ipbeam],TRfpz[ipbeam], cosdirparent);
						pax = TMath::ACos(cosdirparent[0]);
						pay = TMath::ACos(cosdirparent[1]);
						paz = TMath::ACos(cosdirparent[2]);

						/*if (TMath::Abs(kink-az)>0.1) {  //controllo regione in cui kink != az

							cout << " --- Event ID " << ev_ID << " --- " << endl;
							cout << " TRfpx [beam] " << TRfpx[ipbeam] << " TRfpy [beam] " << TRfpy[ipbeam] << " TRfpz [beam] " <<	TRfpz[ipbeam] << endl;
							cout << " TRipx " << TRipx[ip] << " TRipy " << TRipy[ip] << " TRipz " << TRipz[ip] << endl;	
							
							cout << " Kink " << kink << endl;
							cout << " Az [fragment] " << az << endl;
							cout << " Direction Cosine X Y Z [Fragment] " << cosdir[0] << " " << cosdir[1] << " " << cosdir[2] << endl;
							cout << " Direction Cosine X Y Z [Parent] " << cosdirparent[0] << " " << cosdirparent[1] << " " << cosdirparent[2] << endl;
							cout << endl;
										

						}*/

						ix = TRix[ip];
						iy = TRiy[ip];
						iz = TRiz[ip];

						fx = TRfx[ip];
						fy = TRfy[ip];
						fz = TRfz[ip];

                        if(TRfz[ip]>ZLAST) { //trk che escono dal rivelatore
                            EXIT_LASTPLATE++;
                            sqrt1suL=200;
                        }


                        int iplot=0;
                        if(TRcha[ip]<N_CHARGE-1) iplot=TRcha[ip]-1;
                        else iplot=N_CHARGE-1;
                                                
                        if(length>L_THRESH_LONG) ch_reco=1;
                        else ch_reco=0;
                        Z_reco+=ch_reco;
                        
                        if(TRfid[ip]==1) p++;
                        else if(TRfid[ip]==-2) HeavIo++;
                        else if(TRfid[ip]==-3) D++;
                        else if(TRfid[ip]==-4) Tr++;
                        else if(TRfid[ip]==-5) He3++;
                        else if(TRfid[ip]==-6) He4++;
                        else others++;
                        
                        //trigger study
                        
                        if(length>L_THRESH_LONG) longtrack++;
                        if(kink>TH_THRESH) largekinktrack++;
                        if(length>L_THRESH_LONG && kink>TH_THRESH) longandkink++;
                        if(TRreg[ip]>=nregFirstEm2 && TRreg[ip]<=nregLastEm2) {startinemu++; emutrk=1;}
                        else if(TRreg[ip]>=nregFirstPl2 && TRreg[ip]<=nregLastPl2) startinplastic++;
                        //if(startinemu>0 && startinplastic>0) cout << "Event " << EventNumber << endl;

                        
                        treevar->Fill();

                    }
                }
            }

        }


        if(TRn>1 && MULT>0){

			h_mult_all->Fill(MULT);
			h_mult_visible->Fill(MULT_visible);
			h_mult_visible_long->Fill(MULT_visible_long);
            
            if(Z_atom-1>30) {
                beam_interagenti_AgBr++;
                if(justone_v>0) beam_interagenti_AgBr_visibili++;
                
            }

            //if(longtrack>0 || largekinktrack>0 || startinemu>1) RECOVTX++;
            if(longandkink>0 ||  (startinemu>1 && longtrack>0)) RECOVTX++;
            if(startinemu>0 && startinemu>=startinplastic) TOT_vtxemu++;
            if(startinplastic>0 && startinplastic>=startinemu) TOT_vtxpl++;
            
        }
    }
    
    fout->Write();
    fout->Close();

    
    //PLOT
    
    /*float max=0, temp_max=0;
    
    
    c_mult->cd();
    h_mult_visible->SetLineColor(kBlue);
    h_mult_visible->SetLineWidth(2);
    h_mult_all->SetLineColor(kRed);
    h_mult_all->SetLineWidth(2);
    h_mult_all->SetLineStyle(2);
    h_mult_AgBr->SetLineColor(kGreen);
    h_mult_AgBr->SetLineWidth(2);
    h_mult_visible->Draw("");
    h_mult_AgBr->Draw("sames");
    h_mult_all->Draw("sames");
    TLegend *legend_visible = new TLegend(.65,.6,.80,.8);
    legend_visible->AddEntry(h_mult_all,"All");
    legend_visible->AddEntry(h_mult_visible,"Visible");
    legend_visible->AddEntry(h_mult_AgBr,"interactions on AgBr");
    legend_visible->Draw();
    sprintf(title, "Plot2/h_mult_%s_%d.pdf", fascio, geo);
    c_mult->SaveAs(title);
    sprintf(title, "Plot2/h_mult_%s_%d.C", fascio, geo);
    c_mult->SaveAs(title);
    
    c_flukaid->cd();
    c_flukaid->SetLogy();
    h_flukaid_visible->SetLineColor(kBlue);
    h_flukaid_visible->SetLineWidth(2);
    h_flukaid_all->SetLineColor(kRed);
    h_flukaid_all->SetLineWidth(2);
    h_flukaid_all->SetLineStyle(2);
    h_flukaid_all->Draw("");
    h_flukaid_visible->Draw("sames");
    legend_visible->Draw();
    sprintf(title, "Plot2/h_flukaid_%s_%d.pdf", fascio, geo);
    c_flukaid->SaveAs(title);
    sprintf(title, "Plot2/h_flukaid_%s_%d.C", fascio, geo);
    c_flukaid->SaveAs(title);
    
    c_atom->cd();
    h_atom_reco->SetLineColor(kBlue);
    h_atom_reco->SetLineWidth(2);
   // h_atom->SetLineColor(kRed);
    h_atom->SetLineWidth(2);
   // h_atom->SetLineStyle(2);
    h_atom->Draw("");
    //h_atom_reco->Draw("sames");
    TLegend *legend_reco = new TLegend(.65,.6,.80,.8);
    legend_reco->AddEntry(h_atom,"All");
    legend_reco->AddEntry(h_atom_reco,"Reco");
   // legend_reco->Draw();
    sprintf(title, "Plot2/h_atom_%s_%d.pdf", fascio, geo);
    c_atom->SaveAs(title);
    sprintf(title, "Plot2/h_atom_%s_%d.C", fascio, geo);
    c_atom->SaveAs(title);
    
    c_charge->cd();
    c_charge->SetLogy();
    h_charge->SetLineColor(kBlue);
    h_charge->SetLineWidth(2);
    h_charge->Draw("");
    sprintf(title, "Plot2/h_charge_%s_%d.pdf", fascio, geo);
    c_charge->SaveAs(title);
    sprintf(title, "Plot2/h_charge_%s_%d.C", fascio, geo);
    c_charge->SaveAs(title);
    
    c_length->Divide(2,1);
    c_length->cd(1);
    c_length->cd(1)->SetLogy();
    h_length->Draw("");

//    c_length->cd(2);
//    c_length->cd(2)->SetLogy();
//    h_length2->Draw("");
//
    c_length->cd(2);
    c_length->cd(2)->SetLogy();
    h_length3->Draw("");

    sprintf(title, "Plot2/h_length_%s_%d.pdf", fascio, geo);
    c_length->SaveAs(title);
    sprintf(title, "Plot2/h_length_%s_%d.C", fascio, geo);
    c_length->SaveAs(title);
    
    c_length2->cd();
    c_length2->SetLogy();
    h_length2->Draw("");
    sprintf(title, "Plot2/h_length2_%s_%d.pdf", fascio, geo);
    c_length2->SaveAs(title);
    sprintf(title, "Plot2/h_length2_%s_%d.C", fascio, geo);
    c_length2->SaveAs(title);
    
    c_length3->cd();
    c_length3->SetLogy();
    h_length3->Draw("");
    sprintf(title, "Plot2/h_length3_%s_%d.pdf", fascio, geo);
    c_length3->SaveAs(title);
    sprintf(title, "Plot2/h_length3_%s_%d.C", fascio, geo);
    c_length3->SaveAs(title);
    
    c_1sulength->cd();
    c_1sulength->SetLogy();
    char titleZ[10];
    TLegend *legend_Z = new TLegend(.65,.6,.80,.8);
    temp_max=0;
    for(int i=0; i<N_CHARGE; i++){
        h_1sulength[i]->SetLineColor(kColors[i]);
        max = h_1sulength[i]->GetMaximum();
        //if(max>temp_max) {h_1sulength[0]->SetMaximum(max+1000); temp_max=max;}
        if(i==0) h_1sulength[i]->Draw("");
        else h_1sulength[i]->Draw("sames");
        if(i<N_CHARGE-1) sprintf(titleZ, "Z = %d", i+1);
        else sprintf(titleZ, "Z >= %d", N_CHARGE);
        legend_Z->AddEntry(h_1sulength[i],titleZ);
    }
    legend_Z->Draw();
    c_1sulength->Update();
    sprintf(title, "Plot2/h_1sulength_%s_%d.pdf", fascio, geo);
    c_1sulength->SaveAs(title);
    sprintf(title, "Plot2/h_1sulength_%s_%d.C", fascio, geo);
    c_1sulength->SaveAs(title);
    
    c_kink->cd();
    h_kink->Draw("");
    sprintf(title, "Plot2/h_kink_%s_%d.pdf", fascio, geo);
    c_kink->SaveAs(title);
    sprintf(title, "Plot2/h_kink_%s_%d.C", fascio, geo);
    c_kink->SaveAs(title);
    
    c_kink_Z->cd();
    temp_max=0;
    for(int i=0; i<N_CHARGE; i++){
        h_kink_Z[i]->SetLineColor(kColors[i]);
        max = h_kink_Z[i]->GetMaximum();
        if(max>temp_max) {h_kink_Z[0]->SetMaximum(max+50); temp_max=max;}
        if(i==0) h_kink_Z[i]->Draw("");
        else h_kink_Z[i]->Draw("sames");
    }
    legend_Z->Draw();
    sprintf(title, "Plot2/h_kink_Z_%s_%d.pdf", fascio, geo);
    c_kink_Z->SaveAs(title);
    sprintf(title, "Plot2/h_kink_Z_%s_%d.C", fascio, geo);
    c_kink_Z->SaveAs(title);
    
    c_kinkdeg_Z->cd();
    temp_max=0;
    for(int i=0; i<N_CHARGE; i++){
        h_kinkdeg_Z[i]->SetLineColor(kColors[i]);
        max = h_kinkdeg_Z[i]->GetMaximum();
        if(max>temp_max) {h_kinkdeg_Z[0]->SetMaximum(max+50); temp_max=max;}
        if(i==0) h_kinkdeg_Z[i]->Draw("");
        else h_kinkdeg_Z[i]->Draw("sames");
    }
    legend_Z->Draw();
    sprintf(title, "Plot2/h_kinkdeg_Z_%s_%d.pdf", fascio, geo);
    c_kinkdeg_Z->SaveAs(title);
    sprintf(title, "Plot2/h_kinkdeg_Z_%s_%d.C", fascio, geo);
    c_kinkdeg_Z->SaveAs(title);
    
    c_1sulengthVskink->cd();
    temp_max=0;
    for(int i=0; i<N_CHARGE; i++){
        h_1sulengthVskink[i]->SetMarkerStyle(6);
        h_1sulengthVskink[i]->SetMarkerColor(kColors[i]);
        h_1sulengthVskink[0]->GetYaxis()->SetRange(0,max+50);
        h_1sulengthVskink[0]->GetXaxis()->SetRange(0,max+100);
        c_1sulengthVskink->cd();
        if(i==0) h_1sulengthVskink[i]->Draw("");
        else h_1sulengthVskink[i]->Draw("sames");
        fout->cd();
        sprintf(titleZ, "%d", i+1);
        h_1sulengthVskink[i]->Write(titleZ);

    }
    legend_Z->Draw();
    sprintf(title, "Plot2/h_1sulengthVskink_%s_%d.pdf", fascio, geo);
    c_1sulengthVskink->SaveAs(title);
    sprintf(title, "Plot2/h_1sulengthVskink_%s_%d.C", fascio, geo);
    c_1sulengthVskink->SaveAs(title);
    
    c_energy->cd();
    c_energy->SetLogy();
    c_energy->Divide(1,2);
    c_energy->cd(1);
    c_energy->cd(1)->SetLogy();
    h_energy->Draw("");
    c_energy->cd(2);
    h_energy_exiting->Draw("");
    sprintf(title, "Plot2/h_energy_%s_%d.pdf", fascio, geo);
    c_energy->SaveAs(title);
    sprintf(title, "Plot2/h_energy_%s_%d.C", fascio, geo);
    c_energy->SaveAs(title);
    
    c_energy_Z->cd();
    c_energy_Z->SetLogy();
    temp_max=0;
    for(int i=0; i<N_CHARGE; i++){
        h_energy_Z[i]->SetLineWidth(2);
        h_energy_Z[i]->SetLineColor(kColors[i]);
        c_energy_Z->cd();
        max = h_energy_Z[i]->GetMaximum();
        if(max>temp_max) {h_energy_Z[0]->SetMaximum(max+100); temp_max=max;}
        if(i==0) h_energy_Z[i]->Draw("");
        else h_energy_Z[i]->Draw("sames");
//        fout->cd();
//        sprintf(titleZ, "%d", i+1);
//        h_energy_Z[i]->Write(titleZ);
    }
    legend_Z->Draw();
    sprintf(title, "Plot2/h_energy_Z_%s_%d.pdf", fascio, geo);
    c_energy_Z->SaveAs(title);
    sprintf(title, "Plot2/h_energy_Z_%s_%d.C", fascio, geo);
    c_energy_Z->SaveAs(title);
    
    c_trigger->cd();
    c_trigger->Divide(2,2);
    c_trigger->cd(1);
    h_longtrack->Draw("");
    c_trigger->cd(2);
    h_largekinktrack->Draw("");
    c_trigger->cd(3);
    h_startinemu->Draw("");
    c_trigger->cd(4);
    h_startinplastic->Draw("");
    sprintf(title, "Plot2/h_trigger_%s_%d.pdf", fascio, geo);
    c_trigger->SaveAs(title);
    sprintf(title, "Plot2/h_trigger_%s_%d.C", fascio, geo);
    c_trigger->SaveAs(title);
    
    fout->Close();*/

                            
    //PRINT
    cout << endl;
    
    cout << endl;
    
    cout << "N beam: " << events << "\t arrived: " << nbeam << endl;
    
    cout << "Frac assorbiti prima del brick: " << (float)(events-nbeam)/events*100 << "\%\t("<< events-nbeam << ")" << endl;
    
    cout << "Frac interagenti: " << (float)beam_interagenti/events*100 << "\%\t ("<< beam_interagenti << ")" << endl;
    cout << "Frac interagenti su AgBr: " << (float)(beam_interagenti_AgBr)/events*100 << "\%\t ("<< beam_interagenti_AgBr << ")" << "\t frac su interagenti "<< (float)(beam_interagenti_AgBr)/beam_interagenti*100 << endl;
    cout << "Frac interagenti esclusi AgBr: " << (float)(beam_interagenti-beam_interagenti_AgBr)/events*100 << "\%\t ("<< beam_interagenti-beam_interagenti_AgBr << ")" << "\t frac su interagenti "<< (float)(beam_interagenti-beam_interagenti_AgBr)/beam_interagenti*100 << endl;
    cout << "Frac interagenti (tutti) visibili: " << (float)beam_interagenti_visibili/events*100 << "\%\t ("<< beam_interagenti_visibili << ")" << endl;


    
    
    cout << "Beam passanti: " << (float)nbeam_passanti/events*100 << "\%\t (" << nbeam_passanti << ")" << endl << endl;
    
    //cout << "prodotti EXIT_LATERAL: " << EXIT_LATERAL << "(" << (float)EXIT_LATERAL/(prodotti)*100 << "\%)\t";
    cout << std::setprecision(3) << "Secondari (Percentuale su prodotti visibili (" << prodotti << ")): --> " << endl << "\tp=" << p << "("<< (float)p/(prodotti)*100 << "%)" << endl << "\tD=" << D << "(" << (float)D/(prodotti)*100 << "\%)" << endl << "\tTr=" << Tr << "(" << (float)Tr/(prodotti)*100 << "\%)" << endl << "\tHe3=" << He3 << "(" << (float)He3/(prodotti)*100 << "\%)" << endl << "\tHe4=" << He4 << "(" << (float)He4/(prodotti)*100 << "\%)" << endl << "\tHeavy Ions=" << HeavIo << "(" << (float)HeavIo/(prodotti)*100 << "\%)" << endl << "\tothers=" << others << "(" << (float)others/(prodotti)*100 << "\%)" << endl;
    
    cout << "EXIT_LASTPLATE: " << EXIT_LASTPLATE << " (" << (float)EXIT_LASTPLATE/(prodotti)*100 << "\%)"<< endl ;

    
    /*cout << "total tracks: " << h_length->Integral() << "\t" << h_length2->Integral() << endl;
    cout << "Trks below " << L_THRESH << " cm: " << h_length3->Integral(h_length3->FindFixBin(0), h_length3->FindFixBin(L_THRESH*10000)) << "\t(" << h_length3->Integral(h_length3->FindFixBin(0), h_length3->FindFixBin(L_THRESH*10000))/h_length->Integral()*100 << "\%)" << endl;
    cout << "Trks between " << L_THRESH << "cm and " << L_THRESH*5 << " cm: " << h_length3->Integral(h_length3->FindFixBin(L_THRESH*10000), h_length3->FindFixBin(L_THRESH*10000*5)) << "\t(" << h_length3->Integral(h_length3->FindFixBin(L_THRESH*10000), h_length3->FindFixBin(L_THRESH*10000*5))/h_length->Integral()*100 << "\%)" << endl;
    cout << "Trks between " << 0 << "cm and " << L_THRESH_LONG << " cm: " << h_length2->Integral(h_length2->FindFixBin(0), h_length2->FindFixBin(L_THRESH_LONG*10000)) << "\t(" << h_length2->Integral(h_length2->FindFixBin(0), h_length2->FindFixBin(L_THRESH_LONG*10000))/h_length->Integral()*100 << "\%)" << endl;*/


    
    cout << endl << "Potentially Visible vtx: " << beam_interagenti_visibili << "\t(" << (float)beam_interagenti_visibili/beam_interagenti*100 << "\%)" << endl;
    cout << "\tTOT_vtxemu: " << TOT_vtxemu << "\t(" << (float)TOT_vtxemu/(beam_interagenti_visibili)*100 << "\%)" <<  endl;
    cout << "\tTOT_vtxpl: " << TOT_vtxpl << "\t(" << (float)TOT_vtxpl/(beam_interagenti_visibili)*100 << "\%)" <<  endl;
    cout << "Reco* vtx: " << RECOVTX << "\t(" << (float)RECOVTX/beam_interagenti_visibili*100 << "\%)" << endl;
    cout << "\t *Una trk è reco se: " << endl << "- Z!=0" << endl << "- ha una lunghezza > " << L_THRESH*10000 << " um" << endl;
    cout << "\t *Un vtx è reco se vale una delle seguenti: " << endl << "- almeno una traccia più lunga di " << L_THRESH_LONG*10000 << " um e con kink > " << TH_THRESH << " rad" << endl << "- almeno 2 tracce che iniziano in emulsione e una di lunghezza > " << L_THRESH_LONG*10000 << " um" << endl;
    
    
    
    cout << endl << endl;
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


TFile * outputfile;

void merge_files(){

 //outputfile = new TFile(Form("/eos/user/v/viboccia/PRIN22/P018_bot/mt.merged.all.root"),"RECREATE"); 
 //outputfile = new TFile(Form("./mt.merged.all.root"),"RECREATE"); 
 TChain ch("mtracks");

 for (int n=0; n<17493; n = n+200) {//17493
    TString filename = Form("/eos/user/v/viboccia/PRIN22/P018_bot/mt.merged.temp.%i.root",n);
    //TString filename = Form("./mt.merged.temp.%i.root",n);

    if (gSystem->AccessPathName(filename.Data())){//returns False if file exists, True if it does not (yes, I find it weird too)
        cout<<"File does not exist:"<<filename.Data()<<" moving to next file "<<endl;
        continue;
    }

    ch.Add(filename);

 }

 ch.Merge("/eos/user/v/viboccia/PRIN22/P018_bot/mt.merged.all.root");
 //outputfile->cd();
 //ch.Write();
 //outputfile->Close();

}

int main() {
    merge_files();
    return 0;
}
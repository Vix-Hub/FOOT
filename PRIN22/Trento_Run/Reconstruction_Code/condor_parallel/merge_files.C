
TFile * outputfile;
TString path = "/eos/user/v/viboccia/PRIN22/P018_bot/";

void merge_files(){

 //outputfile = new TFile(Form("/eos/user/v/viboccia/PRIN22/P018_bot/mt.merged.all.root"),"RECREATE"); 
 //outputfile = new TFile(Form("./mt.merged.all.root"),"RECREATE"); 
 TChain ch("mtracks");

 for (int n=0; n<17493; n = n+200) {//17493
    TString filename = Form("%smt.merged.temp.%i.root", path.Data(), n);
    //TString filename = Form("./mt.merged.temp.%i.root",n);

    if (gSystem->AccessPathName(filename.Data())){//returns False if file exists, True if it does not (yes, I find it weird too)
        cout<<"File does not exist:"<<filename.Data()<<" moving to next file "<<endl;
        continue;
    }

    ch.Add(filename);

 }

 TFile *outputfile = new TFile(Form("%smt.merged.all.root", path.Data()), "RECREATE");
    
    // Make sure the TFile is open before merging
    if (outputfile->IsOpen()) {
        ch.Merge(outputfile);
        outputfile->Close();
    } else {
        cout << "Error opening output file!" << endl;
    }

}

int main() {
    merge_files();
    return 0;
}
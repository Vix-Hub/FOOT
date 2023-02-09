//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar 17 16:42:42 2022 by ROOT version 6.12/04
// from TTree EventTree/gsimay
// found on file: foot_ecc.root
//////////////////////////////////////////////////////////

#ifndef EventTree_h
#define EventTree_h
#define RAD2DEG  57.2958

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class EventTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

    /*const int nregFirstEm1 = 47;
    const int nregLastEm1 = 48;
    const int nregFirstEm2 = 49;
    const int nregLastEm2 = 248;
    const int nregFirstPl1 = 357; //259;
    const int nregFirstPl2 = 358; //260;
    const int nregLastPl2 = 458; //360;
    const int nregLastPB = 565; //no
    const float ZLAST = 8; //1.96;*/


   // Declaration of leaf types
   Int_t           EventNumber;
   Int_t           TRn;
   Int_t           TRpaid[97];   //[TRn]
   Int_t           TRgen[97];   //[TRn]
   Int_t           TRcha[97];   //[TRn]
   Int_t           TRreg[97];   //[TRn]
   Int_t           TRbar[97];   //[TRn]
   Int_t           TRdead[97];   //[TRn]
   Int_t           TRfid[97];   //[TRn]
   Double_t        TRix[97];   //[TRn]
   Double_t        TRiy[97];   //[TRn]
   Double_t        TRiz[97];   //[TRn]
   Double_t        TRfx[97];   //[TRn]
   Double_t        TRfy[97];   //[TRn]
   Double_t        TRfz[97];   //[TRn]
   Double_t        TRipx[97];   //[TRn]
   Double_t        TRipy[97];   //[TRn]
   Double_t        TRipz[97];   //[TRn]
   Double_t        TRfpx[97];   //[TRn]
   Double_t        TRfpy[97];   //[TRn]
   Double_t        TRfpz[97];   //[TRn]
   Double_t        TRmass[97];   //[TRn]
   Double_t        TRtime[97];   //[TRn]
   Double_t        TRtof[97];   //[TRn]
   Double_t        TRtrlen[97];   //[TRn]
   Int_t           STCn;
   Int_t           STCid[7];   //[STCn]
   Double_t        STCxin[7];   //[STCn]
   Double_t        STCyin[7];   //[STCn]
   Double_t        STCzin[7];   //[STCn]
   Double_t        STCpxin[7];   //[STCn]
   Double_t        STCpyin[7];   //[STCn]
   Double_t        STCpzin[7];   //[STCn]
   Double_t        STCxout[7];   //[STCn]
   Double_t        STCyout[7];   //[STCn]
   Double_t        STCzout[7];   //[STCn]
   Double_t        STCpxout[7];   //[STCn]
   Double_t        STCpyout[7];   //[STCn]
   Double_t        STCpzout[7];   //[STCn]
   Double_t        STCde[7];   //[STCn]
   Double_t        STCal[7];   //[STCn]
   Double_t        STCtim[7];   //[STCn]
   Int_t           BMNn;
   Int_t           BMNid[29];   //[BMNn]
   Int_t           BMNilay[29];   //[BMNn]
   Int_t           BMNiview[29];   //[BMNn]
   Int_t           BMNicell[29];   //[BMNn]
   Double_t        BMNxin[29];   //[BMNn]
   Double_t        BMNyin[29];   //[BMNn]
   Double_t        BMNzin[29];   //[BMNn]
   Double_t        BMNpxin[29];   //[BMNn]
   Double_t        BMNpyin[29];   //[BMNn]
   Double_t        BMNpzin[29];   //[BMNn]
   Double_t        BMNxout[29];   //[BMNn]
   Double_t        BMNyout[29];   //[BMNn]
   Double_t        BMNzout[29];   //[BMNn]
   Double_t        BMNpxout[29];   //[BMNn]
   Double_t        BMNpyout[29];   //[BMNn]
   Double_t        BMNpzout[29];   //[BMNn]
   Double_t        BMNde[29];   //[BMNn]
   Double_t        BMNal[29];   //[BMNn]
   Double_t        BMNtim[29];   //[BMNn]
   Int_t           CROSSn;
   Int_t           CROSSid[7361];   //[CROSSn]
   Int_t           CROSSnreg[7361];   //[CROSSn]
   Int_t           CROSSnregold[7361];   //[CROSSn]
   Double_t        CROSSx[7361];   //[CROSSn]
   Double_t        CROSSy[7361];   //[CROSSn]
   Double_t        CROSSz[7361];   //[CROSSn]
   Double_t        CROSSpx[7361];   //[CROSSn]
   Double_t        CROSSpy[7361];   //[CROSSn]
   Double_t        CROSSpz[7361];   //[CROSSn]
   Double_t        CROSSm[7361];   //[CROSSn]
   Double_t        CROSSch[7361];   //[CROSSn]
   Double_t        CROSSt[7361];   //[CROSSn]

   // List of branches
   TBranch        *b_EventNumber;   //!
   TBranch        *b_TRn;   //!
   TBranch        *b_TRpaid;   //!
   TBranch        *b_TRgen;   //!
   TBranch        *b_TRcha;   //!
   TBranch        *b_TRreg;   //!
   TBranch        *b_TRbar;   //!
   TBranch        *b_TRdead;   //!
   TBranch        *b_TRfid;   //!
   TBranch        *b_TRix;   //!
   TBranch        *b_TRiy;   //!
   TBranch        *b_TRiz;   //!
   TBranch        *b_TRfx;   //!
   TBranch        *b_TRfy;   //!
   TBranch        *b_TRfz;   //!
   TBranch        *b_TRipx;   //!
   TBranch        *b_TRipy;   //!
   TBranch        *b_TRipz;   //!
   TBranch        *b_TRfpx;   //!
   TBranch        *b_TRfpy;   //!
   TBranch        *b_TRfpz;   //!
   TBranch        *b_TRmass;   //!
   TBranch        *b_TRtime;   //!
   TBranch        *b_TRtof;   //!
   TBranch        *b_TRtrlen;   //!
   TBranch        *b_STCn;   //!
   TBranch        *b_STCid;   //!
   TBranch        *b_STCxin;   //!
   TBranch        *b_STCyin;   //!
   TBranch        *b_STCzin;   //!
   TBranch        *b_STCpxin;   //!
   TBranch        *b_STCpyin;   //!
   TBranch        *b_STCpzin;   //!
   TBranch        *b_STCxout;   //!
   TBranch        *b_STCyout;   //!
   TBranch        *b_STCzout;   //!
   TBranch        *b_STCpxout;   //!
   TBranch        *b_STCpyout;   //!
   TBranch        *b_STCpzout;   //!
   TBranch        *b_STCde;   //!
   TBranch        *b_STCal;   //!
   TBranch        *b_STCtim;   //!
   TBranch        *b_BMNn;   //!
   TBranch        *b_BMNid;   //!
   TBranch        *b_BMNilay;   //!
   TBranch        *b_BMNiview;   //!
   TBranch        *b_BMNicell;   //!
   TBranch        *b_BMNxin;   //!
   TBranch        *b_BMNyin;   //!
   TBranch        *b_BMNzin;   //!
   TBranch        *b_BMNpxin;   //!
   TBranch        *b_BMNpyin;   //!
   TBranch        *b_BMNpzin;   //!
   TBranch        *b_BMNxout;   //!
   TBranch        *b_BMNyout;   //!
   TBranch        *b_BMNzout;   //!
   TBranch        *b_BMNpxout;   //!
   TBranch        *b_BMNpyout;   //!
   TBranch        *b_BMNpzout;   //!
   TBranch        *b_BMNde;   //!
   TBranch        *b_BMNal;   //!
   TBranch        *b_BMNtim;   //!
   TBranch        *b_CROSSn;   //!
   TBranch        *b_CROSSid;   //!
   TBranch        *b_CROSSnreg;   //!
   TBranch        *b_CROSSnregold;   //!
   TBranch        *b_CROSSx;   //!
   TBranch        *b_CROSSy;   //!
   TBranch        *b_CROSSz;   //!
   TBranch        *b_CROSSpx;   //!
   TBranch        *b_CROSSpy;   //!
   TBranch        *b_CROSSpz;   //!
   TBranch        *b_CROSSm;   //!
   TBranch        *b_CROSSch;   //!
   TBranch        *b_CROSSt;   //!

   EventTree(TTree *tree=0);
   virtual ~EventTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
    
    
    int geo;
    char fascio[50];
};

#endif

#ifdef EventTree_cxx
EventTree::EventTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ECC_emuasada.root"); //foot_ecc
      if (!f || !f->IsOpen()) {
         f = new TFile("ECC_emuasada.root");
      }
      f->GetObject("EventTree",tree);

   }
   Init(tree);
}

EventTree::~EventTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EventTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EventTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EventTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("TRn", &TRn, &b_TRn);
   fChain->SetBranchAddress("TRpaid", TRpaid, &b_TRpaid);
   fChain->SetBranchAddress("TRgen", TRgen, &b_TRgen);
   fChain->SetBranchAddress("TRcha", TRcha, &b_TRcha);
   fChain->SetBranchAddress("TRreg", TRreg, &b_TRreg);
   fChain->SetBranchAddress("TRbar", TRbar, &b_TRbar);
   fChain->SetBranchAddress("TRdead", TRdead, &b_TRdead);
   fChain->SetBranchAddress("TRfid", TRfid, &b_TRfid);
   fChain->SetBranchAddress("TRix", TRix, &b_TRix);
   fChain->SetBranchAddress("TRiy", TRiy, &b_TRiy);
   fChain->SetBranchAddress("TRiz", TRiz, &b_TRiz);
   fChain->SetBranchAddress("TRfx", TRfx, &b_TRfx);
   fChain->SetBranchAddress("TRfy", TRfy, &b_TRfy);
   fChain->SetBranchAddress("TRfz", TRfz, &b_TRfz);
   fChain->SetBranchAddress("TRipx", TRipx, &b_TRipx);
   fChain->SetBranchAddress("TRipy", TRipy, &b_TRipy);
   fChain->SetBranchAddress("TRipz", TRipz, &b_TRipz);
   fChain->SetBranchAddress("TRfpx", TRfpx, &b_TRfpx);
   fChain->SetBranchAddress("TRfpy", TRfpy, &b_TRfpy);
   fChain->SetBranchAddress("TRfpz", TRfpz, &b_TRfpz);
   fChain->SetBranchAddress("TRmass", TRmass, &b_TRmass);
   fChain->SetBranchAddress("TRtime", TRtime, &b_TRtime);
   fChain->SetBranchAddress("TRtof", TRtof, &b_TRtof);
   fChain->SetBranchAddress("TRtrlen", TRtrlen, &b_TRtrlen);
   fChain->SetBranchAddress("STCn", &STCn, &b_STCn);
   fChain->SetBranchAddress("STCid", STCid, &b_STCid);
   fChain->SetBranchAddress("STCxin", STCxin, &b_STCxin);
   fChain->SetBranchAddress("STCyin", STCyin, &b_STCyin);
   fChain->SetBranchAddress("STCzin", STCzin, &b_STCzin);
   fChain->SetBranchAddress("STCpxin", STCpxin, &b_STCpxin);
   fChain->SetBranchAddress("STCpyin", STCpyin, &b_STCpyin);
   fChain->SetBranchAddress("STCpzin", STCpzin, &b_STCpzin);
   fChain->SetBranchAddress("STCxout", STCxout, &b_STCxout);
   fChain->SetBranchAddress("STCyout", STCyout, &b_STCyout);
   fChain->SetBranchAddress("STCzout", STCzout, &b_STCzout);
   fChain->SetBranchAddress("STCpxout", STCpxout, &b_STCpxout);
   fChain->SetBranchAddress("STCpyout", STCpyout, &b_STCpyout);
   fChain->SetBranchAddress("STCpzout", STCpzout, &b_STCpzout);
   fChain->SetBranchAddress("STCde", STCde, &b_STCde);
   fChain->SetBranchAddress("STCal", STCal, &b_STCal);
   fChain->SetBranchAddress("STCtim", STCtim, &b_STCtim);
   fChain->SetBranchAddress("BMNn", &BMNn, &b_BMNn);
   fChain->SetBranchAddress("BMNid", BMNid, &b_BMNid);
   fChain->SetBranchAddress("BMNilay", BMNilay, &b_BMNilay);
   fChain->SetBranchAddress("BMNiview", BMNiview, &b_BMNiview);
   fChain->SetBranchAddress("BMNicell", BMNicell, &b_BMNicell);
   fChain->SetBranchAddress("BMNxin", BMNxin, &b_BMNxin);
   fChain->SetBranchAddress("BMNyin", BMNyin, &b_BMNyin);
   fChain->SetBranchAddress("BMNzin", BMNzin, &b_BMNzin);
   fChain->SetBranchAddress("BMNpxin", BMNpxin, &b_BMNpxin);
   fChain->SetBranchAddress("BMNpyin", BMNpyin, &b_BMNpyin);
   fChain->SetBranchAddress("BMNpzin", BMNpzin, &b_BMNpzin);
   fChain->SetBranchAddress("BMNxout", BMNxout, &b_BMNxout);
   fChain->SetBranchAddress("BMNyout", BMNyout, &b_BMNyout);
   fChain->SetBranchAddress("BMNzout", BMNzout, &b_BMNzout);
   fChain->SetBranchAddress("BMNpxout", BMNpxout, &b_BMNpxout);
   fChain->SetBranchAddress("BMNpyout", BMNpyout, &b_BMNpyout);
   fChain->SetBranchAddress("BMNpzout", BMNpzout, &b_BMNpzout);
   fChain->SetBranchAddress("BMNde", BMNde, &b_BMNde);
   fChain->SetBranchAddress("BMNal", BMNal, &b_BMNal);
   fChain->SetBranchAddress("BMNtim", BMNtim, &b_BMNtim);
   fChain->SetBranchAddress("CROSSn", &CROSSn, &b_CROSSn);
   fChain->SetBranchAddress("CROSSid", CROSSid, &b_CROSSid);
   fChain->SetBranchAddress("CROSSnreg", CROSSnreg, &b_CROSSnreg);
   fChain->SetBranchAddress("CROSSnregold", CROSSnregold, &b_CROSSnregold);
   fChain->SetBranchAddress("CROSSx", CROSSx, &b_CROSSx);
   fChain->SetBranchAddress("CROSSy", CROSSy, &b_CROSSy);
   fChain->SetBranchAddress("CROSSz", CROSSz, &b_CROSSz);
   fChain->SetBranchAddress("CROSSpx", CROSSpx, &b_CROSSpx);
   fChain->SetBranchAddress("CROSSpy", CROSSpy, &b_CROSSpy);
   fChain->SetBranchAddress("CROSSpz", CROSSpz, &b_CROSSpz);
   fChain->SetBranchAddress("CROSSm", CROSSm, &b_CROSSm);
   fChain->SetBranchAddress("CROSSch", CROSSch, &b_CROSSch);
   fChain->SetBranchAddress("CROSSt", CROSSt, &b_CROSSt);
   Notify();
}

Bool_t EventTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EventTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EventTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

double modulo(double x, double y, double z){
    return sqrt(pow(x,2)+pow(y,2)+pow(z,2));
}

double distanza2D(double x1, double y1, double x2, double y2){
    return sqrt(pow(x1-x2,2)+pow(y1-y2,2));
}

double distanza3D(double x1, double y1, double z1, double x2, double y2, double z2){
    return sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2));
}

void GetDirectionCosine(double px,double py,double pz, Double_t cosdir[3]){
        
    double p = modulo(px,py,pz);
    
    cosdir[0]=px/p;
    cosdir[1]=py/p;
    cosdir[2]=pz/p;
}

Int_t GetDirectionCosineSegment(double sx, double sy, Double_t cosdir[3]){
    
    double rad = sqrt(sx*sx+sy*sy+1);
    cosdir[0]=sx/rad;
    cosdir[1]=sy/rad;
    cosdir[2]=1/rad;
    
    return 0;
}

double Eval_Kink(double px1, double py1, double pz1, double px2, double py2, double pz2){
    
    Double_t kink=0.;
    
    Double_t cosdir1[3];
    Double_t cosdir2[3];
    
    GetDirectionCosine(px1, py1, pz1, cosdir1);
    GetDirectionCosine(px2, py2, pz2, cosdir2);
    
    kink = TMath::ACos(cosdir1[0]*cosdir2[0]+cosdir1[1]*cosdir2[1]+cosdir1[2]*cosdir2[2]);
    
    kink = kink>=0. ? floor(kink*1000.)/1000. : ceil(kink*1000.)/1000.; //ARROTONDO AL MILLIRADIANTE
    
    return kink;
}

#endif // #ifdef EventTree_cxx

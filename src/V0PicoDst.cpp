
#include "V0PicoDst.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void V0PicoDst::Loop()
{
//   In a ROOT session, you can do:
//      root> .L V0PicoDst.C
//      root> V0PicoDst t
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
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}


V0PicoDst::V0PicoDst(TTree *tree) : fChain(0) 
{

   Init(tree);
}

V0PicoDst::~V0PicoDst()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t V0PicoDst::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t V0PicoDst::LoadTree(Long64_t entry)
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

void V0PicoDst::Init(TTree *tree)
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

   fChain->SetBranchAddress("runId", &runId, &b_runId);
   fChain->SetBranchAddress("cent9", &cent9, &b_cent9);
   fChain->SetBranchAddress("cent16", &cent16, &b_cent16);
   fChain->SetBranchAddress("eventWeight", &eventWeight, &b_eventWeight);
   fChain->SetBranchAddress("corrRefMult", &corrRefMult, &b_corrRefMult);
   fChain->SetBranchAddress("tracks", &tracks_, &b_tracks_);
   fChain->SetBranchAddress("tracks.fUniqueID", tracks_fUniqueID, &b_tracks_fUniqueID);
   fChain->SetBranchAddress("tracks.fBits", tracks_fBits, &b_tracks_fBits);
   fChain->SetBranchAddress("tracks.mMomentum.fUniqueID", tracks_mMomentum_fUniqueID, &b_tracks_mMomentum_fUniqueID);
   fChain->SetBranchAddress("tracks.mMomentum.fBits", tracks_mMomentum_fBits, &b_tracks_mMomentum_fBits);
   fChain->SetBranchAddress("tracks.mMomentum.fX", tracks_mMomentum_fX, &b_tracks_mMomentum_fX);
   fChain->SetBranchAddress("tracks.mMomentum.fY", tracks_mMomentum_fY, &b_tracks_mMomentum_fY);
   fChain->SetBranchAddress("tracks.mMomentum.fZ", tracks_mMomentum_fZ, &b_tracks_mMomentum_fZ);
   fChain->SetBranchAddress("tracks.mDedx", tracks_mDedx, &b_tracks_mDedx);
   fChain->SetBranchAddress("tracks.mBeta", tracks_mBeta, &b_tracks_mBeta);
   fChain->SetBranchAddress("tracks.mPathLength", tracks_mPathLength, &b_tracks_mPathLength);
   fChain->SetBranchAddress("tracks.mYLocal", tracks_mYLocal, &b_tracks_mYLocal);
   fChain->SetBranchAddress("tracks.mZLocal", tracks_mZLocal, &b_tracks_mZLocal);
   fChain->SetBranchAddress("tracks.mNHitsFit", tracks_mNHitsFit, &b_tracks_mNHitsFit);
   fChain->SetBranchAddress("tracks.mNHitsPoss", tracks_mNHitsPoss, &b_tracks_mNHitsPoss);
   fChain->SetBranchAddress("tracks.mNHitsDedx", tracks_mNHitsDedx, &b_tracks_mNHitsDedx);
   fChain->SetBranchAddress("tracks.mDcaParams[6]", tracks_mDcaParams, &b_tracks_mDcaParams);
   fChain->SetBranchAddress("tracks.mDcaMatrix[15]", tracks_mDcaMatrix, &b_tracks_mDcaMatrix);
   fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
   Notify();
}

Bool_t V0PicoDst::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void V0PicoDst::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t V0PicoDst::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

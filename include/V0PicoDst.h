//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan  5 17:10:22 2016 by ROOT version 6.04/00
// from TTree V0PicoDst/V0 Pico Dst
// found on file: qa_1CFE6C9E154D39944B9C3A3E1566EBFE_99.root
//////////////////////////////////////////////////////////

#ifndef V0PicoDst_h
#define V0PicoDst_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TObject.h"
#include "TVector3.h"

class V0PicoDst {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static const Int_t kMaxtracks = 10000;

   // Declaration of leaf types
 //StPicoEvent     *event;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Int_t           bField;
   Int_t           runId;
   Char_t          cent9;
   Char_t          cent16;
   Float_t         corrRefMult;
   Float_t         eventWeight;
   Double_t        vtx_mX1;
   Double_t        vtx_mX2;
   Double_t        vtx_mX3;
   Int_t           tracks_;
   UInt_t          tracks_fUniqueID[kMaxtracks];   //[tracks_]
   UInt_t          tracks_fBits[kMaxtracks];   //[tracks_]
   UInt_t          tracks_mMomentum_fUniqueID[kMaxtracks];   //[tracks_]
   UInt_t          tracks_mMomentum_fBits[kMaxtracks];   //[tracks_]
   Double_t        tracks_mMomentum_fX[kMaxtracks];   //[tracks_]
   Double_t        tracks_mMomentum_fY[kMaxtracks];   //[tracks_]
   Double_t        tracks_mMomentum_fZ[kMaxtracks];   //[tracks_]
   Float_t         tracks_mDedx[kMaxtracks];   //[tracks_]
   Float_t         tracks_mBeta[kMaxtracks];   //[tracks_]
   Float_t         tracks_mPathLength[kMaxtracks];   //[tracks_]
   Float_t         tracks_mYLocal[kMaxtracks];   //[tracks_]
   Float_t         tracks_mZLocal[kMaxtracks];   //[tracks_]
   Char_t          tracks_mNHitsFit[kMaxtracks];   //[tracks_]
   Char_t          tracks_mNHitsPoss[kMaxtracks];   //[tracks_]
   UChar_t         tracks_mNHitsDedx[kMaxtracks];   //[tracks_]
   Float_t         tracks_mDcaParams[kMaxtracks][6];   //[tracks_]
   Float_t         tracks_mDcaMatrix[kMaxtracks][15];   //[tracks_]
   Int_t           nTracks;

   // List of branches
   TBranch        *b_event_fUniqueID;   //!
   TBranch        *b_event_fBits;   //!
   TBranch        *b_event_bField;   //!
   TBranch        *b_event_runId;   //!
   TBranch        *b_event_cent9;   //!
   TBranch        *b_event_cent16;   //!
   TBranch        *b_event_corrRefMult;   //!
   TBranch        *b_event_eventWeight;   //!
   TBranch        *b_event_vtx_mX1;   //!
   TBranch        *b_event_vtx_mX2;   //!
   TBranch        *b_event_vtx_mX3;   //!
   TBranch        *b_tracks_;   //!
   TBranch        *b_tracks_fUniqueID;   //!
   TBranch        *b_tracks_fBits;   //!
   TBranch        *b_tracks_mMomentum_fUniqueID;   //!
   TBranch        *b_tracks_mMomentum_fBits;   //!
   TBranch        *b_tracks_mMomentum_fX;   //!
   TBranch        *b_tracks_mMomentum_fY;   //!
   TBranch        *b_tracks_mMomentum_fZ;   //!
   TBranch        *b_tracks_mDedx;   //!
   TBranch        *b_tracks_mBeta;   //!
   TBranch        *b_tracks_mPathLength;   //!
   TBranch        *b_tracks_mYLocal;   //!
   TBranch        *b_tracks_mZLocal;   //!
   TBranch        *b_tracks_mNHitsFit;   //!
   TBranch        *b_tracks_mNHitsPoss;   //!
   TBranch        *b_tracks_mNHitsDedx;   //!
   TBranch        *b_tracks_mDcaParams;   //!
   TBranch        *b_tracks_mDcaMatrix;   //!
   TBranch        *b_nTracks;   //!

   V0PicoDst(TTree *tree=0);
   virtual ~V0PicoDst();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif
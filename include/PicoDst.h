//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 21 21:05:46 2015 by ROOT version 5.26/00
// from TTree PicoDst/StPicoDst
// found on file: st_mtd_15089023_raw_4500029.picoDst.root
//////////////////////////////////////////////////////////

#ifndef PicoDst_h
#define PicoDst_h

#include <iostream>
using namespace std;

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
const Int_t kMaxEvent = 1;
const Int_t kMaxEventPlane = 1;
const Int_t kMaxTracks = 1198;
const Int_t kMaxEmcTrigger = 30;
const Int_t kMaxMtdTrigger = 1;
const Int_t kMaxBTOWHit = 1;
const Int_t kMaxBTofHit = 1;
const Int_t kMaxMtdHit = 19;
const Int_t kMaxEmcPidTraits = 369;
const Int_t kMaxBTofPidTraits = 610;
const Int_t kMaxMtdPidTraits = 8;
const Int_t kMaxV0Ks = 1;
const Int_t kMaxV0L = 1;
const Int_t kMaxV0Lbar = 1;

class PicoDst {
  
private:
  TFile *fileH;
  Char_t *outputFile;
  
  public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Declaration of leaf types
  Int_t           Event_;
  Int_t           Event_mRunId[kMaxEvent];   //[Event_]
  Int_t           Event_mEventId[kMaxEvent];   //[Event_]
  UShort_t        Event_mFillId[kMaxEvent];   //[Event_]
  Float_t         Event_mBField[kMaxEvent];   //[Event_]
  Float_t         Event_mPrimaryVertex_mX1[kMaxEvent];   //[Event_]
  Float_t         Event_mPrimaryVertex_mX2[kMaxEvent];   //[Event_]
  Float_t         Event_mPrimaryVertex_mX3[kMaxEvent];   //[Event_]
  Float_t         Event_mSecondVertex_mX1[kMaxEvent];   //[Event_]
  Float_t         Event_mSecondVertex_mX2[kMaxEvent];   //[Event_]
  Float_t         Event_mSecondVertex_mX3[kMaxEvent];   //[Event_]
  Float_t         Event_mPrimaryVertexError_mX1[kMaxEvent];   //[Event_]
  Float_t         Event_mPrimaryVertexError_mX2[kMaxEvent];   //[Event_]
  Float_t         Event_mPrimaryVertexError_mX3[kMaxEvent];   //[Event_]
  UInt_t          Event_mTriggerWord[kMaxEvent];   //[Event_]
  UInt_t          Event_mTriggerWordMtd[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMultFtpcEast[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMultFtpcWest[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMultNeg[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMultPos[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMult2NegEast[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMult2PosEast[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMult2NegWest[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMult2PosWest[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMultHalfNegEast[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMultHalfPosEast[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMultHalfNegWest[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMultHalfPosWest[kMaxEvent];   //[Event_]
  UShort_t        Event_mGRefMult[kMaxEvent];   //[Event_]
  UShort_t        Event_mNHitsHFT[kMaxEvent][4];   //[Event_]
  UChar_t         Event_mNVpdHitsEast[kMaxEvent];   //[Event_]
  UChar_t         Event_mNVpdHitsWest[kMaxEvent];   //[Event_]
  UShort_t        Event_mNT0[kMaxEvent];   //[Event_]
  Short_t         Event_mVzVpd[kMaxEvent];   //[Event_]
  UInt_t          Event_mZDCx[kMaxEvent];   //[Event_]
  UInt_t          Event_mBBCx[kMaxEvent];   //[Event_]
  Float_t         Event_mBackgroundRate[kMaxEvent];   //[Event_]
  Float_t         Event_mBbcBlueBackgroundRate[kMaxEvent];   //[Event_]
  Float_t         Event_mBbcYellowBackgroundRate[kMaxEvent];   //[Event_]
  Float_t         Event_mBbcEastRate[kMaxEvent];   //[Event_]
  Float_t         Event_mBbcWestRate[kMaxEvent];   //[Event_]
  Float_t         Event_mZdcEastRate[kMaxEvent];   //[Event_]
  Float_t         Event_mZdcWestRate[kMaxEvent];   //[Event_]
  UShort_t        Event_mVpd[kMaxEvent][64];   //[Event_]
  UShort_t        Event_mZdcSumAdcEast[kMaxEvent];   //[Event_]
  UShort_t        Event_mZdcSumAdcWest[kMaxEvent];   //[Event_]
  UShort_t        Event_mZdcSmdEastHorizontal[kMaxEvent][8];   //[Event_]
  UShort_t        Event_mZdcSmdEastVertical[kMaxEvent][8];   //[Event_]
  UShort_t        Event_mZdcSmdWestHorizontal[kMaxEvent][8];   //[Event_]
  UShort_t        Event_mZdcSmdWestVertical[kMaxEvent][8];   //[Event_]
  Float_t         Event_mSpaceCharge[kMaxEvent];   //[Event_]
  UShort_t        Event_mbTofTrayMultiplicity[kMaxEvent];   //[Event_]
  UShort_t        Event_mNumberOfGlobalTracks[kMaxEvent];   //[Event_]
  Float_t         Event_mRanking[kMaxEvent];   //[Event_]
  UShort_t        Event_mNBEMCMatch[kMaxEvent];   //[Event_]
  UShort_t        Event_mNBTOFMatch[kMaxEvent];   //[Event_]
  UShort_t        Event_mBbcAdcEast[kMaxEvent][24];   //[Event_]
  UShort_t        Event_mBbcAdcWest[kMaxEvent][24];   //[Event_]
  UChar_t         Event_mHT_Th[kMaxEvent][4];   //[Event_]
  UChar_t         Event_mJP_Th[kMaxEvent][4];   //[Event_]
  Int_t           EventPlane_;
  UInt_t          EventPlane_fUniqueID[kMaxEventPlane];   //[EventPlane_]
  UInt_t          EventPlane_fBits[kMaxEventPlane];   //[EventPlane_]
  Float_t         EventPlane_mQx_eta_pos[kMaxEventPlane];   //[EventPlane_]
  Float_t         EventPlane_mQy_eta_pos[kMaxEventPlane];   //[EventPlane_]
  Float_t         EventPlane_mQx_eta_neg[kMaxEventPlane];   //[EventPlane_]
  Float_t         EventPlane_mQy_eta_neg[kMaxEventPlane];   //[EventPlane_]
  Short_t         EventPlane_mNtrk_eta_pos[kMaxEventPlane];   //[EventPlane_]
  Short_t         EventPlane_mNtrk_eta_neg[kMaxEventPlane];   //[EventPlane_]
  Float_t         EventPlane_mWeight_eta_pos[kMaxEventPlane];   //[EventPlane_]
  Float_t         EventPlane_mWeight_eta_neg[kMaxEventPlane];   //[EventPlane_]
  Float_t         EventPlane_mQx_chg_pos[kMaxEventPlane];   //[EventPlane_]
  Float_t         EventPlane_mQy_chg_pos[kMaxEventPlane];   //[EventPlane_]
  Float_t         EventPlane_mQx_chg_neg[kMaxEventPlane];   //[EventPlane_]
  Float_t         EventPlane_mQy_chg_neg[kMaxEventPlane];   //[EventPlane_]
  Short_t         EventPlane_mNtrk_chg_pos_eta_pos[kMaxEventPlane];   //[EventPlane_]
  Short_t         EventPlane_mNtrk_chg_pos_eta_neg[kMaxEventPlane];   //[EventPlane_]
  Short_t         EventPlane_mNtrk_chg_neg_eta_pos[kMaxEventPlane];   //[EventPlane_]
  Short_t         EventPlane_mNtrk_chg_neg_eta_neg[kMaxEventPlane];   //[EventPlane_]
  Float_t         EventPlane_mWeight_chg_pos[kMaxEventPlane];   //[EventPlane_]
  Float_t         EventPlane_mWeight_chg_neg[kMaxEventPlane];   //[EventPlane_]
  Float_t         EventPlane_mQx_ran_1[kMaxEventPlane];   //[EventPlane_]
  Float_t         EventPlane_mQy_ran_1[kMaxEventPlane];   //[EventPlane_]
  Float_t         EventPlane_mQx_ran_2[kMaxEventPlane];   //[EventPlane_]
  Float_t         EventPlane_mQy_ran_2[kMaxEventPlane];   //[EventPlane_]
  Short_t         EventPlane_mNtrk_ran_1_eta_pos[kMaxEventPlane];   //[EventPlane_]
  Short_t         EventPlane_mNtrk_ran_1_eta_neg[kMaxEventPlane];   //[EventPlane_]
  Short_t         EventPlane_mNtrk_ran_2_eta_pos[kMaxEventPlane];   //[EventPlane_]
  Short_t         EventPlane_mNtrk_ran_2_eta_neg[kMaxEventPlane];   //[EventPlane_]
  Float_t         EventPlane_mWeight_ran_1[kMaxEventPlane];   //[EventPlane_]
  Float_t         EventPlane_mWeight_ran_2[kMaxEventPlane];   //[EventPlane_]
  Int_t           Tracks_;
  UShort_t        Tracks_mId[kMaxTracks];   //[Tracks_]
  UShort_t        Tracks_mChi2[kMaxTracks];   //[Tracks_]
  Float_t         Tracks_mPMomentum_mX1[kMaxTracks];   //[Tracks_]
  Float_t         Tracks_mPMomentum_mX2[kMaxTracks];   //[Tracks_]
  Float_t         Tracks_mPMomentum_mX3[kMaxTracks];   //[Tracks_]
  UShort_t        Tracks_mDedx[kMaxTracks];   //[Tracks_]
  Char_t          Tracks_mNHitsFit[kMaxTracks];   //[Tracks_]
  Char_t          Tracks_mNHitsMax[kMaxTracks];   //[Tracks_]
  UChar_t         Tracks_mNHitsDedx[kMaxTracks];   //[Tracks_]
  Short_t         Tracks_mNSigmaPion[kMaxTracks];   //[Tracks_]
  Short_t         Tracks_mNSigmaKaon[kMaxTracks];   //[Tracks_]
  Short_t         Tracks_mNSigmaProton[kMaxTracks];   //[Tracks_]
  Short_t         Tracks_mNSigmaElectron[kMaxTracks];   //[Tracks_]
  UInt_t          Tracks_mMap0[kMaxTracks];   //[Tracks_]
  UInt_t          Tracks_mMap1[kMaxTracks];   //[Tracks_]
  Float_t         Tracks_mPar[kMaxTracks][6];   //[Tracks_]
  Float_t         Tracks_mErrMatrix[kMaxTracks][15];   //[Tracks_]
  Short_t         Tracks_mEmcPidTraitsIndex[kMaxTracks];   //[Tracks_]
  Short_t         Tracks_mBTofPidTraitsIndex[kMaxTracks];   //[Tracks_]
  Short_t         Tracks_mMtdPidTraitsIndex[kMaxTracks];   //[Tracks_]
  Int_t           EmcTrigger_;
  UInt_t          EmcTrigger_fUniqueID[kMaxEmcTrigger];   //[EmcTrigger_]
  UInt_t          EmcTrigger_fBits[kMaxEmcTrigger];   //[EmcTrigger_]
  UChar_t         EmcTrigger_mFlag[kMaxEmcTrigger];   //[EmcTrigger_]
  UShort_t        EmcTrigger_mId[kMaxEmcTrigger];   //[EmcTrigger_]
  UShort_t        EmcTrigger_mAdc[kMaxEmcTrigger];   //[EmcTrigger_]
  Int_t           MtdTrigger_;
  UInt_t          MtdTrigger_fUniqueID[kMaxMtdTrigger];   //[MtdTrigger_]
  UInt_t          MtdTrigger_fBits[kMaxMtdTrigger];   //[MtdTrigger_]
  UShort_t        MtdTrigger_mQTtacSum[kMaxMtdTrigger][4][8];   //[MtdTrigger_]
  UShort_t        MtdTrigger_mMT101Tac[kMaxMtdTrigger][4][2];   //[MtdTrigger_]
  UChar_t         MtdTrigger_mMT101Id[kMaxMtdTrigger][4][2];   //[MtdTrigger_]
  UChar_t         MtdTrigger_mTF201TriggerBit[kMaxMtdTrigger];   //[MtdTrigger_]
  Char_t          MtdTrigger_mShouldHaveRejectEvent[kMaxMtdTrigger];   //[MtdTrigger_]
  Int_t           BTOWHit_;
  UInt_t          BTOWHit_fUniqueID[kMaxBTOWHit];   //[BTOWHit_]
  UInt_t          BTOWHit_fBits[kMaxBTOWHit];   //[BTOWHit_]
  UShort_t        BTOWHit_mId[kMaxBTOWHit];   //[BTOWHit_]
  UShort_t        BTOWHit_mAdc[kMaxBTOWHit];   //[BTOWHit_]
  Short_t         BTOWHit_mE[kMaxBTOWHit];   //[BTOWHit_]
  Int_t           BTofHit_;
  UInt_t          BTofHit_fUniqueID[kMaxBTofHit];   //[BTofHit_]
  UInt_t          BTofHit_fBits[kMaxBTofHit];   //[BTofHit_]
  Short_t         BTofHit_mId[kMaxBTofHit];   //[BTofHit_]
  Int_t           MtdHit_;
  UInt_t          MtdHit_fUniqueID[kMaxMtdHit];   //[MtdHit_]
  UInt_t          MtdHit_fBits[kMaxMtdHit];   //[MtdHit_]
  Short_t         MtdHit_mgChannel[kMaxMtdHit];   //[MtdHit_]
  UChar_t         MtdHit_mTriggerFlag[kMaxMtdHit];   //[MtdHit_]
  Float_t         MtdHit_mLeadingEdgeTime_first[kMaxMtdHit];   //[MtdHit_]
  Float_t         MtdHit_mLeadingEdgeTime_second[kMaxMtdHit];   //[MtdHit_]
  Float_t         MtdHit_mTrailingEdgeTime_first[kMaxMtdHit];   //[MtdHit_]
  Float_t         MtdHit_mTrailingEdgeTime_second[kMaxMtdHit];   //[MtdHit_]
  Int_t           EmcPidTraits_;
  UInt_t          EmcPidTraits_fUniqueID[kMaxEmcPidTraits];   //[EmcPidTraits_]
  UInt_t          EmcPidTraits_fBits[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mTrackIndex[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBEMCId[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBTOWADC0[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBTOWE0[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBTOWE[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBEMCDistZ[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBEMCDistPhi[kMaxEmcPidTraits];   //[EmcPidTraits_]
  UChar_t         EmcPidTraits_mBSMDNEta[kMaxEmcPidTraits];   //[EmcPidTraits_]
  UChar_t         EmcPidTraits_mBSMDNPhi[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBTOWId[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Char_t          EmcPidTraits_mBTOWId23[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBTOWE1[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBTOWE2[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBTOWE3[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBTOWDistEta[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBTOWDistPhi[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Int_t           BTofPidTraits_;
  UInt_t          BTofPidTraits_fUniqueID[kMaxBTofPidTraits];   //[BTofPidTraits_]
  UInt_t          BTofPidTraits_fBits[kMaxBTofPidTraits];   //[BTofPidTraits_]
  Short_t         BTofPidTraits_mTrackIndex[kMaxBTofPidTraits];   //[BTofPidTraits_]
  Short_t         BTofPidTraits_mBTofCellId[kMaxBTofPidTraits];   //[BTofPidTraits_]
  UChar_t         BTofPidTraits_mBTofMatchFlag[kMaxBTofPidTraits];   //[BTofPidTraits_]
  UShort_t        BTofPidTraits_mBTof[kMaxBTofPidTraits];   //[BTofPidTraits_]
  UShort_t        BTofPidTraits_mBTofBeta[kMaxBTofPidTraits];   //[BTofPidTraits_]
  Short_t         BTofPidTraits_mBTofYLocal[kMaxBTofPidTraits];   //[BTofPidTraits_]
  Short_t         BTofPidTraits_mBTofZLocal[kMaxBTofPidTraits];   //[BTofPidTraits_]
  Short_t         BTofPidTraits_mBTofHitPosX[kMaxBTofPidTraits];   //[BTofPidTraits_]
  Short_t         BTofPidTraits_mBTofHitPosY[kMaxBTofPidTraits];   //[BTofPidTraits_]
  Short_t         BTofPidTraits_mBTofHitPosZ[kMaxBTofPidTraits];   //[BTofPidTraits_]
  Int_t           MtdPidTraits_;
  UInt_t          MtdPidTraits_fUniqueID[kMaxMtdPidTraits];   //[MtdPidTraits_]
  UInt_t          MtdPidTraits_fBits[kMaxMtdPidTraits];   //[MtdPidTraits_]
  Short_t         MtdPidTraits_mMtdHitChan[kMaxMtdPidTraits];   //[MtdPidTraits_]
  Short_t         MtdPidTraits_mTrackIndex[kMaxMtdPidTraits];   //[MtdPidTraits_]
  Char_t          MtdPidTraits_mMatchFlag[kMaxMtdPidTraits];   //[MtdPidTraits_]
  Float_t         MtdPidTraits_mDeltaY[kMaxMtdPidTraits];   //[MtdPidTraits_]
  Float_t         MtdPidTraits_mDeltaZ[kMaxMtdPidTraits];   //[MtdPidTraits_]
  Float_t         MtdPidTraits_mDeltaTimeOfFlight[kMaxMtdPidTraits];   //[MtdPidTraits_]
  Float_t         MtdPidTraits_mBeta[kMaxMtdPidTraits];   //[MtdPidTraits_]
  Int_t           V0Ks_;
  Short_t         V0Ks_mIndex2Track[kMaxV0Ks][2];   //[V0Ks_]
  //StThreeVector<float> V0Ks_mMomentum[2][kMaxV0Ks];
  Float_t         V0Ks_mV0Pos_mX1[kMaxV0Ks];   //[V0Ks_]
  Float_t         V0Ks_mV0Pos_mX2[kMaxV0Ks];   //[V0Ks_]
  Float_t         V0Ks_mV0Pos_mX3[kMaxV0Ks];   //[V0Ks_]
  UShort_t        V0Ks_mDcaDaughters[kMaxV0Ks];   //[V0Ks_]
  UShort_t        V0Ks_mDca2Vtx[kMaxV0Ks];   //[V0Ks_]
  Float_t         V0Ks_mM[kMaxV0Ks];   //[V0Ks_]
  Int_t           V0L_;
  Short_t         V0L_mIndex2Track[kMaxV0L][2];   //[V0L_]
  //StThreeVector<float> V0L_mMomentum[2][kMaxV0L];
  Float_t         V0L_mV0Pos_mX1[kMaxV0L];   //[V0L_]
  Float_t         V0L_mV0Pos_mX2[kMaxV0L];   //[V0L_]
  Float_t         V0L_mV0Pos_mX3[kMaxV0L];   //[V0L_]
  UShort_t        V0L_mDcaDaughters[kMaxV0L];   //[V0L_]
  UShort_t        V0L_mDca2Vtx[kMaxV0L];   //[V0L_]
  Float_t         V0L_mM[kMaxV0L];   //[V0L_]
  Int_t           V0Lbar_;
  Short_t         V0Lbar_mIndex2Track[kMaxV0Lbar][2];   //[V0Lbar_]
  //StThreeVector<float> V0Lbar_mMomentum[2][kMaxV0Lbar];
  Float_t         V0Lbar_mV0Pos_mX1[kMaxV0Lbar];   //[V0Lbar_]
  Float_t         V0Lbar_mV0Pos_mX2[kMaxV0Lbar];   //[V0Lbar_]
  Float_t         V0Lbar_mV0Pos_mX3[kMaxV0Lbar];   //[V0Lbar_]
  UShort_t        V0Lbar_mDcaDaughters[kMaxV0Lbar];   //[V0Lbar_]
  UShort_t        V0Lbar_mDca2Vtx[kMaxV0Lbar];   //[V0Lbar_]
  Float_t         V0Lbar_mM[kMaxV0Lbar];   //[V0Lbar_]
  
  // List of branches
  TBranch        *b_Event_;   //!
  TBranch        *b_Event_mRunId;   //!
  TBranch        *b_Event_mEventId;   //!
  TBranch        *b_Event_mFillId;   //!
  TBranch        *b_Event_mBField;   //!
  TBranch        *b_Event_mPrimaryVertex_mX1;   //!
  TBranch        *b_Event_mPrimaryVertex_mX2;   //!
  TBranch        *b_Event_mPrimaryVertex_mX3;   //!
  TBranch        *b_Event_mSecondVertex_mX1;   //!
  TBranch        *b_Event_mSecondVertex_mX2;   //!
  TBranch        *b_Event_mSecondVertex_mX3;   //!
  TBranch        *b_Event_mPrimaryVertexError_mX1;   //!
  TBranch        *b_Event_mPrimaryVertexError_mX2;   //!
  TBranch        *b_Event_mPrimaryVertexError_mX3;   //!
  TBranch        *b_Event_mTriggerWord;   //!
  TBranch        *b_Event_mTriggerWordMtd;   //!
  TBranch        *b_Event_mRefMultFtpcEast;   //!
  TBranch        *b_Event_mRefMultFtpcWest;   //!
  TBranch        *b_Event_mRefMultNeg;   //!
  TBranch        *b_Event_mRefMultPos;   //!
  TBranch        *b_Event_mRefMult2NegEast;   //!
  TBranch        *b_Event_mRefMult2PosEast;   //!
  TBranch        *b_Event_mRefMult2NegWest;   //!
  TBranch        *b_Event_mRefMult2PosWest;   //!
  TBranch        *b_Event_mRefMultHalfNegEast;   //!
  TBranch        *b_Event_mRefMultHalfPosEast;   //!
  TBranch        *b_Event_mRefMultHalfNegWest;   //!
  TBranch        *b_Event_mRefMultHalfPosWest;   //!
  TBranch        *b_Event_mGRefMult;   //!
  TBranch        *b_Event_mNHitsHFT;   //!
  TBranch        *b_Event_mNVpdHitsEast;   //!
  TBranch        *b_Event_mNVpdHitsWest;   //!
  TBranch        *b_Event_mNT0;   //!
  TBranch        *b_Event_mVzVpd;   //!
  TBranch        *b_Event_mZDCx;   //!
  TBranch        *b_Event_mBBCx;   //!
  TBranch        *b_Event_mBackgroundRate;   //!
  TBranch        *b_Event_mBbcBlueBackgroundRate;   //!
  TBranch        *b_Event_mBbcYellowBackgroundRate;   //!
  TBranch        *b_Event_mBbcEastRate;   //!
  TBranch        *b_Event_mBbcWestRate;   //!
  TBranch        *b_Event_mZdcEastRate;   //!
  TBranch        *b_Event_mZdcWestRate;   //!
  TBranch        *b_Event_mVpd;   //!
  TBranch        *b_Event_mZdcSumAdcEast;   //!
  TBranch        *b_Event_mZdcSumAdcWest;   //!
  TBranch        *b_Event_mZdcSmdEastHorizontal;   //!
  TBranch        *b_Event_mZdcSmdEastVertical;   //!
  TBranch        *b_Event_mZdcSmdWestHorizontal;   //!
  TBranch        *b_Event_mZdcSmdWestVertical;   //!
  TBranch        *b_Event_mSpaceCharge;   //!
  TBranch        *b_Event_mbTofTrayMultiplicity;   //!
  TBranch        *b_Event_mNumberOfGlobalTracks;   //!
  TBranch        *b_Event_mRanking;   //!
  TBranch        *b_Event_mNBEMCMatch;   //!
  TBranch        *b_Event_mNBTOFMatch;   //!
  TBranch        *b_Event_mBbcAdcEast;   //!
  TBranch        *b_Event_mBbcAdcWest;   //!
  TBranch        *b_Event_mHT_Th;   //!
  TBranch        *b_Event_mJP_Th;   //!
  TBranch        *b_EventPlane_;   //!
  TBranch        *b_EventPlane_fUniqueID;   //!
  TBranch        *b_EventPlane_fBits;   //!
  TBranch        *b_EventPlane_mQx_eta_pos;   //!
  TBranch        *b_EventPlane_mQy_eta_pos;   //!
  TBranch        *b_EventPlane_mQx_eta_neg;   //!
  TBranch        *b_EventPlane_mQy_eta_neg;   //!
  TBranch        *b_EventPlane_mNtrk_eta_pos;   //!
  TBranch        *b_EventPlane_mNtrk_eta_neg;   //!
  TBranch        *b_EventPlane_mWeight_eta_pos;   //!
  TBranch        *b_EventPlane_mWeight_eta_neg;   //!
  TBranch        *b_EventPlane_mQx_chg_pos;   //!
  TBranch        *b_EventPlane_mQy_chg_pos;   //!
  TBranch        *b_EventPlane_mQx_chg_neg;   //!
  TBranch        *b_EventPlane_mQy_chg_neg;   //!
  TBranch        *b_EventPlane_mNtrk_chg_pos_eta_pos;   //!
  TBranch        *b_EventPlane_mNtrk_chg_pos_eta_neg;   //!
  TBranch        *b_EventPlane_mNtrk_chg_neg_eta_pos;   //!
  TBranch        *b_EventPlane_mNtrk_chg_neg_eta_neg;   //!
  TBranch        *b_EventPlane_mWeight_chg_pos;   //!
  TBranch        *b_EventPlane_mWeight_chg_neg;   //!
  TBranch        *b_EventPlane_mQx_ran_1;   //!
  TBranch        *b_EventPlane_mQy_ran_1;   //!
  TBranch        *b_EventPlane_mQx_ran_2;   //!
  TBranch        *b_EventPlane_mQy_ran_2;   //!
  TBranch        *b_EventPlane_mNtrk_ran_1_eta_pos;   //!
  TBranch        *b_EventPlane_mNtrk_ran_1_eta_neg;   //!
  TBranch        *b_EventPlane_mNtrk_ran_2_eta_pos;   //!
  TBranch        *b_EventPlane_mNtrk_ran_2_eta_neg;   //!
  TBranch        *b_EventPlane_mWeight_ran_1;   //!
  TBranch        *b_EventPlane_mWeight_ran_2;   //!
  TBranch        *b_Tracks_;   //!
  TBranch        *b_Tracks_mId;   //!
  TBranch        *b_Tracks_mChi2;   //!
  TBranch        *b_Tracks_mPMomentum_mX1;   //!
  TBranch        *b_Tracks_mPMomentum_mX2;   //!
  TBranch        *b_Tracks_mPMomentum_mX3;   //!
  TBranch        *b_Tracks_mDedx;   //!
  TBranch        *b_Tracks_mNHitsFit;   //!
  TBranch        *b_Tracks_mNHitsMax;   //!
  TBranch        *b_Tracks_mNHitsDedx;   //!
  TBranch        *b_Tracks_mNSigmaPion;   //!
  TBranch        *b_Tracks_mNSigmaKaon;   //!
  TBranch        *b_Tracks_mNSigmaProton;   //!
  TBranch        *b_Tracks_mNSigmaElectron;   //!
  TBranch        *b_Tracks_mMap0;   //!
  TBranch        *b_Tracks_mMap1;   //!
  TBranch        *b_Tracks_mPar;   //!
  TBranch        *b_Tracks_mErrMatrix;   //!
  TBranch        *b_Tracks_mEmcPidTraitsIndex;   //!
  TBranch        *b_Tracks_mBTofPidTraitsIndex;   //!
  TBranch        *b_Tracks_mMtdPidTraitsIndex;   //!
  TBranch        *b_EmcTrigger_;   //!
  TBranch        *b_EmcTrigger_fUniqueID;   //!
  TBranch        *b_EmcTrigger_fBits;   //!
  TBranch        *b_EmcTrigger_mFlag;   //!
  TBranch        *b_EmcTrigger_mId;   //!
  TBranch        *b_EmcTrigger_mAdc;   //!
  TBranch        *b_MtdTrigger_;   //!
  TBranch        *b_MtdTrigger_fUniqueID;   //!
  TBranch        *b_MtdTrigger_fBits;   //!
  TBranch        *b_MtdTrigger_mQTtacSum;   //!
  TBranch        *b_MtdTrigger_mMT101Tac;   //!
  TBranch        *b_MtdTrigger_mMT101Id;   //!
  TBranch        *b_MtdTrigger_mTF201TriggerBit;   //!
  TBranch        *b_MtdTrigger_mShouldHaveRejectEvent;   //!
  TBranch        *b_BTOWHit_;   //!
  TBranch        *b_BTOWHit_fUniqueID;   //!
  TBranch        *b_BTOWHit_fBits;   //!
  TBranch        *b_BTOWHit_mId;   //!
  TBranch        *b_BTOWHit_mAdc;   //!
  TBranch        *b_BTOWHit_mE;   //!
  TBranch        *b_BTofHit_;   //!
  TBranch        *b_BTofHit_fUniqueID;   //!
  TBranch        *b_BTofHit_fBits;   //!
  TBranch        *b_BTofHit_mId;   //!
  TBranch        *b_MtdHit_;   //!
  TBranch        *b_MtdHit_fUniqueID;   //!
  TBranch        *b_MtdHit_fBits;   //!
  TBranch        *b_MtdHit_mgChannel;   //!
  TBranch        *b_MtdHit_mTriggerFlag;   //!
  TBranch        *b_MtdHit_mLeadingEdgeTime_first;   //!
  TBranch        *b_MtdHit_mLeadingEdgeTime_second;   //!
  TBranch        *b_MtdHit_mTrailingEdgeTime_first;   //!
  TBranch        *b_MtdHit_mTrailingEdgeTime_second;   //!
  TBranch        *b_EmcPidTraits_;   //!
  TBranch        *b_EmcPidTraits_fUniqueID;   //!
  TBranch        *b_EmcPidTraits_fBits;   //!
  TBranch        *b_EmcPidTraits_mTrackIndex;   //!
  TBranch        *b_EmcPidTraits_mBEMCId;   //!
  TBranch        *b_EmcPidTraits_mBTOWADC0;   //!
  TBranch        *b_EmcPidTraits_mBTOWE0;   //!
  TBranch        *b_EmcPidTraits_mBTOWE;   //!
  TBranch        *b_EmcPidTraits_mBEMCDistZ;   //!
  TBranch        *b_EmcPidTraits_mBEMCDistPhi;   //!
  TBranch        *b_EmcPidTraits_mBSMDNEta;   //!
  TBranch        *b_EmcPidTraits_mBSMDNPhi;   //!
  TBranch        *b_EmcPidTraits_mBTOWId;   //!
  TBranch        *b_EmcPidTraits_mBTOWId23;   //!
  TBranch        *b_EmcPidTraits_mBTOWE1;   //!
  TBranch        *b_EmcPidTraits_mBTOWE2;   //!
  TBranch        *b_EmcPidTraits_mBTOWE3;   //!
  TBranch        *b_EmcPidTraits_mBTOWDistEta;   //!
  TBranch        *b_EmcPidTraits_mBTOWDistPhi;   //!
  TBranch        *b_BTofPidTraits_;   //!
  TBranch        *b_BTofPidTraits_fUniqueID;   //!
  TBranch        *b_BTofPidTraits_fBits;   //!
  TBranch        *b_BTofPidTraits_mTrackIndex;   //!
  TBranch        *b_BTofPidTraits_mBTofCellId;   //!
  TBranch        *b_BTofPidTraits_mBTofMatchFlag;   //!
  TBranch        *b_BTofPidTraits_mBTof;   //!
  TBranch        *b_BTofPidTraits_mBTofBeta;   //!
  TBranch        *b_BTofPidTraits_mBTofYLocal;   //!
  TBranch        *b_BTofPidTraits_mBTofZLocal;   //!
  TBranch        *b_BTofPidTraits_mBTofHitPosX;   //!
  TBranch        *b_BTofPidTraits_mBTofHitPosY;   //!
  TBranch        *b_BTofPidTraits_mBTofHitPosZ;   //!
  TBranch        *b_MtdPidTraits_;   //!
  TBranch        *b_MtdPidTraits_fUniqueID;   //!
  TBranch        *b_MtdPidTraits_fBits;   //!
  TBranch        *b_MtdPidTraits_mMtdHitChan;   //!
  TBranch        *b_MtdPidTraits_mTrackIndex;   //!
  TBranch        *b_MtdPidTraits_mMatchFlag;   //!
  TBranch        *b_MtdPidTraits_mDeltaY;   //!
  TBranch        *b_MtdPidTraits_mDeltaZ;   //!
  TBranch        *b_MtdPidTraits_mDeltaTimeOfFlight;   //!
  TBranch        *b_MtdPidTraits_mBeta;   //!
  TBranch        *b_V0Ks_;   //!
  TBranch        *b_V0Ks_mIndex2Track;   //!
  TBranch        *b_V0Ks_mV0Pos_mX1;   //!
  TBranch        *b_V0Ks_mV0Pos_mX2;   //!
  TBranch        *b_V0Ks_mV0Pos_mX3;   //!
  TBranch        *b_V0Ks_mDcaDaughters;   //!
  TBranch        *b_V0Ks_mDca2Vtx;   //!
  TBranch        *b_V0Ks_mM;   //!
  TBranch        *b_V0L_;   //!
  TBranch        *b_V0L_mIndex2Track;   //!
  TBranch        *b_V0L_mV0Pos_mX1;   //!
  TBranch        *b_V0L_mV0Pos_mX2;   //!
  TBranch        *b_V0L_mV0Pos_mX3;   //!
  TBranch        *b_V0L_mDcaDaughters;   //!
  TBranch        *b_V0L_mDca2Vtx;   //!
  TBranch        *b_V0L_mM;   //!
  TBranch        *b_V0Lbar_;   //!
  TBranch        *b_V0Lbar_mIndex2Track;   //!
  TBranch        *b_V0Lbar_mV0Pos_mX1;   //!
  TBranch        *b_V0Lbar_mV0Pos_mX2;   //!
  TBranch        *b_V0Lbar_mV0Pos_mX3;   //!
  TBranch        *b_V0Lbar_mDcaDaughters;   //!
  TBranch        *b_V0Lbar_mDca2Vtx;   //!
  TBranch        *b_V0Lbar_mM;   //!
  
  PicoDst(TTree *tree=0, Char_t *output="");
  virtual ~PicoDst();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree, Char_t *outputFile);
  virtual void     Loop( int nEvents);
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PicoDst_cxx

#endif // #ifdef PicoDst_cxx

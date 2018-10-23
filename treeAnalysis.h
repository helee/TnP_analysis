//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 16 11:39:39 2018 by ROOT version 6.06/00
// from TTree fitter_tree/fitter_tree
// found on file: TnPTree_data_63.root
//////////////////////////////////////////////////////////

#ifndef treeAnalysis_h
#define treeAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class treeAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TString         inputfile;
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         el_1overEminus1overP;
   Float_t         el_3charge;
   Float_t         el_5x5_e1x5;
   Float_t         el_5x5_e2x5;
   Float_t         el_5x5_e5x5;
   Float_t         el_5x5_r9;
   Float_t         el_5x5_sieie;
   Float_t         el_5x5_sieip;
   Float_t         el_abseta;
   Float_t         el_chIso;
   Float_t         el_dEtaIn;
   Float_t         el_dEtaSeed;
   Float_t         el_dPhiIn;
   Float_t         el_dr03TkSumPt;
   Float_t         el_e;
   Float_t         el_e1x5;
   Float_t         el_e2x5;
   Float_t         el_e5x5;
   Float_t         el_ecalDriven;
   Float_t         el_ecalEnergy;
   Float_t         el_ecalIso;
   Float_t         el_eoverp_wES;
   Float_t         el_et;
   Float_t         el_eta;
   Float_t         el_etaW;
   Float_t         el_fbrem;
   Float_t         el_hcalIso;
   Float_t         el_hoe;
   Float_t         el_hoe_bc;
   Float_t         el_neuIso;
   Float_t         el_phi;
   Float_t         el_phiW;
   Float_t         el_phoIso;
   Float_t         el_pt;
   Float_t         el_q;
   Float_t         el_r9;
   Float_t         el_reliso03;
   Float_t         el_sc_abseta;
   Float_t         el_sc_e;
   Float_t         el_sc_esE;
   Float_t         el_sc_et;
   Float_t         el_sc_eta;
   Float_t         el_sc_phi;
   Float_t         el_sc_rawE;
   Float_t         el_seed_e;
   Float_t         el_sieie;
   Float_t         el_sumPUPt;
   Float_t         el_tk_eta;
   Float_t         el_tk_phi;
   Float_t         el_tk_pt;
   Float_t         el_trkIso;
   Float_t         el_IsoMVA94X;
   Float_t         el_chisq;
   Float_t         el_dxy;
   Float_t         el_dz;
   Float_t         el_hoe_barrel;
   Float_t         el_hoe_endcap;
   Float_t         el_hzzMVA80X;
   Float_t         el_mHits;
   Float_t         el_miniIso;
   Float_t         el_noIsoMVA94X;
   Float_t         el_nonTrigMVA80X;
   Float_t         el_passConversionVeto;
   Float_t         el_pf_pt;
   Float_t         el_sip3d;
   Int_t           passingHLTsafe;
   Int_t           passingLoose80X;
   Int_t           passingLoose94X;
   Int_t           passingMVA80Xwp80;
   Int_t           passingMVA80Xwp90;
   Int_t           passingMVA94Xwp80iso;
   Int_t           passingMVA94Xwp80noiso;
   Int_t           passingMVA94Xwp90iso;
   Int_t           passingMVA94Xwp90noiso;
   Int_t           passingMVA94XwpLiso;
   Int_t           passingMVA94XwpLnoiso;
   Int_t           passingMedium80X;
   Int_t           passingMedium94X;
   Int_t           passingTight80X;
   Int_t           passingTight94X;
   Int_t           passingVeto80X;
   Int_t           passingVeto94X;
   UInt_t          run;
   UInt_t          lumi;
   ULong64_t       event;
   Int_t           truePU;
   Int_t           event_nPV;
   Float_t         event_met_pfmet;
   Float_t         event_met_pfphi;
   Float_t         event_met_pfsumet;
   Float_t         event_met_pfmetsignificance;
   Float_t         event_PrimaryVertex_x;
   Float_t         event_PrimaryVertex_y;
   Float_t         event_PrimaryVertex_z;
   Float_t         event_BeamSpot_x;
   Float_t         event_BeamSpot_y;
   Float_t         event_BeamSpot_z;
   Float_t         event_rho;
   Float_t         mass;
   Float_t         tag_Ele_1overEminus1overP;
   Float_t         tag_Ele_5x5_sieie;
   Float_t         tag_Ele_abseta;
   Float_t         tag_Ele_dEtaSeed;
   Float_t         tag_Ele_dPhiIn;
   Float_t         tag_Ele_e;
   Float_t         tag_Ele_et;
   Float_t         tag_Ele_eta;
   Float_t         tag_Ele_hoe;
   Float_t         tag_Ele_phi;
   Float_t         tag_Ele_pt;
   Float_t         tag_Ele_q;
   Float_t         tag_sc_abseta;
   Float_t         tag_sc_e;
   Float_t         tag_sc_et;
   Float_t         tag_sc_eta;
   Float_t         tag_Ele_IsoMVA94X;
   Float_t         tag_Ele_dxy;
   Float_t         tag_Ele_dz;
   Float_t         tag_Ele_hoe_barrel;
   Float_t         tag_Ele_hoe_endcap;
   Float_t         tag_Ele_hzzMVA80X;
   Float_t         tag_Ele_mHits;
   Float_t         tag_Ele_miniIso;
   Float_t         tag_Ele_noIsoMVA94X;
   Float_t         tag_Ele_nonTrigMVA80X;
   Float_t         tag_Ele_passConversionVeto;
   Float_t         tag_Ele_sip3d;
   Float_t         tag_Ele_trigMVA;
   Float_t         mc_probe_e;
   Float_t         mc_probe_et;
   Float_t         mc_probe_eta;
   Float_t         mc_probe_phi;
   Int_t           mc_probe_flag;
   Float_t         pair_abseta;
   Float_t         pair_eta;
   Float_t         pair_mass;
   Float_t         pair_pt;

   // List of branches
   TBranch        *b_el_1overEminus1overP;   //!
   TBranch        *b_el_3charge;   //!
   TBranch        *b_el_5x5_e1x5;   //!
   TBranch        *b_el_5x5_e2x5;   //!
   TBranch        *b_el_5x5_e5x5;   //!
   TBranch        *b_el_5x5_r9;   //!
   TBranch        *b_el_5x5_sieie;   //!
   TBranch        *b_el_5x5_sieip;   //!
   TBranch        *b_el_abseta;   //!
   TBranch        *b_el_chIso;   //!
   TBranch        *b_el_dEtaIn;   //!
   TBranch        *b_el_dEtaSeed;   //!
   TBranch        *b_el_dPhiIn;   //!
   TBranch        *b_el_dr03TkSumPt;   //!
   TBranch        *b_el_e;   //!
   TBranch        *b_el_e1x5;   //!
   TBranch        *b_el_e2x5;   //!
   TBranch        *b_el_e5x5;   //!
   TBranch        *b_el_ecalDriven;   //!
   TBranch        *b_el_ecalEnergy;   //!
   TBranch        *b_el_ecalIso;   //!
   TBranch        *b_el_eoverp_wES;   //!
   TBranch        *b_el_et;   //!
   TBranch        *b_el_eta;   //!
   TBranch        *b_el_etaW;   //!
   TBranch        *b_el_fbrem;   //!
   TBranch        *b_el_hcalIso;   //!
   TBranch        *b_el_hoe;   //!
   TBranch        *b_el_hoe_bc;   //!
   TBranch        *b_el_neuIso;   //!
   TBranch        *b_el_phi;   //!
   TBranch        *b_el_phiW;   //!
   TBranch        *b_el_phoIso;   //!
   TBranch        *b_el_pt;   //!
   TBranch        *b_el_q;   //!
   TBranch        *b_el_r9;   //!
   TBranch        *b_el_reliso03;   //!
   TBranch        *b_el_sc_abseta;   //!
   TBranch        *b_el_sc_e;   //!
   TBranch        *b_el_sc_esE;   //!
   TBranch        *b_el_sc_et;   //!
   TBranch        *b_el_sc_eta;   //!
   TBranch        *b_el_sc_phi;   //!
   TBranch        *b_el_sc_rawE;   //!
   TBranch        *b_el_seed_e;   //!
   TBranch        *b_el_sieie;   //!
   TBranch        *b_el_sumPUPt;   //!
   TBranch        *b_el_tk_eta;   //!
   TBranch        *b_el_tk_phi;   //!
   TBranch        *b_el_tk_pt;   //!
   TBranch        *b_el_trkIso;   //!
   TBranch        *b_el_IsoMVA94X;   //!
   TBranch        *b_el_chisq;   //!
   TBranch        *b_el_dxy;   //!
   TBranch        *b_el_dz;   //!
   TBranch        *b_el_hoe_barrel;   //!
   TBranch        *b_el_hoe_endcap;   //!
   TBranch        *b_el_hzzMVA80X;   //!
   TBranch        *b_el_mHits;   //!
   TBranch        *b_el_miniIso;   //!
   TBranch        *b_el_noIsoMVA94X;   //!
   TBranch        *b_el_nonTrigMVA80X;   //!
   TBranch        *b_el_passConversionVeto;   //!
   TBranch        *b_el_pf_pt;   //!
   TBranch        *b_el_sip3d;   //!
   TBranch        *b_passingHLTsafe;   //!
   TBranch        *b_passingLoose80X;   //!
   TBranch        *b_passingLoose94X;   //!
   TBranch        *b_passingMVA80Xwp80;   //!
   TBranch        *b_passingMVA80Xwp90;   //!
   TBranch        *b_passingMVA94Xwp80iso;   //!
   TBranch        *b_passingMVA94Xwp80noiso;   //!
   TBranch        *b_passingMVA94Xwp90iso;   //!
   TBranch        *b_passingMVA94Xwp90noiso;   //!
   TBranch        *b_passingMVA94XwpLiso;   //!
   TBranch        *b_passingMVA94XwpLnoiso;   //!
   TBranch        *b_passingMedium80X;   //!
   TBranch        *b_passingMedium94X;   //!
   TBranch        *b_passingTight80X;   //!
   TBranch        *b_passingTight94X;   //!
   TBranch        *b_passingVeto80X;   //!
   TBranch        *b_passingVeto94X;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_truePU;   //!
   TBranch        *b_mNPV;   //!
   TBranch        *b_mpfMET;   //!
   TBranch        *b_mpfPhi;   //!
   TBranch        *b_mpfSumET;   //!
   TBranch        *b_mpfMETSign;   //!
   TBranch        *b_mPVx;   //!
   TBranch        *b_mPVy;   //!
   TBranch        *b_mPVz;   //!
   TBranch        *b_mBSx;   //!
   TBranch        *b_mBSy;   //!
   TBranch        *b_mBSz;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_mass;   //!
   TBranch        *b_tag_Ele_1overEminus1overP;   //!
   TBranch        *b_tag_Ele_5x5_sieie;   //!
   TBranch        *b_tag_Ele_abseta;   //!
   TBranch        *b_tag_Ele_dEtaSeed;   //!
   TBranch        *b_tag_Ele_dPhiIn;   //!
   TBranch        *b_tag_Ele_e;   //!
   TBranch        *b_tag_Ele_et;   //!
   TBranch        *b_tag_Ele_eta;   //!
   TBranch        *b_tag_Ele_hoe;   //!
   TBranch        *b_tag_Ele_phi;   //!
   TBranch        *b_tag_Ele_pt;   //!
   TBranch        *b_tag_Ele_q;   //!
   TBranch        *b_tag_sc_abseta;   //!
   TBranch        *b_tag_sc_e;   //!
   TBranch        *b_tag_sc_et;   //!
   TBranch        *b_tag_sc_eta;   //!
   TBranch        *b_tag_Ele_IsoMVA94X;   //!
   TBranch        *b_tag_Ele_dxy;   //!
   TBranch        *b_tag_Ele_dz;   //!
   TBranch        *b_tag_Ele_hoe_barrel;   //!
   TBranch        *b_tag_Ele_hoe_endcap;   //!
   TBranch        *b_tag_Ele_hzzMVA80X;   //!
   TBranch        *b_tag_Ele_mHits;   //!
   TBranch        *b_tag_Ele_miniIso;   //!
   TBranch        *b_tag_Ele_noIsoMVA94X;   //!
   TBranch        *b_tag_Ele_nonTrigMVA80X;   //!
   TBranch        *b_tag_Ele_passConversionVeto;   //!
   TBranch        *b_tag_Ele_sip3d;   //!
   TBranch        *b_tag_Ele_trigMVA;   //!
   TBranch        *b_mc_probe_e;   //!
   TBranch        *b_mc_probe_et;   //!
   TBranch        *b_mc_probe_eta;   //!
   TBranch        *b_mc_probe_phi;   //!
   TBranch        *b_mc_probe_flag;   //!
   TBranch        *b_pair_abseta;   //!
   TBranch        *b_pair_eta;   //!
   TBranch        *b_pair_mass;   //!
   TBranch        *b_pair_pt;   //!

//   treeAnalysis(TTree *tree=0);
   treeAnalysis(TString filename = " ", TTree *tree=0);
   virtual ~treeAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef treeAnalysis_cxx
//treeAnalysis::treeAnalysis(TTree *tree) : fChain(0) 
treeAnalysis::treeAnalysis(TString filename, TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
//      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("TnPTree_data_63.root");
      inputfile = filename;
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
      if (!f || !f->IsOpen()) {
//         f = new TFile("TnPTree_data_63.root");
         f = new TFile(filename);
      }
//      TDirectory * dir = (TDirectory*)f->Get("TnPTree_data_63.root:/tnpEleIDs");
      TDirectory * dir = (TDirectory*)f->Get(inputfile+":/tnpEleIDs");
      dir->GetObject("fitter_tree",tree);
   }
   Init(tree);
}

treeAnalysis::~treeAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t treeAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t treeAnalysis::LoadTree(Long64_t entry)
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

void treeAnalysis::Init(TTree *tree)
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

   fChain->SetBranchAddress("el_1overEminus1overP", &el_1overEminus1overP, &b_el_1overEminus1overP);
   fChain->SetBranchAddress("el_3charge", &el_3charge, &b_el_3charge);
   fChain->SetBranchAddress("el_5x5_e1x5", &el_5x5_e1x5, &b_el_5x5_e1x5);
   fChain->SetBranchAddress("el_5x5_e2x5", &el_5x5_e2x5, &b_el_5x5_e2x5);
   fChain->SetBranchAddress("el_5x5_e5x5", &el_5x5_e5x5, &b_el_5x5_e5x5);
   fChain->SetBranchAddress("el_5x5_r9", &el_5x5_r9, &b_el_5x5_r9);
   fChain->SetBranchAddress("el_5x5_sieie", &el_5x5_sieie, &b_el_5x5_sieie);
   fChain->SetBranchAddress("el_5x5_sieip", &el_5x5_sieip, &b_el_5x5_sieip);
   fChain->SetBranchAddress("el_abseta", &el_abseta, &b_el_abseta);
   fChain->SetBranchAddress("el_chIso", &el_chIso, &b_el_chIso);
   fChain->SetBranchAddress("el_dEtaIn", &el_dEtaIn, &b_el_dEtaIn);
   fChain->SetBranchAddress("el_dEtaSeed", &el_dEtaSeed, &b_el_dEtaSeed);
   fChain->SetBranchAddress("el_dPhiIn", &el_dPhiIn, &b_el_dPhiIn);
   fChain->SetBranchAddress("el_dr03TkSumPt", &el_dr03TkSumPt, &b_el_dr03TkSumPt);
   fChain->SetBranchAddress("el_e", &el_e, &b_el_e);
   fChain->SetBranchAddress("el_e1x5", &el_e1x5, &b_el_e1x5);
   fChain->SetBranchAddress("el_e2x5", &el_e2x5, &b_el_e2x5);
   fChain->SetBranchAddress("el_e5x5", &el_e5x5, &b_el_e5x5);
   fChain->SetBranchAddress("el_ecalDriven", &el_ecalDriven, &b_el_ecalDriven);
   fChain->SetBranchAddress("el_ecalEnergy", &el_ecalEnergy, &b_el_ecalEnergy);
   fChain->SetBranchAddress("el_ecalIso", &el_ecalIso, &b_el_ecalIso);
   fChain->SetBranchAddress("el_eoverp_wES", &el_eoverp_wES, &b_el_eoverp_wES);
   fChain->SetBranchAddress("el_et", &el_et, &b_el_et);
   fChain->SetBranchAddress("el_eta", &el_eta, &b_el_eta);
   fChain->SetBranchAddress("el_etaW", &el_etaW, &b_el_etaW);
   fChain->SetBranchAddress("el_fbrem", &el_fbrem, &b_el_fbrem);
   fChain->SetBranchAddress("el_hcalIso", &el_hcalIso, &b_el_hcalIso);
   fChain->SetBranchAddress("el_hoe", &el_hoe, &b_el_hoe);
   fChain->SetBranchAddress("el_hoe_bc", &el_hoe_bc, &b_el_hoe_bc);
   fChain->SetBranchAddress("el_neuIso", &el_neuIso, &b_el_neuIso);
   fChain->SetBranchAddress("el_phi", &el_phi, &b_el_phi);
   fChain->SetBranchAddress("el_phiW", &el_phiW, &b_el_phiW);
   fChain->SetBranchAddress("el_phoIso", &el_phoIso, &b_el_phoIso);
   fChain->SetBranchAddress("el_pt", &el_pt, &b_el_pt);
   fChain->SetBranchAddress("el_q", &el_q, &b_el_q);
   fChain->SetBranchAddress("el_r9", &el_r9, &b_el_r9);
   fChain->SetBranchAddress("el_reliso03", &el_reliso03, &b_el_reliso03);
   fChain->SetBranchAddress("el_sc_abseta", &el_sc_abseta, &b_el_sc_abseta);
   fChain->SetBranchAddress("el_sc_e", &el_sc_e, &b_el_sc_e);
   fChain->SetBranchAddress("el_sc_esE", &el_sc_esE, &b_el_sc_esE);
   fChain->SetBranchAddress("el_sc_et", &el_sc_et, &b_el_sc_et);
   fChain->SetBranchAddress("el_sc_eta", &el_sc_eta, &b_el_sc_eta);
   fChain->SetBranchAddress("el_sc_phi", &el_sc_phi, &b_el_sc_phi);
   fChain->SetBranchAddress("el_sc_rawE", &el_sc_rawE, &b_el_sc_rawE);
   fChain->SetBranchAddress("el_seed_e", &el_seed_e, &b_el_seed_e);
   fChain->SetBranchAddress("el_sieie", &el_sieie, &b_el_sieie);
   fChain->SetBranchAddress("el_sumPUPt", &el_sumPUPt, &b_el_sumPUPt);
   fChain->SetBranchAddress("el_tk_eta", &el_tk_eta, &b_el_tk_eta);
   fChain->SetBranchAddress("el_tk_phi", &el_tk_phi, &b_el_tk_phi);
   fChain->SetBranchAddress("el_tk_pt", &el_tk_pt, &b_el_tk_pt);
   fChain->SetBranchAddress("el_trkIso", &el_trkIso, &b_el_trkIso);
   fChain->SetBranchAddress("el_IsoMVA94X", &el_IsoMVA94X, &b_el_IsoMVA94X);
   fChain->SetBranchAddress("el_chisq", &el_chisq, &b_el_chisq);
   fChain->SetBranchAddress("el_dxy", &el_dxy, &b_el_dxy);
   fChain->SetBranchAddress("el_dz", &el_dz, &b_el_dz);
   fChain->SetBranchAddress("el_hoe_barrel", &el_hoe_barrel, &b_el_hoe_barrel);
   fChain->SetBranchAddress("el_hoe_endcap", &el_hoe_endcap, &b_el_hoe_endcap);
   fChain->SetBranchAddress("el_hzzMVA80X", &el_hzzMVA80X, &b_el_hzzMVA80X);
   fChain->SetBranchAddress("el_mHits", &el_mHits, &b_el_mHits);
   fChain->SetBranchAddress("el_miniIso", &el_miniIso, &b_el_miniIso);
   fChain->SetBranchAddress("el_noIsoMVA94X", &el_noIsoMVA94X, &b_el_noIsoMVA94X);
   fChain->SetBranchAddress("el_nonTrigMVA80X", &el_nonTrigMVA80X, &b_el_nonTrigMVA80X);
   fChain->SetBranchAddress("el_passConversionVeto", &el_passConversionVeto, &b_el_passConversionVeto);
   fChain->SetBranchAddress("el_pf_pt", &el_pf_pt, &b_el_pf_pt);
   fChain->SetBranchAddress("el_sip3d", &el_sip3d, &b_el_sip3d);
   fChain->SetBranchAddress("passingHLTsafe", &passingHLTsafe, &b_passingHLTsafe);
   fChain->SetBranchAddress("passingLoose80X", &passingLoose80X, &b_passingLoose80X);
   fChain->SetBranchAddress("passingLoose94X", &passingLoose94X, &b_passingLoose94X);
   fChain->SetBranchAddress("passingMVA80Xwp80", &passingMVA80Xwp80, &b_passingMVA80Xwp80);
   fChain->SetBranchAddress("passingMVA80Xwp90", &passingMVA80Xwp90, &b_passingMVA80Xwp90);
   fChain->SetBranchAddress("passingMVA94Xwp80iso", &passingMVA94Xwp80iso, &b_passingMVA94Xwp80iso);
   fChain->SetBranchAddress("passingMVA94Xwp80noiso", &passingMVA94Xwp80noiso, &b_passingMVA94Xwp80noiso);
   fChain->SetBranchAddress("passingMVA94Xwp90iso", &passingMVA94Xwp90iso, &b_passingMVA94Xwp90iso);
   fChain->SetBranchAddress("passingMVA94Xwp90noiso", &passingMVA94Xwp90noiso, &b_passingMVA94Xwp90noiso);
   fChain->SetBranchAddress("passingMVA94XwpLiso", &passingMVA94XwpLiso, &b_passingMVA94XwpLiso);
   fChain->SetBranchAddress("passingMVA94XwpLnoiso", &passingMVA94XwpLnoiso, &b_passingMVA94XwpLnoiso);
   fChain->SetBranchAddress("passingMedium80X", &passingMedium80X, &b_passingMedium80X);
   fChain->SetBranchAddress("passingMedium94X", &passingMedium94X, &b_passingMedium94X);
   fChain->SetBranchAddress("passingTight80X", &passingTight80X, &b_passingTight80X);
   fChain->SetBranchAddress("passingTight94X", &passingTight94X, &b_passingTight94X);
   fChain->SetBranchAddress("passingVeto80X", &passingVeto80X, &b_passingVeto80X);
   fChain->SetBranchAddress("passingVeto94X", &passingVeto94X, &b_passingVeto94X);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("truePU", &truePU, &b_truePU);
   fChain->SetBranchAddress("event_nPV", &event_nPV, &b_mNPV);
   fChain->SetBranchAddress("event_met_pfmet", &event_met_pfmet, &b_mpfMET);
   fChain->SetBranchAddress("event_met_pfphi", &event_met_pfphi, &b_mpfPhi);
   fChain->SetBranchAddress("event_met_pfsumet", &event_met_pfsumet, &b_mpfSumET);
   fChain->SetBranchAddress("event_met_pfmetsignificance", &event_met_pfmetsignificance, &b_mpfMETSign);
   fChain->SetBranchAddress("event_PrimaryVertex_x", &event_PrimaryVertex_x, &b_mPVx);
   fChain->SetBranchAddress("event_PrimaryVertex_y", &event_PrimaryVertex_y, &b_mPVy);
   fChain->SetBranchAddress("event_PrimaryVertex_z", &event_PrimaryVertex_z, &b_mPVz);
   fChain->SetBranchAddress("event_BeamSpot_x", &event_BeamSpot_x, &b_mBSx);
   fChain->SetBranchAddress("event_BeamSpot_y", &event_BeamSpot_y, &b_mBSy);
   fChain->SetBranchAddress("event_BeamSpot_z", &event_BeamSpot_z, &b_mBSz);
   fChain->SetBranchAddress("event_rho", &event_rho, &b_rho);
   fChain->SetBranchAddress("mass", &mass, &b_mass);
   fChain->SetBranchAddress("tag_Ele_1overEminus1overP", &tag_Ele_1overEminus1overP, &b_tag_Ele_1overEminus1overP);
   fChain->SetBranchAddress("tag_Ele_5x5_sieie", &tag_Ele_5x5_sieie, &b_tag_Ele_5x5_sieie);
   fChain->SetBranchAddress("tag_Ele_abseta", &tag_Ele_abseta, &b_tag_Ele_abseta);
   fChain->SetBranchAddress("tag_Ele_dEtaSeed", &tag_Ele_dEtaSeed, &b_tag_Ele_dEtaSeed);
   fChain->SetBranchAddress("tag_Ele_dPhiIn", &tag_Ele_dPhiIn, &b_tag_Ele_dPhiIn);
   fChain->SetBranchAddress("tag_Ele_e", &tag_Ele_e, &b_tag_Ele_e);
   fChain->SetBranchAddress("tag_Ele_et", &tag_Ele_et, &b_tag_Ele_et);
   fChain->SetBranchAddress("tag_Ele_eta", &tag_Ele_eta, &b_tag_Ele_eta);
   fChain->SetBranchAddress("tag_Ele_hoe", &tag_Ele_hoe, &b_tag_Ele_hoe);
   fChain->SetBranchAddress("tag_Ele_phi", &tag_Ele_phi, &b_tag_Ele_phi);
   fChain->SetBranchAddress("tag_Ele_pt", &tag_Ele_pt, &b_tag_Ele_pt);
   fChain->SetBranchAddress("tag_Ele_q", &tag_Ele_q, &b_tag_Ele_q);
   fChain->SetBranchAddress("tag_sc_abseta", &tag_sc_abseta, &b_tag_sc_abseta);
   fChain->SetBranchAddress("tag_sc_e", &tag_sc_e, &b_tag_sc_e);
   fChain->SetBranchAddress("tag_sc_et", &tag_sc_et, &b_tag_sc_et);
   fChain->SetBranchAddress("tag_sc_eta", &tag_sc_eta, &b_tag_sc_eta);
   fChain->SetBranchAddress("tag_Ele_IsoMVA94X", &tag_Ele_IsoMVA94X, &b_tag_Ele_IsoMVA94X);
   fChain->SetBranchAddress("tag_Ele_dxy", &tag_Ele_dxy, &b_tag_Ele_dxy);
   fChain->SetBranchAddress("tag_Ele_dz", &tag_Ele_dz, &b_tag_Ele_dz);
   fChain->SetBranchAddress("tag_Ele_hoe_barrel", &tag_Ele_hoe_barrel, &b_tag_Ele_hoe_barrel);
   fChain->SetBranchAddress("tag_Ele_hoe_endcap", &tag_Ele_hoe_endcap, &b_tag_Ele_hoe_endcap);
   fChain->SetBranchAddress("tag_Ele_hzzMVA80X", &tag_Ele_hzzMVA80X, &b_tag_Ele_hzzMVA80X);
   fChain->SetBranchAddress("tag_Ele_mHits", &tag_Ele_mHits, &b_tag_Ele_mHits);
   fChain->SetBranchAddress("tag_Ele_miniIso", &tag_Ele_miniIso, &b_tag_Ele_miniIso);
   fChain->SetBranchAddress("tag_Ele_noIsoMVA94X", &tag_Ele_noIsoMVA94X, &b_tag_Ele_noIsoMVA94X);
   fChain->SetBranchAddress("tag_Ele_nonTrigMVA80X", &tag_Ele_nonTrigMVA80X, &b_tag_Ele_nonTrigMVA80X);
   fChain->SetBranchAddress("tag_Ele_passConversionVeto", &tag_Ele_passConversionVeto, &b_tag_Ele_passConversionVeto);
   fChain->SetBranchAddress("tag_Ele_sip3d", &tag_Ele_sip3d, &b_tag_Ele_sip3d);
   fChain->SetBranchAddress("tag_Ele_trigMVA", &tag_Ele_trigMVA, &b_tag_Ele_trigMVA);
   fChain->SetBranchAddress("mc_probe_e", &mc_probe_e, &b_mc_probe_e);
   fChain->SetBranchAddress("mc_probe_et", &mc_probe_et, &b_mc_probe_et);
   fChain->SetBranchAddress("mc_probe_eta", &mc_probe_eta, &b_mc_probe_eta);
   fChain->SetBranchAddress("mc_probe_phi", &mc_probe_phi, &b_mc_probe_phi);
   fChain->SetBranchAddress("mc_probe_flag", &mc_probe_flag, &b_mc_probe_flag);
   fChain->SetBranchAddress("pair_abseta", &pair_abseta, &b_pair_abseta);
   fChain->SetBranchAddress("pair_eta", &pair_eta, &b_pair_eta);
   fChain->SetBranchAddress("pair_mass", &pair_mass, &b_pair_mass);
   fChain->SetBranchAddress("pair_pt", &pair_pt, &b_pair_pt);
   Notify();
}

Bool_t treeAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void treeAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t treeAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef treeAnalysis_cxx

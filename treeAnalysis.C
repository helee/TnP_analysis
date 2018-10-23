#define treeAnalysis_cxx
#include "treeAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

Float_t weight = 1.0;

void treeAnalysis::Loop()
{
//   In a ROOT session, you can do:
//      root> .L treeAnalysis.C
//      root> treeAnalysis t
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

   const int ptbin = 5;
   float bin1[ptbin + 1] = {75., 100., 200., 300., 400., 1000.};

   bool tag_condi = false, SScharge = false, zmass = false; 
   bool IDcut0 = false, IPcut = false; // IDcut1 = false, IDcut2 = false, IDcut3 = false, IDcut4 = false, IDcut5 = false;  

   TFile *fnew = new TFile("output/output_"+inputfile, "RECREATE");
   TH1F *h_pass_barrel_SS = new TH1F("h_pass_barrel_SS", "", ptbin, bin1); 
   TH1F *h_pass_endcap_SS = new TH1F("h_pass_endcap_SS", "", ptbin, bin1);
   TH1F *h_fail_barrel_SS = new TH1F("h_fail_barrel_SS", "", ptbin, bin1);
   TH1F *h_fail_endcap_SS = new TH1F("h_fail_endcap_SS", "", ptbin, bin1);;
   TH1F *h_pass_barrel_OS = new TH1F("h_pass_barrel_OS", "", ptbin, bin1);
   TH1F *h_pass_endcap_OS = new TH1F("h_pass_endcap_OS", "", ptbin, bin1);
   TH1F *h_fail_barrel_OS = new TH1F("h_fail_barrel_OS", "", ptbin, bin1);
   TH1F *h_fail_endcap_OS = new TH1F("h_fail_endcap_OS", "", ptbin, bin1);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;

     if(tag_Ele_pt > 35 && tag_sc_abseta < 1.4442) tag_condi = true;
     if(tag_Ele_q*el_q > 0) SScharge = true;
     if((pair_mass > 70) && (pair_mass < 110)) zmass = true;

     if(abs(el_sc_eta) < 2.4 && el_miniIso < 0.1 && el_passConversionVeto == 1 && el_mHits < 2) IDcut0 = true;
     if(( abs(el_sc_eta) < 1.479 && el_dxy < 0.05 && el_dz < 0.1 ) || ( abs(el_sc_eta) > 1.479 && el_dxy < 0.1 && el_dz < 0.2 )) IPcut = true;    

     if(tag_condi && SScharge && zmass){
       if(el_pt > 75){
         if(passingLoose94X == 1){
           if(el_sc_abseta < 1.4442) h_pass_barrel_SS->Fill(el_pt, weight);
           if(el_sc_abseta > 1.566 && el_sc_abseta < 2.5) h_pass_endcap_SS->Fill(el_pt, weight);
         }
         else{
           if(el_sc_abseta < 1.4442) h_fail_barrel_SS->Fill(el_pt, weight);
           if(el_sc_abseta > 1.566 && el_sc_abseta < 2.5) h_fail_endcap_SS->Fill(el_pt, weight);
         }
       }
     }
     if(tag_condi && !SScharge && zmass){
       if(el_pt > 75){
         if(passingLoose94X == 1){
           if(el_sc_abseta < 1.4442) h_pass_barrel_OS->Fill(el_pt, weight);
           if(el_sc_abseta > 1.566 && el_sc_abseta < 2.5) h_pass_endcap_OS->Fill(el_pt, weight);
         }
         else{
           if(el_sc_abseta < 1.4442) h_fail_barrel_OS->Fill(el_pt, weight);
           if(el_sc_abseta > 1.566 && el_sc_abseta < 2.5) h_fail_endcap_OS->Fill(el_pt, weight);
         }
       }
     }
     

     SScharge = false; zmass = false; tag_condi = false; IDcut0 = false; IPcut = false;
   }

   fnew->Write();
   fnew->Close();
}

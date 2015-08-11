#include "TChain.h"
#include "TBrowser.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TPrincipal.h"
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;
typedef vector<TLorentzVector> LVs;
typedef TLorentzVector LV;

TLorentzVector getLV(float pt,float eta,float phi,float e){
  TLorentzVector v;
  v.SetPtEtaPhiE(pt,eta,phi,e);
  return v;
}

TLorentzVector getLV(float pt,float eta,float phi){
  TLorentzVector v;
  v.SetPtEtaPhiM(pt,eta,phi,0);
  return v;
}

LVs getLVs(uint n, float* pt,float* eta,float* phi,float* e){
  vector<TLorentzVector> vs;
  for(uint i=0; i<n; i++){
    TLorentzVector v= getLV(pt[i],eta[i],phi[i],e[i]);
    vs.push_back(v);
  }

  return vs;
}

void GetNuVecs(const TLorentzVector & lepvec, const TVector2 & metvec, TLorentzVector & nu1, TLorentzVector & nu2){
  double metvec2 = metvec.Px()*metvec.Px() + metvec.Py()*metvec.Py();
  double mu = (80.4*80.4)/2 + metvec.Px()*lepvec.Px() + metvec.Py()*lepvec.Py();
  double a = (mu*lepvec.Pz())/(lepvec.E()*lepvec.E() - lepvec.Pz()*lepvec.Pz());
  double a2 = TMath::Power(a, 2);
  double b = (TMath::Power(lepvec.E(), 2.)*metvec2 - TMath::Power(mu, 2.)) / (TMath::Power(lepvec.E(), 2)- TMath::Power(lepvec.Pz(), 2));
  float pz1,pz2;
  if (a2-b < 0) { 
    pz1 = a;
    pz2 = a;
  } else {
    double root = sqrt(a2-b);
    pz1 = a + root;
    pz2 = a - root;
  }
  nu1.SetPxPyPzE(metvec.Px(),metvec.Py(),pz1,sqrt(metvec.Mod2()+pz1*pz1));
  nu2.SetPxPyPzE(metvec.Px(),metvec.Py(),pz2,sqrt(metvec.Mod2()+pz2*pz2));
}

void pca(float mw, float mt, float& mw_new , float& mt_new){
  mt_new=(mw-83.7)/21.+(mt-172.5)/33.2;
  mw_new=(mw-83.7)/21.-(mt-172.5)/33.2;

}

void plot(){
  TH1F::SetDefaultSumw2();
  
  //  gStyle->SetOptStat(0);
  TChain* chain = new TChain("MVATree");
  char* filenames = getenv ("FILENAMES");
  char* outfilename = getenv ("OUTFILENAME");
  //  int maxevents = atoi(getenv ("MAXEVENTS"));
  string buf;
  stringstream ss(filenames); 
  while (ss >> buf){
    chain->Add(buf.c_str());
  }
  chain->SetBranchStatus("*",0);

  float Weight;
  chain->SetBranchAddress("Weight",&Weight);
  int N_Jets;
  chain->SetBranchAddress("N_Jets",&N_Jets);
  int N_BTagsM;
  chain->SetBranchAddress("N_BTagsM",&N_BTagsM);
  int N_BTagsT;
  chain->SetBranchAddress("N_BTagsT",&N_BTagsT);
  int N_BTagsL;
  chain->SetBranchAddress("N_BTagsL",&N_BTagsL);
  int N_TightLeptons;
  chain->SetBranchAddress("N_TightLeptons",&N_TightLeptons);
  float* Jet_Pt = new float[20];
  chain->SetBranchAddress("Jet_Pt",Jet_Pt);
  float* Jet_Phi = new float[20];
  chain->SetBranchAddress("Jet_Phi",Jet_Phi);
  float* Jet_Eta = new float[20];
  chain->SetBranchAddress("Jet_Eta",Jet_Eta);
  float* Jet_E = new float[20];
  chain->SetBranchAddress("Jet_E",Jet_E);
  float* Jet_CSV = new float[20];
  chain->SetBranchAddress("Jet_CSV",Jet_CSV);
  float* Jet_Flav = new float[20];
  chain->SetBranchAddress("Jet_Flav",Jet_Flav);
  float Evt_Pt_PrimaryLepton;
  chain->SetBranchAddress("Evt_Pt_PrimaryLepton",&Evt_Pt_PrimaryLepton);
  float Evt_Eta_PrimaryLepton;
  chain->SetBranchAddress("Evt_Eta_PrimaryLepton",&Evt_Eta_PrimaryLepton);
  float Evt_Phi_PrimaryLepton;
  chain->SetBranchAddress("Evt_Phi_PrimaryLepton",&Evt_Phi_PrimaryLepton);
  float Evt_E_PrimaryLepton;
  chain->SetBranchAddress("Evt_E_PrimaryLepton",&Evt_E_PrimaryLepton);
  float Evt_Pt_MET;
  chain->SetBranchAddress("Evt_Pt_MET",&Evt_Pt_MET);
  float Evt_Phi_MET;
  chain->SetBranchAddress("Evt_Phi_MET",&Evt_Phi_MET);
  int GenEvt_I_TTPlusCC;
  chain->SetBranchAddress("GenEvt_I_TTPlusCC",&GenEvt_I_TTPlusCC);
  int GenEvt_I_TTPlusBB;
  chain->SetBranchAddress("GenEvt_I_TTPlusBB",&GenEvt_I_TTPlusBB);
  int N_GenTopHad;
  chain->SetBranchAddress("N_GenTopHad",&N_GenTopHad);
  float* GenTopHad_Pt = new float[20];
  chain->SetBranchAddress("GenTopHad_Pt",GenTopHad_Pt);
  float* GenTopHad_Eta = new float[20];
  chain->SetBranchAddress("GenTopHad_Eta",GenTopHad_Eta);
  float* GenTopHad_Phi = new float[20];
  chain->SetBranchAddress("GenTopHad_Phi",GenTopHad_Phi);
  float* GenTopHad_B_Pt = new float[20];
  chain->SetBranchAddress("GenTopHad_B_Pt",GenTopHad_B_Pt);
  float* GenTopHad_Q1_Pt = new float[20];
  chain->SetBranchAddress("GenTopHad_Q1_Pt",GenTopHad_Q1_Pt);
  float* GenTopHad_Q2_Pt = new float[20];
  chain->SetBranchAddress("GenTopHad_Q2_Pt",GenTopHad_Q2_Pt);
  float* GenTopHad_B_Eta = new float[20];
  chain->SetBranchAddress("GenTopHad_B_Eta",GenTopHad_B_Eta);
  float* GenTopHad_Q1_Eta = new float[20];
  chain->SetBranchAddress("GenTopHad_Q1_Eta",GenTopHad_Q1_Eta);
  float* GenTopHad_Q2_Eta = new float[20];
  chain->SetBranchAddress("GenTopHad_Q2_Eta",GenTopHad_Q2_Eta);
  float* GenTopHad_B_Phi = new float[20];
  chain->SetBranchAddress("GenTopHad_B_Phi",GenTopHad_B_Phi);
  float* GenTopHad_Q1_Phi = new float[20];
  chain->SetBranchAddress("GenTopHad_Q1_Phi",GenTopHad_Q1_Phi);
  float* GenTopHad_Q2_Phi = new float[20];
  chain->SetBranchAddress("GenTopHad_Q2_Phi",GenTopHad_Q2_Phi);
  int N_GenTopLep;
  chain->SetBranchAddress("N_GenTopLep",&N_GenTopLep);
  float* GenTopLep_Pt = new float[20];
  chain->SetBranchAddress("GenTopLep_Pt",GenTopLep_Pt);
  float* GenTopLep_Eta = new float[20];
  chain->SetBranchAddress("GenTopLep_Eta",GenTopLep_Eta);
  float* GenTopLep_Phi = new float[20];
  chain->SetBranchAddress("GenTopLep_Phi",GenTopLep_Phi);
  float* GenTopLep_B_Pt = new float[20];
  chain->SetBranchAddress("GenTopLep_B_Pt",GenTopLep_B_Pt);
  float* GenTopLep_Lep_Pt = new float[20];
  chain->SetBranchAddress("GenTopLep_Lep_Pt",GenTopLep_Lep_Pt);
  float* GenTopLep_Nu_Pt = new float[20];
  chain->SetBranchAddress("GenTopLep_Nu_Pt",GenTopLep_Nu_Pt);
  float* GenTopLep_B_Eta = new float[20];
  chain->SetBranchAddress("GenTopLep_B_Eta",GenTopLep_B_Eta);
  float* GenTopLep_Lep_Eta = new float[20];
  chain->SetBranchAddress("GenTopLep_Lep_Eta",GenTopLep_Lep_Eta);
  float* GenTopLep_Nu_Eta = new float[20];
  chain->SetBranchAddress("GenTopLep_Nu_Eta",GenTopLep_Nu_Eta);
  float* GenTopLep_B_Phi = new float[20];
  chain->SetBranchAddress("GenTopLep_B_Phi",GenTopLep_B_Phi);
  float* GenTopLep_Lep_Phi = new float[20];
  chain->SetBranchAddress("GenTopLep_Lep_Phi",GenTopLep_Lep_Phi);
  float* GenTopLep_Nu_Phi = new float[20];
  chain->SetBranchAddress("GenTopLep_Nu_Phi",GenTopLep_Nu_Phi);
  float GenHiggs_Pt;
  chain->SetBranchAddress("GenHiggs_Pt",&GenHiggs_Pt);
  float GenHiggs_Eta;
  chain->SetBranchAddress("GenHiggs_Eta",&GenHiggs_Eta);
  float GenHiggs_Phi;
  chain->SetBranchAddress("GenHiggs_Phi",&GenHiggs_Phi);
  float GenHiggs_B1_Pt;
  chain->SetBranchAddress("GenHiggs_B1_Pt",&GenHiggs_B1_Pt);
  float GenHiggs_B2_Pt;
  chain->SetBranchAddress("GenHiggs_B2_Pt",&GenHiggs_B2_Pt);
  float GenHiggs_B1_Eta;
  chain->SetBranchAddress("GenHiggs_B1_Eta",&GenHiggs_B1_Eta);
  float GenHiggs_B2_Eta;
  chain->SetBranchAddress("GenHiggs_B2_Eta",&GenHiggs_B2_Eta);
  float GenHiggs_B1_Phi;
  chain->SetBranchAddress("GenHiggs_B1_Phi",&GenHiggs_B1_Phi);
  float GenHiggs_B2_Phi;
  chain->SetBranchAddress("GenHiggs_B2_Phi",&GenHiggs_B2_Phi);
  int N_AdditionalGenBJets;
  chain->SetBranchAddress("N_AdditionalGenBJets",&N_AdditionalGenBJets);
  float* AdditionalGenBJet_Pt = new float[20];
  chain->SetBranchAddress("AdditionalGenBJet_Pt",AdditionalGenBJet_Pt);
  float* AdditionalGenBJet_Eta = new float[20];
  chain->SetBranchAddress("AdditionalGenBJet_Eta",AdditionalGenBJet_Eta);
  float* AdditionalGenBJet_Phi = new float[20];
  chain->SetBranchAddress("AdditionalGenBJet_Phi",AdditionalGenBJet_Phi);
  float* AdditionalGenBJet_E = new float[20];
  chain->SetBranchAddress("AdditionalGenBJet_E",AdditionalGenBJet_E);

  
  // histos

  TH1F* h_M_Higgs_all=new TH1F("M_Higgs_all","M_Higgs_all",60,0,300);
  TH1F* h_M_TopLep_all=new TH1F("M_TopLep_all","M_TopLep_all",100,0,500);
  TH1F* h_M_WHad_all=new TH1F("M_WHad_all","M_WHad_all",60,0,300);
  TH1F* h_M_TopHad_all=new TH1F("M_TopHad_all","M_TopHad_all",100,0,500);
  TH2F* h_M_WHad_vs_M_TopHad_all=new TH2F("h_M_WHad_vs_M_TopHad_all","h_M_WHad_vs_M_TopHad_all",30,0,300,50,0,500);
  TH1F* h_M_WLep_all=new TH1F("M_WLep_all","M_WLep_all",60,0,300);
  TH1F* h_M_WHad_pca_all=new TH1F("M_WHad_pca_all","M_WHad_pca_all",100,-4,4);
  TH1F* h_M_TopHad_pca_all=new TH1F("M_TopHad_pca_all","M_TopHad_pca_all",100,-10,10);
  TH2F* h_M_WHad_vs_M_TopHad_pca_all=new TH2F("h_M_WHad_vs_M_TopHad_pca_all","h_M_WHad_vs_M_TopHad_pca_all",100,-4,4,100,-4,4);

  TH1F* h_M_Higgs_reco=new TH1F("M_Higgs_reco","M_Higgs_reco",60,0,300);
  TH1F* h_M_TopHad_reco=new TH1F("M_TopHad_reco","M_TopHad_reco",100,0,500);
  TH1F* h_M_TopLep_reco=new TH1F("M_TopLep_reco","M_TopLep_reco",100,0,500);
  TH1F* h_M_WHad_reco=new TH1F("M_WHad_reco","M_WHad_reco",60,0,300);
  TH1F* h_M_WHad_pca_reco=new TH1F("M_WHad_pca_reco","M_WHad_pca_reco",100,-4,4);
  TH1F* h_M_TopHad_pca_reco=new TH1F("M_TopHad_pca_reco","M_TopHad_pca_reco",100,-10,10);
  TH2F* h_M_WHad_vs_M_TopHad_reco=new TH2F("h_M_WHad_vs_M_TopHad_reco","h_M_WHad_vs_M_TopHad_reco",30,0,300,50,0,500);
  TH2F* h_M_WHad_vs_M_TopHad_pca_reco=new TH2F("h_M_WHad_vs_M_TopHad_pca_reco","h_M_WHad_vs_M_TopHad_pca_reco",100,-4,4,100,-4,4);
  TPrincipal* prince = new TPrincipal(2);
  double* data = new double[2];
  TH1F* h_M_WLep_reco=new TH1F("M_WLep_reco","M_WLep_reco",60,0,300);

  TH1F* h_M_Higgs_parton=new TH1F("M_Higgs_parton","M_Higgs_parton",60,0,300);
  TH1F* h_M_TopHad_parton=new TH1F("M_TopHad_parton","M_TopHad_parton",100,0,500);
  TH1F* h_M_TopLep_parton=new TH1F("M_TopLep_parton","M_TopLep_parton",100,0,500);
  TH1F* h_M_WHad_parton=new TH1F("M_WHad_parton","M_WHad_parton",60,0,300);
  TH2F* h_M_WHad_vs_M_TopHad_parton=new TH2F("h_M_WHad_vs_M_TopHad_parton","h_M_WHad_vs_M_TopHad_parton",30,0,300,50,0,500);
  TH1F* h_M_WLep_parton=new TH1F("M_WLep_parton","M_WLep_parton",60,0,300);

  TH1F* h_M_Higgs_recoall=new TH1F("M_Higgs_recoall","M_Higgs_recoall",60,0,300);
  TH1F* h_M_TopHad_recoall=new TH1F("M_TopHad_recoall","M_TopHad_recoall",100,0,500);
  TH1F* h_M_TopLep_recoall=new TH1F("M_TopLep_recoall","M_TopLep_recoall",100,0,500);
  TH1F* h_M_WHad_recoall=new TH1F("M_WHad_recoall","M_WHad_recoall",60,0,300);
  TH2F* h_M_WHad_vs_M_TopHad_recoall=new TH2F("h_M_WHad_vs_M_TopHad_recoall","h_M_WHad_vs_M_TopHad_recoall",30,0,300,50,0,500);
  TH1F* h_M_WLep_recoall=new TH1F("M_WLep_recoall","M_WLep_recoall",60,0,300);

  TH1F* h_M_Higgs_best=new TH1F("M_Higgs_best","M_Higgs_best",60,0,300);
  TH1F* h_M_TopHad_best=new TH1F("M_TopHad_best","M_TopHad_best",100,0,500);
  TH1F* h_M_TopHad_best_special=new TH1F("M_TopHad_best_special","M_TopHad_best_special",200,-500,500);
  TH1F* h_M_TopLep_best=new TH1F("M_TopLep_best","M_TopLep_best",100,0,500);
  TH1F* h_M_WHad_best=new TH1F("M_WHad_best","M_WHad_best",60,0,300);
  TH2F* h_M_WHad_vs_M_TopHad_best=new TH2F("h_M_WHad_vs_M_TopHad_best","h_M_WHad_vs_M_TopHad_best",30,0,300,50,0,500);
  TH1F* h_M_WLep_best=new TH1F("M_WLep_best","M_WLep_best",60,0,300);

  TH1F* h_CSV_Higgs_best=new TH1F("CSV_Higgs_best","CSV_Higgs_best",22,-0.1,1);
  TH1F* h_CSV_WHad_best=new TH1F("CSV_WHad_best","CSV_WHad_best",22,-0.1,1);
  TH1F* h_CSV_TopHad_best=new TH1F("CSV_TopHad_best","CSV_TopHad_best",22,-0.1,1);
  TH1F* h_CSV_TopLep_best=new TH1F("CSV_TopLep_best","CSV_TopLep_best",22,-0.1,1);

  TH1F* h_M_WHad_80=new TH1F("M_WHad_80","M_WHad_80",60,0,300);
  TH1F* h_M_WHad_80_matchable=new TH1F("M_WHad_80_matchable","M_WHad_80_matchable",60,0,300);
  TH1F* h_M_WHad_80_unmatchable=new TH1F("M_WHad_80_unmatchable","M_WHad_80_unmatchable",60,0,300);

  TH1F* h_dr2=new TH1F("dr2_reco","dr2_reco",30,0,1.5);
  TH1F* h_B1_dr=new TH1F("B1_dr_reco","B1_dr_reco",30,0,1.5);
  TH1F* h_B2_dr=new TH1F("B2_dr_reco","B2_dr_reco",30,0,1.5);
  TH1F* h_BLep_dr=new TH1F("BLep_dr_reco","BLep_dr_reco",30,0,1.5);
  TH1F* h_BHad_dr=new TH1F("BHad_dr_reco","BHad_dr_reco",30,0,1.5);
  TH1F* h_Q1_dr=new TH1F("Q1_dr_reco","Q1_dr_reco",30,0,1.5);
  TH1F* h_Q2_dr=new TH1F("Q2_dr_reco","Q2_dr_reco",30,0,1.5);
  TH1F* h_Lep_dr=new TH1F("Lep_dr_reco","Lep_dr_reco",30,0,1.5);
  TH1F* h_Nu_dr=new TH1F("Nu_dr_reco","Nu_dr_reco",30,0,1.5);

  TH1F* h_CSV_all=new TH1F("CSV_all","CSV_all",44,-0.1,1);
  TH1F* h_CSV_b=new TH1F("CSV_b","CSV_b",44,-0.1,1);
  TH1F* h_CSV_l_w_c=new TH1F("CSV_l_w_c","CSV_l_w_c",44,-0.1,1);
  TH1F* h_CSV_l_wo_c=new TH1F("CSV_l_wo_c","CSV_l_wo_c",44,-0.1,1);
  TH1F* h_CSV_c=new TH1F("CSV_c","CSV_c",44,-0.1,1);

  TH1F* h_nb=new TH1F("nb","nb",6,-0.5,5.5);
  TH1F* h_N_BTagsL=new TH1F("N_BTagsL","N_BTagsL",6,-0.5,5.5);
  TH1F* h_N_BTagsM=new TH1F("N_BTagsM","N_BTagsM",6,-0.5,5.5);
  TH1F* h_N_BTagsT=new TH1F("N_BTagsT","N_BTagsT",6,-0.5,5.5);
  TH1F* h_N_BTagsL_4b=new TH1F("N_BTagsL_4b","N_BTagsL_4b",6,-0.5,5.5);
  TH1F* h_N_BTagsM_4b=new TH1F("N_BTagsM_4b","N_BTagsM_4b",6,-0.5,5.5);
  TH1F* h_N_BTagsT_4b=new TH1F("N_BTagsT_4b","N_BTagsT_4b",6,-0.5,5.5);
  TH1F* h_N_BTagsL_3b=new TH1F("N_BTagsL_3b","N_BTagsL_3b",6,-0.5,5.5);
  TH1F* h_N_BTagsM_3b=new TH1F("N_BTagsM_3b","N_BTagsM_3b",6,-0.5,5.5);
  TH1F* h_N_BTagsT_3b=new TH1F("N_BTagsT_3b","N_BTagsT_3b",6,-0.5,5.5);
  TH1F* h_N_BTagsL_2b=new TH1F("N_BTagsL_2b","N_BTagsL_2b",6,-0.5,5.5);
  TH1F* h_N_BTagsM_2b=new TH1F("N_BTagsM_2b","N_BTagsM_2b",6,-0.5,5.5);
  TH1F* h_N_BTagsT_2b=new TH1F("N_BTagsT_2b","N_BTagsT_2b",6,-0.5,5.5);

  TH1F* h_M_Total_recoall=new TH1F("M_Total_recoall","M_Total_recoall",60,0,3000);
  TH1F* h_Pt_Total_recoall=new TH1F("Pt_Total_recoall","Pt_Total_recoall",50,0,500);

  TH1F* h_Pt_Total_recoall_6j=new TH1F("Pt_Total_recoall_6j","Pt_Total_recoall_6j",50,0,500);
  TH1F* h_Pt_Total_recoall_7j=new TH1F("Pt_Total_recoall_7j","Pt_Total_recoall_7j",50,0,500);
  TH1F* h_Pt_Total_recoall_8j=new TH1F("Pt_Total_recoall_8j","Pt_Total_recoall_8j",50,0,500);
  TH1F* h_Pt_Total_recoall_9j=new TH1F("Pt_Total_recoall_9j","Pt_Total_recoall_9j",50,0,500);
  TH1F* h_Pt_Total_recoall_10j=new TH1F("Pt_Total_recoall_10j","Pt_Total_recoall_10j",50,0,500);


  TH1F* h_B1_pt_reco=new TH1F("B1_pt_reco","B1_pt_reco",30,0,300);
  TH1F* h_B2_pt_reco=new TH1F("B2_pt_reco","B2_pt_reco",30,0,300);
  TH1F* h_BLep_pt_reco=new TH1F("BLep_pt_reco","BLep_pt_reco",30,0,300);
  TH1F* h_BHad_pt_reco=new TH1F("BHad_pt_reco","BHad_pt_reco",30,0,300);
  TH1F* h_Q1_pt_reco=new TH1F("Q1_pt_reco","Q1_pt_reco",30,0,300);
  TH1F* h_Q2_pt_reco=new TH1F("Q2_pt_reco","Q2_pt_reco",30,0,300);



  // loop
  long nentries = chain->GetEntries(); 
  cout << "total number of events: " << nentries << endl;
  int nMatchable=0;
  int nAll=0;
  for (long iEntry=0;iEntry<nentries;iEntry++) {
    //    if(iEntry==20000) break;
    if(iEntry%10000==0) cout << "analyzing event " << iEntry << endl;

    chain->GetEntry(iEntry); 
    // event selection
    // recect non-semileptonic
    if(N_Jets<6||N_GenTopHad!=1||N_GenTopLep!=1) continue;

    TLorentzVector vHiggs_true=getLV(GenHiggs_Pt,GenHiggs_Eta,GenHiggs_Phi);
    TLorentzVector vTopHad_true=getLV(GenTopHad_Pt[0],GenTopHad_Eta[0],GenTopHad_Phi[0]);
    TLorentzVector vTopLep_true=getLV(GenTopLep_Pt[0],GenTopLep_Eta[0],GenTopLep_Phi[0]);
    TLorentzVector vQ1_true=getLV(GenTopHad_Q1_Pt[0],GenTopHad_Q1_Eta[0],GenTopHad_Q1_Phi[0]);
    TLorentzVector vQ2_true=getLV(GenTopHad_Q2_Pt[0],GenTopHad_Q2_Eta[0],GenTopHad_Q2_Phi[0]);
    TLorentzVector vBHad_true=getLV(GenTopHad_B_Pt[0],GenTopHad_B_Eta[0],GenTopHad_B_Phi[0]);
    TLorentzVector vNu_true=getLV(GenTopLep_Nu_Pt[0],GenTopLep_Nu_Eta[0],GenTopLep_Nu_Phi[0]);
    TLorentzVector vLep_true=getLV(GenTopLep_Lep_Pt[0],GenTopLep_Lep_Eta[0],GenTopLep_Lep_Phi[0]);
    TLorentzVector vBLep_true=getLV(GenTopLep_B_Pt[0],GenTopLep_B_Eta[0],GenTopLep_B_Phi[0]);
    TLorentzVector vB1_true;
    TLorentzVector vB2_true;
    if(GenHiggs_B1_Pt>0.1){
      vB1_true=getLV(GenHiggs_B1_Pt,GenHiggs_B1_Eta,GenHiggs_B1_Phi);
      vB2_true=getLV(GenHiggs_B2_Pt,GenHiggs_B2_Eta,GenHiggs_B2_Phi);
    }
    else if(N_AdditionalGenBJets>=2){
      vB1_true=getLV(AdditionalGenBJet_Pt[0],AdditionalGenBJet_Eta[0],AdditionalGenBJet_Phi[0],AdditionalGenBJet_E[0]);
      vB2_true=getLV(AdditionalGenBJet_Pt[1],AdditionalGenBJet_Eta[1],AdditionalGenBJet_Phi[1],AdditionalGenBJet_E[1]);
    }
    else continue;

    h_M_Higgs_parton->Fill((vB1_true+vB2_true).M(),Weight);
    h_M_TopHad_parton->Fill((vQ2_true+vQ1_true+vBHad_true).M(),Weight);
    h_M_TopLep_parton->Fill((vBLep_true+vNu_true+vLep_true).M(),Weight);
    h_M_WHad_parton->Fill((vQ2_true+vQ1_true).M(),Weight);
    h_M_WHad_vs_M_TopHad_parton->Fill((vQ2_true+vQ1_true).M(),(vQ2_true+vQ1_true+vBHad_true).M(),Weight);
    h_M_WLep_parton->Fill((vNu_true+vLep_true).M(),Weight);

    
    LVs jetvecs = getLVs(N_Jets,Jet_Pt,Jet_Eta,Jet_Phi,Jet_E);
    LV lepvec = getLV(Evt_Pt_PrimaryLepton,Evt_Eta_PrimaryLepton,Evt_Phi_PrimaryLepton,Evt_E_PrimaryLepton);
    TVector2 metvec;
    metvec.SetMagPhi(Evt_Pt_MET,Evt_Phi_MET);
    LV nuvec1;
    LV nuvec2;
    GetNuVecs(lepvec,metvec,nuvec1,nuvec2);
    uint nNu=2;
    if(fabs(nuvec1.Pz()-nuvec2.Pz()<1)) nNu=1;
    bool matchable=false;

    float minDr2=99999;
    float minDrQ1=99999;
    float minDrQ2=99999;
    float minDrBHad=99999;
    float minDrBLep=99999;
    float minDrB1=99999;
    float minDrB2=99999;
    float minDrNu=99999;

    int bestQ1=-1;
    int bestQ2=-1;
    int bestBHad=-1;
    int bestBLep=-1;
    int bestB1=-1;
    int bestB2=-1;
    int bestNu=-1;

    int bestQ1all=-1;
    int bestQ2all=-1;
    int bestBHadall=-1;
    int bestBLepall=-1;
    int bestB1all=-1;
    int bestB2all=-1;
    int bestNuall=-1;

    int ncombis=0;
    float minChi2=99999;
    for(uint iQ1=0; iQ1<jetvecs.size();iQ1++){
      for(uint iQ2=0; iQ2<jetvecs.size();iQ2++){
	if(iQ2!=iQ1){
	  for(uint iBHad=0; iBHad<jetvecs.size();iBHad++){    
	    if(iBHad!=iQ1 && iBHad!=iQ2){
	      for(uint iBLep=0; iBLep<jetvecs.size();iBLep++){
		if(iBLep!=iQ1 && iBLep!=iQ2 && iBLep!=iBHad) {
		  for(uint iB1=0; iB1<jetvecs.size();iB1++){
		    if(iB1!=iQ1 && iB1!=iQ2 && iB1!=iBHad && iB1!=iBLep) {
		      for(uint iB2=0; iB2<jetvecs.size();iB2++){
			if(iB2!=iQ1 && iB2!=iQ2 && iB2!=iBHad && iB2!=iBLep && iB2!=iB1){
			  for(uint iNu=1; iNu<=nNu;iNu++){
			    
			    if(bestQ1!=iQ1&&minDrQ1<0.2) continue;
			    if(bestQ2!=iQ2&&minDrQ2<0.2) continue;
			    if(bestBHad!=iBHad&&minDrBHad<0.2) continue;
			    if(bestBLep!=iBLep&&minDrBLep<0.2) continue;			    
			    if(bestB1!=iB1&&minDrB1<0.2) continue;
			    if(bestB2!=iB2&&minDrB2<0.2) continue;

    /*
    for(uint iQ1=0; iQ1<jetvecs.size()&&minDrQ1>0.2;iQ1++){
      for(uint iQ2=0; iQ2<jetvecs.size()&&minDrQ2>0.2;iQ2++){
	if(iQ2==iQ1) continue;
	for(uint iBHad=0; iBHad<jetvecs.size()&&minDrBHad>0.2;iBHad++){    
	  if(iBHad==iQ1 || iBHad==iQ2) continue;
	  for(uint iBLep=0; iBLep<jetvecs.size()&&minDrBLep>0.2;iBLep++){
	    if(iBLep==iQ1 || iBLep==iQ2 || iBLep==iBHad) continue;
	    for(uint iB1=0; iB1<jetvecs.size()&&minDrB1>0.2;iB1++){
	      if(iB1==iQ1 || iB1==iQ2 || iB1==iBHad || iB1==iBLep) continue;
	      for(uint iB2=0; iB2<jetvecs.size()&&minDrB2>0.2;iB2++){
		if(iB2==iQ1 || iB2==iQ2 || iB2==iBHad || iB2==iBLep || iB2==iB1) continue;
		for(uint iNu=1; iNu<=nNu;iNu++){
    */
			    /*			    cout << "+++++++++++++++++" << endl;
			    cout << minDrQ1 << endl;
			    cout << minDrQ2 << endl;
			    cout << minDrBHad << endl;
			    cout << minDrBLep << endl;
			    cout << minDrB1 << endl;
			    cout << minDrB2 << endl;
			    */
			    // mc matching
			    ncombis++;
			    
			    
			    const float maxdr=0.2;
			    float drQ1=jetvecs[iQ1].DeltaR(vQ1_true);
			    bool matchQ1=drQ1<maxdr;
			    float drQ2=jetvecs[iQ2].DeltaR(vQ2_true);
			    bool matchQ2=drQ2<maxdr;
			    float drBHad=jetvecs[iBHad].DeltaR(vBHad_true);
			    bool matchBHad=drBHad<maxdr;
			    float drBLep=jetvecs[iBLep].DeltaR(vBLep_true);
			    bool matchBLep=drBLep<maxdr;
			    float drB1=jetvecs[iB1].DeltaR(vB1_true);
			    bool matchB1=drB1<maxdr;
			    float drB2=jetvecs[iB2].DeltaR(vB2_true);
			    bool matchB2=drB2<maxdr;
			    float drNu=999;
			    bool matchNu=false;
			    if(iNu==1){
			      drNu=nuvec1.DeltaR(vNu_true);
			      float drNu_alt=nuvec2.DeltaR(vNu_true);
			      matchNu=nNu==1 || drNu<drNu_alt;
			    }
			    else{
			      drNu=nuvec2.DeltaR(vNu_true);
			      float drNu_alt=nuvec1.DeltaR(vNu_true);
			      matchNu=drNu<drNu_alt;
			    }
			    
			    /*			    cout <<"========================== " << endl;
			    cout <<"indexes ";
			    cout << iQ1 << " ";
			    cout << iQ2 << " ";
			    cout << iBHad << " ";
			    cout << iBLep << " ";
			    cout << iB1 << " ";
			    cout << iB2 << " ";
			    cout << iNu << " ";
			    cout << endl;
			    
			    cout <<"dr      ";
			    cout << drQ1 << " ";
			    cout << drQ2 << " ";
			    cout << drBHad << " ";
			    cout << drBLep << " ";
			    cout << drB1 << " ";
			    cout << drB2 << " ";
			    cout << drNu << " ";
			    cout << endl;
			    
			    cout <<"matches ";
			    cout << matchQ1 << " ";
			    cout << matchQ2 << " ";
			    cout << matchBHad << " ";
			    cout << matchBLep << " ";
			    cout << matchB1 << " ";
			    cout << matchB2 << " ";
			    cout << matchNu << " ";
			    cout << endl;*/
			    
			    if(matchQ1 && matchQ2 && matchBHad && matchBLep &&matchB1 && matchB2 && matchNu){
			      matchable=true;
			    }
			    float dr2=drQ1*drQ1+drQ2*drQ2+drBHad*drBHad+drBLep*drBLep+drB1*drB1+drB2*drB2;
			    if(dr2<minDr2&&matchNu){
			      minDr2=dr2;
			      bestQ1all=iQ1;
			      bestQ2all=iQ2;
			      bestBHadall=iBHad;
			      bestBLepall=iBLep;
			      bestB1all=iB1;
			      bestB2all=iB2;
			      bestNuall=iNu;
			    }
			    if(drQ1<minDrQ1){
			      minDrQ1=drQ1;
			      bestQ1=iQ1;
			    }
			    if(drQ2<minDrQ2){
			      minDrQ2=drQ2;
			      bestQ2=iQ2;
			    }
			    if(drBHad<minDrBHad){
			      minDrBHad=drBHad;
			      bestBHad=iBHad;
			    }
			    if(drBLep<minDrBLep){
			      minDrBLep=drBLep;
			      bestBLep=iBLep;
			    }
			    if(drB1<minDrB1){
			      minDrB1=drB1;
			      bestB1=iB1;
			    }
			    if(drB2<minDrB2){
			      minDrB2=drB2;
			      bestB2=iB2;
			    }
			    if(matchNu){
			      minDrNu=drNu;
			      bestNu=iNu;
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  
    /*    cout << "=========" <<endl;
    cout <<"bestQ1 " << bestQ1 << endl;
    cout <<"bestQ2 " << bestQ2 << endl;
    cout <<"bestBHad " << bestBHad << endl;
    cout <<"bestBLep " << bestBLep << endl;
    cout <<"bestB1 " << bestB1 << endl;
    cout <<"bestB2 " << bestB2 << endl;*/

    for(int i=0; i<N_Jets;i++){
      if(fabs(fabs(Jet_Flav[i])-5.)<0.1){
	h_CSV_b->Fill(fmax(-0.0999,Jet_CSV[i]),Weight);
      }
      else if(fabs(fabs(Jet_Flav[i])-4.)<0.1){
	h_CSV_c->Fill(fmax(-0.0999,Jet_CSV[i]),Weight);
	h_CSV_l_w_c->Fill(fmax(-0.0999,Jet_CSV[i]),Weight);
      }
      else{
	h_CSV_l_w_c->Fill(fmax(-0.0999,Jet_CSV[i]),Weight);
	h_CSV_l_wo_c->Fill(fmax(-0.0999,Jet_CSV[i]),Weight);
      }
    }    

    for(uint iB1=0; iB1<jetvecs.size();iB1++){
      for(uint iB2=0; iB2<jetvecs.size();iB2++){
	if(iB1<=iB2) continue;
	if(Jet_CSV[iB1]<.89||Jet_CSV[iB2]<.89) continue;
	h_M_Higgs_all->Fill((jetvecs[iB1]+jetvecs[iB2]).M(),Weight);
      }
    }
    for(uint iQ1=0; iQ1<jetvecs.size();iQ1++){
      for(uint iQ2=0; iQ2<jetvecs.size();iQ2++){
	if(iQ1<=iQ2) continue;
	h_M_WHad_all->Fill((jetvecs[iQ1]+jetvecs[iQ2]).M(),Weight);
	for(uint iBHad=0; iBHad<jetvecs.size();iBHad++){
	  if(iBHad==iQ2||iBHad==iQ1) continue;
	  if(Jet_CSV[iBHad]<.89) continue;
	  float M_TopHad_all=(jetvecs[iBHad]+jetvecs[iQ1]+jetvecs[iQ2]).M();
	  float M_WHad_all=(jetvecs[iQ1]+jetvecs[iQ2]).M();
	  h_M_TopHad_all->Fill(M_TopHad_all,Weight);
	  h_M_WHad_vs_M_TopHad_all->Fill(M_WHad_all,M_TopHad_all,Weight);
	  float pca0;
	  float pca1;
	  pca(M_WHad_all,M_TopHad_all,pca0,pca1);
	  h_M_WHad_vs_M_TopHad_pca_all->Fill(pca0,pca1,Weight);
	  h_M_TopHad_pca_all->Fill(pca1);
	  h_M_WHad_pca_all->Fill(pca0);

	}
      }
    }
    for(uint iNu=1; iNu<=nNu;iNu++){
      if(iNu==1){
	h_M_WLep_all->Fill((lepvec+nuvec1).M(),Weight);
      }
      else {
	h_M_WLep_all->Fill((lepvec+nuvec2).M(),Weight);		    
      }     
      for(uint iBLep=0; iBLep<jetvecs.size();iBLep++){
	if(Jet_CSV[iBLep]<.89) continue;
	if(iNu==1){
	  h_M_TopLep_all->Fill((jetvecs[iBLep]+lepvec+nuvec1).M(),Weight);
	}
	else {
	  h_M_TopLep_all->Fill((jetvecs[iBLep]+lepvec+nuvec2).M(),Weight);
	}     
      }
    }

    
    nAll++;
    if(matchable)
      nMatchable++;
    /*
    cout << N_Jets << endl;
    cout << nNu << endl;
    cout << ncombis << endl;
    if(matchable){
      cout << "++++++++++++++++++++++++++++" << endl;
      cout << TopHad_Q2_Pt[0] << endl;
    }
    
    else{
      cout << "----------------------------" << endl;
      cout << TopHad_Q2_Pt[0] << endl;
    }
    */
    float M_Higgs_reco = -999;
    float M_TopHad_reco = -999;
    float M_TopLep_reco = -999;
    float M_WHad_reco = -999;
    float M_WLep_reco = -999;
    if(bestB1>-1 && bestB2>-1){
      M_Higgs_reco=(jetvecs[bestB1]+jetvecs[bestB2]).M();
    }
    if(bestBHad>-1 && bestQ1>-1 && bestQ2>-1){
      M_TopHad_reco=(jetvecs[bestBHad]+jetvecs[bestQ1]+jetvecs[bestQ2]).M();
    }
    if(bestQ1>-1 && bestQ2>-1){
      M_WHad_reco=(jetvecs[bestQ1]+jetvecs[bestQ2]).M();
    }
    if(bestBLep>-1 && bestNu>-1){
      if(bestNu==1)
	M_TopLep_reco=(jetvecs[bestBLep]+lepvec+nuvec1).M();
      else
	M_TopLep_reco=(jetvecs[bestBLep]+lepvec+nuvec2).M();
    }
    if(bestNu>-1){
      if(bestNu==1)
	M_WLep_reco=(lepvec+nuvec1).M();
      else
	M_WLep_reco=(lepvec+nuvec2).M();
    }
    float maxdr=0.2;
    if(minDrB1<maxdr&&minDrB2<maxdr)
      h_M_Higgs_reco->Fill(M_Higgs_reco,Weight);
    if(minDrQ1<maxdr&&minDrQ2<maxdr&&minDrBHad<maxdr){
      h_M_TopHad_reco->Fill(M_TopHad_reco,Weight);
      h_M_WHad_vs_M_TopHad_reco->Fill(M_WHad_reco,M_TopHad_reco,Weight);
      data[0]=M_WHad_reco;
      data[1]=M_TopHad_reco;
      prince->AddRow(data);
      float pca0;
      float pca1;
      pca(M_WHad_reco,M_TopHad_reco,pca0,pca1);
      h_M_WHad_vs_M_TopHad_pca_reco->Fill(pca0,pca1,Weight);
      h_M_TopHad_pca_reco->Fill(pca1,Weight);
      h_M_WHad_pca_reco->Fill(pca0,Weight);
      
    }
    if(minDrBHad<maxdr)
      h_M_TopLep_reco->Fill(M_TopLep_reco,Weight);
    if(minDrQ1<maxdr&&minDrQ2<maxdr){
      h_M_WHad_reco->Fill(M_WHad_reco,Weight);
    }
    h_M_WLep_reco->Fill(M_WLep_reco,Weight);

    
    float M_Higgs_recoall = -999;
    float M_TopHad_recoall = -999;
    float M_TopLep_recoall = -999;
    float M_WHad_recoall = -999;
    float M_WLep_recoall = -999;
    if(minDr2<0.3){
      if(bestB1all>-1 && bestB2all>-1){
	M_Higgs_recoall=(jetvecs[bestB1all]+jetvecs[bestB2all]).M();
	h_M_Higgs_recoall->Fill(M_Higgs_recoall,Weight);
      }
      if(bestBHadall>-1 && bestQ1all>-1 && bestQ2all>-1){
	M_TopHad_recoall=(jetvecs[bestBHadall]+jetvecs[bestQ1all]+jetvecs[bestQ2all]).M();
	h_M_TopHad_recoall->Fill(M_TopHad_recoall,Weight);
	h_M_WHad_recoall->Fill(M_WHad_recoall,Weight);
	h_M_WHad_vs_M_TopHad_recoall->Fill(M_WHad_recoall,M_TopHad_recoall,Weight);
		
      }
      if(bestQ1all>-1 && bestQ2all>-1){
    	M_WHad_recoall=(jetvecs[bestQ1all]+jetvecs[bestQ2all]).M();
	h_M_WHad_recoall->Fill(M_WHad_recoall,Weight);
      }
      if(bestBLepall>-1 && bestNuall>-1){
	if(bestNuall==1){
	  M_TopLep_recoall=(jetvecs[bestBLepall]+lepvec+nuvec1).M();
	}
	else{
	  M_TopLep_recoall=(jetvecs[bestBLepall]+lepvec+nuvec2).M();
	}
	h_M_TopLep_recoall->Fill(M_TopLep_recoall,Weight);
      }
      if(bestNuall>-1){
	if(bestNuall==1)
	  M_WLep_recoall=(lepvec+nuvec1).M();
	else
	  M_WLep_recoall=(lepvec+nuvec2).M();
	h_M_WLep_recoall->Fill(M_WLep_recoall,Weight);
      }

      TLorentzVector p4_total;
      p4_total+=jetvecs[bestB1all];
      p4_total+=jetvecs[bestB2all];
      p4_total+=jetvecs[bestQ1all];
      p4_total+=jetvecs[bestQ2all];
      p4_total+=jetvecs[bestBHadall];
      p4_total+=jetvecs[bestBLepall];
      p4_total+=lepvec;
      if(bestNuall==1)
	p4_total+=nuvec1;
      else
	p4_total+=nuvec2;

      h_M_Total_recoall->Fill(p4_total.M());
      h_Pt_Total_recoall->Fill(p4_total.Pt());


    }


    float M_Higgs_best = -999;
    float M_TopHad_best = -999;
    float M_TopHad_best_special = -999;
    float M_TopLep_best = -999;
    float M_WHad_best = -999;
    float M_WLep_best = -999;

    float CSV2_Higgs_best = -999;
    float CSV1_Higgs_best = -999;
    float CSV2_WHad_best = -999;
    float CSV1_WHad_best = -999;
    float CSV_TopHad_best = -999;
    float CSV_TopLep_best = -999;

    if(bestB1all>-1 && bestB2all>-1){
      M_Higgs_best=(jetvecs[bestB1all]+jetvecs[bestB2all]).M();
      CSV1_Higgs_best=Jet_CSV[bestB1all];
      CSV2_Higgs_best=Jet_CSV[bestB2all];
    }
    if(bestBHadall>-1 && bestQ1all>-1 && bestQ2all>-1){
      M_TopHad_best=(jetvecs[bestBHadall]+jetvecs[bestQ1all]+jetvecs[bestQ2all]).M();
      if((jetvecs[bestQ1all]+jetvecs[bestQ2all]).M()<110.){
	M_TopHad_best_special=(jetvecs[bestBHadall]+jetvecs[bestQ1all]+jetvecs[bestQ2all]).M();
      }
      else{
	if(minDrQ1<minDrQ2){
	  M_TopHad_best_special=-(jetvecs[bestBHadall]+jetvecs[bestQ1all]).M();
	}
	else{
	  M_TopHad_best_special=-(jetvecs[bestBHadall]+jetvecs[bestQ2all]).M();
	}
      }
      CSV_TopHad_best=Jet_CSV[bestBHadall];
    }
    if(bestQ1all>-1 && bestQ2all>-1){
      M_WHad_best=(jetvecs[bestQ1all]+jetvecs[bestQ2all]).M();
      CSV1_WHad_best=Jet_CSV[bestQ1all];
      CSV2_WHad_best=Jet_CSV[bestQ2all];
    }
    if(bestBLepall>-1 && bestNuall>-1){
      if(bestNuall==1)
	M_TopLep_best=(jetvecs[bestBLepall]+lepvec+nuvec1).M();
      else
	M_TopLep_best=(jetvecs[bestBLepall]+lepvec+nuvec2).M();
      CSV_TopLep_best=Jet_CSV[bestBLepall];
    }
    if(bestNuall>-1){
      if(bestNuall==1)
	M_WLep_best=(lepvec+nuvec1).M();
      else
	M_WLep_best=(lepvec+nuvec2).M();
    }
    
    h_M_Higgs_best->Fill(M_Higgs_best,Weight);
    h_M_TopHad_best->Fill(M_TopHad_best,Weight);
    h_M_TopHad_best_special->Fill(M_TopHad_best_special,Weight);
    h_M_TopLep_best->Fill(M_TopLep_best,Weight);
    h_M_WHad_best->Fill(M_WHad_best,Weight);
    h_M_WHad_vs_M_TopHad_best->Fill(M_WHad_best,M_TopHad_best,Weight);
    h_M_WLep_best->Fill(M_WLep_best,Weight);
    h_CSV_Higgs_best->Fill(CSV2_Higgs_best,Weight);
    h_CSV_Higgs_best->Fill(CSV1_Higgs_best,Weight);
    h_CSV_WHad_best->Fill(CSV2_WHad_best,Weight);
    h_CSV_WHad_best->Fill(CSV1_WHad_best,Weight);
    h_CSV_TopHad_best->Fill(CSV_TopHad_best,Weight);
    h_CSV_TopLep_best->Fill(CSV_TopLep_best,Weight);

    
    h_dr2->Fill(minDr2,Weight);
    h_B1_dr->Fill(minDrB1,Weight);
    h_B2_dr->Fill(minDrB2,Weight);
    h_BLep_dr->Fill(minDrBLep,Weight);
    h_BHad_dr->Fill(minDrBHad,Weight);
    h_Q1_dr->Fill(minDrQ1,Weight);
    h_Q2_dr->Fill(minDrQ2,Weight);
    h_Nu_dr->Fill(minDrNu,Weight);
    h_Lep_dr->Fill(lepvec.DeltaR(vLep_true),Weight);
    int nb=0;
    for(size_t i=0; i<N_Jets; i++){
      if(fabs(fabs(Jet_Flav[i])-5)<0.1){
	nb++;
      }
    }
    h_nb->Fill(nb,Weight);
    h_N_BTagsM->Fill(N_BTagsM,Weight);
    h_N_BTagsL->Fill(N_BTagsL,Weight);
    h_N_BTagsT->Fill(N_BTagsT,Weight);
    if(nb==2){
      h_N_BTagsM_2b->Fill(N_BTagsM,Weight);
      h_N_BTagsL_2b->Fill(N_BTagsL,Weight);
      h_N_BTagsT_2b->Fill(N_BTagsT,Weight);     
    }
    if(nb==3){
      h_N_BTagsM_3b->Fill(N_BTagsM,Weight);
      h_N_BTagsL_3b->Fill(N_BTagsL,Weight);
      h_N_BTagsT_3b->Fill(N_BTagsT,Weight);     
    }
    if(nb==4){
      h_N_BTagsM_4b->Fill(N_BTagsM,Weight);
      h_N_BTagsL_4b->Fill(N_BTagsL,Weight);
      h_N_BTagsT_4b->Fill(N_BTagsT,Weight);     
    }
    
  }
  // Do the actual analysis
  prince->MakePrincipals();
  
  // Print out the result on
  prince->Print();

  // Test the PCA 
  prince->Test();
  prince->MakeCode();
  TBrowser* b = new TBrowser("principalBrowser", prince);
  // write
  TFile* outfile=new TFile(outfilename,"RECREATE");    
  outfile->cd();
  h_M_Higgs_reco->Scale(1./h_M_Higgs_reco->Integral("width"));
  h_M_Higgs_reco->Write();
  h_M_TopHad_reco->Scale(1./h_M_TopHad_reco->Integral("width"));
  h_M_TopHad_reco->Write();
  h_M_TopLep_reco->Scale(1./h_M_TopLep_reco->Integral("width"));
  h_M_TopLep_reco->Write();
  h_M_WHad_reco->Scale(1./h_M_WHad_reco->Integral("width"));
  h_M_WHad_reco->Write();
  h_M_WHad_vs_M_TopHad_reco->Scale(1./h_M_WHad_vs_M_TopHad_reco->Integral("width"));
  h_M_WHad_vs_M_TopHad_reco->Write();
  h_M_WHad_vs_M_TopHad_pca_reco->Scale(1./h_M_WHad_vs_M_TopHad_pca_reco->Integral("width"));
  h_M_WHad_vs_M_TopHad_pca_reco->Write();
  h_M_TopHad_pca_reco->Scale(1./h_M_TopHad_pca_reco->Integral("width"));
  h_M_TopHad_pca_reco->Write();
  h_M_WHad_pca_reco->Scale(1./h_M_WHad_pca_reco->Integral("width"));
  h_M_WHad_pca_reco->Write();

  h_M_WLep_reco->Scale(1./h_M_WLep_reco->Integral("width"));
  h_M_WLep_reco->Write();
  
  h_M_Higgs_recoall->Scale(1./h_M_Higgs_recoall->Integral("width"));
  h_M_Higgs_recoall->Write();
  h_M_TopHad_recoall->Scale(1./h_M_TopHad_recoall->Integral("width"));
  h_M_TopHad_recoall->Write();
  h_M_TopLep_recoall->Scale(1./h_M_TopLep_recoall->Integral("width"));
  h_M_TopLep_recoall->Write();
  h_M_WHad_recoall->Scale(1./h_M_WHad_recoall->Integral("width"));
  h_M_WHad_recoall->Write();
  h_M_WHad_vs_M_TopHad_recoall->Scale(1./h_M_WHad_vs_M_TopHad_recoall->Integral("width"));
  h_M_WHad_vs_M_TopHad_recoall->Write();
  h_M_WLep_recoall->Scale(1./h_M_WLep_recoall->Integral("width"));
  h_M_WLep_recoall->Write();

  h_M_Total_recoall->Scale(1./h_M_Total_recoall->Integral("width"));
  h_M_Total_recoall->Write();
  h_Pt_Total_recoall->Scale(1./h_Pt_Total_recoall->Integral("width"));
  h_Pt_Total_recoall->Write();

  h_M_Higgs_best->Scale(1./h_M_Higgs_best->Integral("width"));
  h_M_Higgs_best->Write();
  h_M_TopHad_best->Scale(1./h_M_TopHad_best->Integral("width"));
  h_M_TopHad_best->Write();
  h_M_TopHad_best_special->Scale(1./h_M_TopHad_best_special->Integral("width"));
  h_M_TopHad_best_special->Write();
  h_M_TopLep_best->Scale(1./h_M_TopLep_best->Integral("width"));
  h_M_TopLep_best->Write();
  h_M_WHad_best->Scale(1./h_M_WHad_best->Integral("width"));
  h_M_WHad_best->Write();
  h_M_WHad_vs_M_TopHad_best->Scale(1./h_M_WHad_vs_M_TopHad_best->Integral("width"));
  h_M_WHad_vs_M_TopHad_best->Write();
  h_M_WLep_best->Scale(1./h_M_WLep_best->Integral("width"));
  h_M_WLep_best->Write();
  h_CSV_Higgs_best->Scale(1./h_CSV_Higgs_best->Integral("width"));
  h_CSV_Higgs_best->Write();
  h_CSV_WHad_best->Scale(1./h_CSV_WHad_best->Integral("width"));
  h_CSV_WHad_best->Write();
  h_CSV_TopHad_best->Scale(1./h_CSV_TopHad_best->Integral("width"));
  h_CSV_TopHad_best->Write();
  h_CSV_TopLep_best->Scale(1./h_CSV_TopLep_best->Integral("width"));
  h_CSV_TopLep_best->Write();


  h_M_Higgs_all->Scale(1./h_M_Higgs_all->Integral("width"));
  h_M_Higgs_all->Write();
  h_M_TopHad_all->Scale(1./h_M_TopHad_all->Integral("width"));
  h_M_TopHad_all->Write();
  h_M_TopLep_all->Scale(1./h_M_TopLep_all->Integral("width"));
  h_M_TopLep_all->Write();
  h_M_WHad_all->Scale(1./h_M_WHad_all->Integral("width"));
  h_M_WHad_all->Write();
  h_M_WHad_vs_M_TopHad_all->Scale(1./h_M_WHad_vs_M_TopHad_all->Integral("width"));
  h_M_WHad_vs_M_TopHad_all->Write();
  h_M_WLep_all->Scale(1./h_M_WLep_all->Integral("width"));
  h_M_WLep_all->Write();
  h_M_TopHad_pca_all->Scale(1./h_M_TopHad_pca_all->Integral("width"));
  h_M_TopHad_pca_all->Write();
  h_M_WHad_pca_all->Scale(1./h_M_WHad_pca_all->Integral("width"));
  h_M_WHad_pca_all->Write();
  h_M_WHad_vs_M_TopHad_pca_all->Scale(1./h_M_WHad_vs_M_TopHad_pca_all->Integral("width"));
  h_M_WHad_vs_M_TopHad_pca_all->Write();


  h_M_Higgs_parton->Scale(1./h_M_Higgs_parton->Integral("width"));
  h_M_Higgs_parton->Write();
  h_M_TopHad_parton->Scale(1./h_M_TopHad_parton->Integral("width"));
  h_M_TopHad_parton->Write();
  h_M_TopLep_parton->Scale(1./h_M_TopLep_parton->Integral("width"));
  h_M_TopLep_parton->Write();
  h_M_WHad_parton->Scale(1./h_M_WHad_parton->Integral("width")); 
  h_M_WHad_parton->Write();
  h_M_WHad_vs_M_TopHad_parton->Scale(1./h_M_WHad_vs_M_TopHad_parton->Integral("width"));
  h_M_WHad_vs_M_TopHad_parton->Write();
  h_M_WLep_parton->Scale(1./h_M_WLep_parton->Integral("width"));
  h_M_WLep_parton->Write();

  h_CSV_b->Scale(1./h_CSV_b->Integral("width"));
  h_CSV_b->Write();
  h_CSV_c->Scale(1./h_CSV_c->Integral("width"));
  h_CSV_c->Write();
  h_CSV_l_w_c->Scale(1./h_CSV_l_w_c->Integral("width"));
  h_CSV_l_w_c->Write();
  h_CSV_l_wo_c->Scale(1./h_CSV_l_wo_c->Integral("width"));
  h_CSV_l_wo_c->Write();

  h_N_BTagsM->Scale(1./h_N_BTagsM->Integral("width"));
  h_N_BTagsM->Write();
  h_N_BTagsL->Scale(1./h_N_BTagsL->Integral("width"));
  h_N_BTagsL->Write();
  h_N_BTagsT->Scale(1./h_N_BTagsT->Integral("width"));
  h_N_BTagsT->Write();
  h_N_BTagsM_2b->Scale(1./h_N_BTagsM_2b->Integral("width"));
  h_N_BTagsM_2b->Write();
  h_N_BTagsL_2b->Scale(1./h_N_BTagsL_2b->Integral("width"));
  h_N_BTagsL_2b->Write();
  h_N_BTagsT_2b->Scale(1./h_N_BTagsT_2b->Integral("width"));
  h_N_BTagsT_2b->Write();
  h_N_BTagsM_3b->Scale(1./h_N_BTagsM_3b->Integral("width"));
  h_N_BTagsM_3b->Write();
  h_N_BTagsL_3b->Scale(1./h_N_BTagsL_3b->Integral("width"));
  h_N_BTagsL_3b->Write();
  h_N_BTagsT_3b->Scale(1./h_N_BTagsT_3b->Integral("width"));
  h_N_BTagsT_3b->Write();
  h_N_BTagsM_4b->Scale(1./h_N_BTagsM_4b->Integral("width"));
  h_N_BTagsM_4b->Write();
  h_N_BTagsL_4b->Scale(1./h_N_BTagsL_4b->Integral("width"));
  h_N_BTagsL_4b->Write();
  h_N_BTagsT_4b->Scale(1./h_N_BTagsT_4b->Integral("width"));
  h_N_BTagsT_4b->Write();

  h_nb->Scale(1./h_nb->Integral("width"));
  h_nb->Write();
  
  h_dr2->Scale(1./h_dr2->Integral("width"));
  h_dr2->Write();
  h_B1_dr->Scale(1./h_B1_dr->Integral("width"));
  h_B1_dr->Write();
  h_B2_dr->Scale(1./h_B2_dr->Integral("width"));
  h_B2_dr->Write();
  h_BLep_dr->Scale(1./h_BLep_dr->Integral("width"));
  h_BLep_dr->Write();
  h_BHad_dr->Scale(1./h_BHad_dr->Integral("width"));
  h_BHad_dr->Write();
  h_Q1_dr->Scale(1./h_Q1_dr->Integral("width"));
  h_Q1_dr->Write();
  h_Q2_dr->Scale(1./h_Q2_dr->Integral("width"));
  h_Q2_dr->Write();
  h_Nu_dr->Scale(1./h_Nu_dr->Integral("width"));
  h_Nu_dr->Write();
  h_Lep_dr->Scale(1./h_Lep_dr->Integral("width"));
  h_Lep_dr->Write();

  cout << "all " << nAll << ", matchable " << nMatchable << " (" << nMatchable/(float)nAll << ")"<<endl;

}

int main(){
  plot();
}

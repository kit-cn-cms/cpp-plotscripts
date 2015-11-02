#include "TChain.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;
// alias for long object names
typedef vector<TLorentzVector> LVs;
typedef TLorentzVector LV;

// helper functions to get LorentzVectors of jets
TLorentzVector getLV(float pt,float eta,float phi,float e){
  TLorentzVector v;
  v.SetPtEtaPhiE(pt,eta,phi,e);
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


void plot(){
  TH1F::SetDefaultSumw2();
  
  // open files
  TChain* chain = new TChain("MVATree");
  char* filenames = getenv ("FILENAMES");
  char* outfilename = getenv ("OUTFILENAME");
  string buf;
  stringstream ss(filenames); 
  while (ss >> buf){
    chain->Add(buf.c_str());
  }
  chain->SetBranchStatus("*",0);

  TFile* outfile=new TFile(outfilename,"RECREATE");    

  // initialize variables from tree
  float Weight;
  chain->SetBranchAddress("Weight",&Weight);
  int N_Jets;
  chain->SetBranchAddress("N_Jets",&N_Jets);
  int N_BTagsM;
  chain->SetBranchAddress("N_BTagsM",&N_BTagsM);
  int N_LooseLeptons;
  chain->SetBranchAddress("N_LooseLeptons",&N_LooseLeptons);
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

  // define histos
  TH1F* h_N_Jets_4j2t=new TH1F("N_Jets_4j2t","Number of jets (4j2t)",15,-0.5,14.5);
  TH1F* h_Pt_Jets_4j2t=new TH1F("Pt_Jets_4j2t","Pt of all jets (4j2t)",60,0,600);
  TH1F* h_Pt_Jet4_4j2t=new TH1F("Pt_Jet4_4j2t","Pt of fourth jet (4j2t)",60,0,300);
  TH1F* h_Dr_Jet1_Jet2_4j2t=new TH1F("Dr_Jet1_Jet2_4j2t","DeltaR(jet1,jet2) (4j2t)",60,0,6);

  // loop over all events
  long nentries = chain->GetEntries(); 
  cout << "total number of events: " << nentries << endl;
  for (long iEntry=0;iEntry<nentries;iEntry++) {
    if(iEntry%10000==0) cout << "analyzing event " << iEntry << endl;
    chain->GetEntry(iEntry); 
    // start lepton + 4j2t selection
    if( (N_BTagsM>1)&&(N_Jets>=4)&&(N_TightLeptons==1)&&(N_LooseLeptons==1)) {


      // Fill histograms
      h_N_Jets_4j2t->Fill(N_Jets,Weight);
      for(int i=0; i<N_Jets;i++){
	h_Pt_Jets_4j2t->Fill(Jet_Pt[i]);
      }
      if(N_Jets>3){
	h_Pt_Jet4_4j2t->Fill(Jet_Pt[3]);
      }
      LV v_jet1;
      LV v_jet2;
      if(N_Jets>0){
	v_jet1=getLV(Jet_Pt[0],Jet_Eta[0],Jet_Phi[0],Jet_E[0]);
      }
      if(N_Jets>1){
	v_jet2=getLV(Jet_Pt[1],Jet_Eta[1],Jet_Phi[1],Jet_E[1]);
      }
      float dr_j1_j2 = v_jet1.DeltaR(v_jet2);
      h_Dr_Jet1_Jet2_4j2t->Fill(dr_j1_j2);
    }
      // end of lepton + 4j2t selection

  } // end of event loop
  // write histos
  outfile->Write();
}

int main(){
  plot();
}

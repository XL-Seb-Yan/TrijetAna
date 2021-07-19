//MC selection
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector2.h>               // 2D vector class
#include <TMath.h>                  // ROOT math library
#include <vector>                   // STL vector class
#include <utility>
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <bitset>
#include "TLorentzVector.h"         // 4-vector class
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TF1.h>
#include <TStopwatch.h>
#include "TH1D.h"
#include "TRandom3.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
// C++ tool
#include <algorithm>
#include <map>
#endif

#define Run 2017

using namespace std;

void select_ML_Fullmatch(const TString samplename="ZprimeTo3Gluon_TuneCP5_13TeV_pythia8", 
								const TString outputDir="ntuples", 
								const int nEvents = -1, 
								const int year = 2017) {
	
	gBenchmark->Start("selectTrijet");
	gROOT->SetBatch(1);
	UInt_t count0=0, count1=0, count2=0, count3=0, count4=0, count5=0, count6=0, count7=0, count8=0, count9=0;

	//--------------------------------------------------------------------------------------------------------------
	// Main analysis code 
	//==============================================================================================================  
	
	// Create output directory
	gSystem->mkdir(outputDir,kTRUE);
	
	TH1F *h_njet = new TH1F("njet","njet",15,0,15);
	TH1F *h_dR_gq0 = new TH1F("dR0","dR gq0",100,0,1);
	TH1F *h_dR_gq1 = new TH1F("dR1","dR gq1",100,0,1);
	TH1F *h_dR_gq2 = new TH1F("dR2","dR gq2",100,0,1);
	
	// Data Structure for output skimmed files
	// Single jet
	float jet_pt_0, jet_eta_0, jet_phi_0, jet_m_0;
	float jet_pt_1, jet_eta_1, jet_phi_1, jet_m_1;
	float jet_pt_2, jet_eta_2, jet_phi_2, jet_m_2;
	// Dijet
	float jj_m_01, jj_m_02, jj_m_12;
	float jj_pt_01, jj_pt_02, jj_pt_12;
	float jj_eta_01, jj_eta_02, jj_eta_12;
	float jj_dR_01, jj_dR_02, jj_dR_12, jj_dR_min, jj_dR_max;
	float jj_dEta_01, jj_dEta_02, jj_dEta_12, jj_dEta_min, jj_dEta_max;
	// Dijet-jet
	float jj_j_dR_01_2, jj_j_dR_02_1, jj_j_dR_12_0, jj_j_dR_min;
	float jj_j_dEtaAbs_01_2, jj_j_dEtaAbs_02_1, jj_j_dEtaAbs_12_0, jj_j_dEtaAbs_max;
	float jj_j_dPhi_01_2, jj_j_dPhi_02_1, jj_j_dPhi_12_0, jj_j_dPhi_min;
	// System
	float M_jjj, pt_jjj, eta_jjj, jet_ptoverM_0, jet_ptoverM_1, jet_ptoverM_2, jj_ptoverM_01, jj_ptoverM_02, jj_ptoverM_12, jj_moverM_01, jj_moverM_02, jj_moverM_12;	   
	// GEN matching
	int gen_dijet_matched;
	// evt info
	UInt_t run_num;
	ULong64_t evt_num;
	int lumi_block;
	int flag_mass;
	float Zprime_mass;
	
	TString outfilename = outputDir + TString("/") + samplename + TString("_ML_study.root");
	TFile *outFile = new TFile(outfilename,"RECREATE"); 
	TTree *outTree = new TTree("Events","Events");
	outTree->Branch("jet_pt_0",       &jet_pt_0);
	outTree->Branch("jet_eta_0",      &jet_eta_0);
	outTree->Branch("jet_phi_0",      &jet_phi_0);
	outTree->Branch("jet_m_0",        &jet_m_0);
	outTree->Branch("jet_pt_1",       &jet_pt_1);
	outTree->Branch("jet_eta_1",      &jet_eta_1);
	outTree->Branch("jet_phi_1",      &jet_phi_1);
	outTree->Branch("jet_m_1",        &jet_m_1);
	outTree->Branch("jet_pt_2",       &jet_pt_2);
	outTree->Branch("jet_eta_2",      &jet_eta_2);
	outTree->Branch("jet_phi_2",      &jet_phi_2);
	outTree->Branch("jet_m_2",        &jet_m_2);
	
	outTree->Branch("jj_m_01",        &jj_m_01);
	outTree->Branch("jj_m_02",        &jj_m_02);
	outTree->Branch("jj_m_12",        &jj_m_12);
	outTree->Branch("jj_pt_01",       &jj_pt_01);
	outTree->Branch("jj_pt_02",       &jj_pt_02);
	outTree->Branch("jj_pt_12",       &jj_pt_12);
	outTree->Branch("jj_eta_01",      &jj_eta_01);
	outTree->Branch("jj_eta_02",      &jj_eta_02);
	outTree->Branch("jj_eta_12",      &jj_eta_12);
	outTree->Branch("jj_dR_01",       &jj_dR_01);
	outTree->Branch("jj_dR_02",       &jj_dR_02);
	outTree->Branch("jj_dR_12",       &jj_dR_12);
	outTree->Branch("jj_dR_max",      &jj_dR_max);
	outTree->Branch("jj_dR_min",      &jj_dR_min);
	outTree->Branch("jj_dEta_01",     &jj_dEta_01);
	outTree->Branch("jj_dEta_02",     &jj_dEta_02);
	outTree->Branch("jj_dEta_12",     &jj_dEta_12);
	outTree->Branch("jj_dEta_max",    &jj_dEta_max);
	outTree->Branch("jj_dEta_min",    &jj_dEta_min);
	
	outTree->Branch("jj_j_dR_01_2",   &jj_j_dR_01_2);
	outTree->Branch("jj_j_dR_02_1",   &jj_j_dR_02_1);
	outTree->Branch("jj_j_dR_12_0",   &jj_j_dR_12_0);
	outTree->Branch("jj_j_dR_min",   &jj_j_dR_min);
	outTree->Branch("jj_j_dEtaAbs_01_2",   &jj_j_dEtaAbs_01_2);
	outTree->Branch("jj_j_dEtaAbs_02_1",   &jj_j_dEtaAbs_02_1);
	outTree->Branch("jj_j_dEtaAbs_12_0",   &jj_j_dEtaAbs_12_0);
	outTree->Branch("jj_j_dEtaAbs_max",   &jj_j_dEtaAbs_max);
	outTree->Branch("jj_j_dPhi_01_2",   &jj_j_dPhi_01_2);
	outTree->Branch("jj_j_dPhi_02_1",   &jj_j_dPhi_02_1);
	outTree->Branch("jj_j_dPhi_12_0",   &jj_j_dPhi_12_0);
	outTree->Branch("jj_j_dPhi_min",   &jj_j_dPhi_min);
	
	outTree->Branch("M_jjj",          &M_jjj);
	outTree->Branch("pt_jjj",         &pt_jjj);
	outTree->Branch("eta_jjj",        &eta_jjj);
	outTree->Branch("jet_ptoverM_0",  &jet_ptoverM_0);
	outTree->Branch("jet_ptoverM_1",  &jet_ptoverM_1);
	outTree->Branch("jet_ptoverM_2",  &jet_ptoverM_2);
	outTree->Branch("jj_ptoverM_01",  &jj_ptoverM_01);
	outTree->Branch("jj_ptoverM_02",  &jj_ptoverM_02);
	outTree->Branch("jj_ptoverM_12",  &jj_ptoverM_12);
	outTree->Branch("jj_moverM_01",   &jj_moverM_01);
	outTree->Branch("jj_moverM_02",   &jj_moverM_02);
	outTree->Branch("jj_moverM_12",   &jj_moverM_12);
	
	outTree->Branch("gen_dijet_matched",   &gen_dijet_matched);
	outTree->Branch("run_num",        &run_num);
	outTree->Branch("evt_num",        &evt_num);
	outTree->Branch("lumi_block",     &lumi_block);
	outTree->Branch("flag_mass",      &flag_mass);
	outTree->Branch("Zprime_mass",    &Zprime_mass);
  
	cout<<"begin loop over files"<<endl;
	TStopwatch stopwatch;
	
	// loop through files
	TTree* eventTree = 0;
	ifstream insample(samplename+TString(".txt"));
	std::string line;
	while (std::getline(insample, line)){
		TString file_name(line);
		
		// Read input file and get the TTrees
		cout << "Processing " << file_name <<endl; cout.flush();
		TFile *infile = TFile::Open(file_name,"READ");
		assert(infile);

		// Access Event Tree
		TTreeReader fReader; 		//!the tree reader
		TTreeReaderValue<UInt_t> run = {fReader, "run"};
    TTreeReaderValue<UInt_t> luminosityBlock = {fReader, "luminosityBlock"};
    TTreeReaderValue<ULong64_t> event = {fReader, "event"};
#if Run == 2016
		//--HLT trigger--
		if(samplename.Contains("2016B2")){
			TTreeReaderValue<Bool_t> HLT_PFHT800 = {fReader, "HLT_PFHT800"};
			TTreeReaderValue<Bool_t> HLT_PFHT900 = {fReader, "HLT_PFHT900"};
			TTreeReaderValue<Bool_t> HHLT_PFJet500 = {fReader, "HLT_PFJet500"};
			TTreeReaderValue<Bool_t> HLT_CaloJet500_NoJetID = {fReader, "HLT_CaloJet500_NoJetID"};
		}
		else if(samplename.Contains("2016H")){
			TTreeReaderValue<Bool_t> HLT_PFHT900 = {fReader, "HLT_PFHT900"};
			TTreeReaderValue<Bool_t> HLT_AK8PFJet450 = {fReader, "HLT_AK8PFJet450"};
			TTreeReaderValue<Bool_t> HLT_AK8PFJet500 = {fReader, "HLT_AK8PFJet500"};
			TTreeReaderValue<Bool_t> HLT_PFJet500 = {fReader, "HLT_PFJet500"};
			TTreeReaderValue<Bool_t> HLT_CaloJet500_NoJetID = {fReader, "HLT_CaloJet500_NoJetID"};
		}
		else{
			TTreeReaderValue<Bool_t> HLT_PFHT800 = {fReader, "HLT_PFHT800"};
			TTreeReaderValue<Bool_t> HLT_PFHT900 = {fReader, "HLT_PFHT900"};
			TTreeReaderValue<Bool_t> HLT_AK8PFJet450 = {fReader, "HLT_AK8PFJet450"};
			TTreeReaderValue<Bool_t> HLT_AK8PFJet500 = {fReader, "HLT_AK8PFJet500"};
			TTreeReaderValue<Bool_t> HLT_PFJet500 = {fReader, "HLT_PFJet500"};
			TTreeReaderValue<Bool_t> HLT_CaloJet500_NoJetID = {fReader, "HLT_CaloJet500_NoJetID"};
		}
#else
		cout<<"====== Using 2017/8 HLT set ======="<<endl;
		TTreeReaderValue<Bool_t> HLT_PFHT1050 = {fReader, "HLT_PFHT1050"};
		TTreeReaderValue<Bool_t> HLT_AK8PFJet550 = {fReader, "HLT_AK8PFJet550"};
		TTreeReaderValue<Bool_t> HLT_AK8PFJet500 = {fReader, "HLT_AK8PFJet500"};
		TTreeReaderValue<Bool_t> HLT_PFJet500 = {fReader, "HLT_PFJet500"};
		TTreeReaderValue<Bool_t> HLT_CaloJet500_NoJetID = {fReader, "HLT_CaloJet500_NoJetID"};
		TTreeReaderValue<Bool_t> HLT_CaloJet550_NoJetID = {fReader, "HLT_CaloJet550_NoJetID"};
#endif
		
		TTreeReaderValue<UInt_t> nJet = {fReader, "nJet"};
		TTreeReaderArray<Float_t> Jet_pt = {fReader, "Jet_pt"};
		TTreeReaderArray<Float_t> Jet_eta = {fReader, "Jet_eta"};
		TTreeReaderArray<Float_t> Jet_phi = {fReader, "Jet_phi"};
		TTreeReaderArray<Float_t> Jet_mass = {fReader, "Jet_mass"};
		TTreeReaderArray<Int_t> Jet_jetId = {fReader, "Jet_jetId"};
		
		// MC Truth
		// TTreeReaderValue<UInt_t> nGenJet = {fReader, "nGenJet"};
		// TTreeReaderArray<Float_t> GenJet_eta = {fReader, "GenJet_eta"};
		// TTreeReaderArray<Float_t> GenJet_mass = {fReader, "GenJet_mass"};
		// TTreeReaderArray<Float_t> GenJet_phi = {fReader, "GenJet_phi"};
		// TTreeReaderArray<Float_t> GenJet_pt = {fReader, "GenJet_pt"};
		// TTreeReaderArray<Int_t> Jet_genJetIdx = {fReader, "Jet_genJetIdx"};
		TTreeReaderValue<UInt_t> nGenPart = {fReader, "nGenPart"};
		TTreeReaderArray<Float_t> GenPart_eta = {fReader, "GenPart_eta"};
		TTreeReaderArray<Float_t> GenPart_mass = {fReader, "GenPart_mass"};
		TTreeReaderArray<Float_t> GenPart_phi = {fReader, "GenPart_phi"};
		TTreeReaderArray<Float_t> GenPart_pt = {fReader, "GenPart_pt"};
		TTreeReaderArray<Int_t> GenPart_genPartIdxMother = {fReader, "GenPart_genPartIdxMother"};
		TTreeReaderArray<Int_t> GenPart_pdgId = {fReader, "GenPart_pdgId"};
		TTreeReaderArray<Int_t> GenPart_status = {fReader, "GenPart_status"};
		
		TTreeReader fReader_tmp;
		TTreeReaderValue<bool> GenModel_m500 = {fReader_tmp, "GenModel_m500"};
		TTreeReaderValue<bool> GenModel_m750 = {fReader_tmp, "GenModel_m750"};
		TTreeReaderValue<bool> GenModel_m1000 = {fReader_tmp, "GenModel_m1000"};
		TTreeReaderValue<bool> GenModel_m1500 = {fReader_tmp, "GenModel_m1500"};
		TTreeReaderValue<bool> GenModel_m2000 = {fReader_tmp, "GenModel_m2000"};
		TTreeReaderValue<bool> GenModel_m2500 = {fReader_tmp, "GenModel_m2500"};
		TTreeReaderValue<bool> GenModel_m3000 = {fReader_tmp, "GenModel_m3000"};
		TTreeReaderValue<bool> GenModel_m3500 = {fReader_tmp, "GenModel_m3500"};
		TTreeReaderValue<bool> GenModel_m4000 = {fReader_tmp, "GenModel_m4000"};
		TTreeReaderValue<bool> GenModel_m4500 = {fReader_tmp, "GenModel_m4500"};
		TTreeReaderValue<bool> GenModel_m5000 = {fReader_tmp, "GenModel_m5000"};
		TTreeReaderValue<bool> GenModel_m6000 = {fReader_tmp, "GenModel_m6000"};
		TTreeReaderValue<bool> GenModel_m7000 = {fReader_tmp, "GenModel_m7000"};
		TTreeReaderValue<bool> GenModel_m8000 = {fReader_tmp, "GenModel_m8000"};
		TTreeReaderValue<bool> GenModel_m9000 = {fReader_tmp, "GenModel_m9000"};
		
		eventTree = (TTree*)infile->Get("Events");
		assert(eventTree);
		
		cout<<"OK"<<endl;
		
		if(eventTree->GetListOfBranches()->FindObject("GenModel_m500")){GenModel_m500 = {fReader, "GenModel_m500"};}
		if(eventTree->GetListOfBranches()->FindObject("GenModel_m750")){GenModel_m750 = {fReader, "GenModel_m750"};}
		if(eventTree->GetListOfBranches()->FindObject("GenModel_m1000")){GenModel_m1000 = {fReader, "GenModel_m1000"};}
		if(eventTree->GetListOfBranches()->FindObject("GenModel_m1500")){GenModel_m1500 = {fReader, "GenModel_m1500"};}
		if(eventTree->GetListOfBranches()->FindObject("GenModel_m2000")){GenModel_m2000 = {fReader, "GenModel_m2000"};}
		if(eventTree->GetListOfBranches()->FindObject("GenModel_m2500")){GenModel_m2500 = {fReader, "GenModel_m2500"};}
		if(eventTree->GetListOfBranches()->FindObject("GenModel_m3000")){GenModel_m3000 = {fReader, "GenModel_m3000"};}
		if(eventTree->GetListOfBranches()->FindObject("GenModel_m3500")){GenModel_m3500 = {fReader, "GenModel_m3500"};}
		if(eventTree->GetListOfBranches()->FindObject("GenModel_m4000")){GenModel_m4000 = {fReader, "GenModel_m4000"};}
		if(eventTree->GetListOfBranches()->FindObject("GenModel_m4500")){GenModel_m4500 = {fReader, "GenModel_m4500"};}
		if(eventTree->GetListOfBranches()->FindObject("GenModel_m5000")){GenModel_m5000 = {fReader, "GenModel_m5000"};}
		if(eventTree->GetListOfBranches()->FindObject("GenModel_m6000")){GenModel_m6000 = {fReader, "GenModel_m6000"};}
		if(eventTree->GetListOfBranches()->FindObject("GenModel_m7000")){GenModel_m7000 = {fReader, "GenModel_m7000"};}
		if(eventTree->GetListOfBranches()->FindObject("GenModel_m8000")){GenModel_m8000 = {fReader, "GenModel_m8000"};}
		if(eventTree->GetListOfBranches()->FindObject("GenModel_m9000")){GenModel_m9000 = {fReader, "GenModel_m9000"};}
		
		fReader.SetTree(eventTree);
		
		for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++){
		//for(UInt_t ientry=0; ientry<10; ientry++){
			count0++;
		
			fReader.SetLocalEntry(ientry);
			
			// Evt info
			run_num = *run;
			evt_num = *event;
			lumi_block = *luminosityBlock;
			std::bitset<15> massbits;
			massbits.set(0, eventTree->GetListOfBranches()->FindObject("GenModel_m500")? *GenModel_m500 : 0);
			massbits.set(1, eventTree->GetListOfBranches()->FindObject("GenModel_m750")? *GenModel_m750 : 0);
			massbits.set(2, eventTree->GetListOfBranches()->FindObject("GenModel_m1000")? *GenModel_m1000 : 0);
			massbits.set(3, eventTree->GetListOfBranches()->FindObject("GenModel_m1500")? *GenModel_m1500 : 0);
			massbits.set(4, eventTree->GetListOfBranches()->FindObject("GenModel_m2000")? *GenModel_m2000 : 0);
			massbits.set(5, eventTree->GetListOfBranches()->FindObject("GenModel_m2500")? *GenModel_m2500 : 0);
			massbits.set(6, eventTree->GetListOfBranches()->FindObject("GenModel_m3000")? *GenModel_m3000 : 0);
			massbits.set(7, eventTree->GetListOfBranches()->FindObject("GenModel_m3500")? *GenModel_m3500 : 0);
			massbits.set(8, eventTree->GetListOfBranches()->FindObject("GenModel_m4000")? *GenModel_m4000 : 0);
			massbits.set(9, eventTree->GetListOfBranches()->FindObject("GenModel_m4500")? *GenModel_m4500 : 0);
			massbits.set(10, eventTree->GetListOfBranches()->FindObject("GenModel_m5000")? *GenModel_m5000 : 0);
			massbits.set(11, eventTree->GetListOfBranches()->FindObject("GenModel_m6000")? *GenModel_m6000 : 0);
			massbits.set(12, eventTree->GetListOfBranches()->FindObject("GenModel_m7000")? *GenModel_m7000 : 0);
			massbits.set(13, eventTree->GetListOfBranches()->FindObject("GenModel_m8000")? *GenModel_m8000 : 0);
			massbits.set(14, eventTree->GetListOfBranches()->FindObject("GenModel_m9000")? *GenModel_m9000 : 0);
			if(massbits.count() != 1) cout<<"Massbit is in error!"<<(int)(massbits.to_ulong())<<endl;
			flag_mass = (int)(massbits.to_ulong());
			
			// HLT selection
			bool passHLT = false;
#if Run == 2016
			if(samplename.Contains("2016B2")){
				passHLT = (*HLT_PFHT800 || *HLT_PFHT900 || *HLT_PFJet500 || *HLT_CaloJet500_NoJetID);
			}
			else if(samplename.Contains("2016H")){
				passHLT = (*HLT_PFHT900 || *HLT_AK8PFJet450 || *HLT_AK8PFJet500 || *HLT_PFJet500 || *HLT_CaloJet500_NoJetID);
			}
			else{
					passHLT = (*HLT_PFHT800 || *HLT_PFHT900 || *HLT_AK8PFJet450 || *HLT_AK8PFJet500 || *HLT_PFJet500 || *HLT_CaloJet500_NoJetID);
			}
#else
			passHLT = (*HLT_PFHT1050 || *HLT_AK8PFJet550 || *HLT_AK8PFJet500 || *HLT_PFJet500 || *HLT_CaloJet500_NoJetID || *HLT_CaloJet550_NoJetID);
#endif
			if (!passHLT) continue;
			count1++;
			
			// Do GEN-RECO matching
			std::vector<TLorentzVector> gen_gq_list; //store the 4-vector of gen level q or g
			std::vector<UInt_t> gen_jet_index; //store the index of resulted gen jets of the decay of R2
			TLorentzVector gen_gq_temp;
			TLorentzVector gen_jet_temp;
			for(UInt_t i=0; i<*nGenPart; i++){
				if(GenPart_pdgId[i] == 9900023)
					Zprime_mass = GenPart_mass[i];
				if(GenPart_pdgId[i] != 21 || GenPart_status[i] != 23) continue;// This selects the 3 gluons we need
				//cout<<"Index: "<<i<<" Pt: "<<GenPart_pt[i]<<" Eta: "<<GenPart_eta[i]<<" Phi: "<<GenPart_phi[i]<<" Mass: "<<GenPart_mass[i]<<" ID: "<<GenPart_pdgId[i]<<" Status: "<<GenPart_status[i]<<" Momindex: "<<GenPart_genPartIdxMother[i]<<" "<<GenPart_pdgId[GenPart_genPartIdxMother[i]]<<endl;
				gen_gq_temp.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
				gen_gq_list.push_back(gen_gq_temp); //The last two 4-vector correspond to the decay production of R2
			}
			if(gen_gq_list.size() != 3) 
				cout<<"Warning: gen_gq_list size is not 3!"<<endl;
			
			// Jet skimming
			std::vector<int> j_sel_index_arr; //store the true indices of RECO jets selected
			std::vector<TLorentzVector> j_sel_arr;
			for(UInt_t i=0; i<*nJet; i++){
				if(Jet_pt[i] < 100) continue;
				if(abs(Jet_eta[i]) > 2.5) continue;
				if(Jet_jetId[i] < 6) continue;
				TLorentzVector vj;
				vj.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
				j_sel_index_arr.push_back(i);
				j_sel_arr.push_back(vj);
			}
			h_njet->Fill(j_sel_index_arr.size());
			if(j_sel_index_arr.size() < 3) continue;
			count2++;
			
			// Jet matching
			std::vector<int> matched_index_arr; //store the indices (within j_sel_arr) of RECO jets which are matched to the GEN particle
			matched_index_arr.push_back(-1);
			matched_index_arr.push_back(-1);
			matched_index_arr.push_back(-1);
			std::vector<float> mindR_arr; //store the dR of RECO jets which are matched to the GEN particle
			mindR_arr.push_back(99);
			mindR_arr.push_back(99);
			mindR_arr.push_back(99);
			for(UInt_t ijet = 0; ijet < j_sel_arr.size(); ijet++){
				TLorentzVector vj = j_sel_arr[ijet];
				int jet_index_orig = j_sel_index_arr[ijet];
				std::vector<float> dR_arr;
				dR_arr.push_back(vj.DeltaR(gen_gq_list[0]));
				dR_arr.push_back(vj.DeltaR(gen_gq_list[1]));
				dR_arr.push_back(vj.DeltaR(gen_gq_list[2]));
				int minGlu_index = std::min_element(dR_arr.begin(),dR_arr.end()) - dR_arr.begin(); // tell this RECO jets best match to which gluon
				float min = dR_arr.at(minGlu_index); // tell the matching dR of the current best match
				if(min < mindR_arr[minGlu_index]){
					if(matched_index_arr[minGlu_index] != -1){ //reconsider last best matchee jets to other gluons
						j_sel_arr.push_back(j_sel_arr[matched_index_arr[minGlu_index]]);
						j_sel_index_arr.push_back(j_sel_index_arr[matched_index_arr[minGlu_index]]);
					}
					matched_index_arr[minGlu_index] = ijet;
					mindR_arr[minGlu_index] = min;
				}
			}		
			
			// for(int i=0; i<3 ;i++){
				// if(matching_index_arr[i] > 2){
					// count3++;
					// for(int j=0; j<j_sel_arr.size(); j++){
						// cout<<"Jet "<<j<<": pt: "<<j_sel_arr[j].Pt()<<" eta: "<<j_sel_arr[j].Eta()<<" phi: "<<j_sel_arr[j].Phi()<<endl;
						// if(j > 4)
							// break;
					// }
					
					// cout<<"Parton 0: pt: "<<gen_gq_list[0].Pt()<<" eta: "<<gen_gq_list[0].Eta()<<" phi: "<<gen_gq_list[0].Phi()<<" is matched to "<<matching_index_arr[0]<<" with dR = "<<mindR_arr[0]<<endl;
					// cout<<"Parton 1: pt: "<<gen_gq_list[1].Pt()<<" eta: "<<gen_gq_list[1].Eta()<<" phi: "<<gen_gq_list[1].Phi()<<" is matched to "<<matching_index_arr[1]<<" with dR = "<<mindR_arr[1]<<endl;
					// cout<<"Parton 2: pt: "<<gen_gq_list[2].Pt()<<" eta: "<<gen_gq_list[2].Eta()<<" phi: "<<gen_gq_list[2].Phi()<<" is matched to "<<matching_index_arr[2]<<" with dR = "<<mindR_arr[2]<<endl;
					// break;
				// }
			// }
			// cout<<endl;
			// h_dR_gq0->Fill(mindR_arr[0]);
			// h_dR_gq1->Fill(mindR_arr[1]);
			// h_dR_gq2->Fill(mindR_arr[2]);

			// Get the info on each of the trijet
			jet_pt_0  = j_sel_arr[0].Pt();
			jet_eta_0 = j_sel_arr[0].Eta();
			jet_phi_0 = j_sel_arr[0].Phi();
			jet_m_0   = j_sel_arr[0].M();
			jet_pt_1  = j_sel_arr[1].Pt();
			jet_eta_1 = j_sel_arr[1].Eta();
			jet_phi_1 = j_sel_arr[1].Phi();
			jet_m_1   = j_sel_arr[1].M();
			jet_pt_2  = j_sel_arr[2].Pt();
			jet_eta_2 = j_sel_arr[2].Eta();
			jet_phi_2 = j_sel_arr[2].Phi();
			jet_m_2   = j_sel_arr[2].M();
			
			jj_m_01 = (j_sel_arr[0]+j_sel_arr[1]).M();
			jj_m_02 = (j_sel_arr[0]+j_sel_arr[2]).M();
			jj_m_12 = (j_sel_arr[1]+j_sel_arr[2]).M();
			jj_pt_01 = (j_sel_arr[0]+j_sel_arr[1]).Pt();
			jj_pt_02 = (j_sel_arr[0]+j_sel_arr[2]).Pt();
			jj_pt_12 = (j_sel_arr[1]+j_sel_arr[2]).Pt();
			jj_eta_01 = (j_sel_arr[0]+j_sel_arr[1]).Eta();
			jj_eta_02 = (j_sel_arr[0]+j_sel_arr[2]).Eta();
			jj_eta_12 = (j_sel_arr[1]+j_sel_arr[2]).Eta();
			std::vector<float> jj_dR;
			jj_dR_01 = j_sel_arr[0].DeltaR(j_sel_arr[1]); jj_dR.push_back(jj_dR_01);
			jj_dR_02 = j_sel_arr[0].DeltaR(j_sel_arr[2]); jj_dR.push_back(jj_dR_02);
			jj_dR_12 = j_sel_arr[1].DeltaR(j_sel_arr[2]); jj_dR.push_back(jj_dR_12);
			jj_dR_min = TMath::MinElement(3, &jj_dR[0]);
			jj_dR_max = TMath::MaxElement(3, &jj_dR[0]);
			std::vector<float> jj_dEta;
			jj_dEta_01 = abs(jet_eta_0 - jet_eta_1); jj_dEta.push_back(jj_dEta_01);
			jj_dEta_02 = abs(jet_eta_0 - jet_eta_2); jj_dEta.push_back(jj_dEta_02);
			jj_dEta_12 = abs(jet_eta_1 - jet_eta_2); jj_dEta.push_back(jj_dEta_12);
			jj_dEta_min = TMath::MinElement(3, &jj_dEta[0]);
			jj_dEta_max = TMath::MaxElement(3, &jj_dEta[0]);
			std::vector<float> jj_j_dR;
			jj_j_dR_01_2 = (j_sel_arr[0]+j_sel_arr[1]).DeltaR(j_sel_arr[2]); jj_j_dR.push_back(jj_j_dR_01_2);
			jj_j_dR_02_1 = (j_sel_arr[0]+j_sel_arr[2]).DeltaR(j_sel_arr[1]); jj_j_dR.push_back(jj_j_dR_02_1);
			jj_j_dR_12_0 = (j_sel_arr[1]+j_sel_arr[2]).DeltaR(j_sel_arr[0]); jj_j_dR.push_back(jj_j_dR_12_0);
			jj_j_dR_min = TMath::MinElement(3, &jj_j_dR[0]);
			std::vector<float> jj_j_dEtaAbs;
			jj_j_dEtaAbs_01_2 = abs(jj_eta_01 - jet_eta_2); jj_j_dEtaAbs.push_back(jj_j_dEtaAbs_01_2);
			jj_j_dEtaAbs_02_1 = abs(jj_eta_02 - jet_eta_1); jj_j_dEtaAbs.push_back(jj_j_dEtaAbs_02_1);
			jj_j_dEtaAbs_12_0 = abs(jj_eta_12 - jet_eta_0); jj_j_dEtaAbs.push_back(jj_j_dEtaAbs_12_0);
			jj_j_dEtaAbs_max = TMath::MaxElement(3, &jj_j_dEtaAbs[0]);
			std::vector<float> jj_j_dPhi;
			float jj_phi = (j_sel_arr[0]+j_sel_arr[1]).Phi();
			jj_j_dPhi_01_2 = abs(jj_phi - jet_phi_2) <= TMath::Pi() ? abs(jj_phi - jet_phi_2) : 2*TMath::Pi()-abs(jj_phi - jet_phi_2); jj_j_dPhi.push_back(jj_j_dPhi_01_2);
			jj_phi = (j_sel_arr[0]+j_sel_arr[2]).Phi();
			jj_j_dPhi_02_1 = abs(jj_phi - jet_phi_1) <= TMath::Pi() ? abs(jj_phi - jet_phi_1) : 2*TMath::Pi()-abs(jj_phi - jet_phi_1); jj_j_dPhi.push_back(jj_j_dPhi_02_1);
			jj_phi = (j_sel_arr[1]+j_sel_arr[2]).Phi();
			jj_j_dPhi_12_0 = abs(jj_phi - jet_phi_0) <= TMath::Pi() ? abs(jj_phi - jet_phi_0) : 2*TMath::Pi()-abs(jj_phi - jet_phi_0); jj_j_dPhi.push_back(jj_j_dPhi_12_0);
			jj_j_dPhi_min = TMath::MinElement(3, &jj_j_dPhi[0]);

			M_jjj = (j_sel_arr[0] + j_sel_arr[1] + j_sel_arr[2]).M();
			pt_jjj = (j_sel_arr[0] + j_sel_arr[1] + j_sel_arr[2]).Pt();
			eta_jjj = (j_sel_arr[0] + j_sel_arr[1] + j_sel_arr[2]).Eta();
			jet_ptoverM_0 = jet_pt_0 / M_jjj;
			jet_ptoverM_1 = jet_pt_1 / M_jjj;
			jet_ptoverM_2 = jet_pt_2 / M_jjj;
			jj_ptoverM_01 = jj_pt_01 / M_jjj;
			jj_ptoverM_02 = jj_pt_02 / M_jjj;
			jj_ptoverM_12 = jj_pt_12 / M_jjj;
			jj_moverM_01 = jj_m_01 / M_jjj;
			jj_moverM_02 = jj_m_02 / M_jjj;
			jj_moverM_12 = jj_m_12 / M_jjj;
			
			// GEN-RECO matching
			int isMatch_j0 = 0;
			int isMatch_j1 = 0;
			int isMatch_j2 = 0;
			for(int iglu = 0; iglu < 3; iglu++){
				if(j_sel_index_arr[0] == j_sel_index_arr[matched_index_arr[iglu]] && mindR_arr[iglu] < 0.4) isMatch_j0 = 1;
			}
			for(int iglu = 0; iglu < 3; iglu++){
				if(j_sel_index_arr[1] == j_sel_index_arr[matched_index_arr[iglu]] && mindR_arr[iglu] < 0.4) isMatch_j1 = 2;
			}
			for(int iglu = 0; iglu < 3; iglu++){
				if(j_sel_index_arr[2] == j_sel_index_arr[matched_index_arr[iglu]] && mindR_arr[iglu] < 0.4) isMatch_j2 = 4;
			}
			
			gen_dijet_matched = isMatch_j0+isMatch_j1+isMatch_j2;
			if(gen_dijet_matched == 7){
				h_dR_gq0->Fill(mindR_arr[0]);
				h_dR_gq1->Fill(mindR_arr[1]);
				h_dR_gq2->Fill(mindR_arr[2]);
				count4++;
			}
			outTree->Fill();
			
			//cout<<"------------------------"<<endl;
		}// End of event loop
		eventTree = 0;
		infile->Close();
	}// End of file loop
	cout<<count0<<" "<<count1<<" "<<count2<<" "<<count3<<" "<<count4<<endl;
	insample.close();
	outFile->cd();
	outFile->Write();
	outFile->Close();
	
	TCanvas *c1 = new TCanvas("","",1200,900);
	c1->cd();
	h_njet->Draw();
	c1->Print("njet.png");
	
	std::vector<int> maxentry;
	maxentry.push_back(h_dR_gq0->GetBinContent(1));
	maxentry.push_back(h_dR_gq1->GetBinContent(1));
	maxentry.push_back(h_dR_gq2->GetBinContent(1));
	int max_h_index = std::max_element(maxentry.begin(),maxentry.end()) - maxentry.begin();
	int max = maxentry.at(max_h_index) * 1.1;
	
	TLegend *leg = new TLegend(0.6,0.75,0.75,0.85);
	TCanvas *c2 = new TCanvas("","",1200,900);
	c2->cd();
	h_dR_gq0->SetLineColor(1);
	h_dR_gq1->SetLineColor(2);
	h_dR_gq2->SetLineColor(4);
	h_dR_gq0->GetYaxis()->SetRangeUser(0,max);
	h_dR_gq0->GetXaxis()->SetTitle("min #Delta R");
	h_dR_gq0->Draw();
	h_dR_gq1->Draw("SAME");
	h_dR_gq2->Draw("SAME");
	leg->AddEntry(h_dR_gq0, "gq0");
	leg->AddEntry(h_dR_gq1, "gq1");
	leg->AddEntry(h_dR_gq2, "gq2");
	leg->Draw();
	c2->Print("dR.png");
	
  gBenchmark->Show("selectTrijet");
}
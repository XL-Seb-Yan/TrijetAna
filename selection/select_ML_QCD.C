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

void select_ML_QCD(const TString samplename="Res1ToRes2GluTo3Glu_M1-1000_R-0p5", 
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
	TH1F *h_dijet_mass = new TH1F("dijet_mass","dijet mass",100,0,3000);
	
	// Data Structure for output skimmed files
	// Dijet-system
	float dijet_pt, dijet_eta, dijet_phi, dijet_m, dR_jj, dEta_jj, dPhi_jj, m_jj;
	// Dijet-single jet
	float jet_pt_0, jet_eta_0, jet_phi_0, jet_m_0, jet_ptoverm_0;
	float jet_pt_1, jet_eta_1, jet_phi_1, jet_m_1, jet_ptoverm_1;
	// Another jet
	float jet_pt_2, jet_eta_2, jet_phi_2, jet_m_2, jet_ptoverm_2, dR_jj_j, dEta_jj_j, dPhi_jj_j;
	// System
	float M_jjj, jet_ptoverM_0, jet_ptoverM_1, jet_ptoverM_2, dijet_ptoverM;	   
	// Dijet in rest frame (use third jet to boost back)
	float dijet_res_dPt, dijet_res_dPhi, dijet_res_dEta;
	// GEN matching
	int gen_dijet_matched;
	float gen_dijet_matched_mass;
	// evt info
	UInt_t run_num;
	ULong64_t evt_num;
	int lumi_block;  		  
	
	TString outfilename = outputDir + TString("/") + samplename + TString("_ML_study.root");
	TFile *outFile = new TFile(outfilename,"RECREATE"); 
	TTree *outTree = new TTree("Events","Events");
	outTree->Branch("dijet_pt",       &dijet_pt);
	outTree->Branch("dijet_eta",      &dijet_eta);
	outTree->Branch("dijet_phi",      &dijet_phi);
	outTree->Branch("dR_jj",          &dR_jj);
	outTree->Branch("dEta_jj",        &dEta_jj);
	outTree->Branch("dPhi_jj",        &dPhi_jj);
	outTree->Branch("m_jj",           &m_jj);
	outTree->Branch("jet_pt_0",       &jet_pt_0);
	outTree->Branch("jet_eta_0",      &jet_eta_0);
	outTree->Branch("jet_phi_0",      &jet_phi_0);
	outTree->Branch("jet_m_0",        &jet_m_0);
	outTree->Branch("jet_ptoverm_0",  &jet_ptoverm_0);
	outTree->Branch("jet_pt_1",       &jet_pt_1);
	outTree->Branch("jet_eta_1",      &jet_eta_1);
	outTree->Branch("jet_phi_1",      &jet_phi_1);
	outTree->Branch("jet_m_1",        &jet_m_1);
	outTree->Branch("jet_ptoverm_1",  &jet_ptoverm_1);
	outTree->Branch("jet_pt_2",       &jet_pt_2);
	outTree->Branch("jet_eta_2",      &jet_eta_2);
	outTree->Branch("jet_phi_2",      &jet_phi_2);
	outTree->Branch("jet_m_2",        &jet_m_2);
	outTree->Branch("jet_ptoverm_2",  &jet_ptoverm_2);
	outTree->Branch("dR_jj_j",        &dR_jj_j);
	outTree->Branch("dEta_jj_j",      &dEta_jj_j);
	outTree->Branch("dPhi_jj_j",      &dPhi_jj_j);
	outTree->Branch("M_jjj",          &M_jjj);
	outTree->Branch("jet_ptoverM_0",  &jet_ptoverM_0);
	outTree->Branch("jet_ptoverM_1",  &jet_ptoverM_1);
	outTree->Branch("jet_ptoverM_2",  &jet_ptoverM_2);
	outTree->Branch("dijet_ptoverM",  &dijet_ptoverM);
	outTree->Branch("dijet_res_dPt",  &dijet_res_dPt);
	outTree->Branch("dijet_res_dEta",  &dijet_res_dEta);
	outTree->Branch("dijet_res_dPhi",  &dijet_res_dPhi);
	outTree->Branch("run_num",        &run_num);
	outTree->Branch("evt_num",        &evt_num);
	outTree->Branch("lumi_block",     &lumi_block);
  
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
		TTreeReader fReader;  //!the tree reader
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
		TTreeReaderValue<UInt_t> nGenPart = {fReader, "nGenPart"};
		TTreeReaderArray<Float_t> GenPart_eta = {fReader, "GenPart_eta"};
		TTreeReaderArray<Float_t> GenPart_mass = {fReader, "GenPart_mass"};
		TTreeReaderArray<Float_t> GenPart_phi = {fReader, "GenPart_phi"};
		TTreeReaderArray<Float_t> GenPart_pt = {fReader, "GenPart_pt"};
		TTreeReaderArray<Int_t> GenPart_genPartIdxMother = {fReader, "GenPart_genPartIdxMother"};
		TTreeReaderArray<Int_t> GenPart_pdgId = {fReader, "GenPart_pdgId"};
		TTreeReaderArray<Int_t> GenPart_status = {fReader, "GenPart_status"};
		// TTreeReaderArray<Int_t> Jet_genJetIdx = {fReader, "Jet_genJetIdx"};
		
		eventTree = (TTree*)infile->Get("Events");
		assert(eventTree);
		fReader.SetTree(eventTree);
		
		for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++){
		//for(UInt_t ientry=0; ientry<3; ientry++){
			count0++;
			fReader.SetLocalEntry(ientry);
			
			// Evt info
			run_num = *run;
			evt_num = *event;
			lumi_block = *luminosityBlock;
			
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
			
			// Select trijet
			std::vector<int> j_sel_index_arr; //store the indices of RECO jets selected
			std::vector<TLorentzVector> j_sel_arr;
			for(UInt_t i=0; i<*nJet; i++){
				if(Jet_pt[i] < 100) continue;
				if(abs(Jet_eta[i]) > 2.5) continue;
				if(Jet_jetId[i] < 6) continue;
				TLorentzVector vj;
				vj.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
				j_sel_index_arr.push_back(i);
				j_sel_arr.push_back(vj);
				if(j_sel_arr.size() == 4)
					break;
			}
			h_njet->Fill(j_sel_index_arr.size());
			if(j_sel_index_arr.size() < 4) continue;
			
			// making trijet from the selected jets, 0,1 are assumed to be coming from Res2 when 2 comes from Res1
			int index_j0[12]={0,0,0,0,0,0,1,1,1,1,2,2}; 
			int index_j1[12]={1,1,2,2,3,3,2,2,3,3,3,3}; 
			int index_j2[12]={2,3,1,3,1,2,0,3,0,2,0,1};
			
			for(int i=0; i<12; i++){
				int j0 = index_j0[i]; 
				int j1 = index_j1[i]; 
				int j2 = index_j2[i];
				
				M_jjj = (j_sel_arr[j0] + j_sel_arr[j1] + j_sel_arr[j2]).M();
				
				// Get the info on each of the trijet
				jet_pt_0  = j_sel_arr[j0].Pt();
				jet_eta_0 = j_sel_arr[j0].Eta();
				jet_phi_0 = j_sel_arr[j0].Phi();
				jet_m_0   = j_sel_arr[j0].M();
				jet_pt_1  = j_sel_arr[j1].Pt();
				jet_eta_1 = j_sel_arr[j1].Eta();
				jet_phi_1 = j_sel_arr[j1].Phi();
				jet_m_1   = j_sel_arr[j1].M();
				jet_pt_2  = j_sel_arr[j2].Pt();
				jet_eta_2 = j_sel_arr[j2].Eta();
				jet_phi_2 = j_sel_arr[j2].Phi();
				jet_m_2   = j_sel_arr[j2].M();
				
				jet_ptoverM_0 = jet_pt_0 / M_jjj;
				jet_ptoverM_1 = jet_pt_1 / M_jjj;
				jet_ptoverM_2 = jet_pt_2 / M_jjj;
				
				TLorentzVector dijet; // Here we use the 1st and 2nd jets to make the dijet, NOTE, it is not necessarily the correct paring
				dijet = (j_sel_arr[j0] + j_sel_arr[j1]);
				dijet_pt = dijet.Pt();
				dijet_eta = dijet.Eta();
				dijet_phi = dijet.Phi();
				m_jj = dijet.M();
				jet_ptoverm_0 = jet_pt_0 / m_jj;
				jet_ptoverm_1 = jet_pt_1 / m_jj;
				jet_ptoverm_2 = jet_pt_2 / m_jj;
				dR_jj = j_sel_arr[j0].DeltaR(j_sel_arr[j1]);
				dEta_jj = abs(jet_eta_0 - jet_eta_1);
				dPhi_jj = abs(jet_phi_0 - jet_phi_1) <= TMath::Pi() ? abs(jet_phi_0 - jet_phi_1) : 2*TMath::Pi()-abs(jet_phi_0 - jet_phi_1);
				dR_jj_j = (j_sel_arr[j0] + j_sel_arr[j1]).DeltaR(j_sel_arr[j2]);
				dEta_jj_j = abs((j_sel_arr[j0] + j_sel_arr[j1]).Eta() - j_sel_arr[j2].Eta());
				dPhi_jj_j = abs((j_sel_arr[j0] + j_sel_arr[j1]).Phi() - j_sel_arr[j2].Phi()) <= TMath::Pi() ? abs((j_sel_arr[j0] + j_sel_arr[j1]).Phi() - j_sel_arr[j2].Phi()) : 2*TMath::Pi()-abs((j_sel_arr[j0] + j_sel_arr[j1]).Phi() - j_sel_arr[j2].Phi());
				dijet_ptoverM = dijet_pt / M_jjj;
				
				TLorentzVector dijet_res_j0, dijet_res_j1;
				dijet_res_j0 = j_sel_arr[j0];
				dijet_res_j1 = j_sel_arr[j1];
				dijet_res_j0.Boost(j_sel_arr[j2].BoostVector());
				dijet_res_j1.Boost(j_sel_arr[j2].BoostVector());
				dijet_res_dPt = abs(dijet_res_j0.Pt() - dijet_res_j1.Pt());
				dijet_res_dEta = abs(dijet_res_j0.Eta() - dijet_res_j1.Eta());
				dijet_res_dPhi = abs(dijet_res_j0.Phi() - dijet_res_j1.Phi());

				outTree->Fill();
			}
			//cout<<"------------------------"<<endl;
		}// End of event loop
		eventTree = 0;
		infile->Close();
	}// End of file loop
	cout<<count0<<" "<<count1<<" "<<count2<<endl;
	insample.close();
	outFile->cd();
	outFile->Write();
	outFile->Close();
	
	TCanvas *c0 = new TCanvas("","",1200,900);
	c0->cd();
	h_dijet_mass->Draw();
	c0->Print("dijet_m.png");
	
	TCanvas *c1 = new TCanvas("","",1200,900);
	c1->cd();
	h_njet->Draw();
	c1->Print("njet.png");
	
  gBenchmark->Show("selectTrijet");
}
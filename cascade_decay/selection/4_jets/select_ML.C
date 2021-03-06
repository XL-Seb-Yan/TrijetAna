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

void select_ML(const TString samplename="Res1ToRes2GluTo3Glu_M1-1000_R-0p5", 
								const float R1M = 1000, 
								const float rho = 0.5,
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
	
	TH1F *h_njet = new TH1F("njet","njet",15,0,30);
	TH1F *h_jet_pt = new TH1F("jet_pt","jet pt",100,0,1500);
	
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
	outTree->Branch("dijet_phi",     &dijet_phi);
	outTree->Branch("dijet_m",        &dijet_m);
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
	outTree->Branch("gen_dijet_matched",   &gen_dijet_matched);
	outTree->Branch("gen_dijet_matched_mass",   &gen_dijet_matched_mass);
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
			h_njet->Fill(*nJet);
			std::vector<TLorentzVector> trijet;
			std::vector<int> matched_genjet_index;
			for(UInt_t i=0; i<*nJet; i++){
				h_jet_pt->Fill(Jet_pt[i]);
				if(Jet_pt[i] < 100) continue;
				if(abs(Jet_eta[i]) > 2.5) continue;
				if(Jet_jetId[i] < 6) continue;
				TLorentzVector vj;
				vj.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
				trijet.push_back(vj);
				//matched_genjet_index.push_back(Jet_genJetIdx[i]);
				if(trijet.size() == 3) break;
			}
			
			if(trijet.size() != 3) continue;
			count2++;
			
			M_jjj = (trijet[0] + trijet[1] + trijet[2]).M();
			
			// Get matched GEN jets;
			/*
			std::vector<TLorentzVector> gen_trijet;
			if(matched_genjet_index[0] > -1 && matched_genjet_index[1] > -1 && matched_genjet_index[2] > -1){
				TLorentzVector gen_jet_0, gen_jet_1, gen_jet_2;
				gen_jet_0.SetPtEtaPhiM(GenJet_pt[matched_genjet_index[0]], GenJet_eta[matched_genjet_index[0]], GenJet_phi[matched_genjet_index[0]], GenJet_mass[matched_genjet_index[0]]);
				gen_jet_1.SetPtEtaPhiM(GenJet_pt[matched_genjet_index[1]], GenJet_eta[matched_genjet_index[1]], GenJet_phi[matched_genjet_index[1]], GenJet_mass[matched_genjet_index[1]]);
				gen_jet_2.SetPtEtaPhiM(GenJet_pt[matched_genjet_index[2]], GenJet_eta[matched_genjet_index[2]], GenJet_phi[matched_genjet_index[2]], GenJet_mass[matched_genjet_index[2]]);
				gen_trijet.push_back(gen_jet_0);
				gen_trijet.push_back(gen_jet_1);
				gen_trijet.push_back(gen_jet_2);
			}
			*/
			
			// Do GEN-RECO matching
			std::vector<TLorentzVector> gen_gq_list; //store the 4-vector of gen level q or g
			std::vector<UInt_t> gen_jet_index; //store the index of resulted gen jets of the decay of R2
			TLorentzVector gen_gq_temp;
			TLorentzVector gen_jet_temp;
			for(UInt_t i=0; i<*nGenPart; i++){
				if(GenPart_pdgId[i] != 21 || GenPart_status[i] != 23) continue;// This selects the 3 gluons we need
				// cout<<"Index: "<<i<<" Pt: "<<GenPart_pt[i]<<" Eta: "<<GenPart_eta[i]<<" Phi: "<<GenPart_phi[i]<<" Mass: "<<GenPart_mass[i]<<" ID: "<<GenPart_pdgId[i]<<" Status: "<<GenPart_status[i]<<" Momindex: "<<GenPart_genPartIdxMother[i]<<endl;
				gen_gq_temp.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
				gen_gq_list.push_back(gen_gq_temp); //The last two 4-vector correspond to the decay production of R2
			}
			if(gen_gq_list.size() != 3) 
				cout<<"Warning: gen_gq_list size is not 3!"<<endl;
			
			// making dijet from trijet, make sure j0 is always the leading jet
			int index_j0[3]={0,1,0}; 
			int index_j1[3]={1,2,2}; 
			int index_j2[3]={2,0,1};
			
			for(int i=0; i<3; i++){
				int j0 = index_j0[i]; 
				int j1 = index_j1[i]; 
				int j2 = index_j2[i];
				
				// Get the info on each of the trijet
				jet_pt_0  = trijet[j0].Pt();
				jet_eta_0 = trijet[j0].Eta();
				jet_phi_0 = trijet[j0].Phi();
				jet_m_0   = trijet[j0].M();
				jet_pt_1  = trijet[j1].Pt();
				jet_eta_1 = trijet[j1].Eta();
				jet_phi_1 = trijet[j1].Phi();
				jet_m_1   = trijet[j1].M();
				jet_pt_2  = trijet[j2].Pt();
				jet_eta_2 = trijet[j2].Eta();
				jet_phi_2 = trijet[j2].Phi();
				jet_m_2   = trijet[j2].M();
				
				jet_ptoverM_0 = jet_pt_0 / M_jjj;
				jet_ptoverM_1 = jet_pt_1 / M_jjj;
				jet_ptoverM_2 = jet_pt_2 / M_jjj;
				
				TLorentzVector dijet; // Here we use the 1st and 2nd jets to make the dijet, NOTE, it is not necessarily the correct paring
				dijet = (trijet[j0] + trijet[j1]);
				dijet_pt = dijet.Pt();
				dijet_eta = dijet.Eta();
				dijet_phi = dijet.Phi();
				dijet_m = dijet.M();
				m_jj = dijet.M();
				jet_ptoverm_0 = jet_pt_0 / m_jj;
				jet_ptoverm_1 = jet_pt_1 / m_jj;
				jet_ptoverm_2 = jet_pt_2 / m_jj;
				dR_jj = trijet[j0].DeltaR(trijet[j1]);
				dEta_jj = abs(jet_eta_0 - jet_eta_1);
				dPhi_jj = abs(jet_phi_0 - jet_phi_1);
				dR_jj_j = (trijet[j0] + trijet[j1]).DeltaR(trijet[j2]);
				dEta_jj_j = abs((trijet[j0] + trijet[j1]).Eta() - trijet[j2].Eta());
				dPhi_jj_j = abs((trijet[j0] + trijet[j1]).Phi() - trijet[j2].Phi());
				dijet_ptoverM = dijet_pt / M_jjj;
				
				TLorentzVector dijet_res_j0, dijet_res_j1;
				dijet_res_j0 = trijet[j0];
				dijet_res_j1 = trijet[j1];
				dijet_res_j0.Boost(trijet[j2].BoostVector());
				dijet_res_j1.Boost(trijet[j2].BoostVector());
				dijet_res_dPt = abs(dijet_res_j0.Pt() - dijet_res_j1.Pt());
				dijet_res_dEta = abs(dijet_res_j0.Eta() - dijet_res_j1.Eta());
				dijet_res_dPhi = abs(dijet_res_j0.Phi() - dijet_res_j1.Phi());
				
				// GEN-RECO matching
				int is_Matched_j0 = 0;
				int is_Matched_j1 = 0;
				for(int i=1; i<3; i++){
					if(trijet[j0].DeltaR(gen_gq_list[i]) < 0.4)
						is_Matched_j0+=i;
				}
				for(int i=1; i<3; i++){
					if(trijet[j1].DeltaR(gen_gq_list[i]) < 0.4)
						is_Matched_j1+=i;
				}
				// if(is_Matched_j0 > 2) cout<<"Warning: j0 matched to both GEN g(q)s, j0 is: "<<j0<<endl;
				// if(is_Matched_j1 > 2) cout<<"Warning: j1 matched to both GEN g(q)s: j1 is: "<<j1<<endl;
				gen_dijet_matched = (is_Matched_j0 > 0 && is_Matched_j0 < 3) &&  (is_Matched_j1 > 0 && is_Matched_j1 < 3) && (is_Matched_j0 != is_Matched_j1);
				// if(gen_dijet_matched && (dijet_m < 700 || dijet_m > 1600)){
					// for(UInt_t i=0; i<*nGenPart; i++){
						// if(GenPart_pdgId[i] != 21 || GenPart_status[i] != 23) continue;// This selects the 3 gluons we need
						// cout<<"Index: "<<i<<" Pt: "<<GenPart_pt[i]<<" Eta: "<<GenPart_eta[i]<<" Phi: "<<GenPart_phi[i]<<" Mass: "<<GenPart_mass[i]<<" ID: "<<GenPart_pdgId[i]<<" Status: "<<GenPart_status[i]<<" Momindex: "<<GenPart_genPartIdxMother[i]<<endl;
					// }
					// cout<<"RECO jet 0 Pt: "<<jet_pt_0<<" Eta: "<<jet_eta_0<<" Phi: "<<jet_phi_0<<endl;
					// cout<<"RECO jet 1 Pt: "<<jet_pt_1<<" Eta: "<<jet_eta_1<<" Phi: "<<jet_phi_1<<endl;
					// cout<<"Dijet mass: "<<m_jj<<endl;
					// cout<<"Delta R for j0: "<<trijet[j0].DeltaR(gen_gq_list[1])<<" "<<trijet[j0].DeltaR(gen_gq_list[2])<<endl;
					// cout<<"Delta R for j1: "<<trijet[j1].DeltaR(gen_gq_list[1])<<" "<<trijet[j1].DeltaR(gen_gq_list[2])<<endl;
					// cout<<"---------------------------------------------------------"<<endl;
				// }
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
	h_jet_pt->Draw();
	c0->Print("pt.png");
	
	TCanvas *c1 = new TCanvas("","",1200,900);
	c1->cd();
	h_njet->Draw();
	c1->Print("njet.png");
	
  gBenchmark->Show("selectTrijet");
}
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
#include <TEfficiency.h>
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

float maxdEta(TLorentzVector vj0, TLorentzVector vj1, TLorentzVector vj2){
	float maxdEta = -1;
	float dEta01 = abs(vj0.Eta() - vj1.Eta());
	if(dEta01 > maxdEta) maxdEta = dEta01;
	float dEta02 = abs(vj0.Eta() - vj2.Eta());
	if(dEta02 > maxdEta) maxdEta = dEta02;
	float dEta12 = abs(vj1.Eta() - vj2.Eta());
	if(dEta12 > maxdEta) maxdEta = dEta12;
	return maxdEta;
}

void HLTEff(const TString samplename="SingleMuon_Run2017B",  
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
	
	TH1F *hist_01 = new TH1F("hist_01","pT lead offline-sel",500,0,10000); // fired mu trigger and pass offline selection
	TH1F *hist_02 = new TH1F("hist_02","M_jjj offline-sel",500,0,10000);
	TH1F *hist_11 = new TH1F("hist_11","pT lead offline-sel & HLT",500,0,10000); // fired mu trigger and pass offline selection and fired jet trigger
	TH1F *hist_12 = new TH1F("hist_12","M_jjj offline-sel & HLT",500,0,10000);
	
	TString outfilename = outputDir + TString("/") + samplename + TString("_HLTEff.root");
	TFile *outFile = new TFile(outfilename,"RECREATE"); 
  
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
		cout<<"====== Using 2016 HLT set ======="<<endl;
		//--HLT trigger--
		TTreeReaderValue<Bool_t> HLT_Mu50     = {fReader, "HLT_Mu50"};
		TTreeReaderValue<Bool_t> HLT_PFHT900  = {fReader, "HLT_PFHT900"};
		TTreeReaderValue<Bool_t> HLT_PFJet500 = {fReader, "HLT_PFJet500"};
		TTreeReaderValue<Bool_t> HLT_AK8PFJet360_TrimMass30 = {fReader, "HLT_AK8PFJet360_TrimMass30"};
		TTreeReaderValue<Bool_t> HLT_CaloJet500_NoJetID     = {fReader, "HLT_CaloJet500_NoJetID"};
		//-- Dummy fReader, allow use HLTs that do not exsit in some PDs
		TTreeReader fReader_tmp;
		TTreeReaderValue<Bool_t> HLT_AK8PFJet450 = {fReader_tmp, "HLT_AK8PFJet450"};
#else
		cout<<"====== Using 2017/8 HLT set ======="<<endl;
		TTreeReaderValue<Bool_t> HLT_Mu50     = {fReader, "HLT_Mu50"};
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
		
		eventTree = (TTree*)infile->Get("Events");
		assert(eventTree);
#if Run == 2016
		if(eventTree->GetListOfBranches()->FindObject("HLT_AK8PFJet450")){HLT_AK8PFJet450 = {fReader, "HLT_AK8PFJet450"};}
#endif
		
		fReader.SetTree(eventTree);
		
		for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++){
			if(ientry % 50000 == 0)
				cout<<"Processing: "<<float(ientry) / eventTree->GetEntries()<<endl;
		//for(UInt_t ientry=0; ientry<100; ientry++){
			count0++;
			fReader.SetLocalEntry(ientry);
			
			// HLT reference
			bool passMuTrig = false;
			passMuTrig = *HLT_Mu50;
			if (!passMuTrig) continue;
			count1++;
			
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
			if(j_sel_index_arr.size() < 3) continue;
			if(maxdEta(j_sel_arr[0], j_sel_arr[1], j_sel_arr[2]) > 1.3) continue;
			count2++;
			
			hist_01->Fill(j_sel_arr[0].Pt());
			hist_02->Fill((j_sel_arr[0] + j_sel_arr[1] + j_sel_arr[2]).M());
			
			bool passHLT = false;
#if Run == 2016
			if(eventTree->GetListOfBranches()->FindObject("HLT_AK8PFJet450")){
				passHLT = (*HLT_PFHT900 || *HLT_PFJet500 || *HLT_AK8PFJet360_TrimMass30 || *HLT_CaloJet500_NoJetID || *HLT_AK8PFJet450);
			}
			else
				passHLT = (*HLT_PFHT900 || *HLT_PFJet500 || *HLT_AK8PFJet360_TrimMass30 || *HLT_CaloJet500_NoJetID);
#else
			passHLT = (*HLT_PFHT1050 || *HLT_AK8PFJet550 || *HLT_AK8PFJet500 || *HLT_PFJet500 || *HLT_CaloJet500_NoJetID || *HLT_CaloJet550_NoJetID);
#endif
			if (!passHLT) continue;
			count3++;
			hist_11->Fill(j_sel_arr[0].Pt());
			hist_12->Fill((j_sel_arr[0] + j_sel_arr[1] + j_sel_arr[2]).M());
			//cout<<"------------------------"<<endl;
		}// End of event loop
		eventTree = 0;
		infile->Close();
	}// End of file loop
	cout<<count0<<" "<<count1<<" "<<count2<<" "<<count3<<" "<<count4<<endl;
	insample.close();
	outFile->cd();
	hist_01->Write();
	hist_02->Write();
	hist_11->Write();
	hist_12->Write();
	
	TEfficiency *Eff_1 = new TEfficiency(*hist_11,*hist_01);
	TEfficiency *Eff_2 = new TEfficiency(*hist_12,*hist_02);
	Eff_1->Write();
	Eff_2->Write();
	
	TCanvas *c1 = new TCanvas("","",1200,900);
	c1->cd();
	Eff_1->Draw();
	TString outname = "TrigEff_pT_lj_" + samplename + ".png";
	c1->Print(outname);
	
	TCanvas *c2 = new TCanvas("","",1200,900);
	c2->cd();
	Eff_2->Draw();
	outname = "TrigEff_M_jjj_" + samplename + ".png";
	c2->Print(outname);
	
	outFile->Close();
  gBenchmark->Show("selectTrijet");
}
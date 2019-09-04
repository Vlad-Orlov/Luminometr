#include "TH1D.h"
#include "TChain.h"
#include <iostream>
#include "TApplication.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include <vector>

#define nSec 100

// TH1D* getHist(TString);


using namespace std;

// TH1D getHist(TString);

#define NHIST 1

TRandom3 *rnd = new TRandom3(0);

Int_t NREBIN = 100;


Bool_t isFired(Double_t *eff, Int_t nPhot){

	// Double_t eff[20] ={	0.6, 1.8, 3.4, 5.7, 8.1, 11.5, 16.4, 24., \
	// 					35.5, 50, 64.5, 76, 83.6, 88.5, 91.9, 94.3, \
	// 					96.6, 98.2, 99.4, 100};

	// Double_t eff[40] ={	0.3, 0.6, 1., 1.8, 2.4, 3.4, 4.1, 5.7, 6.5, 8.1, 11.5, 13., \
	// 					16.4, 19.,  24., 30.3, 35.5, 43.5, 50, 56.5, 64.5, 69.7, 76, \
	// 					81., 83.6, 87., 88.5, 91., 91.9, 93.5, 94.3, 95.9, 96.6, 97.6, \
	// 					98.2, 99., 99.4, 99.9, 100};



	Double_t dice = rnd->Uniform() * 100.;

	// Double_t dice = 43;

	if (nPhot == 0) return false;
	if (nPhot > 99) return true;

	if (dice < eff[nPhot-1]) return true;
	else return false;

	// delete rnd;

}


void getHist(TString filename, TH1D* h1Hnum)
{

	cout << filename << endl;
	// TChain *theChain1 = new TChain("T");

	TFile *file = new TFile(filename);
	TTree *theChain1 = (TTree*)file->Get("T");
	// if (!theChain1) { std::
	// theChain1->Add(filename);


	Double_t eff[100] = {0.06, 0.15, 0.3, 0.4, 0.6, 0.6, 0.9, 1.2, 1.4, \
		1.8, 1.9, 2.2, 2.6, 2.9, 3.4, 3.575, 3.95, 4.5, 4.825, 5.7, 5.9, 6.35, \
		7, 7.4, 8.1, 8.4, 9.05, 10, 10.9, 11.5, 12.1, 12.75, 14, 15, 16.4, 17, \
		18.3, 20.2, 22, 24, 25.6, 26.875, 29.75, 33, 35.5, 37.2, 39.175, 42.85, \
		46, 50, 52, 53.625, 57.25, 60, 64.3, 66.2, 67.625, 70.75, 73.8, 75.6, 76.5, \
		78.35, 80.1, 82, 83.6, 84.2, 84.9, 86.3, 87.3, 88.5, 89.1, 89.35, 90.2, 91, \
		91.9, 92.3, 92.5, 93.1, 94, 94.2, 94.6, 94.875, 95.45, 96, 96.6, 96.8, 97, \
		97.4, 97.7, 98.2, 98.35, 98.5, 98.8, 99.1, 99.4, 99.5, 99.55, 99.7, 99.9, 100 };


	Int_t nPart1;
	Int_t *nPhot = new Int_t[nSec];

	theChain1->SetBranchAddress("nPhot", nPhot);

	// TH1D *h1nColl = new TH1D("nColl1", "nColl1", 10, 0, 10);

	Long_t nEv1 = theChain1->GetEntries();

	
	////////// Loop 1 //////////
	for (Long_t j = 0; j <  nEv1; ++j) {
		theChain1->GetEntry(j);


		Bool_t isChecked[nSec] = {false};
		

		const Int_t nRebin = NREBIN;
		// const Int_t nRebin = 100;


		for (Int_t i = 0; i < nSec/nRebin; ++i){
			Double_t numOfHits = 0;			

			for (Int_t j = 0; j < nRebin; ++j){
				bool fired = false;
				if (isFired(eff, nPhot[nRebin*i+j])) {
					fired = true;
					// break;
				}
				if (fired) {
					numOfHits += 1;	
				}
			}

			h1Hnum->Fill( numOfHits );	
		}



		// for (int i = 0; i < nSec; ++i){

		// 	bool thisSec = isFired(eff, nPhot[2*i]);
		// 	bool lSec = isFired(eff, nPhot[2*i+1]);

		// 	if (thisSec && lSec) {
		// 		numOfHits += 1;	
		// 		++i;
		// 	}
		// 	// if (thisSec && oppSec) numOfHits += 2;	
		// }



		// for (int i = 0; i < nSec; ++i){

		// 	bool thisSec = isFired(eff, nPhot[i]);

		// 	if (thisSec) {
		// 		numOfHits += 1;	
		// 		++i;
		// 	}
		// 	// if (thisSec && oppSec) numOfHits += 2;	
		// }
		

		// h1Hnum->Fill( numOfHits );

	}

	//	h1Hnum->Scale(1. / nEv1);
	
	file->Close();

	// delete TrackID1, ParentID1, StationID1, X1, Y1, Z1, Momentum1, Px1, Py1, Pz1, Time1, PdgID1; 
	// h1Hnum->Draw("E1");
	
}


int main(int argc, char** argv){

	if (argc != 1) NREBIN = atoi(argv[1]);

	TH1D *array[NHIST];

	TH1D *hChi2 = new TH1D ("chi2", "chi2", 1000, 0, 100);

	
	TFile *out = new TFile("histos.root","RECREATE");

	// printHello();
	// getHist(filename);


	for (Int_t i = 0; i < NHIST; ++i) {
	  //TString filename = "newData_";
	  TString filename = "newDataMergedOpt.root";
	  //filename += TString (std::to_string(i+1));
	  //filename += TString(".root");
	  
	  //	  array[i] = new TH1D(TString("Hit Number"+std::to_string(i)),TString("HitNumber"+std::to_string(i)), nSec, -0.5, nSec-0.5);
	  array[i] = new TH1D(TString("Hit Number"+std::to_string(i)),TString("HitNumber"+std::to_string(i)), nSec, 0, nSec);

	  
	  
	  getHist(filename, array[i]);

	  out->cd();
	  array[i]->Write();
	  
	  
	  array[i]->SetLineColor(i);
	  // delete htemp;
	}
	
	
	// TApplication *app = new TApplication("name", &argc, argv);
	
	// TCanvas *c = new TCanvas();

	


	// Double_t binContent[NHIST];
	// Double_t binError[NHIST];
	// Int_t maxI;
	// Int_t maxJ;

	// Double_t maxChi2 = 0;

	// const Int_t nBins = 20;

	// Double_t average[nBins] = {0};
	// Double_t rms[nBins] = {0};

	// for (Int_t binK = 0; binK < nBins; ++binK){
	// 	for (Int_t histI = 0; histI < NHIST; ++histI){
	// 		average[binK] += array[histI]->GetBinContent(binK+1);
	// 	}
	// 	average[binK] /= NHIST;
	// 	if (average[binK] == 0) continue;

	// 	for (Int_t histI = 0; histI < NHIST; ++histI){
	// 		rms[binK] += 
	// 		(array[histI]->GetBinContent(binK+1) - average[binK])*
	// 		(array[histI]->GetBinContent(binK+1) - average[binK]);
	// 	}


	// 	rms[binK] = sqrt(rms[binK]/(NHIST)/(NHIST-1));
	// 	cout << "average: " <<  average[binK] << " \t rms: " << rms[binK]
	// 	<< "\t sqrt: "<< sqrt(average[binK]) << endl;
	// 	maxChi2 += 	(rms[binK] / sqrt(average[binK]))*
	// 				(rms[binK] / sqrt(average[binK]));
	// } 


	// out->Close();
	// c->Show();
	// app->Run();


	// std::cin;

	// delete[] array;
	// delete c;

	return 0;
}

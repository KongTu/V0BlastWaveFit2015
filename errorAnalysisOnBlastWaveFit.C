#include <iostream>
#include <fstream>
#include <math.h>
#include "TMath.h"
#include <stdio.h>
#include <iomanip>
#include <vector>
#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "makeMultiPanelCanvas.C"

#include "fitting.h"

using namespace RooFit;

using namespace std;


void errorAnalysisOnBlastWaveFit(){

	TFile* file1 = new TFile("/Users/kongkong/2015Research/Spectra/pPb2013/rpyDependence/files/MidRpy0_0.8_MultEff_v2.root");
	TFile* file2 = new TFile("/Users/kongkong/2015Research/Spectra/pPb2013/rpyDependence/files/MidRpy0_0.8_MultEff_noEffError_v3.root");

	TH1D* ks_wrong[8];
	TH1D* ks_right[8];

	TH1D* la_wrong[8];
	TH1D* la_right[8];

	stringstream ksName;
	stringstream laName;

	for(int mult = 0 ;mult < 8; mult++){

		ksName.str("");
		laName.str("");

		ksName << "ksSpectra_vtx_";
		ksName << mult+1;

		laName << "laSpectra_vtx_";
		laName << mult+1;

		ks_right[mult] = (TH1D*)file1->Get( ksName.str().c_str() );
		la_right[mult] = (TH1D*)file1->Get( laName.str().c_str() );

		ks_wrong[mult] = (TH1D*)file2->Get( ksName.str().c_str() );
		la_wrong[mult] = (TH1D*)file2->Get( laName.str().c_str() );

	}

	double ks_err_right[8][25];
	double ks_err_wrong[8][25];

	double la_err_right[8][25];
	double la_err_wrong[8][25];

	for(mult = 0; mult < 8; mult++){

		for(int pt = 0; pt < 25; pt++){

			ks_err_right[mult][pt] = ks_right[mult]->GetBinError(pt+3);
			ks_err_wrong[mult][pt] = ks_wrong[mult]->GetBinError(pt+3);

		}

		for(pt = 0; pt < 20; pt++){

			la_err_right[mult][pt] = la_right[mult]->GetBinError(pt+2);
			la_err_wrong[mult][pt] = la_wrong[mult]->GetBinError(pt+2);
		}
	}

	TH1D* ksErrorhist[8];
	TH1D* laErrorhist[8];

	double ks_ptbins_histo[28] = {0.0,0.1,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
	double la_ptbins_histo[22] = {0.0,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};

	for(mult = 0; mult < 8; mult++){

		ksName.str("");
		laName.str("");

		ksName << "ksErrorHist_";
		ksName << mult+1;

		laName << "laErrorHist_";
		laName << mult+1;

		ksErrorhist[mult] = new TH1D( ksName.str().c_str(), ksName.str().c_str(), 27, ks_ptbins_histo );
		laErrorhist[mult] = new TH1D( laName.str().c_str(), laName.str().c_str(), 21, la_ptbins_histo );

		for(pt = 0; pt < 25; pt++){

			ksErrorhist[mult]->SetBinContent(pt+3, ks_err_right[mult][pt]/ks_err_wrong[mult][pt] );
			ksErrorhist[mult]->SetBinError(pt+3, 0.0);
		}

		for(pt = 0; pt < 20; pt++){

			laErrorhist[mult]->SetBinContent(pt+2, la_err_right[mult][pt]/la_err_wrong[mult][pt] );
			laErrorhist[mult]->SetBinError(pt+2, 0.0);
		}
	}



	TCanvas* c3 = new TCanvas();
	c3->Divide(4,2,0,0);
	for(mult = 0; mult < 8; mult++){

		c3->cd(mult+1);
		ksErrorhist[mult]->GetYaxis()->SetRangeUser(1.0,2.0);
		ksErrorhist[mult]->SetMarkerStyle(20);
		ksErrorhist[mult]->Draw("P");
	}

	TCanvas* c4 = new TCanvas();
	c4->Divide(4,2,0,0);
	for(mult = 0; mult < 8; mult++){

		c4->cd(mult+1);
		laErrorhist[mult]->GetYaxis()->SetRangeUser(1.0,2.0);
		laErrorhist[mult]->SetMarkerStyle(20);
		laErrorhist[mult]->Draw("P");
	}


/*

	TCanvas* c1 = new TCanvas();
	c1->Divide(4,2,0,0);

	for(mult = 0; mult < 8; mult++){

		c1->cd(mult+1);
		ks_right[mult]->Divide( ks_wrong[mult] );
		ks_right[mult]->SetMarkerStyle(20);
		ks_right[mult]->GetYaxis()->SetRangeUser(0.95,1.05);
		ks_right[mult]->Draw("P");
	}

	TCanvas* c2 = new TCanvas();
	c2->Divide(4,2,0,0);

	for(mult = 0; mult < 8; mult++){

		c2->cd(mult+1);
		la_right[mult]->Divide( la_wrong[mult] );
		la_right[mult]->SetMarkerStyle(20);
		la_right[mult]->GetYaxis()->SetRangeUser(0.95,1.05);
		la_right[mult]->Draw("P");
	}




*/








}
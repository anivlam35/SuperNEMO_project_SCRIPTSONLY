#include "/sps/nemo/scratch/ikovalen/TKEvent_old/TKEvent/include/TKEvent.h" 
#include "config.h"

R__LOAD_LIBRARY(/sps/nemo/scratch/ikovalen/TKEvent_old/TKEvent/lib/libTKEvent.so);

using namespace std;

void All_Sources_Visualization()
{
	gROOT->SetBatch(kTRUE);

	// 0000000000000000000000000000
	// Build the RMS_Y(x) and RMS_Z(x) for selection ellipses from code Analysis

	TFile* file = new TFile(Form("%sTracks_Run-%d.root", PATH, RUN_N));

	int N_Point = X_OBSERV_MAX - X_OBSERV_MIN + 1;
	double iX = 0;
	
	TH2D* h_vert_real;  

	int    YBINS = int( (Y_MAX - Y_MIN) / 2.0 );
	int    ZBINS = int( (Z_MAX - Z_MIN) / 2.0 );

        stringstream h_name;

        h_name << "Hist for " << iX << " mm x < 0 All Sources ";						   
        h_vert_real = new TH2D(h_name.str().c_str(), h_name.str().c_str(), YBINS, Y_MIN, Y_MAX, ZBINS, Z_MIN, Z_MAX);
        h_vert_real->SetMaximum(750);
        h_vert_real->GetXaxis()->SetRangeUser(Y_MIN, Y_MAX);
        h_vert_real->GetYaxis()->SetRangeUser(Z_MIN, Z_MAX);
        h_vert_real->GetXaxis()->SetTitle("y[mm]");
        h_vert_real->GetYaxis()->SetTitle("z[mm]");
        h_vert_real->GetYaxis()->SetTitleOffset(1.2);
//        h_vert_real->GetZaxis()->SetTitle("N of events / 4 mm^{2}");
	h_vert_real->SetStats(0);

	for(int SRC_NO_i=0; SRC_NO_i < N_SRCPLN; SRC_NO_i++)
	{
		stringstream tr_name;

		tr_name << "Tracks for " << X_BasePlane << " mm x < 0 Source " << SRC_NO_i;
		TTree* tr = (TTree*) file->Get(tr_name.str().c_str());

		cout << tr_name.str() << endl;		

		double parA, parB, parC, parD;
		tr->SetBranchAddress("A", &parA);
		tr->SetBranchAddress("B", &parB);
		tr->SetBranchAddress("C", &parC);
		tr->SetBranchAddress("D", &parD);

		double verY, verZ;

		// Fill histogram of undistorted vertices, and calculate coordinates of mode and the sigmas in both directions

		//TFile *New_file = new TFile("Visual_Histos.root","RECREATE");
		
		for(int i = 0; i < tr->GetEntries(); i++)
		{
			tr->GetEntry(i);

		       	verY = parA*iX + parB;
	        	verZ = parC*iX + parD;
	       		h_vert_real->Fill(verY, verZ);
		}
	}
	TCanvas* C0 = new TCanvas("Canvas", "Canvas", 1500, 1000);
	h_vert_real->Draw("COLZ"); 
	//h_vert_real->Write(h_name.str().c_str());

	C0->SetLogz();
//	C0->Update();
//	TPaveStats *st = (TPaveStats*)h_vert_real->FindObject("stats");
//	st->SetX1NDC(0.7);			
//	st->SetX2NDC(0.9);
//	st->SetY1NDC(0.72);
//	st->SetY2NDC(0.9);

	C0->SaveAs(Form("%sAll_Sources_X_%.1lfmm.png", PATH, iX));

}




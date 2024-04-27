#include "/sps/nemo/scratch/ikovalen/TKEvent_old/TKEvent/include/TKEvent.h" 
#include "config.h"

R__LOAD_LIBRARY(/sps/nemo/scratch/ikovalen/TKEvent_old/TKEvent/lib/libTKEvent.so);

using namespace std;

void Primary_Visualization()
{
	gROOT->SetBatch(kTRUE);

	// 0000000000000000000000000000
	// Build the RMS_Y(x) and RMS_Z(x) for selection ellipses from code Analysis

	TFile* file = new TFile(Form("%sTracks_Run-%d.root", PATH, RUN_N));

	int N_Point = X_OBSERV_MAX - X_OBSERV_MIN + 1;
	double Arr_X[N_Point], Arr_Y[N_Point], Arr_Z[N_Point];
	
	TGraph* grY[N_SRCPLN];
	TGraph* grZ[N_SRCPLN];

	//TFile* pr_hist_file = new TFile(Form("Primary_Histos_Run-%d.root", RUN_N), "RECREATE");

	for(int NSOR=0; NSOR<N_SRCPLN; NSOR++)
	{
		stringstream tr_name;

		tr_name << "Tracks for " << X_BasePlane << " mm x < 0 Source " << NSOR;
		TTree* tr = (TTree*) file->Get(tr_name.str().c_str());

		double parA, parB, parC, parD;
		tr->SetBranchAddress("A", &parA);
		tr->SetBranchAddress("B", &parB);
		tr->SetBranchAddress("C", &parC);
		tr->SetBranchAddress("D", &parD);

		double verY, verZ;

		double YMIN  = Y_MIN + (NSOR % N_COLS + 0.5) * Y_RECT_SIZE - 200;
		double YMAX  = Y_MIN + (NSOR % N_COLS + 0.5) * Y_RECT_SIZE + 200;
		int    YBINS = (YMAX - YMIN) / 2.0;

		double ZMIN  = Z_MAX - (NSOR / N_COLS + 0.5) * Z_RECT_SIZE - 140;
		double ZMAX  = Z_MAX - (NSOR / N_COLS + 0.5) * Z_RECT_SIZE + 140;
		int    ZBINS = (ZMAX - ZMIN) / 2.0;

		TH2D* h_vert_real[2*X_OBSERV_MAX+1];  

		for(int iX = X_OBSERV_MIN; iX<X_OBSERV_MAX+1; iX+=10)
		{
			stringstream h_name;

			h_name << "Hist for " << iX << " mm x < 0 Source " << NSOR;						   
			h_vert_real[iX + X_OBSERV_MAX] = new TH2D(h_name.str().c_str(), h_name.str().c_str(), YBINS, YMIN, YMAX, ZBINS, ZMIN, ZMAX);
			h_vert_real[iX + X_OBSERV_MAX]->SetMaximum(200);
			h_vert_real[iX + X_OBSERV_MAX]->GetXaxis()->SetRangeUser(YMIN, YMAX);
			h_vert_real[iX + X_OBSERV_MAX]->GetYaxis()->SetRangeUser(ZMIN, ZMAX);
			h_vert_real[iX + X_OBSERV_MAX]->GetXaxis()->SetTitle("y[mm]");
			h_vert_real[iX + X_OBSERV_MAX]->GetYaxis()->SetTitle("z[mm]");
			h_vert_real[iX + X_OBSERV_MAX]->GetZaxis()->SetTitle("N of events / 4 mm^{2}");
		}

		for(int i = 0; i < tr->GetEntries(); i++)
		{
			tr->GetEntry(i);

			for(int iX = X_OBSERV_MIN; iX < X_OBSERV_MAX+1; iX+=10)
			{
				verY = parA*iX + parB;
				verZ = parC*iX + parD;
				h_vert_real[iX + X_OBSERV_MAX]->Fill(verY, verZ);
			}
		}

		for(int iX = X_OBSERV_MIN; iX < X_OBSERV_MAX+1; iX+=10)
		{
			TCanvas* C0 = new TCanvas("Canvas", "Canvas", 1000, 700);
			h_vert_real[iX + X_OBSERV_MAX]->Draw("COLZ");
			
			C0->SetLogz();
			C0->Update();
			TPaveStats *st = (TPaveStats*)h_vert_real[iX + X_OBSERV_MAX]->FindObject("stats");
			st->SetX1NDC(0.7);			
			st->SetX2NDC(0.9);
                        st->SetY1NDC(0.72);
                        st->SetY2NDC(0.9);

			C0->SaveAs(Form("%sHISTOS/SRC%02dX%02d.png", PATH, NSOR, iX - X_OBSERV_MIN));
			
			//cout << "X = " << iX << "  " << h_vert_real[iX + X_OBSERV_MAX]->ProjectionX()->GetRMS() << "   " << h_vert_real[iX + X_OBSERV_MAX]->ProjectionY()->GetRMS() <<endl;
		}
		//h_vert_real[0]->Write();
	}
	//pr_hist_file->Close();
}



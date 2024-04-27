#include "/sps/nemo/scratch/ikovalen/TKEvent_old/TKEvent/include/TKEvent.h" 
#include "config.h"

R__LOAD_LIBRARY(/sps/nemo/scratch/ikovalen/TKEvent_old/TKEvent/lib/libTKEvent.so);

using namespace std;

void Final_Visualization()
{
	gROOT->SetBatch(kTRUE);


	// 0000000000000000000000000000
	// Build the RMS_Y(x) and RMS_Z(x) for selection ellipses from code Analysis

	TFile* file = new TFile(Form("%sClear_data_Run-%d.root", PATH, RUN_N));

	int N_Point = (double)(X_OBSERV_MAX - X_OBSERV_MIN) / STEP + 1;
	double Arr_X[N_Point], Arr_Y[N_Point], Arr_Z[N_Point];
	
	TGraph* grY[N_SRCPLN];
	TGraph* grZ[N_SRCPLN];

	for(int SRC_NO_i=0; SRC_NO_i < N_SRCPLN; SRC_NO_i++)
	{
		//stringstream tr_name;
		
		TString tr_name = Form("Tracks for %d mm %s Source %d", X_BasePlane, SIDE, SRC_NO_i);
		//tr_name << "Tracks for " << X_BasePlane << " mm x > 0 Source " << SRC_NO_i;
		TTree* tr = (TTree*) file->Get(tr_name);

		cout << tr_name << endl;		

		double parA, parB, parC, parD;
		tr->SetBranchAddress("A", &parA);
		tr->SetBranchAddress("B", &parB);
		tr->SetBranchAddress("C", &parC);
		tr->SetBranchAddress("D", &parD);

		double verY, verZ;

		// Fill histogram of undistorted vertices, and calculate coordinates of mode and the sigmas in both directions

		//TFile *New_file = new TFile("Visual_Histos.root","RECREATE");

//		double YMIN  = tr->GetMinimum("B") - 20.0;
//		double YMAX  = tr->GetMaximum("B") + 20.0;
//		int    YBINS = int( (YMAX - YMIN) / 2.0 );

//		double ZMIN = tr->GetMinimum("D") - 20.0;
//		double ZMAX = tr->GetMaximum("D") + 20.0;
//		int    ZBINS = int( (ZMAX - ZMIN) / 2.0 );

		double YMIN  = Y_MIN + (SRC_NO_i % N_COLS + 0.5) * Y_RECT_SIZE - 100;
		double YMAX  = Y_MIN + (SRC_NO_i % N_COLS + 0.5) * Y_RECT_SIZE + 100;
		int    YBINS = (YMAX - YMIN) / 2.0;

		double ZMIN  = Z_MAX - (SRC_NO_i / N_COLS + 0.5) * Z_RECT_SIZE - 70;
		double ZMAX  = Z_MAX - (SRC_NO_i / N_COLS + 0.5) * Z_RECT_SIZE + 70;
		int    ZBINS = (ZMAX - ZMIN) / 2.0;

		TH2D* h_vert_real[N_Point];  

		for(int point = 0; point < N_Point; point++)
		{
			//stringstream h_name;
			double X = X_OBSERV_MIN + STEP * point;			
			TString h_name = Form("Hist for %.1lf mm %s Source %d", X, SIDE, SRC_NO_i);
			//h_name << "Hist for " << X << " mm x < 0 Source " << SRC_NO_i;						   
			h_vert_real[point] = new TH2D(h_name, h_name, YBINS, YMIN, YMAX, ZBINS, ZMIN, ZMAX);
			h_vert_real[point]->SetMaximum(750);
			h_vert_real[point]->GetXaxis()->SetRangeUser(YMIN, YMAX);
			h_vert_real[point]->GetYaxis()->SetRangeUser(ZMIN, ZMAX);
			h_vert_real[point]->GetXaxis()->SetTitle("y[mm]");
			h_vert_real[point]->GetYaxis()->SetTitle("z[mm]");
		}

		for(int i = 0; i < tr->GetEntries(); i++)
		{
			tr->GetEntry(i);

			for(int point = 0; point < N_Point; point++)
			{
				double X = X_OBSERV_MIN + STEP * point;
				verY = parA * X + parB;
				verZ = parC * X + parD;
				h_vert_real[point]->Fill(verY, verZ);
			}
		}

		for(int point = 0; point < N_Point; point++)
		{
			double X = X_OBSERV_MIN + STEP * point;
			
			TCanvas* C0 = new TCanvas("Canvas", "Canvas");
			h_vert_real[point]->Draw("COLZ"); 
			//h_vert_real[point]->Write(h_name);

			C0->SetLogz();
			C0->Update();
			TPaveStats *st = (TPaveStats*)h_vert_real[point]->FindObject("stats");
			st->SetX1NDC(0.7);			
			st->SetX2NDC(0.9);
                        st->SetY1NDC(0.72);
                        st->SetY2NDC(0.9);

			C0->SaveAs(Form("%sCUTS/SRC%02d_X%3.1fmm.png", PATH, SRC_NO_i, X));
			
			cout << "X = " << X << "  " << h_vert_real[point]->ProjectionX()->GetRMS() << "   " << h_vert_real[point]->ProjectionY()->GetRMS() <<endl;

			Arr_X[point] = X;
			Arr_Y[point] = h_vert_real[point]->ProjectionX()->GetRMS();
			Arr_Z[point] = h_vert_real[point]->ProjectionY()->GetRMS();
			
			delete C0;
			delete st;
		}
		delete tr;
		for(auto& histo : h_vert_real) delete histo;

		grY[SRC_NO_i] = new TGraph(N_Point, Arr_X, Arr_Y);
		grZ[SRC_NO_i] = new TGraph(N_Point, Arr_X, Arr_Z);
	}
	
	TString SD_fname, position_fname;
	if (SIDE_NUM == 1)
	{
		SD_fname = Form("%sCWs_0mm_french.txt", PATH);
		position_fname = Form("%sSource_plane_position_french.txt", PATH);
	}
	else
	{
		SD_fname = Form("%sCWs_0mm_italian.txt", PATH);
		position_fname = Form("%sSource_plane_position_italian.txt", PATH);
	}
	ofstream SD_file;
	SD_file.open(SD_fname);

	ofstream position_file;
	position_file.open(position_fname);

	SD_file << "##########################################################" << "\n";
	SD_file  << "PATH: " << PATH << ", SIDE: " << SIDE << "\n";
	SD_file  << "##########################################################" << "\n";

	position_file << "##########################################################" << "\n";
	position_file  << "PATH: " << PATH << ", SIDE: " << SIDE << "\n";
	position_file << "##########################################################" << "\n";


	TCanvas* CY = new TCanvas("#Delta y(x)", "GraphY", 1200, 800);

	TPad *pad1y = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
	pad1y->Draw();
	pad1y->cd();

	pad1y->SetBottomMargin(0);

        auto legendy = new TLegend(0.85,0.0,1.0,1.0);
        legendy->SetHeader("#Delta y(x)","C"); // option "C" allows to center the header

	grY[0]->SetMinimum(5.1);
	grY[0]->SetMaximum(12.0);
	grY[0]->SetName("Src 0");
	grY[0]->SetMarkerColor(1);
	grY[0]->SetMarkerStyle(20);
	grY[0]->SetMarkerSize(1);
	grY[0]->SetLineColor(1);
	grY[0]->SetLineWidth(3);
	grY[0]->Draw("APL");
	legendy->AddEntry("Src 0","Src 0","pl");

	grY[0]->GetYaxis()->SetNdivisions(510);
	grY[0]->GetYaxis()->SetLabelFont(43);
	grY[0]->GetYaxis()->SetTitleFont(43);
	grY[0]->GetYaxis()->SetTitleOffset(1);
	grY[0]->GetYaxis()->SetLabelSize(20);
	grY[0]->GetYaxis()->SetTitleSize(20);
	grY[0]->GetYaxis()->SetTitle("#sigma_y");

	grY[0]->GetXaxis()->SetLabelFont(43);
	grY[0]->GetXaxis()->SetTitleFont(43);
	grY[0]->GetXaxis()->SetTitleOffset(3);
	grY[0]->GetXaxis()->SetLabelSize(20);	
	grY[0]->GetXaxis()->SetTitleSize(20);
	grY[0]->GetXaxis()->SetTitle("x [mm]");
	grY[0]->SetTitle("#Delta y (only French side tracks)");

	for(int SRC_NO_i=1; SRC_NO_i<N_SRCPLN; SRC_NO_i++)
	{
	    	char nm[7];
	    	sprintf(nm,"Src %i", SRC_NO_i);

		grY[SRC_NO_i]->SetName(nm);
		grY[SRC_NO_i]->SetMarkerColor(1);
		grY[SRC_NO_i]->SetMarkerStyle(SRC_NO_i % 4 + 20);
		grY[SRC_NO_i]->SetMarkerSize(1);
		grY[SRC_NO_i]->SetLineColor((int)SRC_NO_i / 6 + 1);
		grY[SRC_NO_i]->SetLineWidth(3);
		grY[SRC_NO_i]->Draw("PL same");

		legendy->AddEntry(nm, nm, "pl");
	}

//	legendy->SetTextSize(14);	
//	legendy->Draw("same");	

	CY->cd();
	TPad *pad2y = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
	pad2y->Draw();
	pad2y->cd();

	pad2y->SetBottomMargin(0.2);
	pad2y->SetTopMargin(0);

	TH1F *hY = new TH1F("histo", "histo", 41, X_OBSERV_MIN - 0.5, X_OBSERV_MAX + 0.5);

	hY->SetStats(0);
	hY->SetTitle("");

	hY->GetYaxis()->SetTitleSize(20);
	hY->GetYaxis()->SetLabelSize(15);
	hY->GetYaxis()->SetTitleFont(43);
	hY->GetYaxis()->SetLabelFont(43);
	hY->GetYaxis()->SetTitleOffset(1);
	hY->GetYaxis()->SetNdivisions(10);
	hY->GetYaxis()->SetTitle("N");

	hY->GetXaxis()->SetTickLength(0.1);
	hY->GetXaxis()->SetTitleSize(20);
	hY->GetXaxis()->SetLabelSize(15);
	hY->GetXaxis()->SetTitleFont(43);
	hY->GetXaxis()->SetLabelFont(43);
	hY->GetXaxis()->SetTitleOffset(3);
	hY->GetXaxis()->SetTitle("x [mm]");

	for(int SRC_NO_i=0; SRC_NO_i<N_SRCPLN; SRC_NO_i++)
	{
	    // Get the graph's X and Y arrays
	    const double* xValues = grY[SRC_NO_i]->GetX();
	    const double* yValues = grY[SRC_NO_i]->GetY();
	
	    // Find the index of the minimum Y value
	    double minY = yValues[0];
	    int minIndex = 0;
	    for (int i = 1; i < grY[SRC_NO_i]->GetN(); ++i)
	    {
	        if (yValues[i] < minY)
	        {
	            minY = yValues[i];
	            minIndex = i;
	        }
	    }
	    SD_file << "Source " << SRC_NO_i << " SD_Y_0mm = " << yValues[(int)grY[SRC_NO_i]->GetN() / 2 + 1] << "\n";
	
		    // Get the x-position corresponding to the minimum Y value
	    double minXPosition = xValues[minIndex];
	
	    hY->Fill(minXPosition);
	}

	position_file << "Mean_Y: " << hY->GetMean() << ", SD: " << (N_SRCPLN / (N_SRCPLN - 1)) * hY->GetRMS() << "\n";

	
	hY->Draw();

	CY->SetTitle("#delta y(x)");

	CY->Print(Form("%svY.png", PATH));

        TCanvas* CZ = new TCanvas("#Delta z(x)", "GraphZ", 1200, 800);

        TPad *pad1z = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
        pad1z->Draw();
        pad1z->cd();

        pad1z->SetBottomMargin(0);

        auto legendz = new TLegend(0.9,0.1,1.0,0.9);
        legendz->SetHeader("#Delta z(x)","C"); // option "C" allows to center the header

        grZ[0]->SetMinimum(8);
        grZ[0]->SetMaximum(14.0);
        grZ[0]->SetName("Src 0");
        grZ[0]->SetMarkerColor(1);
        grZ[0]->SetMarkerStyle(20);
        grZ[0]->SetMarkerSize(1);
        grZ[0]->SetLineColor(1);
        grZ[0]->SetLineWidth(3);
        grZ[0]->Draw("APL");
        legendz->AddEntry(grZ[0],"Src 0","pl");

        grZ[0]->GetYaxis()->SetNdivisions(510);
        grZ[0]->GetYaxis()->SetLabelFont(43);
        grZ[0]->GetYaxis()->SetTitleFont(43);
        grZ[0]->GetYaxis()->SetTitleOffset(1);
        grZ[0]->GetYaxis()->SetLabelSize(20);
        grZ[0]->GetYaxis()->SetTitleSize(20);
        grZ[0]->GetYaxis()->SetTitle("#sigma_z");

        grZ[0]->GetXaxis()->SetLabelFont(43);
        grZ[0]->GetXaxis()->SetTitleFont(43);
        grZ[0]->GetXaxis()->SetTitleOffset(3);
        grZ[0]->GetXaxis()->SetLabelSize(20);
        grZ[0]->GetXaxis()->SetTitleSize(20);
        grZ[0]->GetXaxis()->SetTitle("x [mm]");
	grZ[0]->SetTitle("#Delta z (only French side tracks)");

        for(int SRC_NO_i=1; SRC_NO_i<N_SRCPLN; SRC_NO_i++)
        {
                char nm[7];
                sprintf(nm,"Src %i", SRC_NO_i);

                grZ[SRC_NO_i]->SetName(nm);
                grZ[SRC_NO_i]->SetMarkerColor(1);
                grZ[SRC_NO_i]->SetMarkerStyle(SRC_NO_i % 4 + 20);
                grZ[SRC_NO_i]->SetMarkerSize(1);
                grZ[SRC_NO_i]->SetLineColor((int)SRC_NO_i / 2 + 1);
                grZ[SRC_NO_i]->SetLineWidth(3);
                grZ[SRC_NO_i]->Draw("PL same");

                legendz->AddEntry(grZ[SRC_NO_i], nm, "pl");
        }

//	legendz->SetTextSize(16);
//        legendz->Draw("same");

        CZ->cd();
        TPad *pad2z = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
        pad2z->Draw();
        pad2z->cd();

        pad2z->SetBottomMargin(0.2);
        pad2z->SetTopMargin(0);

        TH1F *hZ = new TH1F("histo", "histo", 41, X_OBSERV_MIN - 0.5, X_OBSERV_MAX + 0.5);

        hZ->SetStats(0);
        hZ->SetTitle("");
        hZ->GetYaxis()->SetTitleSize(20);
        hZ->GetYaxis()->SetLabelSize(15);
        hZ->GetYaxis()->SetTitleFont(43);
        hZ->GetYaxis()->SetLabelFont(43);
        hZ->GetYaxis()->SetTitleOffset(1);
        hZ->GetYaxis()->SetNdivisions(10);
        hZ->GetYaxis()->SetTitle("N");

        hZ->GetXaxis()->SetTickLength(0.1);
        hZ->GetXaxis()->SetTitleSize(20);
        hZ->GetXaxis()->SetLabelSize(15);
        hZ->GetXaxis()->SetTitleFont(43);
        hZ->GetXaxis()->SetLabelFont(43);
        hZ->GetXaxis()->SetTitleOffset(3);
        hZ->GetXaxis()->SetTitle("x [mm]");

        for(int SRC_NO_i=0; SRC_NO_i<N_SRCPLN; SRC_NO_i++)
        {
            // Get the graph's X and Y arrays
            const double* xValues = grZ[SRC_NO_i]->GetX();
            const double* zValues = grZ[SRC_NO_i]->GetY();

            // Find the index of the minimum Y value
            double minZ = zValues[0];
            int minIndex = 0;
            for (int i = 1; i < grZ[SRC_NO_i]->GetN(); ++i)
            {
                if (zValues[i] < minZ)
                {
                    minZ = zValues[i];
                    minIndex = i;
                }
            }
	    SD_file << "Source " << SRC_NO_i << " SD_Y_0mm = " << zValues[(int)grY[SRC_NO_i]->GetN() / 2 + 1] << "\n";
            // Get the x-position corresponding to the minimum Y value
            double minXPosition = xValues[minIndex];

            hZ->Fill(minXPosition);
        }

	position_file << "Mean_Z: " << hZ->GetMean() << ", SD: " << (N_SRCPLN / (N_SRCPLN - 1)) * hZ->GetRMS() << "\n";

       	SD_file << "##########################################################" << "\n";
       	position_file << "##########################################################" << "\n";

	SD_file.close();
	position_file.close();

	hZ->Draw();

        CZ->SetTitle("#delta z(x)");

        CZ->Print(Form("%svZ.png", PATH));

	for(int SRC_NO_i=0; SRC_NO_i<N_SRCPLN; SRC_NO_i++)
	{	
		TString print_name_RMS_Y = Form("%sRMS_ALL_SOURCE/RMS_Y_Source_%02d.png", PATH, SRC_NO_i);

		TCanvas* CindY = new TCanvas(print_name_RMS_Y, print_name_RMS_Y);
		grY[SRC_NO_i]->Draw();

		CindY->Print(print_name_RMS_Y);

		TString print_name_RMS_Z = Form("%sRMS_ALL_SOURCE/RMS_Z_Source_%02d.png", PATH, SRC_NO_i);

		TCanvas* CindZ = new TCanvas(print_name_RMS_Z, print_name_RMS_Z);
		grZ[SRC_NO_i]->Draw();

		CindZ->Print(print_name_RMS_Z);
	}
}



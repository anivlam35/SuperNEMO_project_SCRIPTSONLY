#include "../../../TKEvent/TKEvent/include/TKEvent.h"
#include "config.h"

using namespace std;

//R__LOAD_LIBRARY(/sps/nemo/scratch/ikovalen/TKEvent_old/TKEvent/lib/libTKEvent.so);

void OM_xyz_swcr(int OM_num);

void Visualization(){
	TFile* f = new TFile(Form("%sOMs_Tracks_Run-%d.root", PATH, RUN_N));
	
	TTree*   Tree[N_OMs];
	double A_Tree[N_OMs];
	double B_Tree[N_OMs];
	double C_Tree[N_OMs];
	double D_Tree[N_OMs];

        int    Y_bins_general = int(Y_MAX - Y_MIN) / 3;
        int    Z_bins_general = int(Z_MAX - Z_MIN) / 3;
		
	for(int OM_num = 0; OM_num < N_OMs; OM_num++)
	{
		TString treename = Form("Tracks of OM %d", OM_num);		

		Tree[OM_num] = (TTree*) f->Get(treename);
		Tree[OM_num]->SetBranchAddress("A", &A_Tree[OM_num]);
		Tree[OM_num]->SetBranchAddress("B", &B_Tree[OM_num]);
		Tree[OM_num]->SetBranchAddress("C", &C_Tree[OM_num]);
		Tree[OM_num]->SetBranchAddress("D", &D_Tree[OM_num]);
	}

	for(int OM_num = 0; OM_num < 260; OM_num++)
	{
		OM_xyz_swcr(OM_num);
		
		double X_zero_plane = xyz[0];
		
		if (X_zero_plane > 0) X_zero_plane -= mw_sizex / 2 ;
		else X_zero_plane += mw_sizex / 2 ; 

		double Y_hmin = (xyz[1] - mw_sizey / 2) - 200; 		
		double Y_hmax = (xyz[1] + mw_sizey / 2) + 200; 		
		int    Y_bins = int(Y_hmax - Y_hmin) / 4;
		
		double Z_hmin = (xyz[2] - mw_sizez / 2) - 200;
		double Z_hmax = (xyz[2] + mw_sizez / 2) + 200;
		int    Z_bins = int(Z_hmax - Z_hmin) / 4;
	
		TCanvas* C = new TCanvas("Canvas", "Canvas", 1150, 1000);

                C->SetLeftMargin(0.15);
                C->SetRightMargin(0.15);

		TString hname = Form("Tracks of OM#%d", OM_num);
		TH2D* h = new TH2D(hname, hname, Y_bins, Y_hmin, Y_hmax, Z_bins, Z_hmin, Z_hmax);
		
		h->GetXaxis()->SetTitle("y[mm]");
		h->GetYaxis()->SetTitle("z[mm]");
		h->GetZaxis()->SetTitle("N / 16mm^2");
                h->SetTitleOffset(1);
		h->SetStats(0);
		h->SetMinimum(0.001);
		
		cout << Tree[OM_num]->GetEntries() << endl;

		for(int entry = 0; entry < Tree[OM_num]->GetEntries(); entry++)
		{
			Tree[OM_num]->GetEntry(entry);
			
			double Y = A_Tree[OM_num] * X_zero_plane + B_Tree[OM_num];
			double Z = C_Tree[OM_num] * X_zero_plane + D_Tree[OM_num];  
			h->Fill(Y, Z);
		}
		h->Draw("COLZ");

		auto h_line1 = new TLine(Y_hmin, xyz[2] - mw_sizey / 2, Y_hmax, xyz[2] - mw_sizey / 2);
		auto h_line2 = new TLine(Y_hmin, xyz[2] + mw_sizey / 2, Y_hmax, xyz[2] + mw_sizey / 2);
		auto v_line1 = new TLine(xyz[1] - mw_sizez / 2, Z_hmin, xyz[1] - mw_sizez / 2, Z_hmax);
		auto v_line2 = new TLine(xyz[1] + mw_sizez / 2, Z_hmin, xyz[1] + mw_sizez / 2, Z_hmax);
		
		h_line1->SetLineWidth(3);
		h_line2->SetLineWidth(3);
		v_line1->SetLineWidth(3);
		v_line2->SetLineWidth(3);
				
		h_line1->SetLineColor(2);
		h_line2->SetLineColor(2);
		v_line1->SetLineColor(2);
		v_line2->SetLineColor(2);

		h_line1->Draw("Same");
		h_line2->Draw("Same");
		v_line1->Draw("Same");
		v_line2->Draw("Same");
		
		hname = Form("%svisu_OMs_tracks/OM%03d_tracks_Run-%d.png", PATH, OM_num, RUN_N);
		C->SaveAs(hname);

///////////////////////////// GENERAL VIEW ///////////////////////////

		TCanvas* C2 = new TCanvas("Canvas", "Canvas", 2000, 1300);
		
		TString hname2 = Form("Tracks of OM#%03d (x < 0). General view.", OM_num);
        	TH2D* h2 = new TH2D(hname2, hname2, Y_bins_general, Y_MIN, Y_MAX, Z_bins_general, Z_MIN, Z_MAX);

		h2->SetStats(0);
		h2->GetXaxis()->SetTitle("y[mm]");
		h2->GetYaxis()->SetTitle("z[mm]");
                h2->GetYaxis()->SetTitleOffset(1.2);

                for(int entry = 0; entry < Tree[OM_num]->GetEntries(); entry++)
                {
                        Tree[OM_num]->GetEntry(entry);

                        double Y = A_Tree[OM_num] * X_zero_plane + B_Tree[OM_num];
                        double Z = C_Tree[OM_num] * X_zero_plane + D_Tree[OM_num];
                        h2->Fill(Y, Z);
                }
                h2->Draw("COLZ");

		for(int i = 0; i <= 13; i++)
		{
			auto h_line = new TLine(Y_MIN, Z_MIN + i * mw_sizez, Y_MAX, Z_MIN + i * mw_sizez);
			h_line->SetLineWidth(1);
			h_line->SetLineColor(2);
			h_line->Draw("Same");
		}
        	for(int	i = 0; i <= 20; i++)
        	{
                	auto v_line = new TLine(Y_MIN + i * mw_sizey, Z_MIN, Y_MIN + i * mw_sizey, Z_MAX);
                	v_line->SetLineWidth(1);
                	v_line->SetLineColor(2);
                	v_line->Draw("Same");
        	}

		hname2 = Form("%svisu_OMs_tracks/OM%03d_tracks_general_Run-%d.png", PATH, OM_num, RUN_N);
                C2->SaveAs(hname2);

		cout << Form("OM %d is done!", OM_num) << endl;
	}


}

void OM_xyz_swcr(int OM_num)
{
        //mainwall IT
        if(OM_num < 260)
        {
                SWCR[0] = 0;
                SWCR[1] = -1;
                SWCR[2] = OM_num / 13;
                SWCR[3] = OM_num % 13;
        }
	//mainwall FR
        else if(OM_num < 520)
        {
                SWCR[0] = 1;
                SWCR[1] = -1;
                SWCR[2] = (OM_num - 260) / 13;
                SWCR[3] = (OM_num - 260) % 13;
        }
	//Xcalo IT
        else if(OM_num < 584)
        {
                SWCR[0] = 0;
                SWCR[1] = (OM_num - 520) / 32;
                SWCR[2] = ((OM_num - 520) / 16) % 2;
                SWCR[3] = (OM_num -520) % 16;
        }
	//Xcalo FR
        else if(OM_num < 648)
        {
                SWCR[0] = 1;
                SWCR[1] = (OM_num - 520 - 64) / 32;
                SWCR[2] = ((OM_num - 520 - 64) / 16) % 2;
                SWCR[3] = (OM_num -520 - 64) % 16;
        }
	//GVeto IT
        else if(OM_num < 680)
        {
                SWCR[0] = 0;
                SWCR[1] = (OM_num - 520 - 128) / 16;
                SWCR[2] = (OM_num - 520 - 128) % 16;
                SWCR[3] = -1;
        }
	//GVeto FR
        else if(OM_num < 712)
        {
                SWCR[0] = 1;
                SWCR[1] = (OM_num - 520 - 128 - 32) / 16;
                SWCR[2] = (OM_num - 520 - 128 - 32) % 16;
                SWCR[3] = -1;
        }

	int OM_type;

        if(OM_num < 520)
        {
                OM_type = 1302;
        }
	else if(OM_num < 648)
        {
                OM_type = 1232;
        }
	else
	{
                OM_type = 1252;
        }

	switch(OM_type)
        {
                case 1302: //MW
                        if(SWCR[0] == 1)
                                xyz[0] = 532.0;
                        else
                            	xyz[0] = -532.0;
                                xyz[1] = ((double)SWCR[2]- 9.5) * 259.0;
                                xyz[2] = ((double)SWCR[3] - 6) * 259.0;

                        break;

                case 1232: //XW
                        if(SWCR[1] == 1)
                                xyz[1] = 2580.5;
                        else
                            	xyz[1] = -2580.5;

                        if(SWCR[0] == 1)
                        {
                                if(SWCR[2] == 1)
                                        xyz[0] = 333.0;
                                else
                                    	xyz[0] = 130.0;
                        }
                        else
                        {
                                if(SWCR[2] == 1)
                                        xyz[0] = -333.0;
                                else
                                    	xyz[0] = -130.0;
                        }

                        xyz[2] = ((double)SWCR[3] - 7.5) * 212.0;

                        break;

                case 1252: //GV
                        if(SWCR[0] == 1)
                                xyz[0] = 213.5;
                        else
                            	xyz[0] = -213.5;
                        if(SWCR[1] == 1)
                                xyz[2] = 1625.0;
                        else
                            	xyz[2] = -1625.0;
                        if(SWCR[2] > 7)
                                xyz[1] = 161.0 + (((double)SWCR[2]-8) * 311.5);
                        else
                            	xyz[1] = -161.0 + (((double)SWCR[2]-7) * 311.5);
                        break;
        }
}


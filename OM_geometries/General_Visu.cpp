#include "config.h"

using namespace std;
using namespace TMath;

// R__LOAD_LIBRARY(/sps/nemo/scratch/ikovalen/TKEvent_old/TKEvent/lib/libTKEvent.so);

void OM_xyz_swcr(int OM_num);

void General_Visu(){
	TFile* f = new TFile(Form("%sOMs_Tracks_Run-%d.root", PATH, RUN_N));
	
	TTree*   Tree[N_OMs];
	double A_Tree[N_OMs];
	double B_Tree[N_OMs];
	double C_Tree[N_OMs];
	double D_Tree[N_OMs];
		
	for(int OM_num = 0; OM_num < N_OMs; OM_num++)
	{
		TString treename = Form("Tracks of OM %d", OM_num);		

		Tree[OM_num] = (TTree*) f->Get(treename);
		Tree[OM_num]->SetBranchAddress("A", &A_Tree[OM_num]);
		Tree[OM_num]->SetBranchAddress("B", &B_Tree[OM_num]);
		Tree[OM_num]->SetBranchAddress("C", &C_Tree[OM_num]);
		Tree[OM_num]->SetBranchAddress("D", &D_Tree[OM_num]);
	}

        int    Y_bins = int(Y_MAX - Y_MIN) / 5;
        int    Z_bins = int(Z_MAX - Z_MIN) / 5;
	double X_zero_plane;

/////////////////////////////////////// ITALIAN SIDE ///////////////////////////////////////////////////////////

        TCanvas* C0 = new TCanvas("Canvas", "Canvas", 2000, 1300);

        TString hname0 = "Tracks of Main Wall OMs (x < 0)";
        TH2D* h0 = new TH2D(hname0, hname0, Y_bins, Y_MIN, Y_MAX, Z_bins, Z_MIN, Z_MAX);
	
	X_zero_plane = -532 + mw_sizex / 2;

	h0->SetStats(0);
	h0->GetXaxis()->SetTitle("y[mm]");
	h0->GetYaxis()->SetTitle("z[mm]");
	
	for(int OM_num = 0; OM_num < 260; OM_num++)
	{
		for(int entry = 0; entry < Tree[OM_num]->GetEntries(); entry++)
		{
			Tree[OM_num]->GetEntry(entry);
			
			double Y = A_Tree[OM_num] * X_zero_plane + B_Tree[OM_num];
			double Z = C_Tree[OM_num] * X_zero_plane + D_Tree[OM_num];  
			h0->Fill(Y, Z);
		}
		h0->Draw("COLZ");
	}

	for(int i = 0; i <= 13; i++)
	{
		auto h_line = new TLine(Y_MIN, Z_MIN + i * mw_sizez, Y_MAX, Z_MIN + i * mw_sizez);
		h_line->SetLineWidth(1);
		h_line->SetLineColor(2);
		// h_line->Draw("Same");
	}
        for(int	i = 0; i <= 20; i++)
        {
                auto v_line = new TLine(Y_MIN + i * mw_sizey, Z_MIN, Y_MIN + i * mw_sizey, Z_MAX);
                v_line->SetLineWidth(1);
                v_line->SetLineColor(2);
                // v_line->Draw("Same");
        }
        C0->SaveAs(Form("%shist_OMs_tracks_ital.png", PATH));

///////////////////////////////////////////// FRENCH SIDE //////////////////////////////////////////////////////

        TCanvas* C1 = new TCanvas("Canvas", "Canvas", 2000, 1300);

        TString hname1 = "Tracks of Main Wall OMs (x > 0)";
        TH2D* h1 = new TH2D(hname1, hname1, Y_bins, Y_MIN, Y_MAX, Z_bins, Z_MIN, Z_MAX);

	X_zero_plane = 532 - mw_sizex / 2;

	h1->SetStats(0);
	h1->GetXaxis()->SetTitle("y[mm]");
	h1->GetYaxis()->SetTitle("z[mm]");

        for(int OM_num = 260; OM_num < 520; OM_num++)
        {
                for(int entry = 0; entry < Tree[OM_num]->GetEntries(); entry++)
                {
                        Tree[OM_num]->GetEntry(entry);

                        double Y = A_Tree[OM_num] * X_zero_plane + B_Tree[OM_num];
                        double Z = C_Tree[OM_num] * X_zero_plane + D_Tree[OM_num];
                        h1->Fill(Y, Z);
                }
                h1->Draw("COLZ");
        }

	for(int i = 0; i <= 13; i++)
        {
                auto h_line = new TLine(Y_MIN, Z_MIN + i * mw_sizez, Y_MAX, Z_MIN + i * mw_sizez);
                h_line->SetLineWidth(1);
                h_line->SetLineColor(2);
                // h_line->Draw("Same");
        }
	for(int i = 0; i <= 20; i++)
        {
                auto v_line = new TLine(Y_MIN + i * mw_sizey, Z_MIN, Y_MIN + i * mw_sizey, Z_MAX);
                v_line->SetLineWidth(1);
                v_line->SetLineColor(2);
                // v_line->Draw("Same");
        }
	C1->SaveAs(Form("%shist_OMs_tracks_fren.png", PATH));

////////////////////////////// VECTOR FIELD ITALIAN SIDE //////////////////////////////////////////////////////
	
	TArrow* arrows0[260];
	X_zero_plane = -532 + mw_sizex / 2;

	for(int OM_num = 0; OM_num < 260; OM_num++)
	{
		OM_xyz_swcr(OM_num);
		
		double sum_Y = 0;
		double sum_Z = 0;

		for(int i = 0; i < Tree[OM_num]->GetEntries(); i++)
		{
			Tree[OM_num]->GetEntry(i);
			sum_Y += A_Tree[OM_num] * X_zero_plane + B_Tree[OM_num];
			sum_Z += C_Tree[OM_num] * X_zero_plane + D_Tree[OM_num];
		}
		
		double mean_Y = sum_Y / Tree[OM_num]->GetEntries();
		double mean_Z = sum_Z / Tree[OM_num]->GetEntries();

		arrows0[OM_num] = new TArrow(xyz[1], xyz[2], mean_Y, mean_Z, 0.005, "|>");
		//float arrow_l = Sqrt(Power(mean_Y - xyz[1], 2) + Power(mean_Z - xyz[2], 2)) / 300;
		//TColor* color;
		//color->SetRGB(arrow_l, 1 - arrow_l, 0);
		arrows0[OM_num]->SetLineWidth(2);
		//arrow->SetLineColor(color);
	}
        TCanvas* C2 = new TCanvas("Canvas", "Canvas", 2000, 1300);

        TString hname2 = "Vector Field of Main Wall OMs (x < 0)";
        TH2D* h2 = new TH2D(hname2, hname2, Y_bins, Y_MIN, Y_MAX, Z_bins, Z_MIN, Z_MAX);

        h2->SetStats(0);
        h2->GetXaxis()->SetTitle("y[mm]");
        h2->GetYaxis()->SetTitle("z[mm]");

        h2->Draw();

	for(int OM_num = 0; OM_num < 260; OM_num++)
	{
		arrows0[OM_num]->Draw();
	}


	C2->SaveAs(Form("%svect_field_OMs_displ_ital.png", PATH));

////////////////////////////// VECTOR FIELD FRENCH SIDE //////////////////////////////////////////////////////

	TArrow* arrows1[260];
	X_zero_plane = 532 - mw_sizex / 2;

	for(int OM_num = 260; OM_num < 520; OM_num++)
	{
		OM_xyz_swcr(OM_num);
		
		double sum_Y = 0;
		double sum_Z = 0;

		for(int i = 0; i < Tree[OM_num]->GetEntries(); i++)
		{
			Tree[OM_num]->GetEntry(i);
			sum_Y += A_Tree[OM_num] * X_zero_plane + B_Tree[OM_num];
			sum_Z += C_Tree[OM_num] * X_zero_plane + D_Tree[OM_num];
		}
		
		double mean_Y = sum_Y / Tree[OM_num]->GetEntries();
		double mean_Z = sum_Z / Tree[OM_num]->GetEntries();

		arrows1[OM_num - 260] = new TArrow(xyz[1], xyz[2], mean_Y, mean_Z, 0.005, "|>");
		//float arrow_l = Sqrt(Power(mean_Y - xyz[1], 2) + Power(mean_Z - xyz[2], 2)) / 300;
		//TColor* color;
		//color->SetRGB(arrow_l, 1 - arrow_l, 0);
		arrows1[OM_num - 260]->SetLineWidth(2);
		//arrow->SetLineColor(color);
	}
        TCanvas* C3 = new TCanvas("Canvas", "Canvas", 2000, 1300);

        TString hname3 = "Vector Field of Main Wall OMs (x < 0)";
        TH2D* h3 = new TH2D(hname3, hname3, Y_bins, Y_MIN, Y_MAX, Z_bins, Z_MIN, Z_MAX);

        h3->SetStats(0);
        h3->GetXaxis()->SetTitle("y[mm]");
        h3->GetYaxis()->SetTitle("z[mm]");

        h3->Draw();

	for(int OM_num = 260; OM_num < 520; OM_num++)
	{
		arrows1[OM_num - 260]->Draw();
	}


	C3->SaveAs(Form("%svect_field_OMs_displ_french.png", PATH));

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


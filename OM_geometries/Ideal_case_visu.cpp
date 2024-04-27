#include <iostream>
#include "config.h"
using namespace std;

void OM_xyz_swcr(int OM_num);

void Ideal_case_visu()
{
    TFile* f = new TFile(Form("cutted_data_2770mm/OMs_Tracks_Run-%d.root", RUN_N));    

    TTree*   Tree;
	double A_Tree;
	double B_Tree;
	double C_Tree;
	double D_Tree;

    int OM_num = 138;
    OM_xyz_swcr(OM_num);

    TString treename = Form("Tracks of OM %d", OM_num);		

    Tree = (TTree*) f->Get(treename);
    cout << Tree->GetEntries() << endl;
    Tree->SetBranchAddress("A", &A_Tree);
    Tree->SetBranchAddress("B", &B_Tree);
    Tree->SetBranchAddress("C", &C_Tree);
    Tree->SetBranchAddress("D", &D_Tree);

    double Y_hmin = (xyz[1] - mw_sizey / 2) - 200; 		
    double Y_hmax = (xyz[1] + mw_sizey / 2) + 200; 		
    int    Y_bins = int(Y_hmax - Y_hmin) / 1;
    

    TCanvas* C = new TCanvas("Canvas_projections", "Canvas_projections", 800, 600);

    TString hname1 = "Real histo";
    TString hname2 = "Ideal histo";
    TH1D* h1 = new TH1D(hname1, hname1, Y_bins, Y_hmin, Y_hmax);
    TH1D* h2 = new TH1D(hname2, hname2, Y_bins, Y_hmin, Y_hmax);

    double X;
    double X_zero_plane = xyz[0];
        
    if (X_zero_plane > 0) X = X_zero_plane - mw_sizex / 2;
    else X = X_zero_plane + mw_sizex / 2; 
    for(int entry = 0; entry < Tree->GetEntries(); entry++)
    {
        Tree->GetEntry(entry);

        double y = A_Tree * X + B_Tree;
        h1->Fill(y);
    }
    h1->GetXaxis()->SetTitle("y[mm]");
    h1->GetYaxis()->SetTitle("Counts");
    h1->SetTitle(Form("Y projection (real case) OM#%03d X=%.1lfmm", OM_num, X));
    h1->SetStats(0);
    h1->Draw();
    C->SaveAs("Y_projection_real_example.png");



    for (int i = Y_hmin + 0.5; i < Y_hmax; i++)
    {
        if ((i > xyz[1] - mw_sizey / 2) && (i < xyz[1] + mw_sizey / 2))
        {
            for (int j=0; j < 50; j++)
            h2->Fill(i);
        }

    }
    h2->GetXaxis()->SetTitle("y[mm]");
    h2->GetYaxis()->SetTitle("Counts");
    h2->SetTitle(Form("Y projection (ideal case) OM#%03d X=%.1lfmm", OM_num, X));
    h2->SetStats(0);
    h2->Draw();
    C->SaveAs("Y_projection_ideal_example.png");

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
#include <iostream>
#include "config.h"
using namespace std;

void OM_xyz_swcr(int OM_num);
TH2D* create_histo(TTree* _t, int _OM_num, double _X_plane_shift);
int mean_counts(TH1D* h_pr);
double* find_YZ_bounds(TH2D* _h, int _OM_num, int _point, int _DrawOption = 0);
void draw_projections(TH1D* _h_pr_y, TH1D* _h_pr_z, int _mean_Y, int _mean_Z, int _OM_num, int _point);
void draw_projections(TGraph* _gr_y, TGraph* _gr_z, int _mean_Y, int _mean_Z, int _OM_num, int _point);
TGraph* running_average(TH1D* _h, int _run_av_hwidth);
TTree* filter_tree(TTree* _t, int _OM_num, double _bounds0[4], double _bounds1[4]);
double find_intersection(TGraph* _gr, double mean_val, TString _option);


void Filter_data(){
    TFile* f = new TFile(Form("%sOMs_Tracks_Run-%d.root", PATH, RUN_N));
    TFile* clear_data_file = new TFile(Form("%sOMs_Clear_Tracks_Run-%d.root", PATH, RUN_N), "RECREATE");

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

    int OM_num;
    TTree* OM_num_tree = new TTree("List of OMs numbers", "List of OMs numbers");
    OM_num_tree->Branch("OM_num", &OM_num);

    for(OM_num = 0; OM_num < 520; OM_num++)
    {   
        if (OM_num < 13             ||
            OM_num > 246            ||
            OM_num % 13 == 0        ||
            (OM_num - 1) % 13 == 0  ||
            (OM_num + 12) % 13 == 0 ||
            (OM_num + 11) % 13 == 0)
                continue;
        if (Tree[OM_num]->GetEntries() > 1000) 
        {
            OM_num_tree->Fill();
            cout << OM_num << endl;
        
            TH2D* h0 = create_histo(Tree[OM_num], OM_num, X_shifts[0]);
            double* bounds0 = find_YZ_bounds(h0, OM_num, 0, 0);
            TH2D* h1 = create_histo(Tree[OM_num], OM_num, X_shifts[1]);
            double* bounds1 = find_YZ_bounds(h1, OM_num, 0, 0);


            TTree* clear_tree = filter_tree(Tree[OM_num], OM_num, bounds0, bounds1);
            clear_tree->Write(Form("Clear tracks of OM %d", OM_num));

            cout << clear_tree->GetEntries() << endl;

            TH2D* h_clear = create_histo(clear_tree, OM_num, 0);
            TCanvas* C_clear = new TCanvas("Canvas", "Canvas", 1000, 1000);
            h_clear->Draw("colz");
            C_clear->SaveAs(Form("%sClear_tracks/Tracks_of_OM_%03d_Run-%d.png", PATH, OM_num, RUN_N));
        }
    }
    OM_num_tree->Write("List of OMs numbers");
    // OM_xyz_swcr(OM_num);

    // double X_zero_plane = xyz[0];
		
    // if (X_zero_plane > 0) X_zero_plane -= mw_sizex / 2  - 20;
    // else X_zero_plane += mw_sizex / 2 + 20; 

    // double Y_hmin = (xyz[1] - mw_sizey / 2) - 200; 		
    // double Y_hmax = (xyz[1] + mw_sizey / 2) + 200; 		
    // int    Y_bins = int(Y_hmax - Y_hmin) / 4;
    
    // double Z_hmin = (xyz[2] - mw_sizez / 2) - 200;
    // double Z_hmax = (xyz[2] + mw_sizez / 2) + 200;
    // int    Z_bins = int(Z_hmax - Z_hmin) / 4;

    // TCanvas* C = new TCanvas("Canvas", "Canvas", 1000, 1000);

    // TString hname = Form("Tracks of OM#%d", OM_num);
    // TH2D* h = new TH2D(hname, hname, Y_bins, Y_hmin, Y_hmax, Z_bins, Z_hmin, Z_hmax);
    
    // h->GetXaxis()->SetTitle("y[mm]");
    // h->GetYaxis()->SetTitle("z[mm]");
    // h->SetStats(0);
    // h->SetMinimum(0.001);

	// for(int entry = 0; entry < Tree[OM_num]->GetEntries(); entry++)
	// {
    //     Tree[OM_num]->GetEntry(entry);
        
    //     double Y = A_Tree[OM_num] * X_zero_plane + B_Tree[OM_num];
    //     double Z = C_Tree[OM_num] * X_zero_plane + D_Tree[OM_num];  
    //     h->Fill(Y, Z);
	// }

    // h->Draw("colz");
    // C->SaveAs("Tracks_of_OM_140.png");


    // TH1D* h_X = h->ProjectionX();
    // TH1D* h_Y = h->ProjectionY();

    // TCanvas* C_X = new TCanvas("Canvas_X", "Canvas_X", 1000, 500);
    // h_X->GetXaxis()->SetTitle("y[mm]");
    // h_X->GetYaxis()->SetTitle("Counts");
    // h_X->Draw();

    // // Calculate the mean value of the X projection histogram
    // double sum_X = 0;
    // int count_X = 0;
    // for (int bin = 1; bin <= h_X->GetNbinsX(); bin++) {
    //     double binContent = h_X->GetBinContent(bin);
    //     double binCenter = h_X->GetBinCenter(bin);
    //     count_X += binContent;
    // }
    // double mean_X = count_X / h_X->GetNbinsX();

    // TLine* line_X = new TLine(h_X->GetXaxis()->GetXmin(), mean_X, h_X->GetXaxis()->GetXmax(), mean_X);
    // line_X->SetLineColor(kRed);
    // line_X->SetLineWidth(2);
    // line_X->Draw();

    // TCanvas* C_Y = new TCanvas("Canvas_Y", "Canvas_Y", 1000, 500);
    // h_Y->GetXaxis()->SetTitle("z[mm]");
    // h_Y->GetYaxis()->SetTitle("Counts");
    // h_Y->Draw();

    // // Calculate the mean value of the Y projection histogram
    // double sum_Y = 0;
    // int count_Y = 0;
    // for (int bin = 1; bin <= h_Y->GetNbinsX(); bin++) {
    //     double binContent = h_Y->GetBinContent(bin);
    //     double binCenter = h_Y->GetBinCenter(bin);
    //     count_Y += binContent;
    // }
    // double mean_Y = count_Y / h_Y->GetNbinsX();
    // cout << mean_Y << endl;  

    // TLine* line_Y = new TLine(h_Y->GetXaxis()->GetXmin(), mean_Y, h_Y->GetXaxis()->GetXmax(), mean_Y);
    // line_Y->SetLineColor(kRed);
    // line_Y->SetLineWidth(2);
    // line_Y->Draw();

    // C_X->SaveAs("Projection_X.png");
    // C_Y->SaveAs("Projection_Y.png");
}

TH2D* create_histo(TTree* _t, int _OM_num, double _X_plane_shift){
    double parA, parB, parC, parD;
	_t->SetBranchAddress("A", &parA);
	_t->SetBranchAddress("B", &parB);
	_t->SetBranchAddress("C", &parC);
	_t->SetBranchAddress("D", &parD);

    OM_xyz_swcr(_OM_num);

    double X;
    double X_zero_plane = xyz[0];
        
    if (X_zero_plane > 0) X = X_zero_plane - mw_sizex / 2  + _X_plane_shift;
    else X = X_zero_plane + mw_sizex / 2 - _X_plane_shift; 

    double Y_hmin = (xyz[1] - mw_sizey / 2) - 200; 		
    double Y_hmax = (xyz[1] + mw_sizey / 2) + 200; 		
    int    Y_bins = int(Y_hmax - Y_hmin) / 4;
    
    double Z_hmin = (xyz[2] - mw_sizez / 2) - 200;
    double Z_hmax = (xyz[2] + mw_sizez / 2) + 200;
    int    Z_bins = int(Z_hmax - Z_hmin) / 4;

    TString hname = Form("Histo for %.2lf mm OM #%d", X, _OM_num);
    TH2D* h = new TH2D(hname, hname, Y_bins, Y_hmin, Y_hmax, Z_bins, Z_hmin, Z_hmax);
    
    h->GetXaxis()->SetTitle("y[mm]");
    h->GetYaxis()->SetTitle("z[mm]");
    h->SetStats(0);
    h->SetMinimum(0.001);

    for(int entry = 0; entry < _t->GetEntries(); entry++)
    {
        _t->GetEntry(entry);
        
        double Y = parA * X + parB;
        double Z = parC * X + parD;  
        h->Fill(Y, Z);
    }

    return h;
}

int mean_counts(TH1D* h_pr){
    int count = 0;
    for (int bin = 1; bin <= h_pr->GetNbinsX(); bin++) {
        count += h_pr->GetBinContent(bin);
    }
    return count / h_pr->GetNbinsX();
}

double* find_YZ_bounds(TH2D* _h, int _OM_num, int _point, int _DrawOption = 0){
        TH1D* h_pr_y = _h->ProjectionX();
        TH1D* h_pr_z = _h->ProjectionY();

        double* bounds = new double[4];

        int mean_Y = mean_counts(h_pr_y);
        int mean_Z = mean_counts(h_pr_z);

        int run_av_hwidth = 30;

        TGraph* y_graph = running_average(h_pr_y, run_av_hwidth); 
        TGraph* z_graph = running_average(h_pr_z, run_av_hwidth); 


        if (_DrawOption) 
        {
                draw_projections(h_pr_y, h_pr_z, mean_Y, mean_Z, _OM_num, _point);
                draw_projections(y_graph, z_graph, mean_Y, mean_Z, _OM_num, _point);
        }

        double bounds_y[2];
        double bounds_z[2];

        bounds_y[0] = find_intersection(y_graph, mean_Y, TString("min"));
        bounds_y[1] = find_intersection(y_graph, mean_Y, TString("max"));

        bounds_z[0] = find_intersection(z_graph, mean_Z, TString("min"));
        bounds_z[1] = find_intersection(z_graph, mean_Z, TString("max"));


        bounds[0] = bounds_y[0];
        bounds[1] = bounds_y[1];
        bounds[2] = bounds_z[0];
        bounds[3] = bounds_z[1];
        return bounds;
}

double find_intersection(TGraph* _gr, double mean_val, TString _option)
{
        double step, var;
        if (_option == "min") 
        {
                step = 10;
                var = _gr->GetXaxis()->GetXmin();
        }       
        else if (_option == "max")
        {
                step = -10;
                var = _gr->GetXaxis()->GetXmax();
        }
        else return 0;
        
        double diff = _gr->Eval(var, 0, "S") - mean_val * k_mean_spaghetti;
        while (TMath::Abs(diff) > 0.001)
        {       
                var += step;
                double diff_new = _gr->Eval(var, 0, "S") - mean_val * k_mean_spaghetti;
                if ((diff_new > 0 && diff < 0) ||
                    (diff_new < 0 && diff > 0)) step /= -2;
                diff = diff_new;
        }
        return var;
}

TGraph* running_average(TH1D* _h, int _run_av_hwidth){
        int points = _h->GetNbinsX();

        double Arr_X[points], Arr_N[points];
        for (int bin = 0; bin < points; bin++) 
        {
                Arr_X[bin] = _h->GetBinCenter(bin);
                Arr_N[bin] = _h->GetBinContent(bin);
        }

        int points_new = points + 2 * _run_av_hwidth;
        double Arr_N_new[points_new], Arr_X_new[points_new];

        for (int bin = 0; bin < points_new; bin++) 
        {
                double sum = 0;
                for (int i = bin - _run_av_hwidth; i <= bin + _run_av_hwidth; i++) 
                {
                        double a;
                        if ((i < _run_av_hwidth) || (i >= points_new - _run_av_hwidth)) a = 0;
                        else a = Arr_N[i - _run_av_hwidth];
                        sum += a;
                }
                Arr_N_new[bin] = sum / (2 * _run_av_hwidth + 1);

                if (bin < _run_av_hwidth) 
                {
                        Arr_X_new[bin] = Arr_X[0] - (_run_av_hwidth - bin) * (Arr_X[1] - Arr_X[0]);
                }
                else if (bin >= points_new - _run_av_hwidth) 
                {
                        Arr_X_new[bin] = Arr_X[points - 1] + (bin - points_new + _run_av_hwidth + 1) * (Arr_X[1] - Arr_X[0]);
                }
                else 
                {
                        Arr_X_new[bin] = Arr_X[bin - _run_av_hwidth];
                }
        }

        TGraph* graph = new TGraph(points_new, Arr_X_new, Arr_N_new);
        return graph;
}

void draw_projections(TH1D* _h_pr_y, TH1D* _h_pr_z, int _mean_Y, int _mean_Z, int _OM_num, int _point){
        OM_xyz_swcr(_OM_num);

        double X;
        double X_zero_plane = xyz[0];

        if (X_zero_plane > 0) X = X_zero_plane - mw_sizex / 2  + (X_OBSERV_MIN + _point * STEP);
        else X = X_zero_plane + mw_sizex / 2 - (X_OBSERV_MIN + _point * STEP); 
        
        TCanvas* C_pr = new TCanvas("Canvas_projections", "Canvas_projections", 800, 600);

        TLine* line_Y = new TLine(_h_pr_y->GetXaxis()->GetXmin(), _mean_Y * k_mean_spaghetti, _h_pr_y->GetXaxis()->GetXmax(), _mean_Y * k_mean_spaghetti);
        line_Y->SetLineColor(kRed);
        line_Y->SetLineWidth(2);

        TLine* line_Z = new TLine(_h_pr_z->GetXaxis()->GetXmin(), _mean_Z * k_mean_spaghetti, _h_pr_z->GetXaxis()->GetXmax(), _mean_Z * k_mean_spaghetti);
        line_Z->SetLineColor(kRed);
        line_Z->SetLineWidth(2);

 

        _h_pr_y->GetXaxis()->SetTitle("y[mm]");
        _h_pr_y->GetYaxis()->SetTitle("Counts");
        _h_pr_y->SetTitle(Form("Y projection of tracks OM#%d X=%.1lfmm", _OM_num, X)); // Added title
        _h_pr_y->SetStats(0);
        _h_pr_y->Draw();
        line_Y->Draw();

        C_pr->SaveAs(Form("%sClear_Projections_threshold/Projection_Y_OM%03d_point_%03d.png", PATH, _OM_num, _point));


        _h_pr_z->GetXaxis()->SetTitle("z[mm]");
        _h_pr_z->GetYaxis()->SetTitle("Counts");
        _h_pr_z->SetTitle(Form("Z projection of tracks OM#%d X=%.1lfmm", _OM_num, X)); // Added title
        _h_pr_z->SetStats(0);
        _h_pr_z->Draw();
        line_Z->Draw();

        C_pr->SaveAs(Form("%sClear_Projections_threshold/Projection_Z_OM%03d_point_%03d.png", PATH, _OM_num, _point));
}

void draw_projections(TGraph* _gr_y, TGraph* _gr_z, int _mean_Y, int _mean_Z, int _OM_num, int _point){
        OM_xyz_swcr(_OM_num);

        double X;
        double X_zero_plane = xyz[0];

        if (X_zero_plane > 0) X = X_zero_plane - mw_sizex / 2  + (X_OBSERV_MIN + _point * STEP);
        else X = X_zero_plane + mw_sizex / 2 - (X_OBSERV_MIN + _point * STEP); 


        TCanvas* C_Y = new TCanvas("Canvas_projections", "Canvas_projections", 800, 600);

        TLine* line_Y = new TLine(_gr_y->GetXaxis()->GetXmin(), _mean_Y * k_mean_spaghetti, _gr_y->GetXaxis()->GetXmax(), _mean_Y * k_mean_spaghetti);
        line_Y->SetLineColor(kRed);
        line_Y->SetLineWidth(2);

        TLine* line_Z = new TLine(_gr_z->GetXaxis()->GetXmin(), _mean_Z * k_mean_spaghetti, _gr_z->GetXaxis()->GetXmax(), _mean_Z * k_mean_spaghetti);
        line_Z->SetLineColor(kRed);
        line_Z->SetLineWidth(2);

        _gr_y->GetXaxis()->SetTitle("y[mm]");
        _gr_y->GetYaxis()->SetTitle("Counts");
        _gr_y->SetTitle(Form("Y projection of tracks OM#%d X=%.1lfmm", _OM_num, X)); // Added title
        _gr_y->Draw();
        C_Y->SaveAs(Form("%sClear_Projections_RA/Projection_Y_OM%03d_point_%03d.png", PATH, _OM_num, _point));
        line_Y->Draw();

        C_Y->SaveAs(Form("%sClear_Projections_RA_threshold/Projection_Y_OM%03d_point_%03d.png", PATH, _OM_num, _point));

        TCanvas* C_Z = new TCanvas("Canvas_projections", "Canvas_projections", 800, 600);

        _gr_z->GetXaxis()->SetTitle("z[mm]");
        _gr_z->GetYaxis()->SetTitle("Counts");
        _gr_z->SetTitle(Form("Z projection of tracks OM#%d X=%.1lfmm", _OM_num, X)); // Added title
        _gr_z->Draw();
        C_Z->SaveAs(Form("%sClear_Projections_RA/Projection_Z_OM%03d_point_%03d.png", PATH, _OM_num, _point));
        line_Z->Draw();

        C_Z->SaveAs(Form("%sClear_Projections_RA_threshold/Projection_Z_OM%03d_point_%03d.png", PATH, _OM_num, _point));
}

TTree* filter_tree(TTree* _t, int _OM_num, double _bounds0[4], double _bounds1[4]){
    double parA, parB, parC, parD;
	_t->SetBranchAddress("A", &parA);
	_t->SetBranchAddress("B", &parB);
	_t->SetBranchAddress("C", &parC);
	_t->SetBranchAddress("D", &parD);

    TString clear_tree_name = Form("Clear tracks of OM %d", _OM_num);
    TTree* clear_tree = new TTree(clear_tree_name, clear_tree_name);

    double parA_clear, parB_clear, parC_clear, parD_clear;
    clear_tree->Branch("A", &parA_clear);
    clear_tree->Branch("B", &parB_clear);
    clear_tree->Branch("C", &parC_clear);
    clear_tree->Branch("D", &parD_clear);

    double X_plane_0, X_plane_1;
        
    if (xyz[0] > 0) 
    {
        X_plane_0 = xyz[0] - mw_sizex / 2 + X_shifts[0];
        X_plane_1 = xyz[0] - mw_sizex / 2 + X_shifts[1];
    }
    else
    {
        X_plane_0 = xyz[0] + mw_sizex / 2 - X_shifts[0];
        X_plane_1 = xyz[0] + mw_sizex / 2 - X_shifts[1];
    }


    for(int entry = 0; entry < _t->GetEntries(); entry++)
    {
        _t->GetEntry(entry);
        
        double Y0 = parA * X_plane_0 + parB;
        double Z0 = parC * X_plane_0 + parD;
        double Y1 = parA * X_plane_1 + parB;
        double Z1 = parC * X_plane_1 + parD;  

        if (Y0 >= _bounds0[0] && 
            Y0 <= _bounds0[1] && 
            Z0 >= _bounds0[2] && 
            Z0 <= _bounds0[3] &&
            Y1 >= _bounds1[0] &&
            Y1 <= _bounds1[1] &&
            Z1 >= _bounds1[2] &&
            Z1 <= _bounds1[3] ) 
        {
            parA_clear = parA;
            parB_clear = parB;
            parC_clear = parC;
            parD_clear = parD;

            clear_tree->Fill();
        }
    }

    return clear_tree;
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
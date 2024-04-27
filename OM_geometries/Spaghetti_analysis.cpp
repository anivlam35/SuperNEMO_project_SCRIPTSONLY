#include <iostream>
#include "config.h"
using namespace std;

void TGraph_formatting(TGraph* gr, int _OM_num);
void TH1F_formatting(TH1F* h);
void OM_xyz_swcr(int OM_num);
TH2D* create_histo(TTree* _t, int _OM_num, double _X_plane_shift);
int mean_counts(TH1D* h_pr);
double* find_YZ_bounds(TH2D* h, int _OM_num, int _point, int DrawOption = 0);
void draw_projections(TH1D* h_pr_y, TH1D* h_pr_z, int mean_Y, int mean_Z, int _OM_num, int);
void draw_projections(TGraph* _gr_y, TGraph* _gr_z, int _mean_Y, int _mean_Z, int _OM_num, int _point);
double find_intersection(TGraph* _gr, double mean_val, TString _option);
TGraph* running_average(TH1D* _h, int _run_av_hwidth);


void Spaghetti_analysis()
{
        gROOT->SetBatch(kTRUE);
        TFile* f = new TFile(Form("%sOMs_Clear_Tracks_Run-%d.root", PATH, RUN_N));

        TTree* OM_num_tree = (TTree*)f->Get("List of OMs numbers");
        int OM_num;
        OM_num_tree->SetBranchAddress("OM_num", &OM_num);

        int N_OMs_mw = 0;
        
        int start_entry = -1; 

        if (SIDE == 0){
                start_entry = 0;
                for(int entry = 0; entry < OM_num_tree->GetEntries(); entry++)
                {
                        OM_num_tree->GetEntry(entry);
                        if(OM_num < 260) N_OMs_mw++;
                }
        }
        else if (SIDE == 1){
                for(int entry = 0; entry < OM_num_tree->GetEntries(); entry++)
                {
                        OM_num_tree->GetEntry(entry);
                        if(OM_num >= 260) 
                        {
                                N_OMs_mw++;
                                if (start_entry == -1) start_entry = entry;
                        }
                }
        }
        
        OM_num_tree->GetEntry(start_entry);
        OM_xyz_swcr(OM_num);
        double X_BasePlane;
                
        if (xyz[0] > 0) X_BasePlane = xyz[0] - mw_sizex / 2;
        else X_BasePlane = xyz[0] + mw_sizex / 2; 

        int N_Points = int(X_OBSERV_MAX - X_OBSERV_MIN) / STEP + 1;
	double Arr_X[N_Points], Arr_Y[N_Points], Arr_Z[N_Points];
	
	TGraph* grY[N_OMs_mw];
	TGraph* grZ[N_OMs_mw];

        for(int entry = 0; entry < N_OMs_mw; entry++)
	{
                OM_num_tree->GetEntry(start_entry + entry);

                OM_xyz_swcr(OM_num);

		TString treename = Form("Clear tracks of OM %d", OM_num);	
		TTree* tr = (TTree*) f->Get(treename);

                double parA, parB, parC, parD;
		tr->SetBranchAddress("A", &parA);
		tr->SetBranchAddress("B", &parB);
		tr->SetBranchAddress("C", &parC);
		tr->SetBranchAddress("D", &parD);

		for(int point = 0; point < N_Points; point++)
		{
                        double X_shift = X_OBSERV_MIN + point * STEP;

                        TH2D* h_vert = create_histo(tr, OM_num, X_shift);
                        double* bounds = find_YZ_bounds(h_vert, OM_num, point, 0);

                        TCanvas* C0 = new TCanvas("Canvas", "Canvas", 1000, 1000);
                        h_vert->Draw("colz");
                        C0->SaveAs(Form("%sCUTS/OM_%03d_point_%03d.png", PATH, OM_num, point));

                        double X;
                        if (xyz[0] > 0) X = xyz[0] - mw_sizex / 2  + X_shift;
                        else X = xyz[0] + mw_sizex / 2 - X_shift;
                        Arr_X[point] = X;
                        // Arr_Y[point] = bounds[1] - bounds[0];
                        Arr_Z[point] = bounds[3] - bounds[2];

                        delete h_vert;
                        delete bounds;
                        delete C0;
                }
                delete tr;


                // grY[entry] = new TGraph(N_Points, Arr_X, Arr_Y);
                grZ[entry] = new TGraph(N_Points, Arr_X, Arr_Z);
                // grY[entry]->Draw();
	}

        // TCanvas* CY = new TCanvas("#Delta y(x)", "GraphY", 1500, 1000);

        // TPad *pad1y = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
        // pad1y->Draw();
        // pad1y->cd();

        // pad1y->SetBottomMargin(0);

        // OM_num_tree->GetEntry(start_entry);

        // auto legendy = new TLegend(0.9,0.1,1.0,0.9);
        // legendy->SetHeader("#Delta y(x)","C"); // option "C" allows to center the header
        // legendy->AddEntry(grY[0], Form("OM %d", OM_num),"pl");

        // TGraph_formatting(grY[0], OM_num);
        // grY[0]->Draw("APL");
        // grY[0]->GetYaxis()->SetTitle("#Delta y(x) [mm]");
        // grY[0]->GetXaxis()->SetTitle("x [mm]");

        // for(int entry = 1; entry < N_OMs_mw; entry++)
        // {
        //         OM_num_tree->GetEntry(start_entry + entry);

        //         TString nm = Form("OM %d", OM_num);


        //         grY[entry]->SetName(nm);
        //         grY[entry]->SetMarkerColor(1);
        //         grY[entry]->SetMarkerStyle((OM_num % 260) % 4 + 20);
        //         grY[entry]->SetMarkerSize(0.5);
        //         grY[entry]->SetLineColor(+ 5 * int((OM_num % 260)) / 26);
        //         grY[entry]->SetLineWidth(1);
        //         grY[entry]->Draw("PL same");

        //         // legendy->AddEntry(grY[entry], nm, "pl");
        // }

        // // legendy->Draw("same");	

        // CY->cd();
        // TPad *pad2y = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
        // pad2y->Draw();
        // pad2y->cd();

        // pad2y->SetBottomMargin(0.2);
        // pad2y->SetTopMargin(0);

        // TH1F *hY = new TH1F("histo", "histo", N_Points, X_BasePlane + X_OBSERV_MIN - STEP / 2, X_BasePlane + X_OBSERV_MAX + STEP / 2);

        // TH1F_formatting(hY);
        // hY->GetYaxis()->SetTitle("N");
        // hY->GetXaxis()->SetTitle("x [mm]");

        // for(int entry = 0; entry < N_OMs_mw; entry++)
        // {
        //         const double* xValues = grY[entry]->GetX();
        //         const double* yValues = grY[entry]->GetY();

        //         // Find the index of the minimum Y value
        //         double minY = yValues[0];
        //         int minIndex = 0;
        //         for (int i = 1; i < grY[entry]->GetN(); ++i)
        //         {
        //                 if (yValues[i] < minY)
        //                 {
        //                 minY = yValues[i];
        //                 minIndex = i;
        //                 }
        //         }
                

        //         // Get the x-position corresponding to the minimum Y value
        //         double minXPosition = xValues[minIndex];

        //         hY->Fill(minXPosition);
        // }

        // hY->Draw();

        // CY->SetTitle("#Delta y(x)");

        // CY->Print(Form("%svY.png", PATH));

        TCanvas* CZ = new TCanvas("#Delta z(x)", "GraphZ", 1500, 1000);

        TPad *pad1z = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
        pad1z->Draw();
        pad1z->cd();
        pad1z->SetBottomMargin(0);

        OM_num_tree->GetEntry(start_entry);

        auto legendz = new TLegend(0.9,0.1,1.0,0.9);
        legendz->SetHeader("#Delta z(x)","C"); // option "C" allows to center the header
        legendz->AddEntry(grZ[0],Form("OM %d", OM_num),"pl");

        TGraph_formatting(grZ[0], OM_num);
        grZ[0]->Draw("APL");
        grZ[0]->GetYaxis()->SetTitle("#Delta z(x) [mm]");
        grZ[0]->GetXaxis()->SetTitle("x [mm]");
        grZ[0]->SetTitle("OMs x-coordinates (\"french\" side)");

        for(int entry=1; entry < N_OMs_mw; entry++)
        {
                OM_num_tree->GetEntry(start_entry + entry);

                TString nm = Form("OM %d", OM_num);

                grZ[entry]->SetName(nm);
                grZ[entry]->SetMarkerColor(1);
                grZ[entry]->SetMarkerStyle((OM_num % 260) % 4 + 20);
                grZ[entry]->SetMarkerSize(0.5);
                grZ[entry]->SetLineColor(51 + 5 * int(OM_num % 260) / 26);
                grZ[entry]->SetLineWidth(1);
                grZ[entry]->Draw("PL same");

                // legendz->AddEntry(grZ[entry], nm, "pl");
        }

        // legendz->Draw("same");

        CZ->cd();
        TPad *pad2z = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
        pad2z->Draw();
        pad2z->cd();

        pad2z->SetBottomMargin(0.2);
        pad2z->SetTopMargin(0);

        TH1F *hZ = new TH1F("histo", "histo", N_Points, X_BasePlane + X_OBSERV_MIN - STEP / 2, X_BasePlane + X_OBSERV_MAX + STEP / 2);

        TH1F_formatting(hZ);
        hZ->GetYaxis()->SetTitle("N");
        hZ->GetXaxis()->SetTitle("x [mm]");

        for(int entry = 0; entry < N_OMs_mw; entry++)
        {
                // Get the graph's X and Y arrays
                const double* xValues = grZ[entry]->GetX();
                const double* zValues = grZ[entry]->GetY();

                // Find the index of the minimum z value
                double minZ = zValues[0];
                int minIndex = 0;
                for (int i = 1; i < grZ[entry]->GetN(); ++i)
                {
                        if (zValues[i] < minZ)
                        {
                        minZ = zValues[i];
                        minIndex = i;
                        }
                }

                // Get the x-position corresponding to the minimum Y value
                double minXPosition = xValues[minIndex];

                hZ->Fill(minXPosition);
        }

	cout << "Mean_Z: " << hZ->GetMean() << ", SD: " << (N_OMs_mw / (N_OMs_mw - 1)) * hZ->GetRMS() << endl;
        hZ->Draw();

        CZ->SetTitle("#delta z(x)");
        CZ->Print(Form("%svZ.png", PATH));

        for(int entry = 0; entry < N_OMs_mw; entry++)
        {	
                OM_num_tree->GetEntry(start_entry + entry);

                // TString print_name_RMS_Y = Form("%sRMS_ALL_OMs/RMS_Y_OM_%03d.png", PATH, OM_num);
                // TCanvas* CindY = new TCanvas(Form("RMS_Y_OM_%d", OM_num), Form("RMS_Y_OM_%d", OM_num));
                // grY[entry]->SetTitle(Form("CW along y-axis. OM#%d.", OM_num));
                // grY[entry]->GetXaxis()->SetTitle("x [mm]");
                // grY[entry]->GetYaxis()->SetTitle("#Delta y");
                // grY[entry]->Draw();
                // CindY->Print(print_name_RMS_Y);

                TString print_name_RMS_Z = Form("%sRMS_ALL_OMs/RMS_Z_OM_%03d.png", PATH, OM_num);
                TCanvas* CindZ = new TCanvas(Form("RMS_Z_OM_%d", OM_num), Form("RMS_Z_OM_%d", OM_num));
                grZ[entry]->SetTitle(Form("CW along z-axis. OM#%d.", OM_num));
                grZ[entry]->GetXaxis()->SetTitle("x [mm]");
                grZ[entry]->GetYaxis()->SetTitle("#Delta z");
                grZ[entry]->Draw();
                CindZ->Print(print_name_RMS_Z);
        }
}

void TGraph_formatting(TGraph* gr, int _OM_num)
{
        gr->SetMinimum(250);
        gr->SetMaximum(350);
        gr->SetName(Form("OM %d", _OM_num));
        gr->SetMarkerColor(1);
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.5);
        gr->SetLineColor(1);
        gr->SetLineWidth(1);

        gr->GetYaxis()->SetNdivisions(510);
        gr->GetYaxis()->SetLabelFont(43);
        gr->GetYaxis()->SetTitleFont(43);
        gr->GetYaxis()->SetTitleOffset(1.5);
        gr->GetYaxis()->SetLabelSize(20);
        gr->GetYaxis()->SetTitleSize(20);

        gr->GetXaxis()->SetLabelFont(43);
        gr->GetXaxis()->SetTitleFont(43);
        gr->GetXaxis()->SetTitleOffset(3);
        gr->GetXaxis()->SetLabelSize(20);
        gr->GetXaxis()->SetTitleSize(20);

}

void TH1F_formatting(TH1F* h)
{
        h->SetStats(0);
        h->SetTitle("");
        h->GetYaxis()->SetTitleSize(20);
        h->GetYaxis()->SetLabelSize(15);
        h->GetYaxis()->SetTitleFont(43);
        h->GetYaxis()->SetLabelFont(43);
        h->GetYaxis()->SetTitleOffset(1.5);
        h->GetYaxis()->SetNdivisions(10);
        h->GetXaxis()->SetTickLength(0.1);
        h->GetXaxis()->SetTitleSize(20);
        h->GetXaxis()->SetLabelSize(15);
        h->GetXaxis()->SetTitleFont(43);
        h->GetXaxis()->SetLabelFont(43);
        h->GetXaxis()->SetTitleOffset(1);
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
    int    Y_bins = int(Y_hmax - Y_hmin) / 1;
    
    double Z_hmin = (xyz[2] - mw_sizez / 2) - 200;
    double Z_hmax = (xyz[2] + mw_sizez / 2) + 200;
    int    Z_bins = int(Z_hmax - Z_hmin) / 1;

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
        graph->GetXaxis()->SetLimits(_h->GetXaxis()->GetXmin(), _h->GetXaxis()->GetXmax());
        return graph;
}

void draw_projections(TH1D* _h_pr_y, TH1D* _h_pr_z, int _mean_Y, int _mean_Z, int _OM_num, int _point){
        OM_xyz_swcr(_OM_num);

        double X;
        double X_zero_plane = xyz[0];

        if (X_zero_plane > 0) X = X_zero_plane - mw_sizex / 2  + (X_OBSERV_MIN + _point * STEP);
        else X = X_zero_plane + mw_sizex / 2 - (X_OBSERV_MIN + _point * STEP); 
        
        TCanvas* C_pr = new TCanvas("Canvas_projections", "Canvas_projections", 800, 600);

        // TLine* line_Y = new TLine(_h_pr_y->GetXaxis()->GetXmin(), _mean_Y * k_mean_spaghetti, _h_pr_y->GetXaxis()->GetXmax(), _mean_Y * k_mean_spaghetti);
        // line_Y->SetLineColor(kRed);
        // line_Y->SetLineWidth(2);

        TLine* line_Z = new TLine(_h_pr_z->GetXaxis()->GetXmin(), _mean_Z * k_mean_spaghetti, _h_pr_z->GetXaxis()->GetXmax(), _mean_Z * k_mean_spaghetti);
        line_Z->SetLineColor(kRed);
        line_Z->SetLineWidth(2);

 

        // _h_pr_y->GetXaxis()->SetTitle("y[mm]");
        // _h_pr_y->GetYaxis()->SetTitle("Counts");
        // _h_pr_y->SetTitle(Form("Y projection of tracks OM#%d X=%.1lfmm", _OM_num, X)); // Added title
        // _h_pr_y->SetStats(0);
        // _h_pr_y->Draw();
        // C_pr->SaveAs(Form("%sClear_Projections/Projection_Y_OM%03d_point_%03d.png", PATH, _OM_num, _point));
        // line_Y->Draw();

        // C_pr->SaveAs(Form("%sClear_Projections_threshold/Projection_Y_OM%03d_point_%03d.png", PATH, _OM_num, _point));


        _h_pr_z->GetXaxis()->SetTitle("z[mm]");
        _h_pr_z->GetYaxis()->SetTitle("Counts");
        _h_pr_z->SetTitle(Form("Z projection of tracks OM#%d X=%.1lfmm", _OM_num, X)); // Added title
        _h_pr_z->SetStats(0);
        _h_pr_z->Draw();
        C_pr->SaveAs(Form("%sClear_Projections/Projection_Z_OM%03d_point_%03d.png", PATH, _OM_num, _point));
        line_Z->Draw();

        C_pr->SaveAs(Form("%sClear_Projections_threshold/Projection_Z_OM%03d_point_%03d.png", PATH, _OM_num, _point));

        delete C_pr;
        delete line_Z;
}

void draw_projections(TGraph* _gr_y, TGraph* _gr_z, int _mean_Y, int _mean_Z, int _OM_num, int _point){
        OM_xyz_swcr(_OM_num);

        double X;
        double X_zero_plane = xyz[0];

        if (X_zero_plane > 0) X = X_zero_plane - mw_sizex / 2  + (X_OBSERV_MIN + _point * STEP);
        else X = X_zero_plane + mw_sizex / 2 - (X_OBSERV_MIN + _point * STEP); 


        // TCanvas* C_Y = new TCanvas("Canvas_projections", "Canvas_projections", 800, 600);

        // TLine* line_Y = new TLine(_gr_y->GetXaxis()->GetXmin(), _mean_Y * k_mean_spaghetti, _gr_y->GetXaxis()->GetXmax(), _mean_Y * k_mean_spaghetti);
        // line_Y->SetLineColor(kRed);
        // line_Y->SetLineWidth(2);

        TLine* line_Z = new TLine(_gr_z->GetXaxis()->GetXmin(), _mean_Z * k_mean_spaghetti, _gr_z->GetXaxis()->GetXmax(), _mean_Z * k_mean_spaghetti);
        line_Z->SetLineColor(kRed);
        line_Z->SetLineWidth(2);

        // _gr_y->GetXaxis()->SetTitle("y[mm]");
        // _gr_y->GetYaxis()->SetTitle("Counts");
        // _gr_y->SetTitle(Form("Y projection of tracks OM#%d X=%.1lfmm", _OM_num, X)); // Added title
        // _gr_y->Draw();
        // C_Y->SaveAs(Form("%sClear_Projections_RA/Projection_Y_OM%03d_point_%03d.png", PATH, _OM_num, _point));
        // line_Y->Draw();

        // C_Y->SaveAs(Form("%sClear_Projections_RA_threshold/Projection_Y_OM%03d_point_%03d.png", PATH, _OM_num, _point));

        TCanvas* C_Z = new TCanvas("Canvas_projections", "Canvas_projections", 800, 600);

        _gr_z->GetXaxis()->SetTitle("z[mm]");
        _gr_z->GetYaxis()->SetTitle("Counts");
        _gr_z->SetTitle(Form("Z projection of tracks OM#%d X=%.1lfmm", _OM_num, X)); // Added title
        _gr_z->Draw();
        C_Z->SaveAs(Form("%sClear_Projections_RA/Projection_Z_OM%03d_point_%03d.png", PATH, _OM_num, _point));
        line_Z->Draw();

        C_Z->SaveAs(Form("%sClear_Projections_RA_threshold/Projection_Z_OM%03d_point_%03d.png", PATH, _OM_num, _point));

        delete line_Z;
        delete C_Z;
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
#include "/sps/nemo/scratch/ikovalen/TKEvent_old/TKEvent/include/TKEvent.h" 
#include "config.h"

R__LOAD_LIBRARY(/sps/nemo/scratch/ikovalen/TKEvent_old/TKEvent/lib/libTKEvent.so);

vector<double> get_guess(TH1D* _hpy, TH1D* _hpz);
vector<double> find_yz(TTree* _t,                       double _ymin, double _ymax, double _zmin, double _zmax, int _iX, int _iter);
vector<double> find_ab(TTree* _t, double _y, double _z, double _amin, double _amax, double _bmin, double _bmax, int _iX, int _sct);
vector<double> best_par(TTree* _t, int _iX, int _SRC_NO_i, int Draw_opt = 1);
TTree* CleanTree(TTree* _t, vector<vector<double>> _Best_yzab, char _Newtr_name[50]);

using namespace std;
// Set Ranges
const int    X_MIN = -20/*-10*/;
const int    X_MAX = -20/*10*/;

int X_Cl[2] = {-10, 10};

// const int SRC_NO = 5;
const double   R = 30.0;

//const    int k = 4;
//const int    B_BINS  = 6;

void Spaghetti_Analysis()
{
	TFile* file = new TFile(Form("%sTracks_Run-%d.root", PATH, RUN_N));
	// TFile* file = new TFile("Clear_data.root");
	
	TFile *New_file = new TFile(Form("%sClear_data_Run-%d.root", PATH, RUN_N),"RECREATE");
	// Load the data for one source (x<0)

	for(int SRC_NO_i=0; SRC_NO_i < N_SRCPLN; SRC_NO_i++)
	{
		//stringstream tr_name;

		//tr_name << "Tracks for " << X_BasePlane << " mm %s Source " << SRC_NO_i;
		TString tr_name = Form("Tracks for %d mm %s Source %d", X_BasePlane, SIDE, SRC_NO_i);
		TTree* tr = (TTree*) file->Get(tr_name);

		char Newtr_name[50];
		sprintf(Newtr_name,"Tracks for %i mm %s Source %i", X_BasePlane, SIDE, SRC_NO_i);

		vector<vector<double>> Clean_yzab(2);
		Clean_yzab[0] = best_par(tr, X_Cl[0], SRC_NO_i, 1);
		Clean_yzab[1] = best_par(tr, X_Cl[1], SRC_NO_i, 1);

		TTree* Newtr = CleanTree(tr, Clean_yzab, Newtr_name);

		Newtr->Write(Newtr_name);

		cout<<endl;
		cout<<"Done for "<<SRC_NO_i<<" Source"<<endl;
		cout<<endl;
		delete tr;
	}
}

vector<double> get_guess(TH1D* _hpy, TH1D* _hpz)
{
	vector<double> g(4);

	g[0] = _hpy->GetBinCenter(_hpy->GetMaximumBin()); 
	g[1] = _hpz->GetBinCenter(_hpz->GetMaximumBin()); 

	Double_t quant[2] = {0.2, 0.8};  // position where to compute the quantiles in [0,1]
	Double_t res[2]; 		// array to contain the quantiles
	_hpy->GetQuantiles(2, res, quant);
	g[2] = (res[1] - res[0]) / 2.0/*_hpy->GetStdDev()*/; 
	_hpz->GetQuantiles(2, res, quant);
	g[3] = (res[1] - res[0]) / 2.0/*_hpz->GetStdDev()*/; 

	//cout << "GUESS: Y: " << g[0] << " +- " << g[2] << ", Z: " << g[1] << " +- " << g[3] << endl;
	
	return g;
}

vector<double> find_yz(TTree* _t, double _ymin, double _ymax, double _zmin, double _zmax, int _iX, int _iter)
{
	const double R    = 30.0;
	const int    NSEC = 4;

	vector<double> yz(2);

	if (_iter <= 5) 
	{	
		//cout << "Iteration " << _iter << ", R = " << R << ", searching in the yz region:" << endl;
		//cout << _ymin << " < y < " << _ymax << endl;
		//cout << _zmin << " < z < " << _zmax << endl;	
		
		double parA, parB, parC, parD;
		_t->SetBranchAddress("A", &parA);
		_t->SetBranchAddress("B", &parB);
		_t->SetBranchAddress("C", &parC);
		_t->SetBranchAddress("D", &parD);

		double verY, verZ;

		const int Y_BINS = 5;
		const int Z_BINS = 5;
		double    ystp   = (_ymax - _ymin) / double(Y_BINS);
		double    zstp   = (_zmax - _zmin) / double(Z_BINS);
		
		// definition and initization of variables to keep sector counts and Xi2 histograms
		int N[Y_BINS][Z_BINS];
		int S[Y_BINS][Z_BINS][NSEC];
		//TH2D* HistXi2[B_BINS];

		for(int iY = 0; iY < Y_BINS; iY++)
		{
			for(int iZ = 0; iZ < Z_BINS; iZ++)
			{
				N[iY][iZ] = 0;

				for(int iS = 0; iS < NSEC; iS++) 
				{				
					S[iY][iZ][iS] = 0;
				}
			}
		}

		// Use Tree->GetEntries() for all entries in Run
		for(UInt_t i = 0; i < _t->GetEntries(); i++)	// Loop over events
		{
			_t->GetEntry(i);
			verY = parA*_iX + parB;
			verZ = parC*_iX + parD;

			for(int iY = 0; iY < Y_BINS; iY++)
			{
				double cY = _ymin + (iY + 0.5) * ystp;

				for(int iZ = 0; iZ < Z_BINS; iZ++)
				{
					double cZ = _zmin + (iZ + 0.5) * zstp;

					double phi = atan2(     verZ - cZ,             verY - cY       );
					double r   = sqrt( pow( verZ - cZ, 2.0) + pow( verY - cY, 2.0) );

					if(r <= R)
					{
						for(int iS = 0; iS < NSEC; iS++)
						{
							if(phi >  -TMath::Pi() + (2.0*TMath::Pi()* iS     ) / NSEC && 
							   phi <= -TMath::Pi() + (2.0*TMath::Pi()*(iS + 1)) / NSEC ) 
							{
								S[iY][iZ][iS]++;
							}
						}

						N[iY][iZ]++;
					}
				}
			}

			//if (i % 10000 == 0) cout <<"Event No. " << i << " done!" <<endl;
		}	
		
		// Calculate Xi2 maps and save the best values
		double Xi2_bst = 1e10;
		double   Y_bst;
		double   Z_bst;
				
		for(int iY = 0; iY < Y_BINS; iY++)
		{
			double cY = _ymin + (iY + 0.5) * ystp;

			for(int iZ = 0; iZ < Z_BINS; iZ++)
			{
				double cZ = _zmin + (iZ + 0.5) * zstp;

				double Xi2 = 0.0;
				double NoK = double(N[iY][iZ]) / double(NSEC);

				for(int iS = 0; iS < NSEC; iS++) 
				{				
					Xi2 += pow( S[iY][iZ][iS] - NoK, 2.0) / NoK;
				}

				if(Xi2 < Xi2_bst)
				{
					  Y_bst = cY;
					  Z_bst = cZ;
					Xi2_bst  = Xi2;
					
					//cout << "Y_bst = " << Y_bst << "  Z_bst = " << Z_bst << " Xi2_bst / NDF = " << Xi2_bst / (NSEC - 1.0) <<endl;
				}
				
				//HistXi2[iB]->Fill(cY, cZ, Xi2);
			}
		}

		yz = find_yz(_t, Y_bst - 0.5*ystp, Y_bst + 0.5*ystp, Z_bst - 0.5*zstp, Z_bst + 0.5*zstp, _iX, _iter + 1);
	}
	else
	{
		yz[0] = (_ymax + _ymin) / 2.0;
		yz[1] = (_zmax + _zmin) / 2.0;

		return yz;
	}

	return yz;
}

vector<double> find_ab(TTree* _t, double _y, double _z, double _amin, double _amax, double _bmin, double _bmax, int _iX, int _sct)
{
	vector<double> ab(2);

	cout << "Center coordinates are (y, z) = (" << _y << ", " << _z << "), searching in the ab region:" << endl;
	cout << _amin << " < a < " << _amax << endl;
	cout << _bmin << " < b < " << _bmax << endl;	
	
	double parA, parB, parC, parD;
	_t->SetBranchAddress("A", &parA);
	_t->SetBranchAddress("B", &parB);
	_t->SetBranchAddress("C", &parC);
	_t->SetBranchAddress("D", &parD);

	double verY, verZ;

	const int A_BINS = 30;
	const int B_BINS = 30;
	double    astp   = (_amax - _amin) / double(A_BINS);
	double    bstp   = (_bmax - _bmin) / double(B_BINS);
	
	stringstream ss;
	ss << "Xi2 / NDF for (y, z) = (" << _y << ", " << _z << ")";
	
	TH2D* hXi2_ndf = new TH2D(ss.str().c_str(), ss.str().c_str(), A_BINS, _amin, _amax, B_BINS, _bmin, _bmax);

	// definition and initization of variables to keep sector counts and Xi2 histograms
	int N[A_BINS][B_BINS];
	int S[A_BINS][B_BINS][_sct];

	for(int iA = 0; iA < A_BINS; iA++)
	{
		for(int iB = 0; iB < B_BINS; iB++)
		{			
			N [iA][iB] = 0;

			for(int iS = 0; iS < _sct; iS++) 
			{				
				S[iA][iB][iS] = 0;
			}
		}
	}

	// Use Tree->GetEntries() for all entries in Run
	for(UInt_t i = 0; i < _t->GetEntries(); i++)	// Loop over events
	{
		_t->GetEntry(i);
		verY = parA*_iX + parB;
		verZ = parC*_iX + parD;

		for(int iA = 0; iA < A_BINS; iA++)
		{
			double cA  = _amin + (iA + 0.5) * astp;

			for(int iB = 0; iB < B_BINS; iB++)
			{
				double cB  = _bmin + (iB + 0.5) * bstp;

				double phi = atan2(     (verZ - _z)/cB/*(cA*cB)*/,             (verY - _y)/cA);
				double r   = sqrt( pow( (verZ - _z)/cB/*(cA*cB)*/, 2.0) + pow( (verY - _y)/cA, 2.0) );

				if(r <= 1.0)
				{
					for(int iS = 0; iS < _sct; iS++)
					{
						if(phi >  -TMath::Pi() + (2.0*TMath::Pi()* iS     ) / _sct && 
						   phi <= -TMath::Pi() + (2.0*TMath::Pi()*(iS + 1)) / _sct ) 
						{
							S[iA][iB][iS]++;
						}
					}

					N[iA][iB]++;
				}
			}
		}

		//if (i % 10000 == 0) cout <<"Event No. " << i << " done!" <<endl;
	}	
	
	// Calculate Xi2 maps and save the best values
	double Xi2_bst = 1e10;
	double   a_bst;
	double   b_bst;

	double pa[A_BINS];
	double pb[A_BINS];

	for(int iA = 0; iA < A_BINS; iA++)
	{
		double Xi2_bstA = 1e10;
		double cA = _amin + (iA + 0.5) * astp;
		pa[iA] = cA;		

		for(int iB = 0; iB < B_BINS; iB++)
		{
			double cB = _bmin + (iB + 0.5) * bstp;	

			double Xi2 = 0.0;
			double NoK = double(N[iA][iB]) / double(_sct);

			for(int iS = 0; iS < _sct; iS++) 
			{				
				Xi2 += pow( S[iA][iB][iS] - NoK, 2.0) / NoK;
			}

			if(Xi2 < Xi2_bstA)
			{
				  pb[iA] = cB;
				Xi2_bstA = Xi2;
			}

			if(Xi2 < Xi2_bst)
			{
				  a_bst = cA;
				  b_bst = cB;
				Xi2_bst = Xi2;
			}
			
			hXi2_ndf->Fill(cA, cB, Xi2 / (_sct - 1.0));
		}
		
		cout  << "pa = " << pa[iA] << " pb = " << pb[iA] << ", ratio: " << pb[iA] / pa[iA] << " Xi2_bst / NDF = " << Xi2_bstA / (_sct - 1.0) <<endl;
	}

	cout  << "a_bst = " << a_bst << " b_bst = " << b_bst << " Xi2_bst / NDF = " << Xi2_bst / (_sct - 1.0) <<endl;

	// TCanvas *cc = new TCanvas();
	// hXi2_ndf->GetXaxis()->SetTitle("a[mm]");
	// hXi2_ndf->GetYaxis()->SetTitle("b[mm]");
	// hXi2_ndf->Draw("COLZ");	
 
	// TMarker* bd[A_BINS];
	// for(int iA = 0; iA < A_BINS; iA++)
	// {
	// 	bd[iA] = new TMarker(pa[iA], pb[iA], 2);
	// 	bd[iA]->SetMarkerColor(2);
	// 	bd[iA]->Draw("same");
	// }

	ab[0] = a_bst;
	ab[1] = b_bst;

	delete hXi2_ndf;
	return ab;
}

vector<double> best_par(TTree* _t, int _iX, int _SRC_NO_i, int Draw_opt = 1)
{
	vector<double> Best_yzab(4);
	//_t->Print();
	double parA, parB, parC, parD;
	_t->SetBranchAddress("A", &parA);
	_t->SetBranchAddress("B", &parB);
	_t->SetBranchAddress("C", &parC);
	_t->SetBranchAddress("D", &parD);

	double verY = 0, verZ = 0;
	
	// Fill histogram of undistorted vertices, and calculate coordinates of mode and the sigmas in both directions
	TString h_name = Form("Hist for %d mm %s Source %d", _iX, SIDE, _SRC_NO_i);

	//h_name << "Hist for " << _iX << " mm x < 0 Source " << _SRC_NO_i;

	TH2D* h_vert_real = new TH2D(h_name, h_name, int(_t->GetMaximum("B") - _t->GetMinimum("B")), _t->GetMinimum("B"), _t->GetMaximum("B"), 
										int(_t->GetMaximum("D") - _t->GetMinimum("D")), _t->GetMinimum("D"), _t->GetMaximum("D")); 						   

	for(int i = 0; i < _t->GetEntries(); i++)
	{
		_t->GetEntry(i);
		//cout << " PARAMETERS: " << parA << "  " << parB << endl;
		verY = parA*_iX + parB;
		verZ = parC*_iX + parD;

		h_vert_real->Fill(verY, verZ);
	}

	vector<double> guess = get_guess(h_vert_real->ProjectionX(), h_vert_real->ProjectionY());
	vector<double> yz = find_yz(_t, guess[0] - guess[2], guess[0] + guess[2], guess[1] - guess[3], guess[1] + guess[3], _iX, 1);
	vector<double> ab; //= find_ab(_t, yz[0], yz[1],guess[2]/2.0/*0.0*/, 2.0*guess[2], guess[3]/2.0/*0.0*/, 2.0*2.0*guess[3], _iX, 50);
	ab.push_back(30.0);
	ab.push_back(30.0);
	//vector<double> ab = find_ab(tr, yz[0], yz[1], 0.0, 100.0, 0.1, 5.0, 50);

	cout << "Best values: " << yz[0] << ", " << yz[1] << ", " << ab[0] << ", " << ab[1] << endl;

	Best_yzab[0] = yz[0];
	Best_yzab[1] = yz[1];
	Best_yzab[2] = ab[0];
	Best_yzab[3] = ab[1];
	
	if(Draw_opt != 0)
	{
		TCanvas* C0 = new TCanvas(h_name, h_name, 800, 800);
		h_vert_real->GetXaxis()->SetRangeUser(yz[0] - 3.0*guess[3], yz[0] + 3.0*guess[3]);
		h_vert_real->GetYaxis()->SetRangeUser(yz[1] - 3.0*guess[3], yz[1] + 3.0*guess[3]);
		h_vert_real->GetXaxis()->SetTitle("y[mm]");
		h_vert_real->GetYaxis()->SetTitle("z[mm]");
		h_vert_real->SetStats(0);
		h_vert_real->Draw("COLZ"); 

		TEllipse* el = new TEllipse(yz[0], yz[1], ab[0], ab[1]);
		TMarker* bod = new TMarker(yz[0], yz[1], 20);
		bod->SetMarkerColor(2);
		bod->Draw("same");
		el->SetFillStyle(0);
		el->SetLineWidth(4);
		el->SetLineColor(2);
		el->Draw("same");

		C0->SetLogz();
		C0->SaveAs(Form("%sCENTERS/Source_%02d_center.png", PATH, _SRC_NO_i));
		delete C0;
		delete el;
		delete bod;
	}

	return Best_yzab;
}

TTree *CleanTree(TTree *_t, vector<vector<double>> _Best_yzab, char _Newtr_name[50])
{	
	TTree* Newtr = new TTree(_Newtr_name, _Newtr_name);

	double A_Cl, B_Cl, C_Cl, D_Cl;
	TBranch* branchA0 = Newtr->Branch("A", &A_Cl, "A/D");
	TBranch* branchB0 = Newtr->Branch("B", &B_Cl, "B/D");
	TBranch* branchC0 = Newtr->Branch("C", &C_Cl, "C/D");
	TBranch* branchD0 = Newtr->Branch("D", &D_Cl, "D/D");

	double parA_New, parB_New, parC_New, parD_New;
	_t->SetBranchAddress("A", &parA_New);
	_t->SetBranchAddress("B", &parB_New);
	_t->SetBranchAddress("C", &parC_New);
	_t->SetBranchAddress("D", &parD_New);

	double Y_cl[2], Z_cl[2];
	double r[2], phi[2], R[2];

	for(int i = 0; i < _t->GetEntries(); i++)
	{
		_t->GetEntry(i);

		for(int j=0; j<2; j++)
		{
			Y_cl[j] = parA_New*X_Cl[j] + parB_New;
			Z_cl[j] = parC_New*X_Cl[j] + parD_New;

			r[j]   = sqrt( pow( Y_cl[j] - _Best_yzab[j][0], 2.0) + pow( Z_cl[j] - _Best_yzab[j][1], 2.0) );
			// if(_Best_yzab[j][0]>_Best_yzab[j][1]) 
			// {
				phi[j] = atan2(Z_cl[j] - _Best_yzab[j][1], Y_cl[j] - _Best_yzab[j][0]);
				R[j]   = _Best_yzab[j][2]*_Best_yzab[j][3]/(sqrt(pow( _Best_yzab[j][2]*sin(phi[j]), 2.0) + pow( _Best_yzab[j][3]*cos(phi[j]), 2.0)));
			// }
			// else 								  
			// {	
			// 	phi[j] = atan2(Y_cl[j] - _Best_yzab[j][0], -Z_cl[j] + _Best_yzab[j][1]);
			// 	R[j]   = _Best_yzab[j][2]*_Best_yzab[j][3]/(sqrt(pow( _Best_yzab[j][3]*sin(phi[j]), 2.0) + pow( _Best_yzab[j][2]*cos(phi[j]), 2.0)));
			// }
		}

		if(r[0]<=R[0] && r[1]<=R[1])
		{
			A_Cl = parA_New;
			B_Cl = parB_New;
			C_Cl = parC_New;
			D_Cl = parD_New;

			Newtr->Fill();
		}
	}

    return Newtr;
}


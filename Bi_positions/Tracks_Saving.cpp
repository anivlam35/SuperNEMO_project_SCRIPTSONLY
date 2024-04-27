#include "/sps/nemo/scratch/ikovalen/TKEvent/TKEvent/include/TKEvent.h" 
#include "config.h"

R__LOAD_LIBRARY(/sps/nemo/scratch/ikovalen/TKEvent/TKEvent/lib/libTKEvent.so);

using namespace std;
using namespace TMath;

const    int ABS_X_MAX = 10;
const    int ABS_X_MIN = 10;

int event_selection(TKEvent* _event, double _low_qlt, double _high_qlt, int _ambg_mode);
int event_selection(TKEvent* _event, double _low_qlt_r, double _high_qlt_r, double _low_qlt_z, double _high_qlt_z, int _ambg_mode); 
int choose_ambg_track(TKEvent* _event);

void Tracks_Saving()
{
	TFile* f = new TFile(Form("/sps/nemo/scratch/ikovalen/TKEvent/runs/Run-%d.root", RUN_N));
	TTree* s = (TTree*) f->Get("Event");

	// MIRO: ADD RUN NUMBER TO FILENAME
	TFile *New_file = new TFile(Form("%sTracks_Run-%d.root", PATH, RUN_N), "RECREATE"); // new root file for received histograms
	

	// Define Tree and Arrays
	TTree*   Tree[2][N_SRCPLN];
	double A_Tree[2][N_SRCPLN];
	double B_Tree[2][N_SRCPLN];
	double C_Tree[2][N_SRCPLN];
	double D_Tree[2][N_SRCPLN];

	
	for(int NSOR = 0; NSOR<N_SRCPLN; NSOR++)
	{
		stringstream ss0, ss1;
		ss0 << "Tracks for " << 0 << " mm x < 0 Source " << NSOR; 
		ss1 << "Tracks for " << 0 << " mm x > 0 Source " << NSOR; 
		
		// Trees (0 for x<0, and 1 for x>0)
		Tree[0][NSOR] = new TTree(ss0.str().c_str(), ss0.str().c_str());

		TBranch* branchA0 = Tree[0][NSOR]->Branch("A", &A_Tree[0][NSOR], "A/D");
		TBranch* branchB0 = Tree[0][NSOR]->Branch("B", &B_Tree[0][NSOR], "B/D");
		TBranch* branchC0 = Tree[0][NSOR]->Branch("C", &C_Tree[0][NSOR], "C/D");
		TBranch* branchD0 = Tree[0][NSOR]->Branch("D", &D_Tree[0][NSOR], "D/D");

		Tree[1][NSOR] = new TTree(ss1.str().c_str(), ss1.str().c_str());

		TBranch* branchA1 = Tree[1][NSOR]->Branch("A", &A_Tree[1][NSOR], "A/D");
		TBranch* branchB1 = Tree[1][NSOR]->Branch("B", &B_Tree[1][NSOR], "B/D");
		TBranch* branchC1 = Tree[1][NSOR]->Branch("C", &C_Tree[1][NSOR], "C/D");
		TBranch* branchD1 = Tree[1][NSOR]->Branch("D", &D_Tree[1][NSOR], "D/D");
	} 

	// Filling Tree
	// Use s->GetEntries() for all entries in Run
	for(UInt_t i=0; i < s->GetEntries(); i++)	// Loop over events
	{
		TKEvent* Eve = new TKEvent(-1,-1);
		s->SetBranchAddress("Eventdata", &Eve);

		s->GetEntry(i);
		Eve->set_r("Manchester", "distance");
		Eve->set_h(2770.0);
		Eve->reconstruct_ML(0);		

		int event_choice = event_selection(Eve, LOW_Q, UP_Q, AMBG_MODE);
		if (event_choice)
		{	
			int track_num = event_choice - 1;
			double Y = Eve->get_track(track_num)->get_b();
			double Z = Eve->get_track(track_num)->get_d();
		
			if(Z!=0 && Z > Z_MIN && Z < Z_MAX && Y > Y_MIN && Y < Y_MAX)
			{
				int NSOR = N_COLS * int((Z_MAX - Z) / Z_RECT_SIZE) + int((Y - Y_MIN) / Y_RECT_SIZE);	
				int side = Eve->get_track(track_num)->get_side();
				
				A_Tree[side][NSOR] = Eve->get_track(track_num)->get_a();
				B_Tree[side][NSOR] = Eve->get_track(track_num)->get_b();
				C_Tree[side][NSOR] = Eve->get_track(track_num)->get_c();
				D_Tree[side][NSOR] = Eve->get_track(track_num)->get_d();
				
				Tree[side][NSOR]->Fill();
			}
		}
		if (i % 10000 == 0) cout <<"Event No. " << i << " done!" <<endl;
		Eve->~TKEvent();
	}	


	// Writing Histograms
	for(int NSOR = 0; NSOR<N_SRCPLN; NSOR++)
	{
		stringstream ss0, ss1;
		ss0<<"Tracks for " << 0 << " mm x < 0 Source " << NSOR; 
		ss1<<"Tracks for " << 0 << " mm x > 0 Source " << NSOR;  
		
		// write histograms into new root file
		Tree[0][NSOR]->Write(ss0.str().c_str());
		Tree[1][NSOR]->Write(ss1.str().c_str());
	} 

	New_file->Close();
	cout << "DONE!" << endl;
}

int event_selection(TKEvent* _event, double _low_qlt, double _high_qlt, int _ambg_mode)
{
	int choice;
	if(_event->get_no_tracks() == 0)
	{
		choice = 0;
	}
	else if(_event->get_no_tracks()             == 1       &&
	        _event->get_track(0)->get_quality() >  _low_qlt && 
		_event->get_track(0)->get_quality() <  _high_qlt)
	{
		choice = 1;
	}
	else if(_event->get_no_tracks()                    == 2       && 
		_event->get_track(0)->get_ambiguity_type() != 0       && 
		_event->get_track(0)->get_quality()        >  _low_qlt && 
		_event->get_track(0)->get_quality()        <  _high_qlt)
	{
		switch(_ambg_mode)
		{
		case 0:
			choice = 0;
			break;
		case 1:
			choice = choose_ambg_track(_event);
			break;
		}
	}
	else
	{
		choice = 0;
	}	
	return choice;
}

int event_selection(TKEvent* _event, double _low_qlt_r, double _high_qlt_r, double _low_qlt_z, double _high_qlt_z, int _ambg_mode)
{
        int choice;
        if(_event->get_no_tracks() == 0)
        {
                choice = 0;
        }
        else if(_event->get_no_tracks()             == 1             &&
                _event->get_track(0)->get_quality_R() >  _low_qlt_r  &&
                _event->get_track(0)->get_quality_R() <  _high_qlt_r &&
                _event->get_track(0)->get_quality_Z() >  _low_qlt_z  &&
                _event->get_track(0)->get_quality_Z() <  _high_qlt_z)
	{
                choice = 1;
        }
        else if(_event->get_no_tracks()                    == 2      &&
                _event->get_track(0)->get_ambiguity_type() != 0      &&
                _event->get_track(0)->get_quality_R() >  _low_qlt_r  &&
                _event->get_track(0)->get_quality_R() <  _high_qlt_r &&
                _event->get_track(0)->get_quality_Z() >  _low_qlt_z  &&
                _event->get_track(0)->get_quality_Z() <  _high_qlt_z)
        {
                switch(_ambg_mode)
                {
                case 0:
                       	choice = 0;
                        break;
                case 1:
                       	choice = choose_ambg_track(_event);
                        break;
                }
        }
	else
	{
                choice = 0;
        }
	return choice;
}


int choose_ambg_track(TKEvent* _event)
{
	int track;
	double Y1 = _event->get_track(0)->get_b();
        double Z1 = _event->get_track(0)->get_d();

        double Y2 = _event->get_track(1)->get_b();
        double Z2 = _event->get_track(1)->get_d();

        if (Z1 != 0     &&
	    Z2 != 0     && 
	    Z1 >  Z_MIN &&
	    Z1 <  Z_MAX && 
	    Y1 >  Y_MIN && 
	    Y1 <  Y_MAX &&
	    Z2 >  Z_MIN && 
	    Z2 <  Z_MAX && 
	    Y2 >  Y_MIN && 
	    Y2 <  Y_MAX)
        {
        	int NSOR1 = N_COLS * ((int)(Z_MAX - Z1) / Z_RECT_SIZE) + (int)(Y1 - Y_MIN) / Y_RECT_SIZE;
                int NSOR2 = N_COLS * ((int)(Z_MAX - Z2) / Z_RECT_SIZE) + (int)(Y2 - Y_MIN) / Y_RECT_SIZE;

                double Y_SOR1 = Y_MIN + (NSOR1 % N_COLS + 0.5) * Y_RECT_SIZE;
                double Z_SOR1 = Z_MAX - (NSOR1 / N_COLS + 0.5) * Z_RECT_SIZE;

                double Y_SOR2 = Y_MIN + (NSOR2 % N_COLS + 0.5) * Y_RECT_SIZE;
                double Z_SOR2 = Z_MAX - (NSOR2 / N_COLS + 0.5) * Z_RECT_SIZE;

                double DIST_1 = Sqrt(Power(Y1 - Y_SOR1, 2) + Power(Z1 - Z_SOR1, 2));
                double DIST_2 = Sqrt(Power(Y2 - Y_SOR2, 2) + Power(Z2 - Z_SOR2, 2));
		

		if ( DIST_1 < DIST_2 ) { track = 1; }
		else { track = 2; }
	}
	else
	{
		track = 0;
	}
	
	return track;
}

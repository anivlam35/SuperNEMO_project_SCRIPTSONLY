#include "/sps/nemo/scratch/ikovalen/TKEvent/TKEvent/include/TKEvent.h" 
#include "config.h"

R__LOAD_LIBRARY(/sps/nemo/scratch/ikovalen/TKEvent/TKEvent/lib/libTKEvent.so);

using namespace std;
using namespace TMath;

int event_selection(TKEvent* _event, double _low_qlt, double _high_qlt, int _ambg_mode); 
int choose_ambg_track(TKEvent* _event);

void Tracks_Saving()
{
	TFile* f = new TFile(Form("/sps/nemo/scratch/ikovalen/TKEvent/runs/Run-%d.root", RUN_N));
	TTree* s = (TTree*) f->Get("Event");
//	s->Print();

	TKEvent* Eve = new TKEvent(-1,-1);
	s->SetBranchAddress("Eventdata", &Eve);

	TFile *New_file = new TFile(Form("%sOMs_Tracks_Run-%d.root", PATH, RUN_N),"RECREATE"); // new root file for received histograms
	
	// Define Tree and Arrays
	TTree*   Tree[N_OMs];
	double A_Tree[N_OMs];
	double B_Tree[N_OMs];
	double C_Tree[N_OMs];
	double D_Tree[N_OMs];
	
	for(int OM_num = 0; OM_num < N_OMs; OM_num++)
	{
		TString treename = Form("Tracks of OM %d", OM_num);		
		
		Tree[OM_num] = new TTree(treename, treename);
		TBranch* branchA0 = Tree[OM_num]->Branch("A", &A_Tree[OM_num], "A/D");
		TBranch* branchB0 = Tree[OM_num]->Branch("B", &B_Tree[OM_num], "B/D");
		TBranch* branchC0 = Tree[OM_num]->Branch("C", &C_Tree[OM_num], "C/D");
		TBranch* branchD0 = Tree[OM_num]->Branch("D", &D_Tree[OM_num], "D/D");
	}

	// Filling Tree
	// Use s->GetEntries() for all entries in Run
	//for(UInt_t i=0; i < 1000000; i++)	// Loop over events
	for(UInt_t i=0; i < s->GetEntries(); i++)	// Loop over events
	{
		s->GetEntry(i);
		Eve->set_r("Manchester", "distance");
		Eve->set_h(eff_len);
		Eve->reconstruct_ML(0);		
	
		int event_choice = event_selection(Eve, 0.6, 0.98, 1);
		if (event_choice)
		{	
			int track_num = event_choice - 1;
			int OM_num = Eve->get_cluster(0)->get_tr_hits()[0]->get_associated_OMhit()->get_OM_num();

			double Y = Eve->get_track(track_num)->get_b();
			double Z = Eve->get_track(track_num)->get_d();

			if(Z!=0 && Z > Z_MIN && Z < Z_MAX && Y > Y_MIN && Y < Y_MAX)
			{
				A_Tree[OM_num] = Eve->get_track(track_num)->get_a();
				B_Tree[OM_num] = Eve->get_track(track_num)->get_b();
				C_Tree[OM_num] = Eve->get_track(track_num)->get_c();
				D_Tree[OM_num] = Eve->get_track(track_num)->get_d();
				Tree[OM_num]->Fill(); 
			}
		}
		if (i % 100000 == 0) cout <<"Event No. " << i << " done!" <<endl;
	}	


        for(int OM_num = 0; OM_num < N_OMs; OM_num++)
        {
		TString treename = Form("Tracks of OM %d", OM_num);		
		
		// write histograms into new root file
		Tree[OM_num]->Write(treename);
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
	else if(_event->get_no_tracks()             == 1        &&
                _event->get_track(0)->get_quality() >  _low_qlt &&
                _event->get_track(0)->get_quality() <  _high_qlt)
        {
                choice = 1;
        }
	else if(_event->get_no_tracks()                    == 2        &&
                _event->get_track(0)->get_ambiguity_type() != 0        &&
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
                int NSOR1 = N_COLS * ((int)(Z_MAX - Z1) / 468) + (int)(Y1 - Y_MIN) / 833;
                int NSOR2 = N_COLS * ((int)(Z_MAX - Z2) / 468) + (int)(Y2 - Y_MIN) / 833;

                double Y_SOR1 = Y_MIN + (NSOR1 % N_COLS + 0.5) * 833;
                double Z_SOR1 = Z_MAX - (NSOR1 / N_COLS + 0.5) * 468;

                double Y_SOR2 = Y_MIN + (NSOR2 % N_COLS + 0.5) * 833;
                double Z_SOR2 = Z_MAX - (NSOR2 / N_COLS + 0.5) * 468;

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


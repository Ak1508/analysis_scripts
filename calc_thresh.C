// Author Eric Pooser, pooser@jlab.org
// Script to calculate pedestal values on crate, slot, channel basis

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <numeric>
#include <vector>
#include <time.h>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TF1.h"

#define NUMSLOTS 21
#define NADCCHAN 16
#define FADCMODE 10
#define NPOINTS  25
#define NSIGMA   5
#define CHANTHRSH 20

using namespace std;

TFile *rif, *rof;
TH2I *h2_pped[NUMSLOTS];
TH1I *h_pped[NUMSLOTS][NADCCHAN];
TF1 *init_gfit[NUMSLOTS][NADCCHAN], *gfit[NUMSLOTS][NADCCHAN];
uint32_t nentries[NUMSLOTS][NADCCHAN], threshhold[NUMSLOTS][NADCCHAN];
Double_t init_max[NUMSLOTS][NADCCHAN], init_mean[NUMSLOTS][NADCCHAN], init_stddev[NUMSLOTS][NADCCHAN];
Double_t iter_max[NUMSLOTS][NADCCHAN], iter_mean[NUMSLOTS][NADCCHAN], iter_stddev[NUMSLOTS][NADCCHAN];
Double_t finl_max[NUMSLOTS][NADCCHAN], finl_mean[NUMSLOTS][NADCCHAN], finl_stddev[NUMSLOTS][NADCCHAN];
Double_t finl_max_err[NUMSLOTS][NADCCHAN], finl_mean_err[NUMSLOTS][NADCCHAN], finl_stddev_err[NUMSLOTS][NADCCHAN];
Double_t fr_low[NUMSLOTS][NADCCHAN], fr_high[NUMSLOTS][NADCCHAN];

void calc_thresh() {

  // Get the input file and declare the output file
  rif = new TFile("tstfadc.root");
  rof = new TFile("ped.root", "RECREATE");

  // Declare the output text file
  ofstream outputfile;
  
  // Acquire the 2d histos and their respective projections
  for (uint32_t islot = 0; islot < NUMSLOTS; islot++) {
    h2_pped[islot] = (TH2I*) (rif->Get(Form("mode_%d_data/slot_%d/h2_pped", FADCMODE, islot+1)));
    if (h2_pped[islot]) {
      for (uint32_t ichan = 0; ichan < NADCCHAN; ichan++) {
	h_pped[islot][ichan] = (TH1I*) h2_pped[islot]->ProjectionY(Form("h_pped_slot_%d_chan_%d", islot+1, ichan), ichan, ichan);
	// Acquire the number of entries as a restriction for a good fit
	nentries[islot][ichan] = h_pped[islot][ichan]->GetEntries();
	if (h_pped[islot][ichan] && (nentries[islot][ichan] >= NPOINTS)) {
	  // Define the fit
	  init_gfit[islot][ichan] = new TF1(Form("init_fit_slot_%d_chan_%d", islot+1, ichan), "gaus(0)");
	  init_gfit[islot][ichan]->SetLineColor(2);
	  // Acquire the initial fit parameters and initialize the fit
	  init_max[islot][ichan]    = h_pped[islot][ichan]->GetMaximum();
	  init_mean[islot][ichan]   = h_pped[islot][ichan]->GetMean();
	  init_stddev[islot][ichan] = h_pped[islot][ichan]->GetStdDev();
	  init_gfit[islot][ichan]->SetParameters(init_max[islot][ichan],
						 init_mean[islot][ichan],
						 init_stddev[islot][ichan]);
	  // Perform the initial fit over the full range of the histo
	  h_pped[islot][ichan]->Fit(Form("init_fit_slot_%d_chan_%d", islot+1, ichan), "Q");
	  init_gfit[islot][ichan]->Draw();
	  // Declare and initialize the second and final fit
	  gfit[islot][ichan] = new TF1(Form("iter_fit_slot_%d_chan_%d", islot+1, ichan), "gaus(0)");
	  gfit[islot][ichan]->SetLineColor(4);
	  // Acquire the fit parameters from the first fit and initialize the second fit
	  iter_max[islot][ichan]    = init_gfit[islot][ichan]->GetParameter(0);
	  iter_mean[islot][ichan]   = init_gfit[islot][ichan]->GetParameter(1);
	  iter_stddev[islot][ichan] = init_gfit[islot][ichan]->GetParameter(2);
	  gfit[islot][ichan]->SetParameters(iter_max[islot][ichan],
					    iter_mean[islot][ichan],
					    iter_stddev[islot][ichan]);
	  gfit[islot][ichan]->SetParNames("Amp.", "#mu (ADC)", "#sigma (ADC)");
	  // Fit over a mean +/- NSIGMA*mean range
	  fr_low[islot][ichan]  = iter_mean[islot][ichan] - NSIGMA*iter_stddev[islot][ichan];
	  fr_high[islot][ichan] = iter_mean[islot][ichan] + NSIGMA*iter_stddev[islot][ichan];
	  // cout << "MEAN = " << iter_mean[islot][ichan]
	  //      << ", FR LOW = " << fr_low[islot][ichan]
	  //      << ", FR HIGH = " << fr_high[islot][ichan] << endl;
	  gfit[islot][ichan]->SetRange(fr_low[islot][ichan], fr_high[islot][ichan]); 
	  // Perform the fit and store the parameters
	  h_pped[islot][ichan]->Fit(Form("iter_fit_slot_%d_chan_%d", islot+1, ichan), "QR");
	  gfit[islot][ichan]->Draw("same");
	  // Acquire the final fit paramters
	  finl_max[islot][ichan]    = gfit[islot][ichan]->GetParameter(0);
	  finl_mean[islot][ichan]   = gfit[islot][ichan]->GetParameter(1);
	  finl_stddev[islot][ichan] = gfit[islot][ichan]->GetParameter(2);
	  threshhold[islot][ichan]  = (UInt_t) finl_mean[islot][ichan] + CHANTHRSH;
	  // Acquire the final fit parameter errors
	  finl_max_err[islot][ichan]    = gfit[islot][ichan]->GetParError(0);
	  finl_mean_err[islot][ichan]   = gfit[islot][ichan]->GetParError(1);
	  finl_stddev_err[islot][ichan] = gfit[islot][ichan]->GetParError(2);
	}  // Projection objection and nentries requirement
      }  // Channel loop
    }  // 2D histo object requirement
  }  // Slot loop

  outputfile.open("thresholds.dat");
  for (uint32_t islot = 0; islot < NUMSLOTS; islot++) {
    outputfile << "slot=" << islot + 1 << endl;
    for (uint32_t ichan = 0; ichan < NADCCHAN; ichan++)
      outputfile << threshhold[islot][ichan] << endl;
  }  // Slot loop
  
  rof->Write();
  
}


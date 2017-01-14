
// Custom analysis of the hodoscope individual channels
// Author: Eric Pooser, pooser@jlab.org, 01/06/2017
// Script adapted from hcana/examples/hhodrawhists.C

#define NPLANES  4
#define NSIDES   2
#define NSIGNALS 3
//#define NSIGNALS 2
#define NADCSIGNALS 1
#define MAXBARS  16

TString  SPECTROMETER = "H";
TString  DETECTOR     = "hod";

TString  plane_names[NPLANES] = {"1x", "1y", "2x", "2y"};
Int_t    nbars[NPLANES]       = {16, 10, 16, 10};
TString  sides[NSIDES]        = {"neg", "pos"};
//TString  signals[NSIGNALS]    = {"adc", "tdc"};
TString  signals[NSIGNALS]    = {"adc", "tdc", "ped"};
//TString  adc_signals          = {"ped"}

static const Double_t ADC_MIN   = 0.0;
static const Double_t ADC_MAX   = 25000.0;
static const Int_t    ADC_NBINS = 2500;

static const Double_t TDC_MIN   = 0.0;
static const Double_t TDC_MAX   = 8000.0;
static const Int_t    TDC_NBINS = 800;

Int_t    nhits[NPLANES][NSIDES][NSIGNALS];
Double_t paddles[NPLANES][NSIDES][NSIGNALS][MAXBARS];
Double_t values[NPLANES][NSIDES][NSIGNALS][MAXBARS];
Double_t ped_values[NPLANES][NSIDES][NSIGNALS][MAXBARS];

TH1F *h[NPLANES*NSIDES*NSIGNALS*MAXBARS];
TFile *rif, *rof;

TTree *T;

TString base_name, ndata_name, padlist_name, vallist_name, pedlist_name;

Int_t nbins, hmin, hmax, hindex, hindex_base;

TString ibarname, title, name;

Int_t ibar;
Double_t val, ped_val;

TDirectory *h1x_dir, *h1y_dir, *h2x_dir, *h2y_dir;

TH1F *h1x_minus_ht1, *h1y_minus_ht1, *h2x_minus_ht1, *h2y_minus_ht1;
TH1F *h1x_minus_ht2, *h1y_minus_ht2, *h2x_minus_ht2, *h2y_minus_ht2;
TH1F *h1t_minus_ht1, *h2t_minus_ht1, *h1t_minus_ht2, *h2t_minus_ht2;

Double_t h1x_tdc, h1y_tdc, h2x_tdc, h2y_tdc;
Double_t h1t_tdc, h2t_tdc, ht1_tdc, ht2_tdc;

TString h1x_tdc_address = "T.hms.h1X_tdc";
TString h1y_tdc_address = "T.hms.h1Y_tdc";
TString h2x_tdc_address = "T.hms.h2X_tdc";
TString h2y_tdc_address = "T.hms.h2Y_tdc";

TString h1t_tdc_address = "T.hms.h1T_tdc";
TString h2t_tdc_address = "T.hms.h2T_tdc";
TString ht1_tdc_address = "T.hms.hT1_tdc";
TString ht2_tdc_address = "T.hms.hT2_tdc";

void hhodo_analysis(Int_t RunNumber=0, Int_t MaxEvent=0) {

  // Get RunNumber and MaxEvent if not provided.
  if(RunNumber == 0) {
    cout << "Enter a Run Number (-1 to exit): ";
    cin >> RunNumber;
    if( RunNumber<=0 ) return;
  }
  if(MaxEvent == 0) {
    cout << "\nNumber of Events to analyze: ";
    cin >> MaxEvent;
    if(MaxEvent == 0) {
      cerr << "...Invalid entry\n";
      exit;
    }
  } 

  // Define root in/out files
  rif  = new TFile(Form("ROOTfiles/hhodo_htrig_replay_%d.root", RunNumber), "READ");
  rof = new TFile(Form("ROOTfiles/hhodo_analysis_%d.root", RunNumber), "RECREATE");
  rif->cd();
  
  // Acquire tree from root in file
  T = new TTree();
  T = (TTree*) rif->Get("T");

  // Trigger class variables
  T->SetBranchAddress(h1x_tdc_address, &h1x_tdc);
  T->SetBranchAddress(h1y_tdc_address, &h1y_tdc);
  T->SetBranchAddress(h2x_tdc_address, &h2x_tdc);
  T->SetBranchAddress(h2y_tdc_address, &h2y_tdc);

  T->SetBranchAddress(h1t_tdc_address, &h1t_tdc);
  T->SetBranchAddress(h2t_tdc_address, &h2t_tdc);
  T->SetBranchAddress(ht1_tdc_address, &ht1_tdc);
  T->SetBranchAddress(ht2_tdc_address, &ht2_tdc);

  rof->cd();
  
  h1x_minus_ht1 = new TH1F("h1x_minus_ht1", "h1x - ht1 TDC; Raw TDC (TDC Units); Counts / 10 TDC Units", 1000, -5000, 5000);
  h1y_minus_ht1 = new TH1F("h1y_minus_ht1", "h1y - ht1 TDC; Raw TDC (TDC Units); Counts / 10 TDC Units", 1000, -5000, 5000);
  h2x_minus_ht1 = new TH1F("h2x_minus_ht1", "h2x - ht1 TDC; Raw TDC (TDC Units); Counts / 10 TDC Units", 1000, -5000, 5000);
  h2y_minus_ht1 = new TH1F("h2y_minus_ht1", "h2y - ht1 TDC; Raw TDC (TDC Units); Counts / 10 TDC Units", 1000, -5000, 5000);
  h1t_minus_ht1 = new TH1F("h1t_minus_ht1", "h1t - ht1 TDC; Raw TDC (TDC Units); Counts / 10 TDC Units", 1000, -5000, 5000);
  h2t_minus_ht1 = new TH1F("h2t_minus_ht1", "h2t - ht1 TDC; Raw TDC (TDC Units); Counts / 10 TDC Units", 1000, -5000, 5000);

  h1x_minus_ht2 = new TH1F("h1x_minus_ht2", "h1x - ht2 TDC; Raw TDC (TDC Units); Counts / 10 TDC Units", 1000, -5000, 5000);
  h1y_minus_ht2 = new TH1F("h1y_minus_ht2", "h1y - ht2 TDC; Raw TDC (TDC Units); Counts / 10 TDC Units", 1000, -5000, 5000);
  h2x_minus_ht2 = new TH1F("h2x_minus_ht2", "h2x - ht2 TDC; Raw TDC (TDC Units); Counts / 10 TDC Units", 1000, -5000, 5000);
  h2y_minus_ht2 = new TH1F("h2y_minus_ht2", "h2y - ht2 TDC; Raw TDC (TDC Units); Counts / 10 TDC Units", 1000, -5000, 5000);
  h1t_minus_ht2 = new TH1F("h1t_minus_ht2", "h1t - ht2 TDC; Raw TDC (TDC Units); Counts / 10 TDC Units", 1000, -5000, 5000);
  h2t_minus_ht2 = new TH1F("h2t_minus_ht2", "h2t - ht2 TDC; Raw TDC (TDC Units); Counts / 10 TDC Units", 1000, -5000, 5000);

  for(Int_t iplane = 0; iplane < NPLANES; iplane++) {

    //rof->cd();
    
    if (plane_names[iplane] == "1x") {
      h1x_dir = dynamic_cast <TDirectory*> (rof->Get("h1x"));
      if (!h1x_dir) {h1x_dir = rof->mkdir("h1x"); h1x_dir->cd();}
      else rof->cd("h1x");
    }
    
    if (plane_names[iplane] == "1y") {
      h1y_dir = dynamic_cast <TDirectory*> (rof->Get("h1y"));
      if (!h1y_dir) {h1y_dir = rof->mkdir("h1y"); h1y_dir->cd();}
      else rof->cd("h1y");
    }

    if (plane_names[iplane] == "2x") {
      h2x_dir = dynamic_cast <TDirectory*> (rof->Get("h2x"));
      if (!h2x_dir) {h2x_dir = rof->mkdir("h2x"); h2x_dir->cd();}
      else rof->cd("h2x");
    }

    if (plane_names[iplane] == "2y") {
      h2y_dir = dynamic_cast <TDirectory*> (rof->Get("h2y"));
      if (!h2y_dir) {h2y_dir = rof->mkdir("h2y"); h2y_dir->cd();}
      else rof->cd("h2y");
    }
						   
    for(Int_t iside = 0; iside < NSIDES; iside++) {
      for(Int_t isignal = 0; isignal < NSIGNALS; isignal++) {

	// if (signals[isignal] == "adc" || signals[isignal] == "tdc") {
	//   base_name = SPECTROMETER + "." + DETECTOR + "." +
	//     plane_names[iplane] + "." + sides[iside] + signals[isignal];
	//   ndata_name   = "Ndata." + base_name + "pad";
	//   padlist_name = base_name + "pad";
	//   vallist_name = base_name + "val";
	//  }
	if (signals[isignal] == "adc" || signals[isignal] == "tdc") {
	  base_name = SPECTROMETER + "." + DETECTOR + "." +
	    plane_names[iplane] + "." + sides[iside] + signals[isignal];
	  ndata_name   = "Ndata." + base_name + "pad";
	  padlist_name = base_name + "pad";
	  vallist_name = base_name + "val";
	}
	
	//if (signals[isignal] == "adc") pedlist_name = base_name + "ped";
	if (signals[isignal] == "ped") {
	  base_name = SPECTROMETER + "." + DETECTOR + "." +
	    plane_names[iplane] + "." + sides[iside] + "adc";
	  ndata_name   = "Ndata." + base_name + "ped";
	  padlist_name = base_name + "pad";
	  pedlist_name = base_name + "ped";
	}

	// Set branch addresses
	T->SetBranchAddress(ndata_name,   &nhits[iplane][iside][isignal]);

	if (signals[isignal] == "adc" || signals[isignal] == "tdc") {
	  cout << padlist_name << endl;
	  cout << isignal << endl;
	  cout << paddles[iplane][iside][isignal][0] << endl;
	  T->SetBranchAddress(padlist_name, &paddles[iplane][iside][isignal][0]);
	  T->SetBranchAddress(vallist_name, &values[iplane][iside][isignal][0]);
	}
	if (signals[isignal] == "ped") {
	  cout << padlist_name << endl;
	  cout << isignal << endl;
	  cout << paddles[iplane][iside][isignal][0] << endl;
	  //T->SetBranchAddress(padlist_name, &paddles[iplane][iside][isignal][0]);
	  T->SetBranchAddress(pedlist_name, &ped_values[iplane][iside][isignal][0]);
	}
	// Create histograms
	// ADC and TDC histogram for each readout channel
	if(signals[isignal] == "adc")
	  {nbins = ADC_NBINS; hmin = ADC_MIN; hmax = ADC_MAX;}
	if(signals[isignal] == "ped")
	  {nbins = ADC_NBINS; hmin = ADC_MIN; hmax = ADC_MAX;}
	if(signals[isignal] == "tdc")
	  {nbins = TDC_NBINS; hmin = TDC_MIN; hmax = TDC_MAX;}
	
	if (plane_names[iplane] == "1x") rof->cd("h1x");
	if (plane_names[iplane] == "1y") rof->cd("h1y");
	if (plane_names[iplane] == "2x") rof->cd("h2x");
	if (plane_names[iplane] == "2y") rof->cd("h2y");
				
	for(Int_t ibar = 0; ibar < nbars[iplane]; ibar++) {
	  hindex = iplane*NSIDES*NSIGNALS*MAXBARS
	    + iside*NSIGNALS*MAXBARS + isignal*MAXBARS + ibar;
	  ibarname = Form("%d", ibar + 1);
	  if (sides[iside] == "pos") title = "h" + plane_names[iplane] + ibarname + "+ ";
	  if (sides[iside] == "neg") title = "h" + plane_names[iplane] + ibarname + "- ";
	  if (signals[isignal] == "adc")
	    title += Form("FADC Pulse Integral; FADC Pulse Integral (ADC Units); Counts / %d ADC Units", Int_t ((ADC_MAX - ADC_MIN) / ADC_NBINS));
	  if (signals[isignal] == "ped")
	    title += Form("FADC Pulse Pedestal; FADC Pulse Pedestal (ADC Units); Counts / %d ADC Units", Int_t ((ADC_MAX - ADC_MIN) / ADC_NBINS));
	  if (signals[isignal] == "tdc")
	    title += Form("Raw TDC; Raw TDC (TDC Units); Counts / %d TDC Units", Int_t ((TDC_MAX - TDC_MIN) / TDC_NBINS));
	  if (sides[iside] == "pos") name = "h" + plane_names[iplane] + ibarname + "+_" + signals[isignal];
	  if (sides[iside] == "neg") name = "h" + plane_names[iplane] + ibarname + "-_" + signals[isignal];

	  h[hindex] = new TH1F(name, title, nbins, hmin, hmax);
	}
	rif->cd();

      }  // Signal loop
    }  // Side loop
  }  // Plane loop

  // Loop over the events, filling the histograms
  for(Int_t ievent = 0, N = T->GetEntries(); ievent < N; ievent++) {
    
    T->GetEntry(ievent);
    
    h1x_minus_ht1->Fill(TMath::Abs(h1x_tdc - ht1_tdc));
    h1y_minus_ht1->Fill(TMath::Abs(h1y_tdc - ht1_tdc));
    h2x_minus_ht1->Fill(TMath::Abs(h2x_tdc - ht1_tdc));
    h2y_minus_ht1->Fill(TMath::Abs(h2y_tdc - ht1_tdc));
    h1t_minus_ht1->Fill(TMath::Abs(h1t_tdc - ht1_tdc));
    h2t_minus_ht1->Fill(TMath::Abs(h2t_tdc - ht1_tdc));

    h1x_minus_ht2->Fill(TMath::Abs(h1x_tdc - ht2_tdc));
    h1y_minus_ht2->Fill(TMath::Abs(h1y_tdc - ht2_tdc));
    h2x_minus_ht2->Fill(TMath::Abs(h2x_tdc - ht2_tdc));
    h2y_minus_ht2->Fill(TMath::Abs(h2y_tdc - ht2_tdc));
    h1t_minus_ht2->Fill(TMath::Abs(h1t_tdc - ht2_tdc));
    h2t_minus_ht2->Fill(TMath::Abs(h2t_tdc - ht2_tdc));
			
    for(Int_t iplane = 0; iplane < NPLANES; iplane++){
      for(Int_t iside = 0; iside < NSIDES; iside++) {
	for(Int_t isignal = 0; isignal < NSIGNALS; isignal++) {

	  //cout << signals[isignal] << "     " << nhits[iplane][iside][isignal] << endl;

	  hindex_base = iplane*NSIDES*NSIGNALS*MAXBARS
	    + iside*NSIGNALS*MAXBARS + isignal*MAXBARS;

	  for(Int_t ihit = 0; ihit < nhits[iplane][iside][isignal]; ihit++) {
	    
	    ibar = TMath::Nint(paddles[iplane][iside][isignal][ihit]) - 1;
	    
	    if (signals[isignal] == "adc" || signals[isignal] == "tdc") {
	      hindex = hindex_base + ibar;
	      val = values[iplane][iside][isignal][ihit];
	      h[hindex]->Fill(val);
	      //cout << "HINDEX WITH ADC & TDC = " << hindex << endl;
	    }
	    
	    if (signals[isignal] == "ped") {
	      //ibar = TMath::Nint(paddles[iplane][iside][isignal][ihit]) - 1;
	      hindex = hindex_base + ibar;
	      ped_val = ped_values[iplane][iside][isignal][ihit];
	      //cout << "PED = " << ped_val << endl;

	      //cout << "HINDEX WITH PED = " << hindex << endl;
	      
	      //h[hindex]->Fill(ped_val);
	    }
	    
	  }  // Hit loop
	}  // Signal loop
      }  // Side loop 
    }  // Plane loop

    // Display or save the histograms

  }

  // Write the root out file
  rof->cd(); rof->Write();
  
}

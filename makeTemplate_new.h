#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <vector>
#include <cmath>

using namespace std;

const int NCH = 1;
const double bins[NCH] = {100};

const int nSamples = 4;

const TString mcNames[nSamples] = {"ttJetsHadronic", "ttJetsSemiLep", "ttJetsFullLep",
				   "ttA_2to5"};
  
const int mcLayerNumbers[nSamples] = {0, 0, 0, 
				      6};

const double epsilon = 1e-10;

const TString work_dir = "./";
const TString datacard_dir = work_dir+"datacards";

const double lumiSystematic = 0.026;

class BackgroundProcess {

 public:
  BackgroundProcess() {;};
  BackgroundProcess(int layer, TString title, TString binTitle,
		    vector<TH1D*> vhist,
		    vector<TH1D*> vhist_btagWeightUp, vector<TH1D*> vhist_btagWeightDown,
		    vector<TH1D*> vhist_puWeightUp, vector<TH1D*> vhist_puWeightDown,
		    vector<TH1D*> vhist_topPtUp, vector<TH1D*> vhist_topPtDown,
		    vector<TH1D*> vhist_JECup, vector<TH1D*> vhist_JECdown,
		    vector<TH1D*> vhist_leptonSFup, vector<TH1D*> vhist_leptonSFdown,
		    vector<TH1D*> vhist_photonSFup, vector<TH1D*> vhist_photonSFdown,
		    vector<TH1D*> vhist_scaleUp, vector<TH1D*> vhist_scaleDown,
		    vector<TH1D*> vhist_pdfUp, vector<TH1D*> vhist_pdfDown) {

    name = title;
    binName = binTitle;

    TH1D *h,
      *h_btagWeightUp, *h_btagWeightDown,
      *h_puWeightUp, *h_puWeightDown,
      *h_topPtUp, *h_topPtDown,
      *h_JECup, *h_JECdown,
      *h_leptonSFup, *h_leptonSFdown,
      *h_photonSFup, *h_photonSFdown,
      *h_scaleUp, *h_scaleDown,
      *h_pdfUp, *h_pdfDown;
      
    for(unsigned int i = 0; i < vhist.size(); i++) {
      if(mcLayerNumbers[i] == layer) {
	if(i == 0 || mcLayerNumbers[i] != mcLayerNumbers[i-1]) {
	  h = (TH1D*)vhist[i]->Clone();
	  h_btagWeightUp = (TH1D*)vhist_btagWeightUp[i]->Clone();
	  h_btagWeightDown = (TH1D*)vhist_btagWeightDown[i]->Clone();
	  h_puWeightUp = (TH1D*)vhist_puWeightUp[i]->Clone();
	  h_puWeightDown = (TH1D*)vhist_puWeightDown[i]->Clone();
	  h_topPtUp = (TH1D*)vhist_topPtUp[i]->Clone();
	  h_topPtDown = (TH1D*)vhist_topPtDown[i]->Clone();
	  h_JECup = (TH1D*)vhist_JECup[i]->Clone();
	  h_JECdown = (TH1D*)vhist_JECdown[i]->Clone();
	  h_leptonSFup = (TH1D*)vhist_leptonSFup[i]->Clone();
	  h_leptonSFdown = (TH1D*)vhist_leptonSFdown[i]->Clone();
	  h_photonSFup = (TH1D*)vhist_photonSFup[i]->Clone();
	  h_photonSFdown = (TH1D*)vhist_photonSFdown[i]->Clone();
	  h_scaleUp = (TH1D*)vhist_scaleUp[i]->Clone();
	  h_scaleDown = (TH1D*)vhist_scaleDown[i]->Clone();
	  h_pdfUp = (TH1D*)vhist_pdfUp[i]->Clone();
	  h_pdfDown = (TH1D*)vhist_pdfDown[i]->Clone();
	}
      }
    }

    for(unsigned int i = 1; i < vhist.size(); i++) {
      if(mcLayerNumbers[i] == layer && mcLayerNumbers[i] == mcLayerNumbers[i-1]) {
	h->Add(vhist[i]);
	h_btagWeightUp->Add(vhist_btagWeightUp[i]);
	h_btagWeightDown->Add(vhist_btagWeightDown[i]);
	h_puWeightUp->Add(vhist_puWeightUp[i]);
	h_puWeightDown->Add(vhist_puWeightDown[i]);
	h_topPtUp->Add(vhist_topPtUp[i]);
	h_topPtDown->Add(vhist_topPtDown[i]);
	h_JECup->Add(vhist_JECup[i]);
	h_JECdown->Add(vhist_JECdown[i]);
	h_leptonSFup->Add(vhist_leptonSFup[i]);
	h_leptonSFdown->Add(vhist_leptonSFdown[i]);
	h_photonSFup->Add(vhist_photonSFup[i]);
	h_photonSFdown->Add(vhist_photonSFdown[i]);
	h_scaleUp->Add(vhist_scaleUp[i]);
	h_scaleDown->Add(vhist_scaleDown[i]);
	h_pdfUp->Add(vhist_pdfUp[i]);
	h_pdfDown->Add(vhist_pdfDown[i]);
      }
    }

    int binNumber = 0;

    int binLo = h->GetXaxis()->FindBin(bins[binNumber]);
    int binHi = (binNumber < NCH - 1) ? h->GetXaxis()->FindBin(bins[binNumber+1]) - 1 : -1;

    value = h->IntegralAndError(binLo, binHi, stat);
    stat = 1. + stat / value;

    double temp;

    double btag_up = h_btagWeightUp->IntegralAndError(binLo, binHi, temp);
    double btag_down = h_btagWeightDown->IntegralAndError(binLo, binHi, temp);
    btag = max(fabs(btag_up - value), fabs(btag_down - value));
    btag = 1. + btag / value;

    double pileup_up = h_puWeightUp->IntegralAndError(binLo, binHi, temp);
    double pileup_down = h_puWeightDown->IntegralAndError(binLo, binHi, temp);
    pileup = max(fabs(pileup_up - value), fabs(pileup_down - value));
    pileup = 1. + pileup / value;

    double topPt_up = h_topPtUp->IntegralAndError(binLo, binHi, temp);
    double topPt_down = h_topPtDown->IntegralAndError(binLo, binHi, temp);
    topPt = max(fabs(topPt_up - value), fabs(topPt_down - value));
    topPt = 1. + topPt / value;

    double JEC_up = h_JECup->IntegralAndError(binLo, binHi, temp);
    double JEC_down = h_JECdown->IntegralAndError(binLo, binHi, temp);
    jec = max(fabs(JEC_up - value), fabs(JEC_down - value));
    jec = 1. + jec / value;

    double lepton_up = h_leptonSFup->IntegralAndError(binLo, binHi, temp);
    double lepton_down = h_leptonSFdown->IntegralAndError(binLo, binHi, temp);
    leptonID = max(fabs(lepton_up - value), fabs(lepton_down - value));
    leptonID = 1. + leptonID / value;

    double photon_up = h_photonSFup->IntegralAndError(binLo, binHi, temp);
    double photon_down = h_photonSFdown->IntegralAndError(binLo, binHi, temp);
    photonID = max(fabs(photon_up - value), fabs(photon_down - value));
    photonID = 1. + photonID / value;

    double scale_up = h_scaleUp->IntegralAndError(binLo, binHi, temp);
    double scale_down = h_scaleDown->IntegralAndError(binLo, binHi, temp);
    scale = max(fabs(scale_up - value), fabs(scale_down - value));
    scale = 1. + scale / value;

    double pdf_up = h_pdfUp->IntegralAndError(binLo, binHi, temp);
    double pdf_down = h_pdfDown->IntegralAndError(binLo, binHi, temp);
    pdf = max(fabs(pdf_up - value), fabs(pdf_down - value));
    pdf = 1. + pdf / value;
  }
    
  virtual ~BackgroundProcess() {;};

  TString name;
  TString binName;

  double value;
  double stat;
  
  double btag;
  double pileup;
  double scale;
  double pdf;
  double topPt;
  double jec;
  double leptonID;
  double photonID;
 
};

class ObservedData {

 public:
  ObservedData() {;};
  ObservedData(TString binTitle, TH1D * h) {

    binName = binTitle;

    int binNumber = 0;

    int binLo = h->GetXaxis()->FindBin(bins[binNumber]);
    int binHi = (binNumber < NCH - 1) ? h->GetXaxis()->FindBin(bins[binNumber+1]) - 1 : -1;

    value = h->IntegralAndError(binLo, binHi, stat);

    stat = 1. + stat / value;
  }

  virtual ~ObservedData() {;};

  TString binName;

  double value;
  double stat;
};

class SignalYield {

 public:
  SignalYield() {;};
  SignalYield(TString binTitle,
	      TH1D * h,
	      TH1D * h_btagWeightUp, TH1D * h_btagWeightDown,
	      TH1D * h_puWeightUp, TH1D * h_puWeightDown,
	      TH1D * h_topPtUp, TH1D * h_topPtDown,
	      TH1D * h_JECup, TH1D * h_JECdown,
	      TH1D * h_leptonSFup, TH1D * h_leptonSFdown,
	      TH1D * h_photonSFup, TH1D * h_photonSFdown,
	      double crossSection, double crossSection_uncertainty) {
    
    binName = binTitle;

    int binNumber = 0;

    int binLo = h->GetXaxis()->FindBin(bins[binNumber]);
    int binHi = (binNumber < NCH - 1) ? h->GetXaxis()->FindBin(bins[binNumber+1]) - 1 : -1;

    value = h->IntegralAndError(binLo, binHi, stat);
    stat = 1. + stat / value;

    double temp;

    double btag_up = h_btagWeightUp->IntegralAndError(binLo, binHi, temp);
    double btag_down = h_btagWeightDown->IntegralAndError(binLo, binHi, temp);
    btag = max(fabs(btag_up - value), fabs(btag_down - value));
    btag = 1. + btag / value;

    double pileup_up = h_puWeightUp->IntegralAndError(binLo, binHi, temp);
    double pileup_down = h_puWeightDown->IntegralAndError(binLo, binHi, temp);
    pileup = max(fabs(pileup_up - value), fabs(pileup_down - value));
    pileup = 1. + pileup / value;

    double topPt_up = h_topPtUp->IntegralAndError(binLo, binHi, temp);
    double topPt_down = h_topPtDown->IntegralAndError(binLo, binHi, temp);
    topPt = max(fabs(topPt_up - value), fabs(topPt_down - value));
    topPt = 1. + topPt / value;

    double JEC_up = h_JECup->IntegralAndError(binLo, binHi, temp);
    double JEC_down = h_JECdown->IntegralAndError(binLo, binHi, temp);
    jec = max(fabs(JEC_up - value), fabs(JEC_down - value));
    jec = 1. + jec / value;

    double lepton_up = h_leptonSFup->IntegralAndError(binLo, binHi, temp);
    double lepton_down = h_leptonSFdown->IntegralAndError(binLo, binHi, temp);
    leptonID = max(fabs(lepton_up - value), fabs(lepton_down - value));
    leptonID = 1. + leptonID / value;

    double photon_up = h_photonSFup->IntegralAndError(binLo, binHi, temp);
    double photon_down = h_photonSFdown->IntegralAndError(binLo, binHi, temp);
    photonID = max(fabs(photon_up - value), fabs(photon_down - value));
    photonID = 1. + photonID / value;

    xsec = crossSection;
    xsecError = 1. + crossSection_uncertainty / 100.;
  }
    
  void FillHistograms(TH2D*& h_yield, TH2D*& h_stat, TH2D*& h_btag, TH2D*& h_btag, TH2D*& h_pileup, TH2D*& h_jec, TH2D*& h_leptonSF, TH2D*& h_photonSF) {
    h_yield->Fill(value);
    h_stat->Fill(stat - 1);
    h_btag->Fill(btag - 1);
    h_pileup->Fill(pileup - 1);
    h_jec->Fill(jec - 1);
    h_leptonSF->Fill(leptonID - 1);
    h_photonSF->Fill(photonID - 1);
  }
    

  virtual ~SignalYield() {;};

  TString binName;

  double value;
  double stat;
  
  double btag;
  double pileup;
  double topPt;
  double jec;
  double leptonID;
  double photonID;

  double xsec;
  double xsecError;

};

class GridPoint {

 public:
  GridPoint();
  virtual ~GridPoint() {
    data.clear();
    ttbar.clear();
    ttgamma.clear();
    signal.clear();
  }

  void Print();

  int mStop;
  int mBino;

  double xsec;
  double xsecError;

  double lumi_sysError;
  
  vector<BackgroundProcess> ttbar;
  vector<BackgroundProcess> ttgamma;
  vector<ObservedData> data;
  vector<SignalYield> signal;

  double limit;           // observed limit
  double explimit;        // expected limit
  double explimit_1L;     // expected limit -1 sigma
  double explimit_1H;     // expected limit +1 sigma
  double explimit_2L;     // expected limit -2 sigma
  double explimit_2H;     // expected limit +2 sigma
};

GridPoint::GridPoint() {

  mStop = mBino = 0;

  lumi_sysError = 1.026;

  data.clear();
  ttbar.clear();
  ttgamma.clear();
  signal.clear();

  limit = explimit = explimit_1L = explimit_1H = explimit_2L = explimit_2H = 0;
}

void GridPoint::Print() {

  stringstream outname;
  outname << datacard_dir.Data() << "/stop-bino_mst_" << mStop << "_m1_" << mBino << ".dat";
  fstream outfile(outname.str().c_str(), ios::out);

  outfile << "# stop = " << mStop << endl;
  outfile << "# bino = " << mBino << endl;
  outfile << "# xsec = " << xsec << endl;
  outfile << "# xsec uncertainty = " << xsecError << " %" << endl;

  vector<unsigned int> sensitive_bins;
  for(unsigned int i = 0; i < signal.size(); i++) {
    if(signal[i].value > epsilon) sensitive_bins.push_back(i);
  }

   // Compute R_firstguess
  double Rmin = 99999;
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) {
    int bin = sensitive_bins[i];
    
    double unc2 = data[bin].value;
    unc2 += pow(lumi_sysError - 1., 2);

    unc2 += pow(signal[bin].stat - 1., 2);
    unc2 += pow(signal[bin].btag - 1., 2);
    unc2 += pow(signal[bin].pileup - 1., 2);
    unc2 += pow(signal[bin].topPt - 1., 2);
    unc2 += pow(signal[bin].jec - 1., 2);
    unc2 += pow(signal[bin].leptonID - 1., 2);
    unc2 += pow(signal[bin].photonID - 1., 2);
    unc2 += pow(signal[bin].xsecError - 1., 2);

    for(unsigned int j = 0; j < ttbar.size(); j++) {
	unc2 += pow(ttbar[j].stat - 1., 2);
	unc2 += pow(ttbar[j].btag - 1., 2);
	unc2 += pow(ttbar[j].pileup - 1., 2);
	unc2 += pow(ttbar[j].topPt - 1., 2);
	unc2 += pow(ttbar[j].jec - 1., 2);
	unc2 += pow(ttbar[j].leptonID - 1., 2);
	unc2 += pow(ttbar[j].photonID - 1., 2);
	unc2 += pow(ttbar[j].scale - 1., 2);
	unc2 += pow(ttbar[j].pdf - 1., 2);
    }

    for(unsigned int j = 0; j < ttgamma.size(); j++) {
	unc2 += pow(ttgamma[j].stat - 1., 2);
	unc2 += pow(ttgamma[j].btag - 1., 2);
	unc2 += pow(ttgamma[j].pileup - 1., 2);
	unc2 += pow(ttgamma[j].topPt - 1., 2);
	unc2 += pow(ttgamma[j].jec - 1., 2);
	unc2 += pow(ttgamma[j].leptonID - 1., 2);
	unc2 += pow(ttgamma[j].photonID - 1., 2);
	unc2 += pow(ttgamma[j].scale - 1., 2);
	unc2 += pow(ttgamma[j].pdf - 1., 2);
    }

    double R = 2. * sqrt(unc2) / signal[bin].value;
    if(R < Rmin) Rmin = R;
  }

  outfile << "# R_firstguess = " << Rmin << endl;

  outfile << endl << "imax " << sensitive_bins.size() << " number of channels" << endl;
  outfile << "jmax 2 number of backgrounds" << endl;
  outfile << "kmax 23 number of nuisance parameters" << endl;
  outfile << "--------------------" << endl;

  outfile << "bin                 ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) outfile << "\t" << data[sensitive_bins[i]].binName.Data();
  outfile << endl;

  outfile << "observation         ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) outfile << "\t" << data[sensitive_bins[i]].value;
  outfile << endl;

  outfile << "--------------------" << endl;
  outfile << "bin                 ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) outfile << "\t" << signal[sensitive_bins[i]].binName.Data() << "\t" << ttbar[sensitive_bins[i]].binName.Data() << "\t" << ttgamma[sensitive_bins[i]].binName.Data();
  outfile << endl;

  outfile << "process             ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) outfile << "\tsusy ttbar ttgamma";
  outfile << endl;

  outfile << "process             ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) outfile << "\t 0 1 2";
  outfile << endl;

  outfile << "rate                ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) outfile << "\t" << signal[sensitive_bins[i]].value << " " << ttbar[sensitive_bins[i]].value << " " << ttgamma[sensitive_bins[i]].value;
  outfile << endl;
  outfile << "--------------------" << endl;

  outfile << "u_lumi lnN          ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) outfile << "\t" << lumi_sysError << " " << lumi_sysError << " " << lumi_sysError;
  outfile << endl;
  
  outfile << "u_btag lnN          ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) outfile << "\t" << signal[sensitive_bins[i]].btag << " " << ttbar[sensitive_bins[i]].btag << " " << ttgamma[sensitive_bins[i]].btag;
  outfile << endl;

  outfile << "u_pileup lnN        ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) outfile << "\t" << signal[sensitive_bins[i]].pileup << " " << ttbar[sensitive_bins[i]].pileup << " " << ttgamma[sensitive_bins[i]].pileup;
  outfile << endl;

  outfile << "u_susy_xsec lnN     ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) outfile << "\t" << signal[sensitive_bins[i]].xsecError << " - -";
  outfile << endl;

  outfile << "u_ttbar_scale lnN   ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) outfile << "\t- " << ttbar[sensitive_bins[i]].scale << " -";
  outfile << endl;

  outfile << "u_ttgamma_scale lnN ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) outfile << "\t- - " << ttgamma[sensitive_bins[i]].scale;
  outfile << endl;

  outfile << "u_ttbar_pdf lnN    ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) outfile << "\t- " << ttbar[sensitive_bins[i]].pdf << " -";
  outfile << endl;

  outfile << "u_ttgamma_pdf lnN  ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) outfile << "\t- - " << ttgamma[sensitive_bins[i]].pdf;
  outfile << endl;

  outfile << "u_topPt lnN         ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) outfile << "\t- " << ttbar[sensitive_bins[i]].topPt << " " << ttgamma[sensitive_bins[i]].topPt;
  outfile << endl;

  outfile << "u_JEC lnN           ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) outfile << "\t" << signal[sensitive_bins[i]].jec << " " << ttbar[sensitive_bins[i]].jec << " " << ttgamma[sensitive_bins[i]].jec;
  outfile << endl;

  outfile << "u_ele_id lnN        ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) {
    if(signal[sensitive_bins[i]].binName.Contains("ele")) outfile << "\t" << signal[sensitive_bins[i]].leptonID << " " << ttbar[sensitive_bins[i]].leptonID << " " << ttgamma[sensitive_bins[i]].leptonID;
    else outfile << "\t- - -";
  }
  outfile << endl;

  outfile << "u_muon_id lnN       ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) {
    if(signal[sensitive_bins[i]].binName.Contains("muon")) outfile << "\t" << signal[sensitive_bins[i]].leptonID << " " << ttbar[sensitive_bins[i]].leptonID << " " << ttgamma[sensitive_bins[i]].leptonID;
    else outfile << "\t- - -";
  }
  outfile << endl;

  outfile << "u_photon_id lnN     ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) outfile << "\t" << signal[sensitive_bins[i]].photonID << " " << ttbar[sensitive_bins[i]].photonID << " " << ttgamma[sensitive_bins[i]].photonID;
  outfile << endl;

  outfile << "stat_ele_susy lnN   ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) {
    if(signal[sensitive_bins[i]].binName.Contains("ele")) outfile << "\t" << signal[sensitive_bins[i]].stat << " - -";
    else outfile << "\t- - -";
  }
  outfile << endl;

  outfile << "stat_muon_susy lnN  ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) {
    if(signal[sensitive_bins[i]].binName.Contains("muon")) outfile << "\t" << signal[sensitive_bins[i]].stat << " - -";
    else outfile << "\t- - -";
  }
  outfile << endl;

  outfile << "stat_ele_ttbar lnN  ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) {
    if(signal[sensitive_bins[i]].binName.Contains("ele")) outfile << "\t- " << ttbar[sensitive_bins[i]].stat << " -";
    else outfile << "\t- - -";
  }
  outfile << endl;

  outfile << "stat_muon_ttbar lnN ";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) {
    if(signal[sensitive_bins[i]].binName.Contains("muon")) outfile << "\t- " << ttbar[sensitive_bins[i]].stat << " -";
    else outfile << "\t- - -";
  }
  outfile << endl;

  outfile << "stat_ele_ttgamma lnN";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) {
    if(signal[sensitive_bins[i]].binName.Contains("ele")) outfile << "\t- - " << ttgamma[sensitive_bins[i]].stat;
    else outfile << "\t- - -";
  }
  outfile << endl;

  outfile << "stat_muon_ttgamma lnN";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) {
    if(signal[sensitive_bins[i]].binName.Contains("muon")) outfile << "\t- - " << ttgamma[sensitive_bins[i]].stat;
    else outfile << "\t- - -";
  }
  outfile << endl;

  outfile << "u_ttjets_fit_ele lnN";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) {
    if(signal[sensitive_bins[i]].binName.Contains("ele")) outfile << "\t- 1.0665 -";
    else outfile << "\t- - -";
  }
  outfile << endl;

  outfile << "u_ttjets_fit_muon lnN";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) {
    if(signal[sensitive_bins[i]].binName.Contains("muon")) outfile << "\t- 1.1411 -";
    else outfile << "\t- - -";
  }
  outfile << endl;

  outfile << "u_ttgamma_fit_ele lnN";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) {
    if(signal[sensitive_bins[i]].binName.Contains("ele")) outfile << "\t- - 1.0814";
    else outfile << "\t- - -";
  }
  outfile << endl;

  outfile << "u_ttgamma_fit_muon lnN";
  for(unsigned int i = 0; i < sensitive_bins.size(); i++) {
    if(signal[sensitive_bins[i]].binName.Contains("muon")) outfile << "\t- 1.1808 -";
    else outfile << "\t- - -";
  }
  outfile << endl;


}

void GetBackgroundHistograms(TFile * f, TString var, TString req,
			     vector<TH1D*>& h,
			     vector<TH1D*>& btagWeightUp, vector<TH1D*>& btagWeightDown,
			     vector<TH1D*>& puWeightUp, vector<TH1D*>& puWeightDown,
			     vector<TH1D*>& topPtUp, vector<TH1D*>& topPtDown,
			     vector<TH1D*>& JECup, vector<TH1D*>& JECdown,
			     vector<TH1D*>& leptonSFup, vector<TH1D*>& leptonSFdown,
			     vector<TH1D*>& photonSFup, vector<TH1D*>& photonSFdown,
			     vector<TH1D*>& scaleUp, vector<TH1D*>& scaleDown,
			     vector<TH1D*>& pdfUp, vector<TH1D*>& pdfDown) {

  for(int i = 0; i < nSamples; i++) {
    h.push_back((TH1D*)f->Get(var+"_"+mcNames[i]+"_"+req));

    btagWeightUp.push_back((TH1D*)f->Get(var+"_"+mcNames[i]+"_"+req+"_btagWeightUp"));
    btagWeightDown.push_back((TH1D*)f->Get(var+"_"+mcNames[i]+"_"+req+"_btagWeightDown"));

    puWeightUp.push_back((TH1D*)f->Get(var+"_"+mcNames[i]+"_"+req+"_puWeightUp"));
    puWeightDown.push_back((TH1D*)f->Get(var+"_"+mcNames[i]+"_"+req+"_puWeightDown"));

    topPtUp.push_back((TH1D*)f->Get(var+"_"+mcNames[i]+"_"+req+"_topPtUp"));
    topPtDown.push_back((TH1D*)f->Get(var+"_"+mcNames[i]+"_"+req+"_topPtDown"));

    JECup.push_back((TH1D*)f->Get(var+"_"+mcNames[i]+"_"+req+"_JECup"));
    JECdown.push_back((TH1D*)f->Get(var+"_"+mcNames[i]+"_"+req+"_JECdown"));

    leptonSFup.push_back((TH1D*)f->Get(var+"_"+mcNames[i]+"_"+req+"_leptonSFup"));
    leptonSFdown.push_back((TH1D*)f->Get(var+"_"+mcNames[i]+"_"+req+"_leptonSFdown"));

    photonSFup.push_back((TH1D*)f->Get(var+"_"+mcNames[i]+"_"+req+"_photonSFup"));
    photonSFdown.push_back((TH1D*)f->Get(var+"_"+mcNames[i]+"_"+req+"_photonSFdown"));

    scaleUp.push_back((TH1D*)f->Get(var+"_"+mcNames[i]+"_"+req+"_scaleUp"));
    scaleDown.push_back((TH1D*)f->Get(var+"_"+mcNames[i]+"_"+req+"_scaleDown"));

    pdfUp.push_back((TH1D*)f->Get(var+"_"+mcNames[i]+"_"+req+"_pdfUp"));
    pdfDown.push_back((TH1D*)f->Get(var+"_"+mcNames[i]+"_"+req+"_pdfDown"));
  }

}

bool GetSignalHistograms(TFile * f, TString var, TString req,
			 int mStop, int mBino,
			 TH1D*& h,
			 TH1D*& h_btagWeightUp, TH1D*& h_btagWeightDown,
			 TH1D*& h_puWeightUp, TH1D*& h_puWeightDown,
			 TH1D*& h_topPtUp, TH1D*& h_topPtDown,
			 TH1D*& h_JECup, TH1D*& h_JECdown,
			 TH1D*& h_leptonSFup, TH1D*& h_leptonSFdown,
			 TH1D*& h_photonSFup, TH1D*& h_photonSFdown) {

  stringstream code;
  code << "_mst_" << mStop << "_m1_" << mBino;

  h = (TH1D*)f->Get(var+"_gg_"+req+code.str().c_str());

  h_btagWeightUp = (TH1D*)f->Get(var+"_gg_"+req+code.str().c_str()+"_btagWeightUp");
  h_btagWeightDown = (TH1D*)f->Get(var+"_gg_"+req+code.str().c_str()+"_btagWeightDown");

  h_puWeightUp = (TH1D*)f->Get(var+"_gg_"+req+code.str().c_str()+"_puWeightUp");
  h_puWeightDown = (TH1D*)f->Get(var+"_gg_"+req+code.str().c_str()+"_puWeightDown");

  h_topPtUp = (TH1D*)f->Get(var+"_gg_"+req+code.str().c_str()+"_topPtUp");
  h_topPtDown = (TH1D*)f->Get(var+"_gg_"+req+code.str().c_str()+"_topPtDown");

  h_JECup = (TH1D*)f->Get(var+"_gg_"+req+code.str().c_str()+"_JECup");
  h_JECdown = (TH1D*)f->Get(var+"_gg_"+req+code.str().c_str()+"_JECdown");

  h_leptonSFup = (TH1D*)f->Get(var+"_gg_"+req+code.str().c_str()+"_leptonSFup");
  h_leptonSFdown = (TH1D*)f->Get(var+"_gg_"+req+code.str().c_str()+"_leptonSFdown");

  h_photonSFup = (TH1D*)f->Get(var+"_gg_"+req+code.str().c_str()+"_photonSFup");
  h_photonSFdown = (TH1D*)f->Get(var+"_gg_"+req+code.str().c_str()+"_photonSFdown");

  if(!h ||
     !h_btagWeightUp || !h_btagWeightDown ||
     !h_puWeightUp || !h_puWeightDown ||
     !h_topPtUp || !h_topPtDown ||
     !h_JECup || !h_JECdown ||
     !h_leptonSFup || !h_leptonSFdown ||
     !h_photonSFup || !h_photonSFdown) return false;

  return true;
}


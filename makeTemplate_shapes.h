#include "TROOT.h"
#include <TStyle.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TCanvas.h>

#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <vector>
#include <cmath>

using namespace std;

const int nBackgrounds = 6;
const TString backgroundNames[nBackgrounds] = {"ttjets", "wjets", "zjets", "zz", "ttZ", "ttgamma"};
const double scaleUp[nBackgrounds] = {1.025, 1.0065, 1.005, 1.036, 1.092, 1.25};
const double scaleDown[nBackgrounds] = {1.034, 1.0032, 1.0031, 1.036, 1.117, 1.25};
const double pdfUp[nBackgrounds] = {1.026, 1.034, 1.033, 1.036, -1., 1.076};
const double pdfDown[nBackgrounds] = {1.026, 1.034, 1.033, 1.036, -1., 1.099};

const int nSystematics = 6;
const TString systematicNames[nSystematics] = {"btagWeight", "puWeight", "topPt", "JEC", "leptonSF", "photonSF"};


const double epsilon = 1e-10;

const TString datacard_dir = "datacards";

class GridPoint {

 public:
  GridPoint();
  virtual ~GridPoint() {
    backgroundYields_ele.clear();
    backgroundYields_muon.clear();
  }

  void Print();
  bool SetBackgroundYields(TFile * f);
  bool SetSignalYields(TFile * f);

  int mStop;
  int mBino;

  double xsec;
  double xsecError;

  double lumi_sysError;

  vector<double> backgroundYields_ele;
  vector<double> backgroundYields_muon;

  double signalYield_ele;
  double signalYield_muon;

  double obs_ele;
  double obs_muon;

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

  backgroundYields_ele.clear();
  backgroundYields_muon.clear();

  signalYield_ele = signalYield_muon = 0;

  obs_ele = obs_muon = 0;

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

  vector<TString> sensitive_channels;
  if(signalYield_ele > epsilon) sensitive_channels.push_back("ele");
  if(signalYield_muon > epsilon) sensitive_channels.push_back("muon");

  int nchan = 0;
  if(signalYield_ele > epsilon) nchan++;
  if(signalYield_muon > epsilon) nchan++;

  outfile << endl << "imax " << nchan << " number of channels" << endl;
  outfile << endl << "jmax " << nBackgrounds << " number of backgrounds" << endl;
  outfile << "kmax * number of nuisance parameters" << endl;
  outfile << "--------------------" << endl;
  outfile << "shapes * * limitInputs.root $CHANNEL/$PROCESS  $CHANNEL/$PROCESS_$SYSTEMATIC" << endl;
  outfile << "--------------------" << endl;

  outfile << "bin                 ";
  if(signalYield_ele > epsilon) outfile << "\tele";
  if(signalYield_muon > epsilon) outfile << "\tmuon";
  outfile << endl;

  outfile << "observation         ";
  if(signalYield_ele > epsilon) outfile << "\t" << obs_ele;
  if(signalYield_muon > epsilon) outfile << "\t" << obs_muon;
  outfile << endl;

  outfile << "--------------------" << endl;
  outfile << "bin                 ";
  if(signalYield_ele > epsilon) {
    for(int i = 0; i < nBackgrounds + 1; i++) outfile << "\tele";
  }
  if(signalYield_muon > epsilon) {
    for(int i = 0; i < nBackgrounds + 1; i++) outfile << "\tmuon";
  }
  outfile << endl;

  stringstream code;
  code << "_mst_" << mStop << "_m1_" << mBino;
    
  outfile << "process             ";
  if(signalYield_ele > epsilon) {
    outfile << "\tsignal" << code.str();
    for(int i = 0; i < nBackgrounds; i++) outfile << "\t" << backgroundNames[i].Data();
  }
  if(signalYield_muon > epsilon) {
    outfile << "\tsignal" << code.str();
    for(int i = 0; i < nBackgrounds; i++) outfile << "\t" << backgroundNames[i].Data();
  }
  outfile << endl;

  outfile << "process             ";
  if(signalYield_ele > epsilon) {
    for(int i = 0; i < nBackgrounds + 1; i++) outfile << "\t" << i;
  }
  if(signalYield_muon > epsilon) {
    for(int i = 0; i < nBackgrounds + 1; i++) outfile << "\t" << i;
  }
  outfile << endl;

  outfile << "rate                ";
  if(signalYield_ele > epsilon) {
    outfile << "\t" << signalYield_ele;
    for(int i = 0; i < nBackgrounds; i++) outfile << "\t" << backgroundYields_ele[i];
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t" << signalYield_muon;
    for(int i = 0; i < nBackgrounds; i++) outfile << "\t" << backgroundYields_muon[i];
  }
  outfile << endl;
  outfile << "--------------------" << endl;

  outfile << "lumi lnN          ";
  if(signalYield_ele > epsilon) {
    for(int i = 0; i < nBackgrounds + 1; i++) outfile << "\t" << lumi_sysError;
  }
  if(signalYield_muon > epsilon) {
    for(int i = 0; i < nBackgrounds + 1; i++) outfile << "\t" << lumi_sysError;
  }
  outfile << endl;

  for(int i = 0; i < nSystematics; i++) {
    outfile << systematicNames[i].Data() << " shape          ";
    if(signalYield_ele > epsilon) {
      for(int j = 0; j < nBackgrounds + 1; j++) outfile << "\t1.0";
    }
    if(signalYield_muon > epsilon) {
      for(int j = 0; j < nBackgrounds + 1; j++) outfile << "\t1.0";
    }
    outfile << endl;
  }

  for(int j = 0; j < nBackgrounds; j++) {

    outfile << "scale_" << backgroundNames[j].Data() << " lnN          ";
    if(signalYield_ele > epsilon) {
      outfile << "\t-";
      for(int k = 0; k < nBackgrounds; k++) {
	if(k == j) {
	  if(scaleUp[k] == scaleDown[k]) outfile << "\t" << scaleUp[k];
	  else outfile << "\t" << 2. - scaleDown[k] << "/" << scaleUp[k];
	}
	else outfile << "\t-";
      }
    }
    if(signalYield_muon > epsilon) {
      outfile << "\t-";
      for(int k = 0; k < nBackgrounds; k++) {
	if(k == j) {
	  if(scaleUp[k] == scaleDown[k]) outfile << "\t" << scaleUp[k];
	  else outfile << "\t" << 2. - scaleDown[k] << "/" << scaleUp[k];
	}
	else outfile << "\t-";
      }
    }
    outfile << endl;

  }

  for(int j = 0; j < nBackgrounds; j++) {

    if(pdfUp[j] < 0. || pdfDown[j] < 0.) continue;

    outfile << "pdf_" << backgroundNames[j].Data() << " lnN          ";
    if(signalYield_ele > epsilon) {
      outfile << "\t-";
      for(int k = 0; k < nBackgrounds; k++) {
	if(k == j) {
	  if(pdfUp[k] == pdfDown[k]) outfile << "\t" << pdfUp[k];
	  else outfile << "\t" << 2. - pdfDown[k] << "/" << pdfUp[k];
	}
	else outfile << "\t-";
      }
    }
    if(signalYield_muon > epsilon) {
      outfile << "\t-";
      for(int k = 0; k < nBackgrounds; k++) {
	if(k == j) {
	  if(pdfUp[k] == pdfDown[k]) outfile << "\t" << pdfUp[k];
	  else outfile << "\t" << 2. - pdfDown[k] << "/" << pdfUp[k];
	}
	else outfile << "\t-";
      }
    }
    outfile << endl;

  }

  outfile << "susy_xsec lnN     ";
  if(signalYield_ele > epsilon) {
    outfile << "\t" << 1. + xsecError/100.;
    for(int i = 0; i < nBackgrounds; i++) outfile << "\t-";
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t" << 1. + xsecError/100.;
    for(int i = 0; i < nBackgrounds; i++) outfile << "\t-";
  }
  outfile << endl;

  outfile << "ttbar_float lnN";
  if(signalYield_ele > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(i == 0 || i == 5) outfile << "\treplaceme";
      else outfile << "\t-";
    }
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(i == 0 || i == 5) outfile << "replaceme";
      else outfile << "\t-";
    }
  }

}

bool GridPoint::SetBackgroundYields(TFile * f) {

  TH1D * h = (TH1D*)f->Get("ele/data_obs");
  if(!h) return false;
  obs_ele = h->Integral();

  h = (TH1D*)f->Get("muon/data_obs");
  if(!h) return false;
  obs_muon = h->Integral();

  for(int i = 0; i < nBackgrounds; i++) {
    h = (TH1D*)f->Get("ele/"+backgroundNames[i]);
    if(!h) return false;
    backgroundYields_ele.push_back(h->Integral());

    h = (TH1D*)f->Get("muon/"+backgroundNames[i]);
    if(!h) return false;
    backgroundYields_muon.push_back(h->Integral());
  }

  return true;

}

bool GridPoint::SetSignalYields(TFile * f) {
  
  stringstream code;
  code << "_mst_" << mStop << "_m1_" << mBino;
  TString code_t = code.str();

  TH1D * h = (TH1D*)f->Get("ele/signal"+code_t);
  if(!h) return false;
  signalYield_ele = h->Integral();

  h = (TH1D*)f->Get("muon/signal"+code_t);
  if(!h) return false;
  signalYield_muon = h->Integral();

  return true;
}

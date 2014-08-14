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

const int nBackgrounds = 10;
const TString backgroundNames[nBackgrounds] = {"ttjets", "wjets", "zjets", "singleTop", "ww", "wz", "zz", "ttW", "ttZ", "ttgamma"};

// don't forget scale and pdf
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

  double qcdYield_ele;
  double qcdYield_muon;

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
  outfile << endl << "jmax " << nBackgrounds + 1 << " number of backgrounds" << endl;
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
    outfile << "\tele\tele";
    for(int i = 0; i < nBackgrounds; i++) outfile << "\tele";
  }
  if(signalYield_muon > epsilon) {
    outfile << "\tmuon\tmuon";
    for(int i = 0; i < nBackgrounds; i++) outfile << "\tmuon";
  }
  outfile << endl;

  stringstream code;
  code << "_mst_" << mStop << "_m1_" << mBino;
    
  outfile << "process             ";
  if(signalYield_ele > epsilon) {
    outfile << "\tsignal" << code.str() << "\tqcd";
    for(int i = 0; i < nBackgrounds; i++) outfile << "\t" << backgroundNames[i].Data();
  }
  if(signalYield_muon > epsilon) {
    outfile << "\tsignal" << code.str() << "\tqcd";
    for(int i = 0; i < nBackgrounds; i++) outfile << "\t" << backgroundNames[i].Data();
  }
  outfile << endl;

  outfile << "process             ";
  if(signalYield_ele > epsilon) {
    for(int i = 0; i < nBackgrounds + 2; i++) outfile << "\t" << i;
  }
  if(signalYield_muon > epsilon) {
    for(int i = 0; i < nBackgrounds + 2; i++) outfile << "\t" << i;
  }
  outfile << endl;

  outfile << "rate                ";
  if(signalYield_ele > epsilon) {
    outfile << "\t" << signalYield_ele;
    outfile << "\t" << qcdYield_ele;
    for(int i = 0; i < nBackgrounds; i++) outfile << "\t" << backgroundYields_ele[i];
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t" << signalYield_muon;
    outfile << "\t" << qcdYield_muon;
    for(int i = 0; i < nBackgrounds; i++) outfile << "\t" << backgroundYields_muon[i];
  }
  outfile << endl;
  outfile << "--------------------" << endl;

  outfile << "lumi lnN          ";
  if(signalYield_ele > epsilon) {
    outfile << "\t" << lumi_sysError;
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) outfile << "\t" << lumi_sysError;
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t" << lumi_sysError;
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) outfile << "\t" << lumi_sysError;
  }
  outfile << endl;

  for(int i = 0; i < nSystematics; i++) {
    outfile << systematicNames[i].Data() << " shapeN2          ";
    if(signalYield_ele > epsilon) {
      for(int j = 0; j < nBackgrounds + 2; j++) {
	if(j == 1) outfile << "\t-";
	else outfile << "\t1.0";
      }
    }
    if(signalYield_muon > epsilon) {
      for(int j = 0; j < nBackgrounds + 2; j++) {
	if(j == 1) outfile << "\t-";
	else outfile << "\t1.0";
      }
    }
    outfile << endl;
  }

  for(int j = 0; j < nBackgrounds; j++) {

    outfile << "scale_" << backgroundNames[j].Data() << " shapeN2          ";
    if(signalYield_ele > epsilon) {
      outfile << "\t-\t-";
      for(int k = 0; k < nBackgrounds; k++) {
	if(k == j) outfile << "\t1.0";
	else outfile << "\t-";
      }
    }
    if(signalYield_muon > epsilon) {
      outfile << "\t-\t-";
      for(int k = 0; k < nBackgrounds; k++) {
	if(k == j) outfile << "\t1.0";
	else outfile << "\t-";
      }
    }
    outfile << endl;

  }

  for(int j = 0; j < nBackgrounds; j++) {

    outfile << "pdf_" << backgroundNames[j].Data() << " shapeN2          ";
    if(signalYield_ele > epsilon) {
      outfile << "\t-\t-";
      for(int k = 0; k < nBackgrounds; k++) {
	if(k == j) outfile << "\t1.0";
	else outfile << "\t-";
      }
    }
    if(signalYield_muon > epsilon) {
      outfile << "\t-\t-";
      for(int k = 0; k < nBackgrounds; k++) {
	if(k == j) outfile << "\t1.0";
	else outfile << "\t-";
      }
    }
    outfile << endl;

  }

  outfile << "susy_xsec lnN     ";
  if(signalYield_ele > epsilon) {
    outfile << "\t" << 1. + xsecError/xsec;
    for(int i = 0; i < nBackgrounds + 1; i++) outfile << "\t-";
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t" << 1. + xsecError/xsec;
    for(int i = 0; i < nBackgrounds + 1; i++) outfile << "\t-";
  }
  outfile << endl;

  outfile << "ttjets_fit_ele lnN";
  if(signalYield_ele > epsilon) {
    outfile << "\t-\t-\t1.0665";
    for(int i = 1; i < nBackgrounds; i++) outfile << "\t-";
  }
  if(signalYield_muon > epsilon) {
    for(int i = 0; i < nBackgrounds + 2; i++) outfile << "\t-";
  }
  outfile << endl;

  outfile << "ttjets_fit_muon lnN";
  if(signalYield_ele > epsilon) {
    for(int i = 0; i < nBackgrounds + 2; i++) outfile << "\t-";
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t-\t-\t1.1411";
    for(int i = 1; i < nBackgrounds; i++) outfile << "\t-";
  }
  outfile << endl;

  outfile << "ttgamma_fit_ele lnN";
  if(signalYield_ele > epsilon) {
    for(int i = 0; i < nBackgrounds + 1; i++) outfile << "\t-";
    outfile << "\t1.0814";
  }
  if(signalYield_muon > epsilon) {
    for(int i = 0; i < nBackgrounds + 2; i++) outfile << "\t-";
  }
  outfile << endl;

  outfile << "ttgamma_fit_muon lnN";
  if(signalYield_ele > epsilon) {
    for(int i = 0; i < nBackgrounds + 2; i++) outfile << "\t-";
  }
  if(signalYield_muon > epsilon) {
    for(int i = 0; i < nBackgrounds + 1; i++) outfile << "\t-";
    outfile << "\t1.1808";
  }
  outfile << endl;

}

bool GridPoint::SetBackgroundYields(TFile * f) {

  TH1D * h = (TH1D*)f->Get("ele/data_obs");
  if(!h) return false;
  obs_ele = h->Integral();

  h = (TH1D*)f->Get("muon/data_obs");
  if(!h) return false;
  obs_muon = h->Integral();

  h = (TH1D*)f->Get("ele/qcd");
  if(!h) return false;
  qcdYield_ele = h->Integral();
  
  h = (TH1D*)f->Get("muon/qcd");
  if(!h) return false;
  qcdYield_muon = h->Integral();

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

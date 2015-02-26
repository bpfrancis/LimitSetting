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

const int nBackgrounds = 8;
const TString backgroundNames[nBackgrounds] = {"ttjets", "wjets", "zjets", "singleTop", "diboson", "ttW", "ttZ", "ttgamma"};

const bool scale_tt[nBackgrounds] = {true, false, false, true, false, true, true, true};
const bool scale_V[nBackgrounds]  = {false, true, true, false, false, false, false, false};
const bool scale_VV[nBackgrounds] = {false, false, false, false, true, false, false, false};

const bool pdf_gg[nBackgrounds]   = {true, false, false, false, false, false, false/*no pdf given for ttZ, but it would be here*/, true};
const bool pdf_qq[nBackgrounds]   = {false, true, true, false, true, true, false, false};
const bool pdf_qg[nBackgrounds]   = {false, false, false, true, false, false, false, false};

const int nSystematics = 7;
const TString systematicNames[nSystematics] = {"btagWeight", "puWeight", "topPt", "JEC", "leptonSF", "photonSF", "extraSystematic"};

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

  outfile << "scale_tt shape ";
  if(signalYield_ele > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(scale_tt[i]) outfile << "\t1.0";
      else outfile << "\t-";
    }
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(scale_tt[i]) outfile << "\t1.0";
      else outfile << "\t-";
    }
  }
  outfile << endl;
  
  outfile << "scale_V shape ";
  if(signalYield_ele > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(scale_V[i]) outfile << "\t1.0";
      else outfile << "\t-";
    }
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(scale_V[i]) outfile << "\t1.0";
      else outfile << "\t-";
    }
  }
  outfile << endl;

  outfile << "scale_VV shape ";
  if(signalYield_ele > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(scale_VV[i]) outfile << "\t1.0";
      else outfile << "\t-";
    }
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(scale_VV[i]) outfile << "\t1.0";
      else outfile << "\t-";
    }
  }
  outfile << endl;

  outfile << "pdf_gg shape ";
  if(signalYield_ele > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(pdf_gg[i]) outfile << "\t1.0";
      else outfile << "\t-";
    }
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(pdf_gg[i]) outfile << "\t1.0";
      else outfile << "\t-";
    }
  }
  outfile << endl;

  outfile << "pdf_qq shape ";
  if(signalYield_ele > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(pdf_qq[i]) outfile << "\t1.0";
      else outfile << "\t-";
    }
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(pdf_qq[i]) outfile << "\t1.0";
      else outfile << "\t-";
    }
  }
  outfile << endl;

  outfile << "pdf_qg shape ";
  if(signalYield_ele > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(pdf_qg[i]) outfile << "\t1.0";
      else outfile << "\t-";
    }
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(pdf_qg[i]) outfile << "\t1.0";
      else outfile << "\t-";
    }
  }
  outfile << endl;

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

  /* no fit or float anymore
  outfile << "ttjets_fit_ele lnN";
  if(signalYield_ele > epsilon) {
    outfile << "\t-\t-\t1.011";
    for(unsigned int i = 1; i < nBackgrounds; i++) outfile << "\t-";
  }
  if(signalYield_muon > epsilon) {
    for(unsigned int i = 0; i < nBackgrounds + 2; i++) outfile << "\t-";
  }
  outfile << endl;
  
  outfile << "ttjets_fit_muon lnN";
  if(signalYield_ele > epsilon) {
    for(unsigned int i = 0; i < nBackgrounds + 2; i++) outfile << "\t-";
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t-\t-\t1.2";
    for(unsigned int i = 1; i < nBackgrounds; i++) outfile << "\t-";
  }
  outfile << endl;

  outfile << "ttgamma_fit_ele lnN";
  if(signalYield_ele > epsilon) {
    for(unsigned int i = 0; i < nBackgrounds + 1; i++) outfile << "\t-";
    outfile << "\t1.041";
  }
  if(signalYield_muon > epsilon) {
    for(unsigned int i = 0; i < nBackgrounds + 2; i++) outfile << "\t-";
  }
  outfile << endl;

  outfile << "ttgamma_fit_muon lnN";
  if(signalYield_ele > epsilon) {
    for(unsigned int i = 0; i < nBackgrounds + 2; i++) outfile << "\t-";
  }
  if(signalYield_muon > epsilon) {
    for(unsigned int i = 0; i < nBackgrounds + 1; i++) outfile << "\t-";
    outfile << "\t1.19";
  }
  outfile << endl;

  outfile << "ttbar_float lnN";
  if(signalYield_ele > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(i == 0 || i == 5) outfile << "\t2.0";
      else outfile << "\t-";
    }
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(i == 0 || i == 5) outfile << "\t2.0";
      else outfile << "\t-";
    }
  }
  outfile << endl;
  */
  //durp
  for(int ibin = 1; ibin <= 6; ibin++) {

    for(int i = 0; i < nBackgrounds; i++) {
      outfile << backgroundNames[i] << "_SR2_stat_bin" << ibin << " shape";
      
      if(signalYield_ele > epsilon) {
	outfile << "\t-";
	for(int j = 0; j < nBackgrounds; j++) {
	  if(j == i) outfile << "\t1.0";
	  else outfile << "\t-";
	}
      }
      if(signalYield_muon > epsilon) {
	outfile << "\t-";
	for(int j = 0; j < nBackgrounds; j++) {
	  if(j == i) outfile << "\t1.0";
	  else outfile << "\t-";
	}
      }

      outfile << endl;

    }

  } // stats block
      

} // Print()

bool GridPoint::SetBackgroundYields(TFile * f) {

  TH1D * h = (TH1D*)f->Get("ele_SR2/data_obs");
  if(!h) return false;
  obs_ele = h->Integral();

  h = (TH1D*)f->Get("muon_SR2/data_obs");
  if(!h) return false;
  obs_muon = h->Integral();

  for(int i = 0; i < nBackgrounds; i++) {
    h = (TH1D*)f->Get("ele_SR2/"+backgroundNames[i]);
    if(!h) return false;
    backgroundYields_ele.push_back(h->Integral());

    h = (TH1D*)f->Get("muon_SR2/"+backgroundNames[i]);
    if(!h) return false;
    backgroundYields_muon.push_back(h->Integral());
  }

  return true;

}

bool GridPoint::SetSignalYields(TFile * f) {
  
  stringstream code;
  code << "_mst_" << mStop << "_m1_" << mBino;
  TString code_t = code.str();

  TH1D * h = (TH1D*)f->Get("ele_SR2/signal"+code_t);
  if(!h) return false;
  signalYield_ele = h->Integral();

  h = (TH1D*)f->Get("muon_SR2/signal"+code_t);
  if(!h) return false;
  signalYield_muon = h->Integral();

  return true;
}

void GridPoint::SetUseStatError(TFile * f) {

  durp;

}

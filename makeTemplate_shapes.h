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

const double scaleUp_tt[nBackgrounds] = {1.025, -1., -1., -1., 1.092, 1.25};
const double scaleDown_tt[nBackgrounds] = {1.034, -1., -1., -1., 1.117, 1.25};
const double pdfUp_gg[nBackgrounds] = {1.026, -1., -1., -1., -1., 1.076};
const double pdfDown_gg[nBackgrounds] = {1.026, -1., -1., -1., -1., 1.099};

const double scaleUp_V[nBackgrounds] = {-1., 1.0065, 1.005, -1., -1., -1.};
const double scaleDown_V[nBackgrounds] = {-1., 1.0032, 1.0031, -1., -1., -1.};
const double pdfUp_qq[nBackgrounds] = {-1., 1.034, 1.033, -1., -1., -1.};
const double pdfDown_qq[nBackgrounds] = {-1., 1.034, 1.033, -1., -1., -1.};

const double scaleUp_VV[nBackgrounds] = {-1., -1., -1., 1.036, -1., -1.};
const double scaleDown_VV[nBackgrounds] = {-1., -1., -1., 1.036, -1., -1.};
const double pdfUp_qg[nBackgrounds] = {-1., -1., -1., 1.036, -1., -1.};
const double pdfDown_qg[nBackgrounds] = {-1., -1., -1., 1.036, -1., -1.};

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

  outfile << "scale_tt lnN ";
  if(signalYield_ele > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(scaleUp_tt[i] < 0) outfile << "\t-";
      else if(scaleUp_tt[i] == scaleDown_tt[i]) outfile << "\t" << scaleUp_tt[i];
      else outfile << "\t" << 2. - scaleDown_tt[i] << "/" << scaleUp_tt[i];
    }
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(scaleUp_tt[i] < 0) outfile << "\t-";
      else if(scaleUp_tt[i] == scaleDown_tt[i]) outfile << "\t" << scaleUp_tt[i];
      else outfile << "\t" << 2. - scaleDown_tt[i] << "/" << scaleUp_tt[i];
    }
  }
  outfile << endl;

  outfile << "scale_V lnN ";
  if(signalYield_ele > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(scaleUp_V[i] < 0) outfile << "\t-";
      else if(scaleUp_V[i] == scaleDown_V[i]) outfile << "\t" << scaleUp_V[i];
      else outfile << "\t" << 2. - scaleDown_V[i] << "/" << scaleUp_V[i];
    }
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(scaleUp_V[i] < 0) outfile << "\t-";
      else if(scaleUp_V[i] == scaleDown_V[i]) outfile << "\t" << scaleUp_V[i];
      else outfile << "\t" << 2. - scaleDown_V[i] << "/" << scaleUp_V[i];
    }
  }
  outfile << endl;

  outfile << "scale_VV lnN ";
  if(signalYield_ele > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(scaleUp_VV[i] < 0) outfile << "\t-";
      else if(scaleUp_VV[i] == scaleDown_VV[i]) outfile << "\t" << scaleUp_VV[i];
      else outfile << "\t" << 2. - scaleDown_VV[i] << "/" << scaleUp_VV[i];
    }
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(scaleUp_VV[i] < 0) outfile << "\t-";
      else if(scaleUp_VV[i] == scaleDown_VV[i]) outfile << "\t" << scaleUp_VV[i];
      else outfile << "\t" << 2. - scaleDown_VV[i] << "/" << scaleUp_VV[i];
    }
  }
  outfile << endl;

  outfile << "pdf_gg lnN ";
  if(signalYield_ele > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(pdfUp_gg[i] < 0) outfile << "\t-";
      else if(pdfUp_gg[i] == pdfDown_gg[i]) outfile << "\t" << pdfUp_gg[i];
      else outfile << "\t" << 2. - pdfDown_gg[i] << "/" << pdfUp_gg[i];
    }
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(pdfUp_gg[i] < 0) outfile << "\t-";
      else if(pdfUp_gg[i] == pdfDown_gg[i]) outfile << "\t" << pdfUp_gg[i];
      else outfile << "\t" << 2. - pdfDown_gg[i] << "/" << pdfUp_gg[i];
    }
  }
  outfile << endl;

  outfile << "pdf_qq lnN ";
  if(signalYield_ele > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(pdfUp_qq[i] < 0) outfile << "\t-";
      else if(pdfUp_qq[i] == pdfDown_qq[i]) outfile << "\t" << pdfUp_qq[i];
      else outfile << "\t" << 2. - pdfDown_qq[i] << "/" << pdfUp_qq[i];
    }
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(pdfUp_qq[i] < 0) outfile << "\t-";
      else if(pdfUp_qq[i] == pdfDown_qq[i]) outfile << "\t" << pdfUp_qq[i];
      else outfile << "\t" << 2. - pdfDown_qq[i] << "/" << pdfUp_qq[i];
    }
  }
  outfile << endl;

  outfile << "pdf_qg lnN ";
  if(signalYield_ele > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(pdfUp_qg[i] < 0) outfile << "\t-";
      else if(pdfUp_qg[i] == pdfDown_qg[i]) outfile << "\t" << pdfUp_qg[i];
      else outfile << "\t" << 2. - pdfDown_qg[i] << "/" << pdfUp_qg[i];
    }
  }
  if(signalYield_muon > epsilon) {
    outfile << "\t-";
    for(int i = 0; i < nBackgrounds; i++) {
      if(pdfUp_qg[i] < 0) outfile << "\t-";
      else if(pdfUp_qg[i] == pdfDown_qg[i]) outfile << "\t" << pdfUp_qg[i];
      else outfile << "\t" << 2. - pdfDown_qg[i] << "/" << pdfUp_qg[i];
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

  for(int ibin = 1; ibin <= 6; ibin++) {
    outfile << "signal_stat_bin" << ibin << " shape";
    if(signalYield_ele > epsilon) {
      outfile << "\t1.0";
      for(int i = 0; i < nBackgrounds; i++) outfile << "\t-";
    }
    if(signalYield_muon > epsilon) {
      outfile << "\t1.0";
      for(int i = 0; i < nBackgrounds; i++) outfile << "\t-";
    }
    outfile << endl;

    for(int i = 0; i < nBackgrounds; i++) {
      outfile << backgroundNames[i] << "_stat_bin" << ibin << " shape";
      
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

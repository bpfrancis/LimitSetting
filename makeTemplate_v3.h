#include "TROOT.h"
#include <TStyle.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TCanvas.h>
#include <TMath.h>

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

template <typename T>
class vector2d {
 public:
 vector2d(size_t _d1 = 0, size_t _d2 = 0, T& t=T()) :
  d1(_d1), d2(_d2), data(_d1*_d2, t)
    {}
  
  vector2d& operator=(vector2d& vec) {
    d1 = vec.sizeX();
    d2 = vec.sizeY();
    data = vec.GetData();
    return *this;
  }

  T & operator()(size_t i, size_t j) {
    return data[i*d2 + j];
  }
  
  void clear() { data.clear(); }
  size_t size() { return d1*d2; }
  size_t sizeX() { return d1; }
  size_t sizeY() { return d2; }
  vector<T> GetData() { return data; }
  
 private:
  size_t d1,d2;
  vector<T> data;
};

template <typename T>
class vector3d {
 public:
 vector3d(size_t _d1 = 0, size_t _d2 = 0, size_t _d3 = 0, T& t=T()) :
  d1(_d1), d2(_d2), d3(_d3), data(_d1*_d2*_d3, t)
    {}
  
  vector3d& operator=(vector3d& vec) {
    d1 = vec.sizeX();
    d2 = vec.sizeY();
    d3 = vec.sizeZ();
    data = vec.GetData();
    return *this;
  }

  T & operator()(size_t i, size_t j, size_t k) {
    return data[i*d2*d3 + j*d3 + k];
  }
  
  void clear() { data.clear(); }
  size_t size() { return d1*d2*d3; }
  size_t sizeX() { return d1; }
  size_t sizeY() { return d2; }
  size_t sizeZ() { return d3; }
  vector<T> GetData() { return data; }
  
 private:
  size_t d1,d2,d3;
  vector<T> data;
};

class GridPoint {

 public:
  GridPoint();
  virtual ~GridPoint() {
    channels.clear();
    backgroundYields.clear();
    
    useStatErrors.clear();
    val.clear();
    val_err.clear();

    bkg.clear();
    bkg_err.clear();
    data_err.clear();
    sig.clear();

    isSensitive.clear();
    signalYields.clear();
    obs.clear();
  }

  void AddChannel(TString chan) { 
    channels.push_back(chan);
  };

  void Init() {
    // chan/bkg
    vector2d<double> bkgY(channels.size(), nBackgrounds, 0.);
    backgroundYields = bkgY;

    // chan/bkg/bin
    vector3d<bool> useSt(channels.size(), nBackgrounds, 6, false);
    useStatErrors = useSt;
    
    vector3d<double> vl(channels.size(), nBackgrounds, 6, 0.);
    val = vl;
    val_err = vl;

    // chan/bin
    vector2d<double> bg(channels.size(), 6, 0.);
    bkg = bg;
    bkg_err = bg;
    data_err = bg;
    sig = bg;

    // chan
    signalYields.resize(channels.size());
    isSensitive.resize(channels.size());
    obs.resize(channels.size());
  }

  void Print();
  bool SetBackgroundYields(TFile * f);
  bool SetSignalYields(TFile * f);
  void SetUseStatError();

  int mStop;
  int mBino;

  double xsec;
  double xsecError;

  double lumi_sysError;

  vector<TString> channels;

  vector2d<double> backgroundYields;
  vector<double> signalYields;
  vector<bool> isSensitive;
  vector<double> obs;

  vector3d<bool> useStatErrors;
  vector3d<double> val;
  vector3d<double> val_err;

  vector2d<double> bkg;
  vector2d<double> bkg_err;
  vector2d<double> data_err;
  vector2d<double> sig;
  
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

  channels.clear();
  backgroundYields.clear();
  
  useStatErrors.clear();
  val.clear();
  val_err.clear();
  
  bkg.clear();
  bkg_err.clear();
  data_err.clear();
  sig.clear();
  
  isSensitive.clear();
  signalYields.clear();
  obs.clear();
    
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

  outfile << endl << "imax * number of channels" << endl;
  outfile << endl << "jmax * number of backgrounds" << endl;
  outfile << "kmax * number of nuisance parameters" << endl;
  outfile << "--------------------" << endl;
  outfile << "shapes * * limitInputs_bjj.root $CHANNEL/$PROCESS  $CHANNEL/$PROCESS_$SYSTEMATIC" << endl;
  outfile << "--------------------" << endl;

  outfile << "bin                 ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) outfile << "\t" << channels[i];
  }
  outfile << endl;

  outfile << "observation         ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) outfile << "\t" << obs[i];
  }
  outfile << endl;

  outfile << "--------------------" << endl;
  outfile << "bin                 ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) {
      for(int j = 0; j < nBackgrounds + 1; j++) outfile << "\t" << channels[i];
    }
  }
  outfile << endl;

  stringstream code;
  code << "_mst_" << mStop << "_m1_" << mBino;
    
  outfile << "process             ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) {
      outfile << "\tsignal" << code.str();
      for(int j = 0; j < nBackgrounds; j++) outfile << "\t" << backgroundNames[j];
    }
  }
  outfile << endl;

  outfile << "process             ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) {
      for(int j = 0; j < nBackgrounds + 1; j++) outfile << "\t" << i;
    }
  }
  outfile << endl;

  outfile << "rate                ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) {
      outfile << "\t" << signalYields[i];
      for(int j = 0; j < nBackgrounds; j++) outfile << "\t" << backgroundYields(i, j);
    }
  }
  outfile << endl;
  outfile << "--------------------" << endl;

  outfile << "lumi lnN          ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) {
      for(int j = 0; j < nBackgrounds + 1; j++) outfile << "\t" << lumi_sysError;
    }
  }
  outfile << endl;

  for(int iS = 0; iS < nSystematics; iS++) {

    outfile << systematicNames[iS].Data() << " shape          ";
    for(unsigned int i = 0; i < channels.size(); i++) {
      if(isSensitive[i]) {
	for(int j = 0; j < nBackgrounds + 1; j++) outfile << "\t1.0";
      }
    }
    outfile << endl;

  }

  outfile << "scale_tt shape ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) {
      outfile << "\t-";
      for(int j = 0; j < nBackgrounds; j++) {
	if(scale_tt[j]) outfile << "\t1.0";
	else outfile << "\t-";
      }
    }
  }
  outfile << endl;

  outfile << "scale_V shape ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) {
      outfile << "\t-";
      for(int j = 0; j < nBackgrounds; j++) {
	if(scale_V[j]) outfile << "\t1.0";
	else outfile << "\t-";
      }
    }
  }
  outfile << endl;

  outfile << "scale_VV shape ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) {
      outfile << "\t-";
      for(int j = 0; j < nBackgrounds; j++) {
	if(scale_VV[j]) outfile << "\t1.0";
	else outfile << "\t-";
      }
    }
  }
  outfile << endl;
  

  outfile << "pdf_gg shape ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) {
      outfile << "\t-";
      for(int j = 0; j < nBackgrounds; j++) {
	if(pdf_gg[j]) outfile << "\t1.0";
	else outfile << "\t-";
      }
    }
  }
  outfile << endl;

  outfile << "pdf_qq shape ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) {
      outfile << "\t-";
      for(int j = 0; j < nBackgrounds; j++) {
	if(pdf_qq[j]) outfile << "\t1.0";
	else outfile << "\t-";
      }
    }
  }
  outfile << endl;

  outfile << "pdf_qg shape ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) {
      outfile << "\t-";
      for(int j = 0; j < nBackgrounds; j++) {
	if(pdf_qg[j]) outfile << "\t1.0";
	else outfile << "\t-";
      }
    }
  }
  outfile << endl;

  outfile << "susy_xsec lnN     ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) {
      outfile << "\t" << 1. + xsecError/100.;
      for(int j = 0; j < nBackgrounds; j++) outfile << "\t-";
    }
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

  for(int ibin = 0; ibin < 6; ibin++) {

    for(int ibkg = 0; ibkg < nBackgrounds; ibkg++) {
      
      for(unsigned int ichan = 0; ichan < channels.size(); ichan++) {

	if(!useStatErrors(ichan, ibkg, ibin)) continue;

	if(isSensitive[ichan]) {

	  outfile << backgroundNames[ibkg] << "_" << channels[ichan] << "_stat_bin" << ibin + 1 << " shape";
	  for(unsigned int jchan = 0; jchan < channels.size(); jchan++) {
	    if(isSensitive[jchan]) {
	      outfile << "\t-";
	      for(int jbkg = 0; jbkg < nBackgrounds; jbkg++) {
		if(jbkg == ibkg && jchan == ichan) outfile << "\t1.0";
		else outfile << "\t-";
	      }
	    }
	  }
	  outfile << endl;

	}

      }
      
      // so for ibin=0 and ibkg=0:
      // ttjets_ele_SR1_stat_bin1 shape  (- 1 -----) (- - -----) (- - -----) (- - -----)
      // ttjets_muon_SR1_stat_bin1 shape (- - -----) (- 1 -----) (- - -----) (- - -----)
      // ttjets_ele_SR2_stat_bin1 shape  (- - -----) (- - -----) (- 1 -----) (- - -----)
      // ttjets_muon_SR2_stat_bin1 shape (- - -----) (- - -----) (- - -----) (- 1 -----)
      // etc...
    }

  } // stats block


} // Print()

bool GridPoint::SetBackgroundYields(TFile * f) {

  for(unsigned int i = 0; i < channels.size(); i++) {

    TH1D * h = (TH1D*)f->Get(channels[i]+"/data_obs");
    if(!h) return false;
    obs[i] = h->Integral();
    for(int ibin = 0; ibin < 6; ibin++) data_err(i, ibin) = h->GetBinError(ibin+1);
    
    for(int j = 0; j < nBackgrounds; j++) {
      h = (TH1D*)f->Get(channels[i]+"/"+backgroundNames[j]);
      if(!h) return false;
      backgroundYields(i, j) = h->Integral();
      for(int ibin = 0; ibin < 6; ibin++) {
	val(i, j, ibin) = h->GetBinContent(ibin+1);
	val_err(i, j, ibin) = h->GetBinError(ibin+1);

	bkg(i, ibin) += h->GetBinContent(ibin+1);
	bkg_err(i, ibin) = TMath::Sqrt(bkg_err(i, ibin)*bkg_err(i, ibin) + h->GetBinError(ibin+1)*h->GetBinError(ibin+1));
      }

    }

  }

  return true;
}

bool GridPoint::SetSignalYields(TFile * f) {
  
  stringstream code;
  code << "_mst_" << mStop << "_m1_" << mBino;
  TString code_t = code.str();

  for(unsigned int i = 0; i < channels.size(); i++) {
    TH1D * h = (TH1D*)f->Get(channels[i]+"/signal"+code_t);
    if(!h) return false;
    signalYields[i] = h->Integral();
    isSensitive[i] = (h->Integral() > epsilon);
    for(int ibin = 0; ibin < 6; ibin++) sig(i, ibin) = h->GetBinContent(ibin+1);
  }

  return true;
}

void GridPoint::SetUseStatError() {
  
  for(unsigned int chan = 0; chan < channels.size(); chan++) {
    for(int ibkg = 0; ibkg < nBackgrounds; ibkg++) {
      for(int bin = 0; bin < 6; bin++) {

	bool negligable = val(chan, ibkg, bin) < 0.01 ||
	  bkg_err(chan, bin) < data_err(chan, bin) / 5. ||
	  TMath::Sqrt(bkg_err(chan, bin)*bkg_err(chan, bin) - val_err(chan, ibkg, bin)*val_err(chan, ibkg, bin)) / bkg_err(chan, bin) > 0.95 ||
	  sig(chan, bin) / bkg(chan, bin) < 0.01;

	useStatErrors(chan, ibkg, bin) = !negligable;

      }
    }
  }

}

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
const TString systematicNames[nSystematics] = {"btagWeight", "puWeight", "topPt", "JEC", "eleSF", "muonSF", "photonSF"};

const double epsilon = 1e-10;

const TString datacard_dir = "datacards";

class GridPoint {

 public:
  GridPoint();
  virtual ~GridPoint() {
    channels.clear();
    backgroundYields.clear();
    
    useStatErrors.clear();
    useQCDStatErrors.clear();
    val.clear();
    val_err.clear();

    bkg.clear();
    bkg_err.clear();
    data_err.clear();
    sig.clear();

    qcd.clear();
    qcd_err.clear();

    isSensitive.clear();
    signalYields.clear();
    obs.clear();

    // chan
    qcdYields.clear();
    // chan/bin
    qcd.clear();
    qcd_err.clear();
  }

  void AddChannel(TString chan) { 
    channels.push_back(chan);
  };

  void Init() {

    // chan/bkg
    vector<double> bkgY(nBackgrounds, 0.);
    for(unsigned int i = 0; i < channels.size(); i++) backgroundYields.push_back(bkgY);

    // chan/bkg/bin
    vector< vector<bool> > useSt(nBackgrounds, vector<bool>(6, false));
    for(unsigned int i = 0; i < channels.size(); i++) useStatErrors.push_back(useSt);

    vector< vector<double> > vl(nBackgrounds, vector<double>(6, 0.));
    for(unsigned int i = 0; i < channels.size(); i++) val.push_back(vl);

    vector< vector<double> > vle(nBackgrounds, vector<double>(6, 0.));
    for(unsigned int i = 0; i < channels.size(); i++) val_err.push_back(vle);

    // chan/bin
    vector<double> bg(6, 0.);
    for(unsigned int i = 0; i < channels.size(); i++) bkg.push_back(bg);

    vector<double> bge(6, 0.);
    for(unsigned int i = 0; i < channels.size(); i++) bkg_err.push_back(bge);

    vector<double> de(6, 0.);
    for(unsigned int i = 0; i < channels.size(); i++) data_err.push_back(de);

    vector<double> sg(6, 0.);
    for(unsigned int i = 0; i < channels.size(); i++) sig.push_back(sg);

    // chan
    signalYields.resize(channels.size());
    isSensitive.resize(channels.size());
    obs.resize(channels.size());

    // chan
    qcdYields.resize(channels.size());

    // chan/bin
    vector<bool> useQCDSt(6, false);
    for(unsigned int i = 0; i < channels.size(); i++) useQCDStatErrors.push_back(useQCDSt);

    vector<double> v_qcd(6, 0.);
    for(unsigned int i = 0; i < channels.size(); i++) qcd.push_back(v_qcd);

    vector<double> v_qcde(6, 0.);
    for(unsigned int i = 0; i < channels.size(); i++) qcd_err.push_back(v_qcde);

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

  vector< vector<double> > backgroundYields;
  vector<double> signalYields;
  vector<bool> isSensitive;
  vector<double> obs;

  vector<double> qcdYields;

  vector< vector< vector<bool> > > useStatErrors;
  vector< vector< vector<double> > > val;
  vector< vector< vector<double> > > val_err;

  vector< vector<double> > bkg;
  vector< vector<double> > bkg_err;
  vector< vector<double> > data_err;
  vector< vector<double> > sig;

  vector< vector<bool> > useQCDStatErrors;
  vector< vector<double> > qcd;
  vector< vector<double> > qcd_err;

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

  qcdYields.clear();
  qcd.clear();
  qcd_err.clear();
  useQCDStatErrors.clear();
    
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
      for(int j = 0; j < nBackgrounds + 2; j++) outfile << "\t" << channels[i];
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
      outfile << "\tqcd";
    }
  }
  outfile << endl;

  outfile << "process             ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) {
      for(int j = 0; j < nBackgrounds + 2; j++) outfile << "\t" << j;
    }
  }
  outfile << endl;

  outfile << "rate                ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) {
      outfile << "\t" << signalYields[i];
      for(int j = 0; j < nBackgrounds; j++) outfile << "\t" << backgroundYields[i][j];
      outfile << "\t" << qcdYields[i];
    }
  }
  outfile << endl;
  outfile << "--------------------" << endl;

  outfile << "lumi lnN          ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) {
      for(int j = 0; j < nBackgrounds + 1; j++) outfile << "\t" << lumi_sysError;
      outfile << "\t-";
    }
  }
  outfile << endl;

  for(int iS = 0; iS < nSystematics; iS++) {

    outfile << systematicNames[iS].Data() << " shape          ";
    for(unsigned int i = 0; i < channels.size(); i++) {
      if(isSensitive[i]) {

	for(int j = 0; j < nBackgrounds + 1; j++) {
	  if(systematicNames[iS] == "eleSF" && channels[i].Contains("ele")) outfile << "\t1.0";
	  else if(systematicNames[iS] == "muonSF" && channels[i].Contains("muon")) outfile << "\t1.0";
	  else outfile << "\t-";
	}
	outfile << "\t-";

      }
    }
    outfile << endl;

  }

  outfile << "extraSystematic shape ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) {
      outfile << "\t-";
      for(int j = 0; j < nBackgrounds + 1; j++) outfile << "\t1.0";
    }
  }
  outfile << endl;

  outfile << "scale_tt shape ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) {
      outfile << "\t-";
      for(int j = 0; j < nBackgrounds; j++) {
	if(scale_tt[j]) outfile << "\t1.0";
	else outfile << "\t-";
      }
      outfile << "\t-";
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
      outfile << "\t-";
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
      outfile << "\t-";
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
      outfile << "\t-";
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
      outfile << "\t-";
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
      outfile << "\t-";
    }
  }
  outfile << endl;

  outfile << "susy_xsec lnN     ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) {
      outfile << "\t" << 1. + xsecError/100.;
      for(int j = 0; j < nBackgrounds + 1; j++) outfile << "\t-";
    }
  }
  outfile << endl;

  outfile << "def shape ";
  for(unsigned int i = 0; i < channels.size(); i++) {
    if(isSensitive[i]) {
      for(int j = 0; j < nBackgrounds + 1; j++) outfile << "\t-";
      outfile << "\t1.0";
    }
  }
  outfile << endl;

  for(int ibin = 0; ibin < 6; ibin++) {

    for(int ibkg = 0; ibkg < nBackgrounds; ibkg++) {
      
      for(unsigned int ichan = 0; ichan < channels.size(); ichan++) {

	if(!useStatErrors[ichan][ibkg][ibin]) continue;

	if(isSensitive[ichan]) {

	  outfile << backgroundNames[ibkg] << "_" << channels[ichan] << "_stat_bin" << ibin + 1 << " shape";
	  for(unsigned int jchan = 0; jchan < channels.size(); jchan++) {
	    if(isSensitive[jchan]) {
	      outfile << "\t-";
	      for(int jbkg = 0; jbkg < nBackgrounds; jbkg++) {
		if(jbkg == ibkg && jchan == ichan) outfile << "\t1.0";
		else outfile << "\t-";
	      }
	      outfile << "\t-";
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



    for(unsigned int ichan = 0; ichan < channels.size(); ichan++) {

      if(!useQCDStatErrors[ichan][ibin]) continue;

      if(isSensitive[ichan]) {

	outfile << "qcd_" << channels[ichan] << "_stat_bin" << ibin + 1 << " shape";
	for(unsigned int jchan = 0; jchan < channels.size(); jchan++) {
	  if(isSensitive[jchan]) {
	    for(int jbkg = 0; jbkg < nBackgrounds + 1; jbkg++) outfile << "\t-";
	    if(jchan == ichan) outfile << "\t1.0";
	    else outfile << "\t-";
	  }
	}
	outfile << endl;

      }

    }

  } // stats block


} // Print()

bool GridPoint::SetBackgroundYields(TFile * f) {

  for(unsigned int i = 0; i < channels.size(); i++) {

    TH1D * h = (TH1D*)f->Get(channels[i]+"/data_obs");
    if(!h) return false;
    obs[i] = h->Integral();
    for(int ibin = 0; ibin < 6; ibin++) data_err[i][ibin] = h->GetBinError(ibin+1);
    
    for(int j = 0; j < nBackgrounds; j++) {
      h = (TH1D*)f->Get(channels[i]+"/"+backgroundNames[j]);
      if(!h) return false;
      backgroundYields[i][j] = h->Integral();
      for(int ibin = 0; ibin < 6; ibin++) {
	val[i][j][ibin] = h->GetBinContent(ibin+1);
	val_err[i][j][ibin] = h->GetBinError(ibin+1);

	bkg[i][ibin] += h->GetBinContent(ibin+1);
	bkg_err[i][ibin] = TMath::Sqrt(bkg_err[i][ibin]*bkg_err[i][ibin] + h->GetBinError(ibin+1)*h->GetBinError(ibin+1));
      }

    }

    h = (TH1D*)f->Get(channels[i]+"/qcd");
    if(!h) return false;
    qcdYields[i] = h->Integral();
    for(int ibin = 0; ibin < 6; ibin++) {
      qcd[i][ibin] = h->GetBinContent(ibin+1);
      qcd_err[i][ibin] = h->GetBinError(ibin+1);

      bkg[i][ibin] += h->GetBinContent(ibin+1);
      bkg_err[i][ibin] = TMath::Sqrt(bkg_err[i][ibin]*bkg_err[i][ibin] + h->GetBinError(ibin+1)*h->GetBinError(ibin+1));
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
    for(int ibin = 0; ibin < 6; ibin++) sig[i][ibin] = h->GetBinContent(ibin+1);
  }

  return true;
}

void GridPoint::SetUseStatError() {
  
  for(unsigned int chan = 0; chan < channels.size(); chan++) {

    for(int ibkg = 0; ibkg < nBackgrounds; ibkg++) {

      for(int bin = 0; bin < 6; bin++) {
	bool negligable = val[chan][ibkg][bin] < 0.01 ||
	  bkg_err[chan][bin] < data_err[chan][bin] / 5. ||
	  TMath::Sqrt(bkg_err[chan][bin]*bkg_err[chan][bin] - val_err[chan][ibkg][bin]*val_err[chan][ibkg][bin]) / bkg_err[chan][bin] > 0.95 ||
	  sig[chan][bin] / bkg[chan][bin] < 0.01;

	useStatErrors[chan][ibkg][bin] = !negligable;
      }

    }

    for(int bin = 0; bin < 6; bin++) {
      bool negligable = qcd[chan][bin] < 0.01 ||
	bkg_err[chan][bin] < data_err[chan][bin] / 5. ||
	TMath::Sqrt(bkg_err[chan][bin]*bkg_err[chan][bin] - qcd_err[chan][bin]*qcd_err[chan][bin]) / bkg_err[chan][bin] > 0.95 ||
	sig[chan][bin] / bkg[chan][bin] < 0.01;

      useQCDStatErrors[chan][bin] = !negligable;
    }

  }

}

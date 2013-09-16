#include <TFile.h>
#include <TH1F.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <vector>
#include <cmath>

using namespace std;

const int NCH = 5;
const double bins[NCH] = {50,60,70,80,100};
const double epsilon = 1e-10;
const double luminosity = 19712;
const double DataMCScale = 1.005; // 1.005 +- 0.001 (stat) +- 3e-7 (fit) +- 0.006 (PU) +- 5e-6 (bias) +- 0.002 (e/g difference)

const double triggerEff = 0.883;  
const double triggerEff_err = 0.001 / triggerEff; 

const TString work_dir = "./";
const TString datacard_dir = work_dir+"datacards";

class BinInfo {

 public:
  BinInfo() { 
    x = y = error = 0;
    name = "";
  }
  ~BinInfo() {}

  int x;
  double y;
  double error;

  TString name;
};

class GridPoint {

 public:
  GridPoint();
  virtual ~GridPoint() {
    ggBins.clear();
    ewBins.clear();
    qcdBins.clear();
    qcdSysErrors.clear();
    sig_ggBins.clear();
    sig_ffBins.clear();
    sig_gfBins.clear();
    sigBins.clear();
  }

  int mStop;
  int mBino;

  int ngen;
  double acc;

  double lumi;
  double xsecValue;
  double xsecError;
  double lumi_sysError;
  double qcd_sysError;
  double ew_sysError;
  double sig_sysError;

  vector<BinInfo> ggBins;
  vector<BinInfo> ewBins;
  vector<BinInfo> qcdBins;
  vector<BinInfo> qcdSysErrors;
  vector<BinInfo> sig_ggBins;
  vector<BinInfo> sig_ffBins;
  vector<BinInfo> sig_gfBins;
  vector<BinInfo> sigBins;

  double limit;           // observed limit
  double explimit;        // expected limit
  double explimit_1L;     // expected limit -1 sigma
  double explimit_1H;     // expected limit +1 sigma
  double explimit_2L;     // expected limit -2 sigma
  double explimit_2H;     // expected limit +2 sigma
};

GridPoint::GridPoint() {
  mStop = mBino = 0;
  ngen = 15000;
  acc = 0;
  lumi = luminosity;
  xsecValue = 0;
  xsecError = 0;

  lumi_sysError = 1.026;
  ew_sysError = 1.0497; // fakerate systematic error
  qcd_sysError = 1.25; // shape difference between ee and ff -- default placeholder of 25%, will calculate bin-by-bin later

  // .00199 for electron/photon difference
  // 2.985e-7 for signal fit over/underestimation
  // .00597 for pileup
  // 4.975e-6 for bias
  // total error of 0.0063 (.63%), and the scale factor is 1.005 +/- .001
  sig_sysError = 0.00637 * sqrt(2);

  ggBins.clear();
  ewBins.clear();
  qcdBins.clear();
  qcdSysErrors.clear();
  sig_ggBins.clear();
  sig_ffBins.clear();
  sigBins.clear();

  limit = explimit = explimit_1L = explimit_1H = explimit_2L = explimit_2H = 0;
}

void getBins(TFile * f_eleJets, TFile * f_muJets, TFile * f_hadronic, TString category, vector<BinInfo>& binInfos) {

  TH1F * ele = (TH1F*)f_eleJets->Get(category+"_eleJets");
  TH1F * mu = (TH1F*)f_muJets->Get(category+"_muJets");
  TH1F * had = (TH1F*)f_hadronic->Get(category+"_hadronic");

  BinInfo bin_ele;
  bin_ele.x = 80;
  bin_ele.y = ele->IntegralAndError(ele->GetXaxis()->FindBin(80.), -1, bin_ele.error);
  bin_ele.name = "eleJets";
  binInfos.push_back(bin_ele);

  BinInfo bin_mu;
  bin_mu.x = 80;
  bin_mu.y = mu->IntegralAndError(mu->GetXaxis()->FindBin(80.), -1, bin_mu.error);
  bin_mu.name = "muJets";
  binInfos.push_back(bin_mu);

  BinInfo bin_had;
  bin_had.x = 100;
  bin_had.y = had->IntegralAndError(had->GetXaxis()->FindBin(100.), -1, bin_had.error);
  bin_had.name = "hadronic";
  binInfos.push_back(bin_had);

}

void getSigBins(TFile * f, TString category, TString code, vector<BinInfo>& binInfos) {

  TH1F * ele = (TH1F*)f->Get(category+"_eleJets"+code);
  TH1F * mu = (TH1F*)f->Get(category+"_muJets"+code);
  TH1F * had = (TH1F*)f->Get(category+"_hadronic"+code);

  if(!ele) {
    cout << "eleJets is unavailable for " << code << " !!" <<endl;
    return;
  }

  if(!mu) {
    cout << "muJets is unavailable for " << code << " !!" <<endl;
    return;
  }

  if(!had) {
    cout << "hadronic is unavailable for " << code << " !!" <<endl;
    return;
  }

  BinInfo bin_ele;
  bin_ele.x = 80;
  bin_ele.y = ele->IntegralAndError(ele->GetXaxis()->FindBin(80.), -1, bin_ele.error);
  bin_ele.name = "eleJets";
  binInfos.push_back(bin_ele);

  BinInfo bin_mu;
  bin_mu.x = 80;
  bin_mu.y = mu->IntegralAndError(mu->GetXaxis()->FindBin(80.), -1, bin_mu.error);
  bin_mu.name = "muJets";
  binInfos.push_back(bin_mu);

  BinInfo bin_had;
  bin_had.x = 100;
  bin_had.y = had->IntegralAndError(had->GetXaxis()->FindBin(100.), -1, bin_had.error);
  bin_had.name = "hadronic";
  binInfos.push_back(bin_had);

}

void printBins(vector<BinInfo>& binInfos){
  for(vector<BinInfo>::iterator it = binInfos.begin(); it != binInfos.end(); it++)
    printf("%s : %3d(%6.2f +/- %5.2f) ", (it->name).Data(), it->x, it->y, it->error);
  printf("\n");
}

void printBins(fstream& of, vector<BinInfo>& binInfos) {
  for(vector<BinInfo>::iterator it = binInfos.begin(); it != binInfos.end(); it++)
    of << "(" << it->name << " : " << it->x << ", " << it->y << " +/- " << it->error << ") ";
  of << endl;
}

void readData(TFile * f_eleJets, TFile * f_muJets, TFile * f_hadronic,
	      vector<BinInfo>& ggBins,
	      vector<BinInfo>& qcdBins,
	      vector<BinInfo>& ewBins,
	      vector<BinInfo>& qcdSysErrors) {

  ggBins.clear();
  getBins(f_eleJets, f_muJets, f_hadronic, "pfMET_gg", ggBins);
  cout << "gg  events ----------" << endl;
  printBins(ggBins);

  qcdBins.clear();
  getBins(f_eleJets, f_muJets, f_hadronic, "pfMET_gf", qcdBins);
  cout << "qcd events ----------" << endl;
  printBins(qcdBins);

  ewBins.clear();
  getBins(f_eleJets, f_muJets, f_hadronic, "pfMET_eg", ewBins);
  cout << "ew  events ----------" << endl;
  printBins(ewBins);

  vector<BinInfo> qcd_ffBins;
  getBins(f_eleJets, f_muJets, f_hadronic, "pfMET_ff", qcd_ffBins);
  cout << "qcd_ff events ----------" << endl;
  printBins(qcd_ffBins);

  for(unsigned int i = 0; i < qcd_ffBins.size(); i++){
    double diff = qcd_ffBins[i].y - qcdBins[i].y;
    BinInfo bin;
    bin.x = qcdBins[i].x;
    bin.y = fabs(diff);
    bin.name = qcdBins[i].name;
    qcdSysErrors.push_back(bin);
  }

  cout << "qcd_sysErrors ---------" << endl;
  printBins(qcdSysErrors);

}

void readSig(TFile * f, int mStop, int mBino,
	     vector<BinInfo>& sig_ggBins,
	     vector<BinInfo>& sig_ffBins,
	     vector<BinInfo>& sig_gfBins,
	     vector<BinInfo>& sigBins) {

  

  stringstream code;
  code << "_mst_" << mStop << "_m1_" << mBino;

  sig_ggBins.clear();
  getSigBins(f, "met_gg", code.str().c_str(), sig_ggBins);

  sig_gfBins.clear();
  getSigBins(f, "met_gf", code.str().c_str(), sig_gfBins);

  sig_ffBins.clear();
  getSigBins(f, "met_ff", code.str().c_str(), sig_ffBins);

  sigBins.clear();
  for(unsigned int k = 0; k < sig_ggBins.size(); k++) {
    BinInfo b;
    b.x = sig_ggBins[k].x;
    b.y = sig_ggBins[k].y - sig_gfBins[k].y;
    b.name = sig_ggBins[k].name;
    if(b.y < 1e-6) b.y = 0.0;
    b.error = sqrt(sig_ggBins[k].error * sig_ggBins[k].error + sig_gfBins[k].error * sig_gfBins[k].error);
    sigBins.push_back(b);
  }

}

void GetXSection(vector<GridPoint>& grid, TString datafile)
{
  ifstream fin;
  fin.open(datafile.Data());
      
  while(1) {
    int mStop;
    double xsec, xsecErr;

    fin >> mStop >> xsec >> xsecErr;
    
    if(!fin.good()) break;
    for(unsigned int ig = 0; ig < grid.size(); ig++) {
      if(mStop == grid[ig].mStop) {
        grid[ig].mStop = mStop;
        grid[ig].xsecValue = abs(xsec);
        grid[ig].xsecError = abs(xsecErr);
      }
    }
  }
  fin.close();
}

void makeSignalGains(vector<GridPoint>& grid) {

  for(vector<GridPoint>::iterator it = grid.begin(); it != grid.end(); it++) {

    double nevt = 0;
    for(vector<BinInfo>::iterator bit = it->sigBins.begin(); bit != it->sigBins.end(); bit++) {
      nevt += bit->y;
      double acc = bit->y * DataMCScale * DataMCScale / it->ngen;
      double nexpected = it->xsecValue * it->lumi * acc * triggerEff;
      double error = bit->error * it->xsecValue * it->lumi * DataMCScale * DataMCScale / it->ngen;
      bit->y = nexpected;
      bit->error = error;
    }
    it->sig_sysError = sqrt(it->sig_sysError*it->sig_sysError + it->xsecError*it->xsecError + triggerEff_err*triggerEff_err);
    it->sig_sysError += 1; // relative error
    it->acc = nevt/it->ngen;
  }

}

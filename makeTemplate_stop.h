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
  BinInfo() { x = y = error = 0; }
  ~BinInfo() {}

  int x;
  double y;
  double error;
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

void getBins(TH1F* h, vector<BinInfo>& binInfos) {

  for(int i=0; i<NCH; i++) {
    int ibin = h->GetXaxis()->FindBin(bins[i]);
    int jbin = (i < NCH - 1) ? h->GetXaxis()->FindBin(bins[i+1]) - 1 : -1;

    BinInfo binInfo;
    binInfo.x = bins[i];
    binInfo.y = h->IntegralAndError(ibin,jbin,binInfo.error);
    binInfos.push_back(binInfo);
  }

}

void printBins(vector<BinInfo>& binInfos){
  for(vector<BinInfo>::iterator it = binInfos.begin(); it != binInfos.end(); it++)
    printf("%3d(%6.2f +/- %5.2f) ",it->x,it->y,it->error);
  printf("\n");
}

void printBins(fstream& of, vector<BinInfo>& binInfos) {
  for(vector<BinInfo>::iterator it = binInfos.begin(); it != binInfos.end(); it++)
    of << "(" << it->x << ", " << it->y << " +/- " << it->error << ") ";
  of << endl;
}

void readData(TFile* f, TString req,
	      vector<BinInfo>& ggBins,
	      vector<BinInfo>& qcdBins,
	      vector<BinInfo>& ewBins,
	      vector<BinInfo>& qcdSysErrors) {

  TH1F * gg = (TH1F*) f->Get("met_gg_"+req);
  ggBins.clear();
  getBins(gg, ggBins);
  cout << "gg  events ----------" << endl;
  printBins(ggBins);

  TH1F * qcd = (TH1F*) f->Get("met_gf_"+req);
  qcdBins.clear();
  getBins(qcd, qcdBins);
  cout << "qcd events ----------" << endl;
  printBins(qcdBins);

  TH1F * ew = (TH1F*) f->Get("met_eg_"+req);
  ewBins.clear();
  getBins(ew, ewBins);
  cout << "ew  events ----------" << endl;
  printBins(ewBins);

  TH1F * qcd_ff = (TH1F*) f->Get("met_ff_"+req);
  vector<BinInfo> qcd_ffBins;
  getBins(qcd_ff, qcd_ffBins);
  cout << "qcd_ff events ----------" << endl;
  printBins(qcd_ffBins);

  for(unsigned int i = 0; i < qcd_ffBins.size(); i++){
    double diff = qcd_ffBins[i].y - qcdBins[i].y;
    BinInfo bin;
    bin.x = qcdBins[i].x;
    bin.y = fabs(diff);
    qcdSysErrors.push_back(bin);
  }

  cout << "qcd_sysErrors ---------" << endl;
  printBins(qcdSysErrors);

}

void readSig(TFile* f, TString req, int mStop, int mBino,
	     vector<BinInfo>& sig_ggBins,
	     vector<BinInfo>& sig_ffBins,
	     vector<BinInfo>& sig_gfBins,
	     vector<BinInfo>& sigBins) {

  stringstream ggname;
  ggname << "met_gg_" << req << "_mst_" << mStop << "_m1_" << mBino;

  TH1F * gg = (TH1F*) f->Get(ggname.str().c_str());
  if(!gg) {
    cout << "histogram " << ggname.str() << " is not available!!!" << endl;
    return;
  }
  sig_ggBins.clear();
  getBins(gg, sig_ggBins);

  stringstream gfname;
  gfname << "met_gf_" << req << "_mst_" << mStop << "_m1_" << mBino;

  TH1F * gf = (TH1F*) f->Get(gfname.str().c_str());
  sig_gfBins.clear();
  getBins(gf, sig_gfBins);

  stringstream ffname;
  ffname << "met_ff_" << req << "_mst_" << mStop << "_m1_" << mBino;

  TH1F* ff = (TH1F*) f->Get(ffname.str().c_str());
  sig_ffBins.clear();
  getBins(ff, sig_ffBins);

  sigBins.clear();
  int n = sig_ggBins.size();
  for(unsigned int k = 0; k < sig_ggBins.size(); k++) {
    BinInfo b;
    b.x = sig_ggBins[k].x;
    b.y = sig_ggBins[k].y - sig_gfBins[k].y;
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

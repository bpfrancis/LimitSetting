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

  int mG;
  int mS;
  int mN;
  int mW; //wino

  int ngen;
  double acc;

  double lumi;
  double xsecValue;       // NLO XS from Prospino in pb
  double xsecPDFError;    // PDF error
  double xsecRSErrorNeg;  // renormalization scale error
  double xsecRSErrorPos;  // renormalization scale error
  double accErrorPDF;     // acceptance errors on XS PDF
  double lumi_sysError;
  double qcd_sysError;
  double ew_sysError;
  double sig_sysError;

  vector<BinInfo> ggBins;
  vector<BinInfo> ewBins;
  vector<BinInfo> qcdBins;
  vector<BinInfo> qcdSysErrors;
  vector<BinInfo> sig_ggBins;
  vector<BinInfo> sig_gfBins;
  vector<BinInfo> sig_ffBins;
  vector<BinInfo> sigBins;

  double limit;           // observed limit
  double explimit;        // expected limit
  double explimit_1L;     // expected limit -1 sigma
  double explimit_1H;     // expected limit +1 sigma
  double explimit_2L;     // expected limit -2 sigma
  double explimit_2H;     // expected limit +2 sigma

};

GridPoint::GridPoint() {
  mG = mS = mN = mW = 0;
  ngen = 10000;
  acc = 0;
  lumi = luminosity; // int. luminosity
  xsecValue = 0;
  xsecPDFError = 0;
  xsecRSErrorNeg = 0;
  xsecRSErrorPos = 0;
  accErrorPDF = 0;

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

void getBins(TH1F* h, vector<BinInfo>& binInfos){

  for(int i=0; i<NCH; i++) {
    int ibin = h->GetXaxis()->FindBin(bins[i]);
    int jbin = -1;
    if(i<NCH-1) jbin = h->GetXaxis()->FindBin(bins[i+1]) - 1;
    BinInfo binInfo;
    binInfo.x = bins[i];
    //    binInfo.y = h->IntegralAndError(ibin,jbin,binInfo.error,"width");
    binInfo.y = h->IntegralAndError(ibin,jbin,binInfo.error);
    binInfos.push_back(binInfo);
  }

}

void printBins(vector<BinInfo>& binInfos){

  for(vector<BinInfo>::iterator it = binInfos.begin();
      it != binInfos.end(); it++) {
    printf("%3d(%6.2f +/- %5.2f) ",it->x,it->y,it->error);
  }
  printf("\n");

}

void printBins(fstream& of, vector<BinInfo>& binInfos) {
  for(vector<BinInfo>::iterator it = binInfos.begin();
      it != binInfos.end(); it++) {
    of << "(" << it->x << ", " << it->y << " +/- " << it->error << ") ";
  }
  of << endl;
}

void readData(TFile* f, TString jet,
	      vector<BinInfo>& ggBins,
	      vector<BinInfo>& qcdBins,
	      vector<BinInfo>& ewBins,
	      vector<BinInfo>& qcdSysErrors) {

  TH1F* gg = (TH1F*) f->Get("pfMET_gg_"+jet);
  ggBins.clear();
  getBins(gg,ggBins);
  cout << "gg  events ----------" << endl;
  printBins(ggBins);

  //  TH1F* qcd = (TH1F*) f->Get("met_qcd_avg_"+jet);
  TH1F* qcd = (TH1F*) f->Get("pfMET_ff_"+jet);
  qcdBins.clear();
  getBins(qcd,qcdBins);
  cout << "qcd events ----------" << endl;
  printBins(qcdBins);

  TH1F* ew = (TH1F*) f->Get("pfMET_eg_"+jet);
  ewBins.clear();
  getBins(ew,ewBins);
  cout << "ew  events ----------" << endl;
  printBins(ewBins);

  TH1F* qcd_gf = (TH1F*) f->Get("pfMET_gf_"+jet);
  vector<BinInfo> qcd_gfBins;
  getBins(qcd_gf,qcd_gfBins);
  cout << "qcd_gf events ----------" << endl;
  printBins(qcd_gfBins);

  int N = int(qcd_gfBins.size());
  for(int i=0; i<N; i++){
    double diff = qcd_gfBins[i].y - qcdBins[i].y;
    BinInfo bin;
    bin.x = qcdBins[i].x;
    bin.y = fabs(diff);
    qcdSysErrors.push_back(bin);
  }

  cout << "qcd_sysErrors ---------" << endl;
  printBins(qcdSysErrors);

}

void readSig(TFile* f, TString bino, TString jet, int mS, int mG, int mN,
	     vector<BinInfo>& sig_ggBins,
	     vector<BinInfo>& sig_ffBins,
	     vector<BinInfo>& sig_gfBins,
	     vector<BinInfo>& sigBins) {

  stringstream ggname;
  if(bino.Contains("SMSScan")) ggname << "h_gg_met_" << jet << "_mG" << mG << "_mN_" << mN;
  else ggname << "h_gg_met_" << jet << "_mS" << mS << "_mG" << mG << "_mN" << mN;

  TH1F* gg = (TH1F*) f->Get(ggname.str().c_str());
  if(!gg) {
    cout << "histogram " << ggname.str() << " is not available!!!" << endl;
    return;
  }
  sig_ggBins.clear();
  getBins(gg,sig_ggBins);

  stringstream ffname;
  if(bino.Contains("SMSScan")) ffname << "h_ff_met_" << jet << "_mG" << mG << "_mN_" << mN;
  else ffname << "h_ff_met_" << jet << "_mS" << mS << "_mG" << mG << "_mN" << mN;

  TH1F* ff = (TH1F*) f->Get(ffname.str().c_str());
  sig_ffBins.clear();
  getBins(ff,sig_ffBins);

  stringstream gfname;
  if(bino.Contains("SMSScan")) gfname << "h_gf_met_" << jet << "_mG" << mG << "_mN_" << mN;
  else gfname << "h_gf_met_" << jet << "_mS" << mS << "_mG" << mG << "_mN" << mN;

  TH1F * gf = (TH1F*) f->Get(gfname.str().c_str());
  sig_gfBins.clear();
  getBins(gf, sig_gfBins);

  sigBins.clear();
  int n = sig_ggBins.size();
  for(int k=0; k<n; k++) {
    BinInfo b;
    b.x = sig_ggBins[k].x;
    b.y = sig_ggBins[k].y - sig_ffBins[k].y;
    if(b.y < 1e-6) b.y = 0.0;
    b.error = sqrt(sig_ggBins[k].error * sig_ggBins[k].error + sig_ffBins[k].error * sig_ffBins[k].error);
    sigBins.push_back(b);
  }//k
  //printBins(sigBins);

}

void GetNGen(vector<GridPoint>& grid, TString datafile) {
  unsigned int ngrid = grid.size();
  ifstream fin;
  fin.open(datafile.Data());
  while(1){
    int mG, mS, mChi, nGen;
    double dummy;
    fin >> mG >> mS >> mChi >> nGen >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
    if(!fin.good()) break;
    for(unsigned int ig=0; ig<ngrid; ig++){
      if( (mG ==grid[ig].mG) && (mS ==grid[ig].mS) ){
	grid[ig].ngen = nGen;
      }
    }
  }
  fin.close();
}

void GetXSection(vector<GridPoint>& grid, TString datafile)
{
  unsigned int ngrid = grid.size();
  ifstream fin;
  // RS errors are absolute.   
  fin.open(datafile.Data());
      
  while(1) {
    int nGen_, mS_, mG_, mB_, mW_;
    double xsec, rsp, rsm, dummy;   
    string dummyS;
    //     nGen     mS     mG     mBino  mWino  LO:       LOxsec   +         LOerr+   -         LOerr-   NLO:      xsec    +       NLOerr+  -         NLOerr-
    fin >> nGen_ >> mS_ >> mG_ >> mB_ >> mW_ >> dummyS >> dummy >> dummyS >> dummy >> dummyS >> dummy >> dummyS >> xsec >> dummyS >> rsp >> dummyS >> rsm;
    if(!fin.good()) break;
    for(unsigned int ig = 0; ig < ngrid; ig++) {
      if((mG_ == grid[ig].mG) && (mS_ == grid[ig].mS)) {
        grid[ig].ngen = nGen_;
        grid[ig].xsecValue = abs(xsec);
        grid[ig].xsecRSErrorPos = abs(rsp/xsec);
        grid[ig].xsecRSErrorNeg = abs(rsm/xsec);
      }
    }
  }
  fin.close();
}

void GetPDFErrorsXSection(vector<GridPoint>& grid, TString datafile)
{
  unsigned int ngrid = grid.size();
  ifstream fin;
  fin.open(datafile.Data());
  if(datafile.Contains("mNScan")) {
    while(1){
      int mG, mS, mN;
      double exs;
      fin >> mG >> mS >> mN >> exs;
      if(!fin.good()) break;
      for(unsigned int ig=0; ig<ngrid; ig++){
	if( (mG == grid[ig].mG) && (mS == grid[ig].mS) && (mN == grid[ig].mN) ) {
	  // errors are in %. convert to relative error
	  grid[ig].xsecPDFError = 0.01*exs;
	}
      } 
    }
  }
  else {
    while(1){
      int mG, mS;
      double exs, dummy;
      //     ngen     mG    mS    mB       mW       pdfer  accPdfErr
      fin >> dummy >> mG >> mS >> dummy >> dummy >> exs >> dummy;
      //fin >> mG >> mS >> exs;
      if(!fin.good()) break;
      for(unsigned int ig=0; ig<ngrid; ig++){
        if( (mG == grid[ig].mG) && (mS == grid[ig].mS) ) {
          // errors are in %. convert to relative error
          grid[ig].xsecPDFError = 0.01*exs;
        }
      }
    }
  }
  fin.close();
}

void GetPDFErrorsAcceptance(vector<GridPoint>& grid, TString datafile)
{
  unsigned int ngrid = grid.size();
  ifstream fin;
  fin.open(datafile.Data());
  if(datafile.Contains("mNScan")) {
    while(1){
      int mG, mS, mN;
      double exs;
      fin >> mG >> mS >> mN >> exs;
      if(!fin.good()) break;
      for(unsigned int ig=0; ig<ngrid; ig++){
	if( (mG == grid[ig].mG) && (mS == grid[ig].mS) && (mN == grid[ig].mN) ) {
	  // errors are in %. convert to relative error
	  grid[ig].accErrorPDF = 0.01*exs;
	}
      } 
    }
  }
  else {
    while(1){
      int mG, mS;
      double exs, dummy;
      //     ngen     mG    mS    mB       mW       pdfer  accPdfErr
      fin >> dummy >> mG >> mS >> dummy >> dummy >> dummy >> exs;
      //fin >> mG >> mS >> exs;
      if(!fin.good()) break;
      for(unsigned int ig=0; ig<ngrid; ig++){
        if( (mG == grid[ig].mG) && (mS == grid[ig].mS) ) {
          // errors are in %. convert to relative error
          grid[ig].accErrorPDF = 0.01*exs;
        }
      }
    }
  }
  fin.close();
}

void GetSMSXSection(vector<GridPoint>& grid, TString datafile) {

  cout << "Read SMS XSection file " << datafile.Data() << " ------------------" << endl;

  unsigned int ngrid = grid.size();

  ifstream fin;
  fin.open(datafile.Data());

  int mG;
  double xsec, rsp, rsm;
  while(1){

    // SMS Scan MC : rsp and rsm contains PDF errors and RS errors
    fin >> mG >> xsec >> rsp;
    rsm = rsp;
    //    cout << "mG " << mG << " xsec " << xsec << " rsp " << rsp << " rsm " << rsm << endl;

    if(!fin.good()) break;
    if(abs(xsec) < 1.0e-20) continue;

    for(unsigned int ig=0; ig<ngrid; ig++){
      if(mG == grid[ig].mG) {
	grid[ig].xsecValue = abs(xsec);
	grid[ig].xsecRSErrorPos = abs(rsp);  // already relative error
	grid[ig].xsecRSErrorNeg = abs(rsm);  // already relative error
	grid[ig].xsecPDFError = 0;
	grid[ig].accErrorPDF = 0;
      }
    }
  }
  fin.close();
}

void makeSignalGains(vector<GridPoint>& grid) {

  for(vector<GridPoint>::iterator it = grid.begin();
      it != grid.end(); it++) {
    double nevt = 0;
    for(vector<BinInfo>::iterator bit = it->sigBins.begin();
	bit != it->sigBins.end(); bit++) {
      nevt += bit->y;
      double acc = bit->y * DataMCScale * DataMCScale / it->ngen;
      double nexpected = it->xsecValue * it->lumi * acc * triggerEff;
      double error = bit->error * it->xsecValue * it->lumi * DataMCScale * DataMCScale / it->ngen;
      bit->y = nexpected;
      bit->error = error;
      //      cout << "y, xsec, acc : " << bit->y << ", " << it->xsecValue << ", " << acc << endl;
    }// bit
    it->sig_sysError = sqrt(it->sig_sysError*it->sig_sysError + it->accErrorPDF*it->accErrorPDF + triggerEff_err*triggerEff_err);
    it->sig_sysError += 1; // relative error
    it->acc = nevt/it->ngen;
    //    printBins(it->sigBins);
  }// it

}



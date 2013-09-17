#include <TFile.h>
#include <TH1F.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <vector>
#include <cmath>

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

  GridPoint() {Init();}
  ~GridPoint() {Init();}

  void Init();

  int mG;
  int mS;
  int mN;
  int mW; //wino
  int mStop;

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

  std::vector<BinInfo> ggBins;
  std::vector<BinInfo> ewBins;
  std::vector<BinInfo> qcdBins;
  std::vector<BinInfo> qcdSysErrors;
  std::vector<BinInfo> sig_ggBins;
  std::vector<BinInfo> sig_ffBins;
  std::vector<BinInfo> sigBins;

  double limit;           // observed limit
  double explimit;        // expected limit
  double explimit_1L;     // expected limit -1 sigma
  double explimit_1H;     // expected limit +1 sigma
  double explimit_2L;     // expected limit -2 sigma
  double explimit_2H;     // expected limit +2 sigma

};

void GridPoint::Init() {
  mG = mS = mN = mW = mStop = 0;
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


void readData(TFile* f, TString jet,
	      std::vector<BinInfo>& ggBins,
	      std::vector<BinInfo>& qcdBins,
	      std::vector<BinInfo>& ewBins,
	      std::vector<BinInfo>& qcdSysErrors);
void readSig(TFile* f, TString bino, TString jet, int mS, int mG, int mN,
	     std::vector<BinInfo>& sig_ggBins,
	     std::vector<BinInfo>& sig_ffBins,
	     std::vector<BinInfo>& sigBins);
void getBins(TH1F* h, std::vector<BinInfo>& binInfos);
void printBins(std::vector<BinInfo>& binInfos);
void printBins(std::fstream& of, std::vector<BinInfo>& binInfos);
void GetNGen(std::vector<GridPoint>& grid, TString datafile);
void GetXSection(std::vector<GridPoint>& grid, TString datafile);
void GetPDFErrorsXSection(std::vector<GridPoint>& grid, TString datafile);
void GetPDFErrorsAcceptance(std::vector<GridPoint>& grid, TString datafile);
void GetSMSXSection(std::vector<GridPoint>& grid, TString datafile);
void GetStopXSection(std::vector<GridPoint>& grid, TString datafile);
void makeSignalGains(std::vector<GridPoint>& grid);
void makeDataCard(std::vector<GridPoint>& grid, TString bino, TString jet);
void testDataCard(std::vector<GridPoint>& grid, TString bino, TString jet);


void makeTemplate(TString bino = "bino", TString jet="1jet") {

  TString hist_dir = work_dir+"inputHists";
  TString dataHist = hist_dir + "/" + "met_reweighted_"+jet+".root";
  TString sigHist = hist_dir + "/" + "signal_contamination_" + bino + "_chi0375.root";
  if(bino.Contains("mNScan")) sigHist = hist_dir + "/" + "signal_contamination_bino_chi0.root";
  if(bino.Contains("SMSScan")) sigHist = hist_dir + "/" + "signal_contamination_T5gg.root";
  if(bino.Contains("stop-bino")) sigHist = hist_dir + "/contamination_stop.root";

  TFile* fData = new TFile(dataHist,"READ");
  std::vector<BinInfo> ggBins;
  std::vector<BinInfo> qcdBins;
  std::vector<BinInfo> ewBins;
  std::vector<BinInfo> qcdSysErrors;
  readData(fData,jet,ggBins,qcdBins,ewBins,qcdSysErrors);

  TFile* fSig = new TFile(sigHist,"READ");

  std::vector<GridPoint> grids;

  int npoint = 0;

  if(bino.Contains("mNScan")) {
    for(int iN=150; iN<=1050; iN+=100) {
      for(int iG=160; iG<=2000; iG+=80) {
	npoint++;
	GridPoint grid;
	grid.mG = iG;
	grid.mS = 2500;
	grid.mN = iN;
	grid.ngen = 10000;
	grid.lumi = luminosity;

	std::vector<BinInfo> sig_ggBins;
	std::vector<BinInfo> sig_ffBins;
	std::vector<BinInfo> sigBins;
	readSig(fSig,bino,jet,2500,iG,iN,sig_ggBins,sig_ffBins,sigBins);

	grid.ggBins = ggBins;
	grid.qcdBins = qcdBins;
	grid.ewBins = ewBins;
	grid.qcdSysErrors = qcdSysErrors;
	grid.sig_ggBins = sig_ggBins;
	grid.sig_ffBins = sig_ffBins;
	grid.sigBins = sigBins;
	grids.push_back(grid);
      }// iG
    }// iN
    GetXSection(grids,"xsecdat/binoNLOxsec_mNScan.dat"); // get XS and RS error
    GetPDFErrorsXSection(grids,"xsecdat/xsectionPDFErrors_mNScan.dat");    // PDF errors on XS
    GetPDFErrorsAcceptance(grids,"xsecdat/acceptancePDFErrors_mNScan.dat");  // PDF errors on acceptance
  }
  else if(bino.Contains("stop-bino")) {

    int mst[29] = {110, 160, 185, 210, 235, 260, 285, 310, 335, 360, 385, 410, 460, 510, 560, 610, 660, 710, 810, 910, 1010, 1110, 1210, 1310, 1410, 1510, 1710, 2010, 5010};
    int mBino[31] = {25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 375, 425, 475, 525, 575, 625, 675, 725, 825, 925, 1025, 1125, 1225, 1325, 1425, 1525, 1725, 2025};

    for(int i = 0; i < 899; i++) {
      npoint++;
      GridPoint grid;
      grid.mStop = mst[int(i)/31];
      grid.mN = mBino[int(i)%31];
      grid.ngen = 15000;
      grid.lumi = luminosity;
      
      std::vector<BinInfo> sig_ggBins;
      std::vector<BinInfo> sig_ffBins;
      std::vector<BinInfo> sigBins;
      readSig(fSig, bino, jet, mst[int(i)/31], 2500, mBino[int(i)%31], sig_ggBins, sig_ffBins, sigBins);
      
      grid.ggBins = ggBins;
      grid.qcdBins = qcdBins;
      grid.ewBins = ewBins;
      grid.qcdSysErrors = qcdSysErrors;
      grid.sig_ggBins = sig_ggBins;
      grid.sig_ffBins = sig_ffBins;
      grid.sigBins = sigBins;
      grids.push_back(grid);
    }
    GetStopXSection(grids,"xsecdat/2012/stop.dat");
  }
  else if(bino.Contains("SMSScan")) {
    std::cout << "Create grids for " << bino.Data() << std::endl;
    for(int iN=25; iN<=2000; iN+=25) {
      for(int iG=400; iG<=2000; iG+=25) {
	npoint++;
	GridPoint grid;
	grid.mG = iG;
	grid.mS = 0;
	grid.mN = iN;
	grid.ngen = 10000;
	grid.lumi = luminosity;

	std::vector<BinInfo> sig_ggBins;
	std::vector<BinInfo> sig_ffBins;
	std::vector<BinInfo> sigBins;
	readSig(fSig,bino,jet,grid.mS,iG,iN,sig_ggBins,sig_ffBins,sigBins);

	grid.ggBins = ggBins;
	grid.qcdBins = qcdBins;
	grid.ewBins = ewBins;
	grid.qcdSysErrors = qcdSysErrors;
	grid.sig_ggBins = sig_ggBins;
	grid.sig_ffBins = sig_ffBins;
	grid.sigBins = sigBins;
	grids.push_back(grid);
      }// iG
    }// iN
    std::cout << "Done with creating grids" << std::endl;
    GetSMSXSection(grids,"xsecdat/2012/SMS/SMS_xsec.dat"); // get XS and RS error
  }
  else {
    for(int iS=400; iS<=2000; iS+=100) {
      for(int iG=420; iG<=2020; iG+=100) {
	npoint++;
	GridPoint grid;
	grid.mG = iG;
	grid.mS = iS;
	grid.mN = 375;
	grid.ngen = 10000;
	grid.lumi = luminosity;

	std::vector<BinInfo> sig_ggBins;
	std::vector<BinInfo> sig_ffBins;
	std::vector<BinInfo> sigBins;
	readSig(fSig,bino,jet,iS,iG,375,sig_ggBins,sig_ffBins,sigBins);

	grid.ggBins = ggBins;
	grid.qcdBins = qcdBins;
	grid.ewBins = ewBins;
	grid.qcdSysErrors = qcdSysErrors;
	grid.sig_ggBins = sig_ggBins;
	grid.sig_ffBins = sig_ffBins;
	grid.sigBins = sigBins;
	grids.push_back(grid);
      }// iG
    }// iS

    if(bino.Contains("bino")) {
      GetXSection(grids,"xsecdat/2012/Spectra_gsq_B_8TeV.xsec"); // get XS and RS error
      GetPDFErrorsXSection(grids,"xsecdat/2012/Spectra_gsq_B_dipho_envpdfuncert.dat");    // PDF errors on XS
      GetPDFErrorsAcceptance(grids,"xsecdat/2012/Spectra_gsq_B_dipho_envpdfuncert.dat");  // PDF errors on acceptance
    }
    else if(bino.Contains("wino")) {
      GetXSection(grids,"xsecdat/2012/Spectra_gsq_W_8TeV.xsec"); // get XS and RS error
      GetPDFErrorsXSection(grids,"xsecdat/2012/Spectra_gsq_W_dipho_envpdfuncert.dat");    // PDF errors on XS
      GetPDFErrorsAcceptance(grids,"xsecdat/2012/Spectra_gsq_W_dipho_envpdfuncert.dat");  // PDF errors on acceptance
    }

  }

  //if(bino.Contains("wino")) GetNGen(grids,"xsecdat/mcAccMap_wino_mN375.dat");

  std::cout << "calculate signal gains" << std::endl;

  makeSignalGains(grids);

  std::cout << "make data cards" << std::endl;

  makeDataCard(grids,bino,jet);
  //  testDataCard(grids,bino,jet);

  std::cout << "npoint : " << npoint << std::endl;

}



void readData(TFile* f, TString jet,
	      std::vector<BinInfo>& ggBins,
	      std::vector<BinInfo>& qcdBins,
	      std::vector<BinInfo>& ewBins,
	      std::vector<BinInfo>& qcdSysErrors) {

  TH1F* gg = (TH1F*) f->Get("pfMET_gg_"+jet);
  ggBins.clear();
  getBins(gg,ggBins);
  std::cout << "gg  events ----------" << std::endl;
  printBins(ggBins);

  //  TH1F* qcd = (TH1F*) f->Get("met_qcd_avg_"+jet);
  TH1F* qcd = (TH1F*) f->Get("pfMET_ff_"+jet);
  qcdBins.clear();
  getBins(qcd,qcdBins);
  std::cout << "qcd events ----------" << std::endl;
  printBins(qcdBins);

  TH1F* ew = (TH1F*) f->Get("pfMET_eg_"+jet);
  ewBins.clear();
  getBins(ew,ewBins);
  std::cout << "ew  events ----------" << std::endl;
  printBins(ewBins);

  TH1F* qcd_gf = (TH1F*) f->Get("pfMET_gf_"+jet);
  std::vector<BinInfo> qcd_gfBins;
  getBins(qcd_gf,qcd_gfBins);
  std::cout << "qcd_gf events ----------" << std::endl;
  printBins(qcd_gfBins);

  int N = int(qcd_gfBins.size());
  for(int i=0; i<N; i++){
    double diff = qcd_gfBins[i].y - qcdBins[i].y;
    BinInfo bin;
    bin.x = qcdBins[i].x;
    bin.y = fabs(diff);
    qcdSysErrors.push_back(bin);
  }

  std::cout << "qcd_sysErrors ---------" << std::endl;
  printBins(qcdSysErrors);

}


void readSig(TFile* f, TString bino, TString jet, int mS, int mG, int mN,
	     std::vector<BinInfo>& sig_ggBins,
	     std::vector<BinInfo>& sig_ffBins,
	     std::vector<BinInfo>& sigBins) {

  std::stringstream ggname;
  if(bino.Contains("SMSScan")) ggname << "h_gg_met_" << jet << "_mG" << mG << "_mN_" << mN;
  else if(bino.Contains("stop-bino")) ggname << "met_gg_" << jet << "_mst_" << mS << "_m1_" << mN;
  else ggname << "h_gg_met_" << jet << "_mS" << mS << "_mG" << mG << "_mN" << mN;

  TH1F* gg = (TH1F*) f->Get(ggname.str().c_str());
  if(!gg) {
    std::cout << "histogram " << ggname.str() << " is not available!!!" << std::endl;
    return;
  }
  sig_ggBins.clear();
  getBins(gg,sig_ggBins);

  std::stringstream ffname;
  if(bino.Contains("SMSScan")) ffname << "h_ff_met_" << jet << "_mG" << mG << "_mN_" << mN;
  else if(bino.Contains("stop-bino")) ffname << "met_ff_" << jet << "_mst_" << mS << "_m1_" << mN;
  else ffname << "h_ff_met_" << jet << "_mS" << mS << "_mG" << mG << "_mN" << mN;

  TH1F* ff = (TH1F*) f->Get(ffname.str().c_str());
  sig_ffBins.clear();
  getBins(ff,sig_ffBins);

  sigBins.clear();
  int n = sig_ggBins.size();
  for(int k=0; k<n; k++) {
    BinInfo b;
    b.x = sig_ggBins[k].x;
    b.y = sig_ggBins[k].y - sig_ffBins[k].y;
    if(b.y < 1e-6) b.y = 0.0;
    b.error = std::sqrt(sig_ggBins[k].error * sig_ggBins[k].error + sig_ffBins[k].error * sig_ffBins[k].error);
    sigBins.push_back(b);
  }//k
  //printBins(sigBins);

}


void getBins(TH1F* h, std::vector<BinInfo>& binInfos){

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

void printBins(std::vector<BinInfo>& binInfos){

  for(std::vector<BinInfo>::iterator it = binInfos.begin();
      it != binInfos.end(); it++) {
    printf("%3d(%6.2f +/- %5.2f) ",it->x,it->y,it->error);
  }
  printf("\n");

}

void printBins(std::fstream& of, std::vector<BinInfo>& binInfos) {
  for(std::vector<BinInfo>::iterator it = binInfos.begin();
      it != binInfos.end(); it++) {
    of << "(" << it->x << ", " << it->y << " +/- " << it->error << ") ";
  }
  of << std::endl;
}


void GetNGen(std::vector<GridPoint>& grid, TString datafile) {
  unsigned int ngrid = grid.size();
  std::ifstream fin;
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


void GetXSection(std::vector<GridPoint>& grid, TString datafile)
{
  unsigned int ngrid = grid.size();
  std::ifstream fin;
  // RS errors are absolute.   
  fin.open(datafile.Data());
      
  while(1) {
    int nGen_, mS_, mG_, mB_, mW_;
    double xsec, rsp, rsm, dummy;   
    std::string dummyS;
    //     nGen     mS     mG     mBino  mWino  LO:       LOxsec   +         LOerr+   -         LOerr-   NLO:      xsec    +       NLOerr+  -         NLOerr-
    fin >> nGen_ >> mS_ >> mG_ >> mB_ >> mW_ >> dummyS >> dummy >> dummyS >> dummy >> dummyS >> dummy >> dummyS >> xsec >> dummyS >> rsp >> dummyS >> rsm;
    if(!fin.good()) break;
    for(unsigned int ig = 0; ig < ngrid; ig++) {
      if((mG_ == grid[ig].mG) && (mS_ == grid[ig].mS)) {
        grid[ig].ngen = nGen_;
        grid[ig].xsecValue = std::abs(xsec);
        grid[ig].xsecRSErrorPos = std::abs(rsp/xsec);
        grid[ig].xsecRSErrorNeg = std::abs(rsm/xsec);
      }
    }
  }
  fin.close();
}


void GetPDFErrorsXSection(std::vector<GridPoint>& grid, TString datafile)
{
  unsigned int ngrid = grid.size();
  std::ifstream fin;
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


void GetPDFErrorsAcceptance(std::vector<GridPoint>& grid, TString datafile)
{
  unsigned int ngrid = grid.size();
  std::ifstream fin;
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


void GetSMSXSection(std::vector<GridPoint>& grid, TString datafile) {

  std::cout << "Read SMS XSection file " << datafile.Data() << " ------------------" << std::endl;

  unsigned int ngrid = grid.size();

  std::ifstream fin;
  fin.open(datafile.Data());

  int mG;
  double xsec, rsp, rsm;
  while(1){

    // SMS Scan MC : rsp and rsm contains PDF errors and RS errors
    fin >> mG >> xsec >> rsp;
    rsm = rsp;
    //    std::cout << "mG " << mG << " xsec " << xsec << " rsp " << rsp << " rsm " << rsm << std::endl;

    if(!fin.good()) break;
    if(std::abs(xsec) < 1.0e-20) continue;

    for(unsigned int ig=0; ig<ngrid; ig++){
      if(mG == grid[ig].mG) {
	grid[ig].xsecValue = std::abs(xsec);
	grid[ig].xsecRSErrorPos = std::abs(rsp);  // already relative error
	grid[ig].xsecRSErrorNeg = std::abs(rsm);  // already relative error
	grid[ig].xsecPDFError = 0;
	grid[ig].accErrorPDF = 0;
      }
    }
  }
  fin.close();
}



void makeSignalGains(std::vector<GridPoint>& grid) {

  for(std::vector<GridPoint>::iterator it = grid.begin();
      it != grid.end(); it++) {
    double nevt = 0;
    for(std::vector<BinInfo>::iterator bit = it->sigBins.begin();
	bit != it->sigBins.end(); bit++) {
      nevt += bit->y;
      double acc = bit->y * DataMCScale * DataMCScale / it->ngen;
      double nexpected = it->xsecValue * it->lumi * acc * triggerEff;
      double error = bit->error * it->xsecValue * it->lumi * DataMCScale * DataMCScale / it->ngen;
      bit->y = nexpected;
      bit->error = error;
      //      std::cout << "y, xsec, acc : " << bit->y << ", " << it->xsecValue << ", " << acc << std::endl;
    }// bit
    it->sig_sysError = std::sqrt(it->sig_sysError*it->sig_sysError + it->accErrorPDF*it->accErrorPDF + triggerEff_err*triggerEff_err);
    it->sig_sysError += 1; // relative error
    it->acc = nevt/it->ngen;
    //    printBins(it->sigBins);
  }// it

}


void makeDataCard(std::vector<GridPoint>& grid, TString bino, TString jet) {

  for(std::vector<GridPoint>::iterator it = grid.begin();
      it != grid.end(); it++) {

    if(it->sigBins.size() == 0) continue;

    std::stringstream outname;
    if(bino.Contains("stop-bino")) outname << datacard_dir.Data() << "/multiChannel/" << bino.Data() << "_mst_" << it->mStop << "_m1_" << it->mN << "_" << jet << ".dat";
    else outname << datacard_dir.Data() << "/multiChannel/" << bino.Data() << "_mS" << it->mS << "_mG" << it->mG << "_mN" << it->mN << "_" << jet << ".dat";
    std::fstream outfile(outname.str().c_str(),std::ios::out);

    outfile << "# gluino = " << it->mG << std::endl;
    outfile << "# squark = " << it->mS << std::endl;
    outfile << "# chi1 = " << it->mN << std::endl;
    outfile << "# cha1 = " << it->mW << std::endl;
    outfile << "# Xsection.NLO = " << it->xsecValue << std::endl;
    outfile << "# Luminosity = " << it->lumi << std::endl;
    outfile << "# signal.scale.uncertainty = " << it->xsecRSErrorPos << std::endl;
    outfile << "# signal.scale.uncertainty.UP = " << it->xsecRSErrorPos << std::endl;
    outfile << "# signal.scale.uncertainty.DN = " << it->xsecRSErrorNeg << std::endl;
    outfile << "# signal.PDF.uncertainty = " << it->xsecPDFError << std::endl;
    outfile << "# signal.PDFacc.uncertainty = " << it->accErrorPDF << std::endl;
    outfile << "# signal.ngen = " << it->ngen << std::endl;
    outfile << "# signal.acceptance = " << it->acc << std::endl;
    outfile << "# lumi_sysError = " << it->lumi_sysError << std::endl;
    outfile << "# qcd_sysError = " << it->qcd_sysError << std::endl;
    outfile << "# ew_sysError = " << it->ew_sysError << std::endl;
    outfile << "# sig_sysError = " << it->sig_sysError << std::endl;

    std::vector<int> sensitive_bins;
    for(int i=0; i<int(it->sigBins.size()); i++){
      if(it->sigBins[i].y > epsilon) sensitive_bins.push_back(i);
    }

    int nch = sensitive_bins.size();

    // for R_firstgues
    double d=0,b=0,s=0,cont=0,R,Rmin=99999;
    for(int i=0; i<nch; i++) {
      int bin = sensitive_bins[i];
      d += it->ggBins[bin].y;
      b += it->qcdBins[bin].y + it->ewBins[bin].y;
      s += it->sigBins[bin].y;
      cont += it->sig_ffBins[bin].y;
      double unc2= it->ggBins[bin].y;
      unc2 += pow(it->lumi_sysError-1.,2);
      unc2 += pow(it->qcd_sysError-1.,2);
      unc2 += pow(it->ew_sysError-1.,2);
      unc2 += pow(it->sig_sysError-1.,2);
      unc2 += pow(0.02,2); // JES
      unc2 += pow((it->ggBins[bin].error/it->ggBins[bin].y),2);
      unc2 += pow((it->qcdBins[bin].error/it->qcdBins[bin].y),2);
      unc2 += pow((it->ewBins[bin].error/it->ewBins[bin].y),2);
      unc2 += pow((it->sigBins[bin].error/it->sigBins[bin].y),2);
      R = 2.0*sqrt(unc2)/(it->sigBins[bin].y);
      if (R < Rmin) Rmin = R;	  
    }
    outfile << "# data = " << d << "\n";
    outfile << "# background = " << b << "\n";
    outfile << "# signal = " << s << "\n";
    outfile << "# R_firstguess = " << Rmin << "\n";
    outfile << "## data = " << d << "\n";
    outfile << "## background = " << b << "\n";
    outfile << "## signal = " << s << "\n";
    outfile << "## signal.contamination = " << cont << "\n";
    outfile << "## acceptance = " << s/(it->lumi * it->xsecValue) << "\n";
    outfile << "## ggBins = "; printBins(outfile,it->ggBins);
    outfile << "## ewBins = "; printBins(outfile,it->ewBins);
    outfile << "## qcdBins = "; printBins(outfile,it->qcdBins);
    outfile << "## qcdSysErrors = "; printBins(outfile,it->qcdSysErrors);
    outfile << "## sig_ggBins = "; printBins(outfile,it->sig_ggBins);
    outfile << "## sig_ffBins = "; printBins(outfile,it->sig_ffBins);
    outfile << "## sigBins = "; printBins(outfile,it->sigBins);
    


    outfile << "imax " << nch << " number of channels" << std::endl;
    outfile << "jmax 2 number of backgrounds" << std::endl;
    outfile << "kmax " << (5 + nch*3) << " number of nuisance parameters" << std::endl;
    outfile << "--------------" << std::endl;

    outfile << "bin         ";
    for(int i=0; i<nch; i++) outfile << "\t" << i;
    outfile << std::endl;

    outfile << "observation ";
    for(int i=0; i<nch; i++) {
      int bin = sensitive_bins[i];
      outfile << "\t" << it->ggBins[bin].y;
    }
    outfile << std::endl;

    outfile << "--------------" << std::endl;
    outfile << "bin         ";
    for(int i=0; i<nch; i++) outfile << "\t" << i << " " << i << " " << i;
    outfile << std::endl;

    outfile << "process       ";
    for(int i=0; i<nch; i++) outfile << "\t susy qcd ew";
    outfile << std::endl;

    outfile << "process       ";
    for(int i=0; i<nch; i++) outfile << "\t 0 1 2";
    outfile << std::endl;

    outfile << "rate          ";
    for(int i=0; i<nch; i++) {
      int bin = sensitive_bins[i];
      outfile << "\t" << it->sigBins[bin].y << " " << it->qcdBins[bin].y << " " << it->ewBins[bin].y;
    }
    outfile << std::endl;
    outfile << "--------------" << std::endl;

    outfile << "u_lumi lnN    ";
    for(int i=0; i<nch; i++) outfile << "\t" << it->lumi_sysError << " - - ";
    outfile << std::endl;

    outfile << "u_id lnN      ";
    for(int i=0; i<nch; i++) outfile << "\t" << it->sig_sysError << " - - ";
    outfile << std::endl;

    outfile << "u_JES lnN     ";
    for(int i=0; i<nch; i++) outfile << "\t" << "1.02 - - ";
    outfile << std::endl;

    outfile << "u_qcd lnN     ";
    for(int i=0; i<nch; i++){
      int bin = sensitive_bins[i];
      outfile << "\t" << "- " << 1+(it->qcdSysErrors[bin].y/it->qcdBins[bin].y) << " - ";
    }
    outfile << std::endl;
    
    outfile << "u_ew lnN      ";
    for(int i=0; i<nch; i++) outfile << "\t" << "- - " << it->ew_sysError;
    outfile << std::endl;
    
    for(int i=0; i<nch; i++) {
      int bin = sensitive_bins[i];
      outfile << "stat_susy_bin" << i << " lnN ";
      for(int j=0; j<nch; j++) {
	if(i==j) outfile << 1+(it->sigBins[bin].error/it->sigBins[bin].y) << " - - ";
	else outfile << " - - - ";
      }//j
      outfile << std::endl;
    }// i

    for(int i=0; i<nch; i++) {
      int bin = sensitive_bins[i];
      outfile << "stat_qcd_bin" << i << " lnN ";
      for(int j=0; j<nch; j++) {
	if(i==j) outfile << " - " << 1 + (it->qcdBins[bin].error/it->qcdBins[bin].y) << " - ";
	else outfile << " - - - ";
      }//j
      outfile << std::endl;
    }// i

    for(int i=0; i<nch; i++) {
      int bin = sensitive_bins[i];
      outfile << "stat_ew_bin" << i << " lnN ";
      for(int j=0; j<nch; j++) {
	if(i==j) outfile << " - - " << 1 + (it->ewBins[bin].error/it->ewBins[bin].y);
	else outfile << " - - - ";
      }//j
      outfile << std::endl;
    }// i

  }// it


}


void testDataCard(std::vector<GridPoint>& grid, TString bino, TString jet) {

  for(std::vector<GridPoint>::iterator it = grid.begin();
      it != grid.end(); it++) {

    if(it->mG != 1425) continue;
    if(it->mN != 1175) continue;

    if(it->sigBins.size() == 0) continue;

    std::stringstream outname;
    outname << bino.Data() << "_mS" << it->mS << "_mG" << it->mG << "_mN" << it->mN << "_" << jet << ".dat";

    std::cout << "-------------- " << outname << " ---------------" << std::endl;

    std::cout << "# gluino = " << it->mG << std::endl;
    std::cout << "# squark = " << it->mS << std::endl;
    std::cout << "# chi1 = " << it->mN << std::endl;
    std::cout << "# Xsection.NLO = " << it->xsecValue << std::endl;
    std::cout << "# Luminosity = " << it->lumi << std::endl;
    std::cout << "# signal.scale.uncertainty = " << it->xsecRSErrorPos << std::endl;
    std::cout << "# signal.scale.uncertainty.UP = " << it->xsecRSErrorPos << std::endl;
    std::cout << "# signal.scale.uncertainty.DN = " << it->xsecRSErrorNeg << std::endl;
    std::cout << "# signal.PDF.uncertainty = " << it->xsecPDFError << std::endl;
    std::cout << "# signal.PDFacc.uncertainty = " << it->accErrorPDF << std::endl;
    std::cout << "# signal.ngen = " << it->ngen << std::endl;
    std::cout << "# signal = " << it->acc << std::endl;
    std::cout << "# lumi_sysError = " << it->lumi_sysError << std::endl;
    std::cout << "# qcd_sysError = " << it->qcd_sysError << std::endl;
    std::cout << "# ew_sysError = " << it->ew_sysError << std::endl;
    std::cout << "# sig_sysError = " << it->sig_sysError << std::endl;

    std::vector<int> sensitive_bins;
    for(int i=0; i<int(it->sigBins.size()); i++){
      if(it->sigBins[i].y > epsilon) sensitive_bins.push_back(i);
    }

    int nch = sensitive_bins.size();

    // for R_firstgues
    double d=0,b=0,s=0,cont=0,R,Rmin=99999;
    for(int i=0; i<nch; i++) {
      int bin = sensitive_bins[i];
      d += it->ggBins[bin].y;
      b += it->qcdBins[bin].y + it->ewBins[bin].y;
      s += it->sigBins[bin].y;
      cont += it->sig_ffBins[bin].y;
      double unc2= it->ggBins[bin].y;
      unc2 += pow(it->lumi_sysError-1.,2);
      unc2 += pow(it->qcd_sysError-1.,2);
      unc2 += pow(it->ew_sysError-1.,2);
      unc2 += pow(it->sig_sysError-1.,2);
      unc2 += pow(0.02,2); // JES
      unc2 += pow((it->ggBins[bin].error/it->ggBins[bin].y),2);
      unc2 += pow((it->qcdBins[bin].error/it->qcdBins[bin].y),2);
      unc2 += pow((it->ewBins[bin].error/it->ewBins[bin].y),2);
      unc2 += pow((it->sigBins[bin].error/it->sigBins[bin].y),2);
      R = 2.0*sqrt(unc2)/(it->sigBins[bin].y);
      if (R < Rmin) Rmin = R;	  
    }
    std::cout << "# R_firstguess = " << Rmin << "\n";
    std::cout << "## data = " << d << "\n";
    std::cout << "## background = " << b << "\n";
    std::cout << "## signal = " << s << "\n";
    std::cout << "## signal.contamination = " << cont << "\n";
    std::cout << "## acceptance = " << s/(it->lumi * it->xsecValue) << "\n";
    
  }// it


}

void GetStopXSection(vector<GridPoint>& grid, TString datafile)
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
        grid[ig].accErrorPDF = abs(xsecErr);
      }
    }
  }
  fin.close();
}

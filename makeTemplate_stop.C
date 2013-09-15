#include "makeTemplate_stop.h"

void makeDataCard(vector<GridPoint>& grid, TString req) {

  for(vector<GridPoint>::iterator it = grid.begin();
      it != grid.end(); it++) {

    if(it->sigBins.size() == 0) continue;

    stringstream outname;
    outname << datacard_dir.Data() << "/multiChannel/stop-bino_" << it->mStop << "_m1_" << it->mBino << "_" << req << ".dat";
    fstream outfile(outname.str().c_str(), ios::out);

    outfile << "# stop = " << it->mStop << endl;
    outfile << "# bino = " << it->mBino << endl;
    outfile << "# Xsection.NLO = " << it->xsecValue << endl;
    outfile << "# Xsection.Error = " << it->xsecError << endl;
    outfile << "# Luminosity = " << it->lumi << endl;
    outfile << "# signal.ngen = " << it->ngen << endl;
    outfile << "# signal.acceptance = " << it->acc << endl;
    outfile << "# lumi_sysError = " << it->lumi_sysError << endl;
    outfile << "# qcd_sysError = " << it->qcd_sysError << endl;
    outfile << "# ew_sysError = " << it->ew_sysError << endl;
    outfile << "# sig_sysError = " << it->sig_sysError << endl;

    vector<int> sensitive_bins;
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
    outfile << "## ggBins = "; printBins(outfile, it->ggBins);
    outfile << "## ewBins = "; printBins(outfile, it->ewBins);
    outfile << "## qcdBins = "; printBins(outfile, it->qcdBins);
    outfile << "## qcdSysErrors = "; printBins(outfile, it->qcdSysErrors);
    outfile << "## sig_ggBins = "; printBins(outfile, it->sig_ggBins);
    outfile << "## sig_gfBins = "; printBins(outfile, it->sig_gfBins);
    outfile << "## sig_ffBins = "; printBins(outfile, it->sig_ffBins);
    outfile << "## sigBins = "; printBins(outfile, it->sigBins);
    
    outfile << "imax " << nch << " number of channels" << endl;
    outfile << "jmax 2 number of backgrounds" << endl;
    outfile << "kmax " << (5 + nch*3) << " number of nuisance parameters" << endl;
    outfile << "--------------" << endl;

    outfile << "bin         ";
    for(int i=0; i<nch; i++) outfile << "\t" << i;
    outfile << endl;

    outfile << "observation ";
    for(int i=0; i<nch; i++) {
      int bin = sensitive_bins[i];
      outfile << "\t" << it->ggBins[bin].y;
    }
    outfile << endl;

    outfile << "--------------" << endl;
    outfile << "bin         ";
    for(int i=0; i<nch; i++) outfile << "\t" << i << " " << i << " " << i;
    outfile << endl;

    outfile << "process       ";
    for(int i=0; i<nch; i++) outfile << "\t susy qcd ew";
    outfile << endl;

    outfile << "process       ";
    for(int i=0; i<nch; i++) outfile << "\t 0 1 2";
    outfile << endl;

    outfile << "rate          ";
    for(int i=0; i<nch; i++) {
      int bin = sensitive_bins[i];
      outfile << "\t" << it->sigBins[bin].y << " " << it->qcdBins[bin].y << " " << it->ewBins[bin].y;
    }
    outfile << endl;
    outfile << "--------------" << endl;

    outfile << "u_lumi lnN    ";
    for(int i=0; i<nch; i++) outfile << "\t" << it->lumi_sysError << " - - ";
    outfile << endl;

    outfile << "u_id lnN      ";
    for(int i=0; i<nch; i++) outfile << "\t" << it->sig_sysError << " - - ";
    outfile << endl;

    outfile << "u_JES lnN     ";
    for(int i=0; i<nch; i++) outfile << "\t" << "1.02 - - ";
    outfile << endl;

    outfile << "u_qcd lnN     ";
    for(int i=0; i<nch; i++){
      int bin = sensitive_bins[i];
      outfile << "\t" << "- " << 1+(it->qcdSysErrors[bin].y/it->qcdBins[bin].y) << " - ";
    }
    outfile << endl;
    
    outfile << "u_ew lnN      ";
    for(int i=0; i<nch; i++) outfile << "\t" << "- - " << it->ew_sysError;
    outfile << endl;
    
    for(int i=0; i<nch; i++) {
      int bin = sensitive_bins[i];
      outfile << "stat_susy_bin" << i << " lnN ";
      for(int j=0; j<nch; j++) {
	if(i==j) outfile << 1+(it->sigBins[bin].error/it->sigBins[bin].y) << " - - ";
	else outfile << " - - - ";
      }//j
      outfile << endl;
    }// i

    for(int i=0; i<nch; i++) {
      int bin = sensitive_bins[i];
      outfile << "stat_qcd_bin" << i << " lnN ";
      for(int j=0; j<nch; j++) {
	if(i==j) outfile << " - " << 1 + (it->qcdBins[bin].error/it->qcdBins[bin].y) << " - ";
	else outfile << " - - - ";
      }//j
      outfile << endl;
    }// i

    for(int i=0; i<nch; i++) {
      int bin = sensitive_bins[i];
      outfile << "stat_ew_bin" << i << " lnN ";
      for(int j=0; j<nch; j++) {
	if(i==j) outfile << " - - " << 1 + (it->ewBins[bin].error/it->ewBins[bin].y);
	else outfile << " - - - ";
      }//j
      outfile << endl;
    }// i

  }// it


}

void makeTemplate_stop(TString req="bj") {

  TString hist_dir = "inputHists";
  TString dataHist = hist_dir+"/"+"met_reweighted_"+req+".root";
  TString sigHist = hist_dir+"/"+"contamination_stop.root";

  TFile* fData = new TFile(dataHist, "READ");
  vector<BinInfo> ggBins;
  vector<BinInfo> qcdBins;
  vector<BinInfo> ewBins;
  vector<BinInfo> qcdSysErrors;
  readData(fData, req, ggBins, qcdBins, ewBins, qcdSysErrors);

  TFile* fSig = new TFile(sigHist, "READ");

  vector<GridPoint> grids;
  int npoint = 0;

  int mst[29] = {110, 160, 185, 210, 235, 260, 285, 310, 335, 360, 385, 410, 460, 510, 560, 610, 660, 710, 810, 910, 1010, 1110, 1210, 1310, 1410, 1510, 1710, 2010, 5010};
  int mBino[31] = {25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 375, 425, 475, 525, 575, 625, 675, 725, 825, 925, 1025, 1125, 1225, 1325, 1425, 1525, 1725, 2025};

  for(int i = 0; i < 899; i++) {
    npoint++;
    GridPoint grid;
    grid.mStop = mst[int(i)/31];
    grid.mBino = mBino[int(i)%31];
    grid.ngen = 15000;
    grid.lumi = luminosity;

    vector<BinInfo> sig_ggBins;
    vector<BinInfo> sig_ffBins;
    vector<BinInfo> sig_gfBins;
    vector<BinInfo> sigBins;
    readSig(fSig, req, mst[int(i)/31], mBino[int(i)%31], sig_ggBins, sig_ffBins, sig_gfBins, sigBins);

    grid.ggBins = ggBins;
    grid.qcdBins = qcdBins;
    grid.ewBins = ewBins;
    grid.qcdSysErrors = qcdSysErrors;
    grid.sig_ggBins = sig_ggBins;
    grid.sig_gfBins = sig_gfBins;
    grid.sigBins = sigBins;
    grids.push_back(grid);
  }
  GetXSection(grids, "xsecdat/2012/stop.dat");

  cout << "calculate signal gains" << endl;

  makeSignalGains(grids);

  cout << "make data cards" << endl;

  makeDataCard(grids, req);

  cout << "npoint : " << npoint << endl;

}


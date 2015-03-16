#include "makeTemplate.h"

void makeDataCard(vector<GridPoint>& grid, TString bino, TString jet) {

  for(vector<GridPoint>::iterator it = grid.begin();
      it != grid.end(); it++) {

    if(it->sigBins.size() == 0) continue;

    stringstream outname;
    outname << datacard_dir.Data() << "/multiChannel/" << bino.Data() << "_mS" << it->mS << "_mG" << it->mG << "_mN" << it->mN << "_" << jet << ".dat";
    fstream outfile(outname.str().c_str(),ios::out);
    
    outfile << "# gluino = " << it->mG << endl;
    outfile << "# squark = " << it->mS << endl;
    outfile << "# chi1 = " << it->mN << endl;
    outfile << "# cha1 = " << it->mW << endl;
    outfile << "# Xsection.NLO = " << it->xsecValue << endl;
    outfile << "# Luminosity = " << it->lumi << endl;
    outfile << "# signal.scale.uncertainty = " << it->xsecRSErrorPos << endl;
    outfile << "# signal.scale.uncertainty.UP = " << it->xsecRSErrorPos << endl;
    outfile << "# signal.scale.uncertainty.DN = " << it->xsecRSErrorNeg << endl;
    outfile << "# signal.PDF.uncertainty = " << it->xsecPDFError << endl;
    outfile << "# signal.PDFacc.uncertainty = " << it->accErrorPDF << endl;
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
    outfile << "## ggBins = "; printBins(outfile,it->ggBins);
    outfile << "## ewBins = "; printBins(outfile,it->ewBins);
    outfile << "## qcdBins = "; printBins(outfile,it->qcdBins);
    outfile << "## qcdSysErrors = "; printBins(outfile,it->qcdSysErrors);
    outfile << "## sig_ggBins = "; printBins(outfile,it->sig_ggBins);
    outfile << "## sig_ffBins = "; printBins(outfile,it->sig_ffBins);
    outfile << "## sig_gfBins = "; printBins(outfile, it->sig_gfBins);
    outfile << "## sigBins = "; printBins(outfile,it->sigBins);
    
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

void makeTemplate(TString bino, TString jet) {

  TString hist_dir = work_dir+"inputHists";
  TString dataHist = hist_dir + "/" + "met_reweighted_"+jet+".root";
  TString sigHist = hist_dir + "/" + "signal_contamination_" + bino + "_chi0375.root";
  if(bino.Contains("mNScan")) sigHist = hist_dir + "/" + "signal_contamination_bino_chi0.root";
  if(bino.Contains("SMSScan")) sigHist = hist_dir + "/" + "signal_contamination_T5gg.root";

  TFile* fData = new TFile(dataHist,"READ");
  vector<BinInfo> ggBins;
  vector<BinInfo> qcdBins;
  vector<BinInfo> ewBins;
  vector<BinInfo> qcdSysErrors;
  readData(fData, jet, ggBins, qcdBins, ewBins, qcdSysErrors);

  TFile* fSig = new TFile(sigHist,"READ");

  vector<GridPoint> grids;

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

	vector<BinInfo> sig_ggBins;
	vector<BinInfo> sig_ffBins;
	 vector<BinInfo> sig_gfBins;
	vector<BinInfo> sigBins;
	readSig(fSig, bino, jet, 2500, iG, iN, sig_ggBins, sig_ffBins, sig_gfBins, sigBins);

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
  else if(bino.Contains("SMSScan")) {
    cout << "Create grids for " << bino.Data() << endl;
    for(int iN=25; iN<=2000; iN+=25) {
      for(int iG=400; iG<=2000; iG+=25) {
	npoint++;
	GridPoint grid;
	grid.mG = iG;
	grid.mS = 0;
	grid.mN = iN;
	grid.ngen = 10000;
	grid.lumi = luminosity;

	vector<BinInfo> sig_ggBins;
	vector<BinInfo> sig_ffBins;
	vector<BinInfo> sig_gfBins;
	vector<BinInfo> sigBins;
	readSig(fSig, bino, jet, grid.mS, iG, iN, sig_ggBins, sig_ffBins, sig_gfBins, sigBins);

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
    cout << "Done with creating grids" << endl;
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

	vector<BinInfo> sig_ggBins;
	vector<BinInfo> sig_gfBins;
	vector<BinInfo> sig_ffBins;
	vector<BinInfo> sigBins;
	readSig(fSig, bino, jet, iS, iG, 375, sig_ggBins, sig_ffBins, sig_gfBins, sigBins);

	grid.ggBins = ggBins;
	grid.qcdBins = qcdBins;
	grid.ewBins = ewBins;
	grid.qcdSysErrors = qcdSysErrors;
	grid.sig_ggBins = sig_ggBins;
	grid.sig_ffBins = sig_ffBins;
	grid.sig_gfBins = sig_gfBins;
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

  cout << "calculate signal gains" << endl;

  makeSignalGains(grids);

  cout << "make data cards" << endl;

  makeDataCard(grids,bino,jet);
  //  testDataCard(grids,bino,jet);

  cout << "npoint : " << npoint << endl;

}

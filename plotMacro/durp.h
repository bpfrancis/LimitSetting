#include <iostream>
#include <fstream>
#include <vector>

#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphSmooth.h>
#include <TLatex.h>
#include <TLegend.h>

#include "util.h"

using namespace std;

const int nX = 16;
const int nY = 16;

Double_t mst[nX] = {235, 260, 285, 310, 335, 360, 385, 410, 460, 510, 560, 610, 660, 710, 810, 910};
Double_t mBino[nY] = {150, 175, 200, 225, 250, 275, 300, 325, 375, 425, 475, 525, 575, 625, 675, 725};

const int nxs = 3;
TString xsname[nxs] = {"xsec", "xsec_obs","xsec_exp"};
  
const int nlimit = 8;
TString limitname[nlimit] = {"obs", "obs_1L", "obs_1H", "exp", "exp_1L","exp_1H", "obs_up3", "obs_down3"};

class PlotMaker {

 public:
  PlotMaker();
  virtual ~PlotMaker() {
    h_xs.clear();
    h_limit.clear();
    curv.clear();
    curvS.clear();
  }

  void InitStyle();
  void FillBin(double mstop, double mbino,
	       double xsec, double xsecError,
	       double obsLimit,
	       double expLimit, double exp_m1s, double exp_m2s, double exp_p1s, double exp_p2s);

  void SetRange(double xlo, double xhi, double ylo, double yhi) {
    xMin = xlo;
    xMax = xhi;
    yMin = ylo;
    yMax = yhi;
  }

  void SetAxisTitles(TString x, TString y) {
    xLabel = x;
    yLabel = y;
  }

  void FillPotHoles(TH2D*& h);
  void FillPotHoles() {
    for(unsigned int i = 0; i < h_xs.size(); i++) FillPotHoles(h_xs[i]);
    for(unsigned int i = 0; i < h_limit.size(); i++) FillPotHoles(h_limit[i]);
  }

  void getMinMaxValues(TH2D *h, double& minVal, double& maxVal);

  void DrawCrossSection();
  void DrawUpperLimit();
  
  void GetContours();
  void SmoothContours();

  void DrawExclusion();
  void DrawExclusionOnLimit();

  void Save(TString output_dir);

 private:

  double xMin;
  double xMax;
  double yMin;
  double yMax;

  TString scan;

  TString xLabel;
  TString yLabel;

  TString option2D;
  int legendFillColor;

  TLatex * lat;
  TLatex * lat2;

  TGraph * diagonalRegion;
  TLatex * lat_diagonal;
  
  vector<TH2D*> h_xs;
  vector<TH2D*> h_limit;
  vector<TGraph*> curv;
  vector<TGraph*> curvS;

  TH2D * h_back;

  TGraph * excludedRegion;
  TGraph * exp1sigma_aroundExp;

};

PlotMaker::PlotMaker() {

  h_xs.clear();
  h_limit.clear();
  curv.clear();
  curvS.clear();

  Double_t xBins[nX+1];
  xBins[0] = 222.5;
  for(int i = 1; i < nX; i++) xBins[i] = (mst[i] + mst[i-1])/2.;
  xBins[nX] = 960;

  Double_t yBins[nY+1];
  yBins[0] = 137.5;
  for(int i = 1; i < nY; i++) yBins[i] = (mBino[i] + mBino[i-1])/2.;
  yBins[nY] = 775;

  for(int i = 0; i < nxs; i++) h_xs.push_back(new TH2D(xsname[i], xsname[i], nX, xBins, nY, yBins));
  for(int i = 0; i < nlimit; i++) h_limit.push_back(new TH2D(limitname[i], limitname[i], nX, xBins, nY, yBins));

  lat = new TLatex(0.18, 0.92, "#bf{CMS Preliminary}  #sqrt{s} = 8 TeV");
  lat->SetNDC(true);
  lat->SetTextFont(43);
  lat->SetTextSize(21/*25*/);

  lat2 = new TLatex(0.49, 0.92, "       L_{int} = 19.7 fb^{-1},   #geq 3 jets, #geq 1 btag, #geq 2 #gamma");
  lat2->SetNDC(true);
  lat2->SetTextFont(43);
  lat2->SetTextSize(25);

  diagonalRegion = new TGraph(3);
  diagonalRegion->SetPoint(0, TMath::Min(xMin, yMin), TMath::Min(xMin, yMin));
  diagonalRegion->SetPoint(1, TMath::Min(xMin, yMin), yBins[nY]);
  diagonalRegion->SetPoint(2, TMath::Max(xBins[nX], yBins[nY]), TMath::Max(xBins[nX], yBins[nY]));
  diagonalRegion->SetFillColor(legendFillColor);

  lat_diagonal = new TLatex(0.6, 0.85, "#tilde{t} NLSP");
  lat_diagonal->SetNDC(true);
  lat_diagonal->SetTextSize(0.04);

  h_back = new TH2D("h_back_"+scan, ";"+xLabel+";"+yLabel, 100, xMin, xMax, 100, yMin, yMax);

}
  
  

void PlotMaker::InitStyle() {

  gROOT->ForceStyle();

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  //For the temperature plots
  gStyle->SetPadRightMargin(0.2);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleOffset(1.4, "xz");
  gStyle->SetTitleOffset(1.9, "y");

  gStyle->SetNdivisions(505);
  gStyle->SetTitleFont(43, "xyz");
  gStyle->SetTitleSize(32, "xyz");
  gStyle->SetLabelFont(42, "xyz");
  gStyle->SetLabelSize(0.04, "xyz");
}

void PlotMaker::FillBin(double mstop, double mbino,
			double xsec, double xsecError,
			double obsLimit,
			double expLimit, double exp_m1s, double exp_m2s, double exp_p1s, double exp_p2s) {

  if(h_xs.size() < 3 || h_limit.size() < 8) {
    cout << "Requested filling something that doesn't exist!" << endl;
    continue;
  }

  double oneSigma_L = xsecError / 100.;
  double oneSigma_H = xsecError / 100.;

  h_xs[0]->Fill(mstop, mbino, xsec);
  h_xs[1]->Fill(mstop,mbino,obsLimit*xsec);
  h_xs[2]->Fill(mstop,mbino,expLimit*xsec);
   
  h_limit[0]->Fill(mstop, mbino, obsLimit);
  h_limit[1]->Fill(mstop, mbino, obsLimit * (1 - oneSigma_L));
  h_limit[2]->Fill(mstop, mbino, obsLimit * (1 + oneSigma_H));
  h_limit[3]->Fill(mstop, mbino, expLimit);
  h_limit[4]->Fill(mstop, mbino, exp_m1s);
  h_limit[5]->Fill(mstop, mbino, exp_p1s);
  h_limit[6]->Fill(mstop, mbino, obsLimit / 3.);
  h_limit[7]->Fill(mstop, mbino, obsLimit * 3.);
}

void PlotMaker::FillPotHoles(TH2D*& h) {

  int nbinsX = h->GetXaxis()->GetNbins();
  int nbinsY = h->GetYaxis()->GetNbins();

  double epsilon = 1e-10;

  for(int ix = 1; ix <= nbinsX; ix++) {
    for(int iy=1; iy <= nbinsY; iy++) {

      double xx = h->GetXaxis()->GetBinCenter(ix);
      double yy = h->GetYaxis()->GetBinCenter(iy) + 100;
      if(xx < yy) continue;

      double val = h->GetBinContent(ix, iy);
      if(TMath::IsNaN(val)) h->SetBinContent(ix, iy, 0); // checking for NAN
    }
  }

  for(int ix = 1; ix <= nbinsX; ix++) {
    for(int iy = 1; iy <= nbinsY; iy++) {

      double xx = h->GetXaxis()->GetBinCenter(ix);
      double yy = h->GetYaxis()->GetBinCenter(iy) + 100;
      if(xx < yy) continue;

      double val = h->GetBinContent(ix, iy);
      if(fabs(val) > epsilon) continue;

      int ncnt = 0;
      double sum = 0;
      double sumErr = 0;

      double up = h->GetBinContent(ix, iy+1);
      if(fabs(up) > epsilon && iy < nbinsY){
	sum += up;
	sumErr += h->GetBinError(ix, iy+1)*h->GetBinError(ix, iy+1);
	ncnt++;
      }

      double down = h->GetBinContent(ix, iy-1);
      if(fabs(down) > epsilon && iy > 1){
	sum += down;
	sumErr += h->GetBinError(ix, iy-1)*h->GetBinError(ix, iy-1);
	ncnt++;
      }

      double left = h->GetBinContent(ix-1,iy);
      if(fabs(left) > epsilon && ix > 1){
	sum += left;
	sumErr += h->GetBinError(ix-1, iy)*h->GetBinError(ix-1, iy);
	ncnt++;
      }

      double right = h->GetBinContent(ix+1,iy);
      if(fabs(right) > epsilon && ix < nbinsX){
	sum += right;
	sumErr += h->GetBinError(ix+1, iy)*h->GetBinError(ix+1, iy);
	ncnt++;
      }

      if(ncnt > 0) {
	h->SetBinContent(ix, iy, sum/ncnt);
	h->SetBinError(ix, iy, sqrt(sumErr)/ncnt);
      }

    } // for iy
  } // for ix

}

void PlotMaker::getMinMaxValues(TH2D *h, double& minVal, double& maxVal) {

  maxVal = 0;
  minVal = 9999999;

  int nbinsX = h->GetXaxis()->GetNbins();
  int nbinsY = h->GetYaxis()->GetNbins();

  for(int ix=1; ix <= nbinsX; ix++) {
    for(int iy=1; iy <= nbinsY; iy++) {
      double xx = h->GetXaxis()->GetBinCenter(ix);
      double yy = h->GetYaxis()->GetBinCenter(iy) + 100;
      if(xx < yy) continue;
    }
    
    double val = h->GetBinContent(ix,iy);
    if(TMath::IsNaN(val)) continue;
    if(val > maxVal) maxVal = val;
    if(val < minVal) minVal = val;
  }

}


void PlotMaker::DrawCrossSection() {

  TCanvas * can_xs = new TCanvas("can_xsec_"+scan, "can_xsec_"+scan, 900, 800);
  can_xs->SetLogz();

  double maxVal = 0;
  double minVal = 9999999;
  getMinMaxValues(h_xs[0], minVal, maxVal);
  h_xs[0]->SetMaximum(1.1 * maxVal);
  h_xs[0]->SetMinimum(0.9 * minVal);

  h_xs[0]->SetTitle(";" + xLabel + ";" + yLabel + ";Cross Section [pb]");
  h_xs[0]->Draw(option2D);

  lat->Draw("same");
  lat2->Draw("same");

  can_xs->RedrawAxis();

  can_xs->Print("", ".gif");
  can_xs->Print("", ".pdf");

  delete can_xs;
}

void PlotMaker::DrawUpperLimit() {

  TCanvas* can_limit = new TCanvas("can_limit_"+scan,"can_limit_"+scan,900,800);
  can_limit->SetLogz();

  double maxVal = 0;
  double minVal = 9999999;
  getMinMaxValues(h_xs[1], minVal, maxVal);
  h_xs[1]->SetMaximum(1.1 * maxVal);
  h_xs[1]->SetMinimum(0.9 * minVal);

  h_xs[1]->SetTitle(";" + xLabel + ";" + yLabel + ";95% CL cross section upper limit [pb]");
  h_xs[1]->Draw(option2D);

  lat->Draw("same");
  lat2->Draw("same");

  can_limit->RedrawAxis();

  can_limit->Print("",".gif");
  can_limit->Print("",".pdf");

  delete can_limit;
}

void PlotMaker::GetContours() {

  cout << "Now get contours..." << endl;

  double contours[2] = {0., 1.};

  TCanvas * can_excl01 = new TCanvas("can_contour_"+scan, "can_contour_"+scan, 1200, 800);

  can_excl01->Divide(nlimit/2,2);

  for(unsigned int i = 0; i < h_limit.size(); i++) {
    can_excl01->cd(i+1);
    h_back->Draw();

    h_limit[i]->SetContour(2, contours);
    h_limit[i]->Draw("SAME CONT LIST");
    gPad->Update();
    
    TObjArray * contsM = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    TList * contLevel = (TList*)contsM->At(0);
    curv.push_back((TGraph*)contLevel->First()->Clone("exclusion_contour_"+limitname[i]));

    //for(int iDurp = 0; iDurp < 16; iDurp++) curv[i]->RemovePoint(curv[i]->GetN() - 1);

  }// for i

  can_excl01->Print("", ".gif");
  can_excl01->Print("", ".pdf");

  delete can_excl01;
}

void PlotMaker::SmoothContours() {

  cout << "Smoothing..." << endl;

  TGraphSmooth* gs[nlimit];

  for(int i = 0; i < nlimit; i++) {
    gs[i] = new TGraphSmooth("normal");
    curvS.push_back(gs[i]->SmoothSuper(curv[i]));
    curvS[i]->SetName("exclusion_smooth_contour_"+limitname[i]);
  }

}

void PlotMaker::SetExclusion() {

  excludedRegion = new TGraph(curvS[0]->GetN()+3);
  int nbins = curvS[0]->GetN();

  for(int i = 0; i < nbins; i++){
    double x, y;
    curvS[0]->GetPoint(i,x,y);
    excludedRegion->SetPoint(i,x,y);
  }

  excludedRegion->SetPoint(nbins, xMin, yMax);
  excludedRegion->SetPoint(nbins+1, xMin, yMin);
  excludedRegion->SetPoint(nbins+2, xMax, yMin);

  excludedRegion->SetFillColor(kBlue-10);
  excludedRegion->SetFillStyle(1001);

  // experimental 1 sigma error around expected limit
  curvS[4]->SetLineStyle(3);
  curvS[4]->SetLineWidth(3);
  curvS[4]->SetLineColor(kOrange+9);

  curvS[5]->SetLineStyle(3);
  curvS[5]->SetLineWidth(3);
  curvS[5]->SetLineColor(kOrange+9);

  // expected limit
  curvS[3]->SetLineStyle(9);
  curvS[3]->SetLineWidth(3);
  curvS[3]->SetLineColor(kOrange+9);

  // theory 1 sigma around observed limit
  curvS[1]->SetLineStyle(3);
  curvS[1]->SetLineWidth(2);
  curvS[1]->SetLineColor(4);

  curvS[2]->SetLineStyle(3);
  curvS[2]->SetLineWidth(2);
  curvS[2]->SetLineColor(4);

  // observed limit
  curvS[0]->SetLineWidth(3);
  curvS[0]->SetLineColor(4);

  // observed limit with xsec scaled up by 3
  curvS[6]->SetLineStyle(9);
  curvS[6]->SetLineWidth(2);
  curvS[6]->SetLineColor(kBlack);

  // observed limit with xsec scaled down by 3
  curvS[7]->SetLineStyle(2);
  curvS[7]->SetLineWidth(2);
  curvS[7]->SetLineColor(kBlack);

  // experimental 1 sigma band around expected limit

  int n1 = curvS[4]->GetN();
  int n2 = curvS[5]->GetN();

  exp1sigma_aroundExp = new TGraph(n1+n2);

  for(int i = 0; i < n1; i++){
    double x, y;
    curvS[4]->GetPoint(i, x, y);
    exp1sigma_aroundExp->SetPoint(i, x, y);
  }
  for(int i = 0; i < n2; i++){
    double x, y;
    curvS[5]->GetPoint(n2-1-i, x, y);
    exp1sigma_aroundExp->SetPoint(n1+i, x, y);
  }
  
  exp1sigma_aroundExp->SetFillColor(kOrange-3);
  exp1sigma_aroundExp->SetFillStyle(1001);

}

void PlotMaker::DrawExclusion() {

  TCanvas * can_excl02 = new TCanvas("can_exclusion_"+scan, "can_exclusion_"+scan,900,800);
  h_back->Draw();
  can_excl02->SetRightMargin(0.08);

  //excludedRegion->Draw("SAME F");

  exp1sigma_aroundExp->Draw("SAME F");
  curvS[3]->Draw("SAME L");
  curvS[1]->Draw("SAME L");
  curvS[2]->Draw("SAME L");
  curvS[0]->Draw("SAME L");
  lat->Draw("same");
  lat2->Draw("same");

  diagonalRegion->Draw("same f");

  double leg_xmin = 0.58 - 0.35;
  double leg_xmax = 0.9 - 0.35;
  double leg_ymin = 0.64;
  double leg_ymax = 0.87;

  TString legendTitle = "pp #rightarrow #tilde{t}#tilde{t},  #tilde{t} #rightarrow t + #tilde{#chi}^{0}_{2},  #tilde{#chi}^{0}_{2} #rightarrow #gamma + #tilde{#chi}^{0}_{1}";

  TLegend * leg = new TLegend(leg_xmin, leg_ymin, leg_xmax, leg_ymax, legendTitle, "brNDC");
  leg->SetFillColor(legendFillColor);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->AddEntry(curvS[0],"Observed","L");
  leg->AddEntry(curvS[1],"Observed #pm1#sigma theory","L");
  TGraph * legGraph = (TGraph*) exp1sigma_aroundExp->Clone();
  legGraph->SetLineColor(kOrange+9);
  legGraph->SetLineStyle(9);
  legGraph->SetLineWidth(3);
  leg->AddEntry(legGraph,"Expected #pm1#sigma exp.","FL");
  leg->Draw("same");

  double xv = 0.25;
  double yv = 0.25;
  
  TLatex * lat4 = new TLatex(xv, yv, "Excluded");
  lat4->SetNDC(true);
  lat4->SetTextFont(43);
  lat4->SetTextSize(30);
  lat4->Draw("same");

  can_excl02->RedrawAxis();

  can_excl02->Print("",".gif");
  can_excl02->Print("",".pdf");

  delete can_excl02;
}

void PlotMaker::DrawExclusionOnLimit() {

  TCanvas * can_exclusionOnLimit = new TCanvas("can_exclusionOnLimit_"+scan,"can_exclusionOnLimit_"+scan,900,800);
  can_exclusionOnLimit->SetLogz();
  h_xs[1]->SetTitle(";M_{gluino} [GeV];M_{Neutralino} [GeV];95% CL cross section upper limit [pb]");
  h_xs[1]->SetMinimum(0);
  //h_xs[1]->GetZaxis()->SetRangeUser(9.e-4, 2.1e-2);
  //h_xs[1]->GetZaxis()->SetNdivisions(210);
  h_xs[1]->Draw(option2D);

  curvS[0]->SetLineColor(kBlack);
  
  curvS[6]->Draw("SAME L");
  curvS[7]->Draw("SAME L");
  curvS[0]->Draw("SAME L");
  
  diagonalRegion->SetFillColor(0);
  diagonalRegion->Draw("same f");
  
  TLegend * leg2 = new TLegend(leg_xmin,leg_ymin,leg_xmax - 0.05,leg_ymax,legendTitle,"brNDC");
  //leg2->SetFillColor(legendFillColor);
  leg2->SetFillColor(0);
  leg2->SetLineColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(62);
  leg2->SetTextSize(0.03);
  //leg2->AddEntry("NULL","NLO+NLL Limits","h");
  leg2->AddEntry("NULL", "m(#tilde{q}) >> m(#tilde{g})", "h");
  leg2->AddEntry(curvS[0],"#sigma^{NLO-QCD}","L");
  leg2->AddEntry(curvS[6],"3#times#sigma^{NLO-QCD}","L");
  leg2->AddEntry(curvS[7],"1/3#times#sigma^{NLO-QCD}","L");
  leg2->Draw("same");
  
  lat->Draw("same");
  lat2->Draw("same");
  
  can_exclusionOnLimit->RedrawAxis();
  
  can_exclusionOnLimit->Print("",".gif");
  can_exclusionOnLimit->Print("",".pdf");

  delete can_exclusionOnLimit;

}

void PlotMaker::Save(TString output_dir) {

  TFile* fOut = new TFile(output_dir + "/hist_sms_gg_1jet_output.root", "RECREATE");
  fOut->cd();

  h_xs[1]->Write("h2_95CL_limit");
  exp1sigma_aroundExp->Write("contour_exp_1s_band");
  curvS[3]->Write("contour_exp");
  curvS[4]->Write("contour_exp_1s_down");
  curvS[5]->Write("contour_exp_1s_up");
  curvS[1]->Write("contour_theory_obs_1s_down");
  curvS[2]->Write("contour_theory_obs_1s_up");
  curvS[0]->Write("contour_obs");
  curvS[6]->Write("contour_obs_up3");
  curvS[7]->Write("contour_obs_down3");
  h_xs[0]->Write("xsec");
  h_xs[1]->Write("observed_limit_in_xsecUnit");
  h_xs[2]->Write("expected_limit_in_xsecUnit");
  h_limit[0]->Write("observed_limit_in_R");
  h_limit[1]->Write("observed_limit_1s_down_in_R");
  h_limit[2]->Write("observed_limit_1s_up_in_R");
  h_limit[3]->Write("expected_limit_in_R");
  h_limit[4]->Write("expected_limit_1s_down_in_R");
  h_limit[5]->Write("expected_limit_1s_up_in_R");

  //can_exclusionOnLimit->Write();

  for(int i = 0; i < nlimit; i++){
    curv[i]->Write();
    curvS[i]->Write();
  }

  fOut->Write();
  fOut->Close();
  
}

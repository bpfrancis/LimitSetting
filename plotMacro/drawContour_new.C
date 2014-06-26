#include <iostream>
#include <fstream>

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

void drawContour_stop(TString scan="stop-bino", bool print=false) {

  TString data_dir = "table";
  TString output_dir = "hist";

  gSystem->mkdir(output_dir,true);

  bool useSmooth = false;
  bool useCustomGetContour = false;

  int diagonal = 1;
  int legendFillColor = 16;

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


  //  TString option2D = "CONT4 Z";
  TString option2D = "COL Z";

  int xMin = 0;
  int xMax = 5010;
  int yMin = 0;
  int yMax = 2025;

  TString xLabel = "m_{Stop} [GeV]";
  TString yLabel = "m_{Bino} [GeV]";

  Double_t mst[29] = {110, 160, 185, 210, 235, 260, 285, 310, 335, 360, 385, 410, 460, 510, 560, 610, 660, 710, 810, 910, 1010, 1110, 1210, 1310, 1410, 1510, 1710, 2010, 5010};
  Double_t mBino[31] = {25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 375, 425, 475, 525, 575, 625, 675, 725, 825, 925, 1025, 1125, 1225, 1325, 1425, 1525, 1725, 2025};

  const int nX = 30;
  const int nY = 32;

  Double_t xBins[nX+1];
  xBins[0] = 0;
  xBins[1] = 55;
  for(int i = 1; i < 29; i++) xBins[i+1] = (mst[i] + mst[i-1])/2.;
  xBins[30] = 6510;

  Double_t yBins[nY+1];
  yBins[0] = 0;
  yBins[1] = 12.5;
  for(int i = 1; i < 31; i++) yBins[i+1] = (mBino[i] + mBino[i-1])/2.;
  yBins[32] = 2175;

  TFile* fout = new TFile(output_dir+ "/hist_exclusion_"+scan+".root","RECREATE");

  const int nxs = 2;
  TString xsname[nxs] = {"xsec_obs","xsec_exp"};
  TH2D* h_xs[nxs];
  for(int i=0; i<nxs; i++){
    h_xs[i] = new TH2D(xsname[i],xsname[i],nX,xBins,nY,yBins);
  }

  const int nlimit = 6;
  TString limitname[nlimit] = {"obs", "obs_1L", "obs_1H", "exp", "exp_1L","exp_1H"};
  TH2D* h_limit[nlimit];
  for(int i=0; i<nlimit; i++){
    h_limit[i] = new TH2D(limitname[i],limitname[i],nX,xBins,nY,yBins);
  }

  TString datafile = data_dir + "/" + scan + ".table";

  std::ifstream fin;
  fin.open(datafile.Data());

  while(1){

    //mStop mBino xsec xsecError obsLimit expLimit exp_m1s exp_m2s exp_p1s exp_p2s

    int ms, mb;
    double xsec, xsecError, obsLimit, expLimit, exp_m1s, exp_m2s, exp_p1s, exp_p2s;
    fin >> ms >> mb >> xsec >> xsecError >> obsLimit >> expLimit >> exp_m1s >> exp_m2s >> exp_p1s >> exp_p2s;
    if(!fin.good()) break;

    double oneSigma_L = xsecError;
    double oneSigma_H = xsecError;

    double xx = ms;
    double yy = mb;

    h_xs[0]->Fill(xx,yy,obsLimit*xsec);
    h_xs[1]->Fill(xx,yy,expLimit*xsec);

    h_limit[0]->Fill(xx,yy,obsLimit);
    h_limit[1]->Fill(xx,yy,obsLimit*(1-oneSigma_L));
    h_limit[2]->Fill(xx,yy,obsLimit*(1+oneSigma_H));
    h_limit[3]->Fill(xx,yy,expLimit);
    h_limit[4]->Fill(xx,yy,exp_m1s);
    h_limit[5]->Fill(xx,yy,exp_p1s);
    h_limit[6]->Fill(xx,yy,obsLimit/3.);
    h_limit[7]->Fill(xx,yy,obsLimit*3.);

  }// while
  fin.close();

  // fix pot holes for all TH2D's
  for(int i=0; i<nxs; i++) fillPotHoles(h_xs[i],diagonal);
  for(int i=0; i<nlimit; i++) fillPotHoles(h_limit[i],diagonal);

  // common labels for every plot
  TLatex* lat = new TLatex(0.18,0.92,"#bf{CMS Preliminary}  #sqrt{s} = 8 TeV");
  lat->SetNDC(true);
  lat->SetTextFont(43);
  lat->SetTextSize(21/*25*/);

  TLatex* lat2 = new TLatex(0.49, 0.92, "       L_{int} = 19.7 fb^{-1},   #geq 2 #gamma's,   #geq "+scan);
  lat2->SetNDC(true);
  lat2->SetTextFont(43);
  lat2->SetTextSize(21/*25*/);

  TLatex* lat3 = new TLatex(0.49, 0.92, "       L_{int} = 19.7 fb^{-1},   #geq 2 #gamma's,   #geq "+scan);
  lat3->SetNDC(true);
  lat3->SetTextFont(43);
  lat3->SetTextSize(25);


  // for diagonal == 1
  TGraph* upperDiagonalRegion = new TGraph(3);
  upperDiagonalRegion->SetPoint(0,TMath::Min(xMin,yMin),TMath::Min(xMin,yMin));
  upperDiagonalRegion->SetPoint(1,TMath::Min(xMin,yMin),yBins[nY]);
  upperDiagonalRegion->SetPoint(2,TMath::Max(xBins[nX],yBins[nY]),TMath::Max(xBins[nX],yBins[nY]));
  upperDiagonalRegion->SetFillColor(16);

  TLatex* lat_upperDiagonal = new TLatex(0.6,0.85,"#tilde{t} NLSP");
  lat_upperDiagonal->SetNDC(true);
  lat_upperDiagonal->SetTextSize(0.04);

  // for diagonal == 2
  TGraph* lowerDiagonalRegion = new TGraph(3);
  lowerDiagonalRegion->SetPoint(0,TMath::Min(xMin,yMin),TMath::Min(xMin,yMin));
  lowerDiagonalRegion->SetPoint(1,xBins[nX],TMath::Min(xMin,yMin));
  lowerDiagonalRegion->SetPoint(2,TMath::Max(xBins[nX],yBins[nY]),TMath::Max(xBins[nX],yBins[nY]));
  lowerDiagonalRegion->SetFillColor(16);

  TLatex* lat_lowerDiagonal;
  if(scan.Contains("WB")) lat_lowerDiagonal = new TLatex(0.6,0.25,"M_{bino} > M_{wino}");
  else lat_lowerDiagonal = new TLatex(0.6,0.25,"#tilde{t} NLSP");
  lat_lowerDiagonal->SetNDC(true);
  lat_lowerDiagonal->SetTextSize(0.04);

  TString title;
  double minVal, maxVal;

  TCanvas* can_limit = new TCanvas("can_limit_"+scan,"can_limit_"+scan,900,800);
  getMinMaxValues(h_xs[0],minVal,maxVal,diagonal);
  std::cout << "limit min=" << minVal << ", max=" << maxVal << std::endl;
  h_xs[0]->SetMaximum(1.1*maxVal);
  h_xs[0]->SetMinimum(0.9*minVal);
  title = ";" + xLabel + ";" + yLabel + ";95% CL cross section upper limit [pb]";
  h_xs[0]->SetTitle(title);
  can_limit->SetLogz();
  h_xs[0]->Draw(option2D);
  lat->Draw("same");
  lat2->Draw("same");
//   if(diagonal==1) {
//     upperDiagonalRegion->Draw("same f");
//     //    lat_upperDiagonal->Draw("same");
//   }
//   else if(diagonal==2) {
//     lowerDiagonalRegion->Draw("same f");
//     //    lat_lowerDiagonal->Draw("same");
//   }
  can_limit->RedrawAxis();
  if(print) {
    can_limit->Print("",".gif");
    can_limit->Print("",".pdf");
  }

  TH2D* h_back;
  if(scan.Contains("gB")) h_back = new TH2D("h_back_"+scan,";"+xLabel+";"+yLabel,100,xMin,xMax,100,yMin,yMax+500);
  else h_back = new TH2D("h_back_"+scan,";"+xLabel+";"+yLabel,100,xMin,xMax,100,yMin,yMax);

  std::cout << "Now get contours..." << std::endl;

  TGraph* curv[nlimit];

  //double contours[2]={ 1.0, 1.5};
  double contours[2]={ 0.0, 1.0};
  TCanvas *can_excl01 = new TCanvas("can_contour_"+scan, "can_contour_"+scan,1200,800);
  can_excl01->Divide(nlimit/2,2);
  for(int i=0; i<nlimit; i++) {
    can_excl01->cd(i + 1);
    h_back->Draw();

    h_limit[i]->SetContour(2,contours);
    h_limit[i]->Draw("SAME CONT LIST");
    gPad->Update();
    
    TObjArray *contsM = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
    TList* contLevel = (TList*)contsM->At(0);
    curv[i] = (TGraph*)contLevel->First()->Clone("exclusion_contour_"+limitname[i]);
  }// for i

  std::cout << "Smoothing..." << std::endl;

  TGraphSmooth* gs[nlimit];
  TGraph* curvS[nlimit];
  for(int i=0; i<nlimit; i++) {
    gs[i] = new TGraphSmooth("normal");
    if(useSmooth) curvS[i] = gs[i]->SmoothSuper(curv[i]);
    else curvS[i] = (TGraph*) curv[i]->Clone();
    curvS[i]->SetName("exclusion_smooth_contour_"+limitname[i]);
  }


  // make excluded region
  TGraph* excludedRegion = new TGraph(curvS[0]->GetN()+3);
  int nbins = curvS[0]->GetN();
  for(int i=0; i<nbins; i++){
    double x,y;
    curvS[0]->GetPoint(i,x,y);
    excludedRegion->SetPoint(i,x,y);
  }

//   excludedRegion->SetPoint(nbins,xMax,yMin);
//   excludedRegion->SetPoint(nbins+1,xMin,yMin);
//   excludedRegion->SetPoint(nbins+2,xMin,yMax);

  excludedRegion->SetPoint(nbins,xMin,yMax);
  excludedRegion->SetPoint(nbins+1,xMin,yMin);
  excludedRegion->SetPoint(nbins+2,xMax,yMin);

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

  TGraph* exp1sigma_aroundExp = makeBandGraph(curvS[4],curvS[5]);

  // experimental 1 sigma band around expected limit
  exp1sigma_aroundExp->SetFillColor(kOrange-3);
  exp1sigma_aroundExp->SetFillStyle(1001);


  TCanvas* can_excl02 = new TCanvas("can_exclusion_"+scan, "can_exclusion_"+scan,900,800);
  h_back->Draw();
  can_excl02->SetRightMargin(0.08);

  // ecluded region
  //  excludedRegion->Draw("SAME F");

  exp1sigma_aroundExp->Draw("SAME F");
  curvS[3]->Draw("SAME L");
  curvS[1]->Draw("SAME L");
  curvS[2]->Draw("SAME L");
  curvS[0]->Draw("SAME L");
  lat->Draw("same");
  lat3->Draw("same");

  //  PrintPoints(curvS[0]);

  if(diagonal==1) {
    upperDiagonalRegion->Draw("same f");
    //    lat_upperDiagonal->Draw("same");
  }
  else if(diagonal==2){
    lowerDiagonalRegion->Draw("same f");
    //    lat_lowerDiagonal->Draw("same");
  }
  //      PrintPoints(curv[i][j]);
  //      RemovePoints(curv[i][j]);



  double leg_xmin = 0.58;
  double leg_xmax = 0.9;
  double leg_ymin = 0.64;
  double leg_ymax = 0.87;
  if(diagonal){
    leg_xmin -= 0.35;
    leg_xmax -= 0.35;
  }
  if(scan.Contains("combined")){
    leg_xmin = 0.2;
    leg_xmax = 0.5;
    legendFillColor = kBlue-10;
  }

  TString legendTitle = "pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow t + #tilde{#chi}^{0}_{2}, #tilde{#chi}^{0}_{2} #rightarrow #gamma + #tilde{#chi}^{0}_{1}";

  TLegend* leg = new TLegend(leg_xmin,leg_ymin,leg_xmax,leg_ymax,legendTitle,"brNDC");
  leg->SetFillColor(legendFillColor);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->AddEntry(curvS[0],"Observed","L");
  leg->AddEntry(curvS[1],"Observed #pm1#sigma theory","L");
  TGraph* legGraph = (TGraph*) exp1sigma_aroundExp->Clone();
  legGraph->SetLineColor(kOrange+9);
  legGraph->SetLineStyle(9);
  legGraph->SetLineWidth(3);
  leg->AddEntry(legGraph,"Expected #pm1#sigma exp.","FL");
  leg->Draw("same");


  double xv = 0.25;
  double yv = 0.25;
  if(scan.Contains("gsq_W")){
    yv = 0.18;
  }
  else if(scan.Contains("gB") || scan.Contains("gW") || scan.Contains("WB")) {
    xv = 0.22;
    yv = 0.3;
  }
  TLatex* lat4 = new TLatex(xv,yv,"Excluded");
  lat4->SetNDC(true);
  lat4->SetTextFont(43);
  lat4->SetTextSize(30);
  lat4->Draw("same");

  can_excl02->RedrawAxis();

  if(print) {
    can_excl02->Print("",".gif");
    can_excl02->Print("",".pdf");
  }


  TCanvas* can_exclusionOnLimit;
  if(scan.Contains("sms_gg") || scan.Contains("stop-bino")){
    can_exclusionOnLimit = new TCanvas("can_exclusionOnLimit_"+scan,"can_exclusionOnLimit_"+scan,900,800);
    can_exclusionOnLimit->SetLogz();
    h_xs[0]->SetTitle(";M_{gluino} [GeV];M_{Neutralino} [GeV];95% CL cross section upper limit [pb]");
    h_xs[0]->SetMinimum(0);
    h_xs[0]->GetZaxis()->SetRangeUser(9.e-4, 2.1e-2); //durp
    //h_xs[0]->GetZaxis()->SetNdivisions(210);
    h_xs[0]->Draw(option2D);

    curvS[0]->SetLineColor(kBlack);

    curvS[6]->Draw("SAME L");
    curvS[7]->Draw("SAME L");
    curvS[0]->Draw("SAME L");

    upperDiagonalRegion->SetFillColor(0);
    upperDiagonalRegion->Draw("same f");

    TLegend* leg2 = new TLegend(leg_xmin,leg_ymin,leg_xmax - 0.05,leg_ymax,legendTitle,"brNDC");
    //    leg2->SetFillColor(legendFillColor);
    leg2->SetFillColor(0);
    leg2->SetLineColor(0);
    leg2->SetBorderSize(0);
    leg2->SetTextFont(62);
    leg2->SetTextSize(0.03);
    //    leg2->AddEntry("NULL","NLO+NLL Limits","h");
    leg2->AddEntry("NULL", "m(#tilde{q}) >> m(#tilde{g})", "h");
    leg2->AddEntry(curvS[0],"#sigma^{NLO-QCD}","L");
    leg2->AddEntry(curvS[6],"3#times#sigma^{NLO-QCD}","L");
    leg2->AddEntry(curvS[7],"1/3#times#sigma^{NLO-QCD}","L");
    leg2->Draw("same");

    lat->Draw("same");
    lat2->Draw("same");

    can_exclusionOnLimit->RedrawAxis();

    if(print) {
      can_exclusionOnLimit->Print("",".gif");
      can_exclusionOnLimit->Print("",".pdf");
    }

    TFile* fsms = new TFile(output_dir+ "/hist_sms_gg_1jet_output.root","RECREATE");
    fsms->cd();
    h_xs[1]->Write("h2_95CL_limit");
    exp1sigma_aroundExp->Write("contour_exp_1s_band");
    curvS[3]->Write("contour_exp");
    curvS[4]->Write("contour_exp_1s_down");
    curvS[5]->Write("contour_exp_1s_up");
    curvS[1]->Write("contour_theory_obs_1s_down");
    curvS[2]->Write("contour_theory_obs_1s_up");
    curvS[0]->Write("contour_obs");
    h_xs[0]->Write("observed_limit_in_xsecUnit");
    h_xs[1]->Write("expected_limit_in_xsecUnit");
    h_limit[0]->Write("observed_limit_in_R");
    h_limit[1]->Write("observed_limit_1s_down_in_R");
    h_limit[2]->Write("observed_limit_1s_up_in_R");
    h_limit[3]->Write("expected_limit_in_R");
    h_limit[4]->Write("expected_limit_1s_down_in_R");
    h_limit[5]->Write("expected_limit_1s_up_in_R");
    can_exclusionOnLimit->Write();
    fsms->Write();
    fsms->Close();
  }


  fout->cd();
  fout->Write();

  can_limit->Write();
  can_excl01->Write();
  can_excl02->Write();
  if(scan.Contains("sms_gg")) can_exclusionOnLimit->Write();

  for(int i=0; i<nlimit; i++){
    curv[i]->Write();
    curvS[i]->Write();
  }


}




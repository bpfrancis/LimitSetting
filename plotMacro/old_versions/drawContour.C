
#include "util.C"

void drawContour(TString bino="bino", TString jet="jet", bool print=false) {

  bool useCustomGetContour = false;

  TString data_dir = "table/multiChannel";

  gStyle->SetOptStat(0);

  TString label = bino + "_mN375_met100" + jet;
  if(bino.Contains("mNScan")) label = bino + "_met100" + jet;

  gStyle->SetPalette(1);
  gStyle->SetPadLeftMargin(0.15);

  //  TString option2D = "CONT4 Z";
  TString option2D = "COL Z";

  // for bino & wino
  const int nG = 17;
  const int nS = 17;
  float mGVals[nG+1] = {400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100};
  float mSVals[nS+1] = {420, 520, 620, 720, 820, 920, 1020, 1120, 1220, 1320, 1420, 1520, 1620, 1720, 1820, 1920, 2020, 2120};

  // to make valuse on the center of each bin
  float* mGBins = new float[nG+1];
  float* mSBins = new float[nS+1];
  for(int i=0; i<nG+1; i++) mGBins[i] = mGVals[i]-50;
  for(int i=0; i<nS+1; i++) mSBins[i] = mSVals[i]-50;


  // for mNScan
  const int nX = 10;
  const int nY = 23;
  float xVals[nX+1] = {150, 250, 350, 450, 550, 650, 750, 850, 950, 1050, 1150};
  float yVals[nY+1] = {240, 320, 400, 480, 560, 640, 720, 800, 880, 960, 1040, 1120, 1200, 1280, 1360, 1440, 1520, 1600, 1680, 1760, 1840, 1920, 2000, 2100};

  float* xBins = new float[nX+1];
  float* yBins = new float[nY+1];

  for(int i=0; i<nX+1; i++) xBins[i] = xVals[i]-50;
  for(int i=0; i<nY+1; i++) yBins[i] = yVals[i]-40;


  TFile* fout = new TFile("hist_exclusion_"+label+".root","RECREATE");

  const int nxs = 6;
  TString xsname[nxs] = {"xsec","xsec_1L","xsec_1H","xsec_2L","xsec_2H","acc"};
  TH2D* h_xs[nxs];
  for(int i=0; i<nxs; i++) {
    if(bino.Contains("mNScan")) h_xs[i] = new TH2D(xsname[i],xsname[i],nX,xBins,nY,yBins);
    else h_xs[i] = new TH2D(xsname[i],xsname[i],nS,mSBins,nG,mGBins);
  }

  const int nlimit = 6;
  TString limitname[nlimit] = {"limit", "exp", "exp_1L","exp_1H","exp_2L","exp_2H"};
  TH2D* h_limit[nlimit];
  for(int i=0; i<nlimit; i++) {
    if(bino.Contains("mNScan")) h_limit[i] = new TH2D(limitname[i],limitname[i],nX,xBins,nY,yBins);
    else h_limit[i] = new TH2D(limitname[i],limitname[i],nS,mSBins,nG,mGBins);
  }

  TString datafile = data_dir + "/" + bino + jet + ".table";

  std::cout << std::endl << "Using file: " << data_dir << "/" << datafile << std::endl << std::endl;

  std::ifstream fin;
  fin.open(datafile.Data());
  while(1){
    // #echo "mS mG mN acc xsec xsecPDFError xsecRSErrorNeg xsecRSErrorPos obsLimit expLimit exp_m1s exp_m2s exp_p1s exp_p2s"
    int mS, mG, mN;
    double acc, xsec, xsecPDFError, xsecRSErrorNeg, xsecRSErrorPos, obsLimit, expLimit, exp_m1s, exp_m2s, exp_p1s, exp_p2s;
    fin >> mS >> mG >> mN >> acc >> xsec >> xsecPDFError >> xsecRSErrorNeg >> xsecRSErrorPos >> obsLimit >> expLimit >> exp_m1s >> exp_m2s >> exp_p1s >> exp_p2s;
    if(!fin.good()) break;
/*
     std::cout << mS << ", " << mG << ", " << mN << ", " << acc << ", " << xsec << ", " << xsecPDFError << ", "
 	      << xsecRSErrorNeg << ", " << xsecRSErrorPos << ", " << obsLimit << ", " << expLimit << ", "
 	      << exp_m1s << ", " << exp_m2s << ", " << exp_p1s << ", " << exp_p2s << std::endl;
*/
    double oneSigma_L = std::sqrt(xsecRSErrorNeg * xsecRSErrorNeg + xsecPDFError * xsecPDFError);
    double oneSigma_H = std::sqrt(xsecRSErrorPos * xsecRSErrorPos + xsecPDFError * xsecPDFError);

    if(bino.Contains("mNScan")) {
      if(mS != 2500) continue;
      h_xs[5]->Fill(mN,mG,acc);
      h_xs[0]->Fill(mN,mG,xsec);
      h_xs[1]->Fill(mN,mG,xsec - xsec*oneSigma_L);
      h_xs[2]->Fill(mN,mG,xsec + xsec*oneSigma_H);
      h_xs[3]->Fill(mN,mG,xsec - xsec*2*oneSigma_L);
      h_xs[4]->Fill(mN,mG,xsec + xsec*2*oneSigma_H);

      h_limit[0]->Fill(mN,mG,obsLimit*xsec);
      h_limit[1]->Fill(mN,mG,expLimit*xsec);
      h_limit[2]->Fill(mN,mG,exp_m1s*xsec);
      h_limit[3]->Fill(mN,mG,exp_p1s*xsec);
      h_limit[4]->Fill(mN,mG,exp_m2s*xsec);
      h_limit[5]->Fill(mN,mG,exp_p2s*xsec);
    }
    else {
      if(mN != 375) continue;
      h_xs[5]->Fill(mS,mG,acc);
      h_xs[0]->Fill(mS,mG,xsec);
      h_xs[1]->Fill(mS,mG,xsec - xsec*oneSigma_L);
      h_xs[2]->Fill(mS,mG,xsec + xsec*oneSigma_H);
      h_xs[3]->Fill(mS,mG,xsec - xsec*2*oneSigma_L);
      h_xs[4]->Fill(mS,mG,xsec + xsec*2*oneSigma_H);

      h_limit[0]->Fill(mS,mG,obsLimit*xsec);
      h_limit[1]->Fill(mS,mG,expLimit*xsec);
      h_limit[2]->Fill(mS,mG,exp_m1s*xsec);
      h_limit[3]->Fill(mS,mG,exp_p1s*xsec);
      h_limit[4]->Fill(mS,mG,exp_m2s*xsec);
      h_limit[5]->Fill(mS,mG,exp_p2s*xsec);
    }// if - else

  }// while
  fin.close();

  for(int i=0; i<nxs; i++) fillPotHoles(h_xs[i]);
  for(int i=0; i<nlimit; i++) fillPotHoles(h_limit[i]);


  TGraph* noRegion2 = new TGraph(3);
  noRegion2->SetPoint(0,200,200);
  noRegion2->SetPoint(1,1100,200);
  noRegion2->SetPoint(2,1100,1100);
  noRegion2->SetFillColor(16);

  TLatex* lat44 = new TLatex(0.7,0.25,"#tilde{g} NLSP");
  lat44->SetNDC(true);
  lat44->SetTextSize(0.04);


  TString title;

  TCanvas* can_acc = new TCanvas("can_acc_"+label,"can_acc_"+label,1000,800);
  can_acc->SetRightMargin(0.19);
  h_xs[5]->GetXaxis()->SetNdivisions(505);
  h_xs[5]->GetYaxis()->SetNdivisions(505);
  h_xs[5]->GetYaxis()->SetTitleOffset(1.3);
  h_xs[5]->GetZaxis()->SetTitleOffset(1.3);
  if(bino.Contains("bino")) h_xs[5]->GetZaxis()->SetRangeUser(0.14, 0.275/*0.19, 0.35*/);
  title = ";m_{#tilde{q}} (GeV/c^{2});m_{#tilde{g}} (GeV/c^{2});Acceptance";
  if(bino.Contains("mNScan")) title = ";m_{#chi^{0}} (GeV/c^{2});m_{#tilde{g}} (GeV/c^{2});Acceptance";
  h_xs[5]->SetTitle(title);
  h_xs[5]->Draw(option2D);
  if(bino.Contains("mNScan")){
    noRegion2->Draw("same f");
    lat44->Draw("same");
  }
  if(print) {
    can_acc->Print("",".gif");
    can_acc->Print("",".pdf");
  }



  TCanvas* can_xs = new TCanvas("can_xsec_"+label,"can_xsec_"+label,1000,800);
  can_xs->SetRightMargin(0.17);
  can_xs->SetLogz();
  h_xs[0]->GetXaxis()->SetNdivisions(505);
  h_xs[0]->GetYaxis()->SetNdivisions(505);
  h_xs[0]->GetYaxis()->SetTitleOffset(1.3);
  h_xs[0]->GetZaxis()->SetTitleOffset(1.2);
  //h_xs[0]->GetZaxis()->SetRangeUser(1e-4, 20);
  title = ";m_{#tilde{q}} (GeV/c^{2});m_{#tilde{g}} (GeV/c^{2});Cross Section (pb)";
  if(bino.Contains("mNScan")) title = ";m_{#chi^{0}} (GeV/c^{2});m_{#tilde{g}} (GeV/c^{2});Cross Section (pb)";
  h_xs[0]->SetTitle(title);
  h_xs[0]->Draw(option2D);
  if(bino.Contains("mNScan")){
    noRegion2->Draw("same f");
    lat44->Draw("same");
  }
  if(print) {
    can_xs->Print("",".gif");
    can_xs->Print("",".pdf");
  }

  TCanvas* can_limit = new TCanvas("can_limit_"+label,"can_limit_"+label,1000,800);
  can_limit->SetRightMargin(0.2);
  h_limit[0]->GetXaxis()->SetNdivisions(505);
  h_limit[0]->GetYaxis()->SetNdivisions(505);
  h_limit[0]->GetYaxis()->SetTitleOffset(1.3);
  h_limit[0]->GetZaxis()->SetTitleOffset(1.6);
  if(bino.Contains("bino")) h_limit[0]->GetZaxis()->SetRangeUser(0.005, 0.015);
  if(bino.Contains("wino")) h_limit[0]->GetZaxis()->SetRangeUser(0, 9);
  if(bino.Contains("wino")){
    can_limit->SetRightMargin(0.17);
    h_limit[0]->GetZaxis()->SetTitleOffset(1.0);
  }

  title = ";m_{#tilde{q}} (GeV/c^{2});m_{#tilde{g}} (GeV/c^{2});95% CL Upper Limit (pb)";
  if(bino.Contains("mNScan")) title = ";m_{#chi^{0}} (GeV/c^{2});m_{#tilde{g}} (GeV/c^{2});95% CL Upper Limit (pb)";
  h_limit[0]->SetTitle(title);
  h_limit[0]->Draw(option2D);
  if(bino.Contains("mNScan")){
    noRegion2->Draw("same f");
    lat44->Draw("same");
  }

  TLatex * lat_temp = new TLatex(0.2, 0.92, "CMS Preliminary #int #font[12]{L}dt = 4.04fb^{-1}, #sqrt{s} = 8 TeV");
  lat_temp->SetNDC(true);
  lat_temp->SetTextFont(43);
  lat_temp->SetTextSize(30);
  lat_temp->Draw("same");

  if(print) {
    can_limit->Print("",".gif");
    can_limit->Print("",".pdf");
  }

  // now find exclusion curves
  TCanvas* can_diff = new TCanvas("can_diff_"+label,"can_diff_"+label,1200,800);
  //  can_diff->Divide(nlimit,nxs);
  TH2D* h_excl[nlimit][nxs-1];
  for(int i=0; i<nlimit; i++) {
    for(int j=0; j<nxs-1; j++) {

      h_excl[i][j] = (TH2D*) h_limit[i]->Clone("exclusion_"+limitname[i]+"_"+xsname[j]);
      int nbinsx = h_excl[i][j]->GetXaxis()->GetNbins();
      int nbinsy = h_excl[i][j]->GetYaxis()->GetNbins();

      for( int ibx=1; ibx<=nbinsx; ++ibx){
	for( int iby=1; iby<=nbinsy; ++iby){
	  double x1 = h_limit[i]->GetBinContent(ibx,iby);
	  double x2 = h_xs[j]->GetBinContent(ibx,iby);
	  h_excl[i][j]->SetBinContent(ibx,iby,x2-x1);
	  x1 = h_limit[i]->GetBinError(ibx,iby);
	  x2 = h_xs[j]->GetBinError(ibx,iby);
	  h_excl[i][j]->SetBinError(ibx,iby,std::sqrt(x1*x1+x2*x2));
	}// for iby
      }// for ibx
      fixBadCells(h_excl[i][j]);
      if(i==0 && j==0) h_excl[i][j]->Draw("TEXT");
    }// for j
  }// for i


  float xmin = 400;
  float ymin = 400;
  float xmax = 2000;
  float ymax = 2000;
  if(bino.Contains("mNScan")){
    xmin = 200;
    xmax = 1500;
    ymin = 300;
    ymax = 2000;
  }


  TGraph* curv[nlimit][nxs-1];
  TGraphSmooth* gsmooth = new TGraphSmooth("normal");

  TH2D* h_back;
  if(bino.Contains("mNScan")) h_back = new TH2D("h_back",";m_{#tilde{#chi^{0}}} (GeV/c^{2});m_{#tilde{g}} (GeV/c^{2})",100,xmin,xmax,100,ymin,ymax);
  else h_back = new TH2D("h_back",";m_{#tilde{q}} (GeV/c^{2});m_{#tilde{g}} (GeV/c^{2})",100,xmin,xmax,100,ymin,ymax);

  h_back->GetXaxis()->SetNdivisions(505);
  h_back->GetYaxis()->SetNdivisions(505);
  h_back->GetYaxis()->SetTitleOffset(1.2);


  double contours[2]={ 0.0, 1.0 };
  TCanvas *can_excl01 = new TCanvas("can_excl01_"+label, "can_excl01_"+label,1200,800);
  can_excl01->Divide(nlimit,nxs-1);
  for(int i=0; i<nlimit; i++) {
    for(int j=0; j<nxs-1; j++) {

      can_excl01->cd(j*nlimit + i + 1);
      h_back->Draw();

      if(useCustomGetContour) {
        curv[i][j] = getContour(h_excl[i][j],"excl_curv_"+limitname[i]+"_"+xsname[j]);
        curv[i][j]->Draw("SAME L");
      }
      else {
	h_excl[i][j]->SetContour(2,contours);
	h_excl[i][j]->Draw("SAME CONT LIST");
	gPad->Update();

	TObjArray *contsM = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
	TList* contLevel = (TList*)contsM->At(0);
	curv[i][j] = (TGraph*)contLevel->First()->Clone("excl_curv_"+limitname[i]+"_"+xsname[j]);
      }
      //      PrintPoints(curv[i][j]);
      //      RemovePoints(curv[i][j]);

    }// for j
  }// for i

  if(bino.Contains("mNScan")) {
    for(int i=0; i<nlimit; i++) {
      for(int j=0; j<nxs-1; j++) {
	RemovePoints(curv[i][j]);
      }
    }

    for(int i=0; i<nlimit; i++) {
      for(int j=0; j<nxs-1; j++) {
	double x,y;
	int whichone = curv[i][j]->GetN()-1;
	curv[i][j]->GetPoint(whichone-1,x,y);
	curv[i][j]->SetPoint(whichone,y,y);
      }
    }
  }

  TGraphSmooth* gs[nlimit][nxs-1];
  TGraph* curvS[nlimit][nxs-1];
  for(int i=0; i<nlimit; i++) {
    for(int j=0; j<nxs-1; j++) {
      //      gs[i][j] = new TGraphSmooth("normal");
      //      curvS[i][j] = gs[i][j]->SmoothSuper(curv[i][j]);
      curvS[i][j] = (TGraph*) curv[i][j]->Clone();
      curvS[i][j]->SetName("excl_curvS_"+limitname[i]+"_"+xsname[j]);
    }
  }
/*
  std::cout << "curv[3][0]----------------- N : " << curv[3][0]->GetN() << std::endl;
  PrintPoints(curv[3][0]);
  std::cout << "curvS[3][0]----------------- N : " << curvS[3][0]->GetN() << std::endl;
  PrintPoints(curvS[3][0]);
  std::cout << "---------------------------------" << std::endl;
*/

  // make excluded region
  TGraph* excludedRegion = new TGraph(curvS[0][0]->GetN()+3);
  int nbins = curvS[0][0]->GetN();
  for(int i=0; i<nbins; i++){
    double x,y;
    curvS[0][0]->GetPoint(i,x,y);
    excludedRegion->SetPoint(i,x,y);
  }

  excludedRegion->SetPoint(nbins,xmax,ymin);
  excludedRegion->SetPoint(nbins+1,xmin,ymin);
  excludedRegion->SetPoint(nbins+2,xmin,ymax);

  // make band graph
  TGraph* exp1sigma_aroundExp = makeBandGraph(curvS[2][0],curvS[3][0]);


  TCanvas* can_excl02 = new TCanvas("can_excl02_"+label, "can_excl02_"+label,1000,800);
  h_back->Draw();
  //  can_excl02->SetGrid(1,1);

  // ecluded region
  excludedRegion->SetFillColor(kBlue-10);
  excludedRegion->SetFillStyle(1001);
  excludedRegion->Draw("SAME F");

  // experimental 1 sigma band around expected limit
  exp1sigma_aroundExp->SetFillColor(kOrange-3);
  exp1sigma_aroundExp->SetFillStyle(1001);
  exp1sigma_aroundExp->Draw("SAME F");

  // expected limit
  curvS[1][0]->SetLineStyle(9);
  curvS[1][0]->SetLineWidth(2);
  curvS[1][0]->SetLineColor(kOrange+9);
  curvS[1][0]->Draw("SAME L");

  // theory 1 sigma around expected limit
  curvS[1][1]->SetLineStyle(3);
  curvS[1][1]->SetLineWidth(1);
  curvS[1][1]->SetLineColor(kOrange+9);
  curvS[1][1]->Draw("SAME L");

  curvS[1][2]->SetLineStyle(3);
  curvS[1][2]->SetLineWidth(1);
  curvS[1][2]->SetLineColor(kOrange+9);
  curvS[1][2]->Draw("SAME L");

  // observed limit
  curvS[0][0]->SetLineWidth(3);
  curvS[0][0]->SetLineColor(4);
  curvS[0][0]->Draw("SAME L");

  // theory 1 sigma around observed limit
  curvS[0][1]->SetLineStyle(3);
  curvS[0][1]->SetLineWidth(2);
  curvS[0][1]->SetLineColor(kBlue);
  curvS[0][1]->Draw("SAME L");

  curvS[0][2]->SetLineStyle(3);
  curvS[0][2]->SetLineWidth(2);
  curvS[0][2]->SetLineColor(kBlue);
  curvS[0][2]->Draw("SAME L");


  PrintPoints(curvS[0][0]);


  float leg_xmin = 0.2;
  float leg_xmax = 0.45;
  float leg_ymin = 0.4;
  float leg_ymax = 0.7;
  if(bino.Contains("mNScan")){
    leg_xmin -= 0.45;
    leg_xmax -= 0.45;
  }
  if(bino.Contains("wino")){
    leg_xmin = 0.55;
    leg_xmax = 0.85;
    leg_ymin = 0.45;
    leg_ymax = 0.75;
  }

  gPad->SetTicks();

  TLegend* leg;
  if(bino.Contains("mNScan")) leg = new TLegend(leg_xmin,leg_ymin,leg_xmax,leg_ymax,"","brNDC");
  else leg = new TLegend(leg_xmin,leg_ymin,leg_xmax,leg_ymax,"GGM "+bino+"-like #tilde{#chi}^{0}","brNDC");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  if(bino.Contains("bino")) leg->SetFillColor(kBlue-10);
  if(bino.Contains("bino")) leg->SetLineColor(kBlue-10);
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.03);
  if(bino.Contains("mNScan")) leg->AddEntry("NULL","m_{#tilde{q}} = 2500 (GeV/c^{2})","h");
  else leg->AddEntry("NULL","m_{#tilde{#chi}^{0}} = 375 (GeV/c^{2})","h");
  TString jetRequirement = "No Jet Requirement";
  if(jet.Contains("_jet")) jetRequirement = "#geq1 Jet Requirement";
  if(jet.Contains("_btag")) jetRequirement = "At least 1 btag requirement";
  leg->AddEntry("NULL",jetRequirement,"h");
  leg->AddEntry("NULL","NLO Limits","h");
  leg->AddEntry(curvS[0][0],"Observed","L");
  leg->AddEntry(curvS[0][1],"#pm1#sigma (theory)","L");
  leg->AddEntry(curvS[1][0],"Expected","L");
  leg->AddEntry(curvS[1][1],"#pm1#sigma (theory)","L");
  leg->AddEntry(exp1sigma_aroundExp,"#pm1#sigma (experimental)","F");
  leg->Draw("same");

  Double_t lat_x = leg_xmin+0.02;
  Double_t lat_y = leg_ymax+0.07; if(bino.Contains("wino")) lat_y = 0.85;

  //TLatex* lat = new TLatex(lat_x,lat_y,"CMS Preliminary 2012");
  TLatex* lat = new TLatex(lat_x,lat_y,"CMS Preliminary");
  lat->SetNDC(true);
  lat->SetTextFont(43);
  lat->SetTextSize(30);
  lat->Draw("same");

  Double_t lat_x2 = leg_xmin;
  Double_t lat_y2 = leg_ymax+0.03; if(bino.Contains("wino")) lat_y2 = 0.81;

  TLatex* lat2 = new TLatex(lat_x2,lat_y2,"#int #font[12]{L}dt = 4.04fb^{-1}, #sqrt{s} = 8 TeV");
  lat2->SetNDC(true);
  lat2->SetTextFont(43);
  lat2->SetTextSize(24);
  lat2->Draw("same");

  float xv = 0.25;
  float yv = 0.25;
  if(bino.Contains("wino")){
    xv = 0.19;
    yv = 0.19;
  }
  if(bino.Contains("mNScan")){
    xv = 0.2;
    yv = 0.3;
  }
  TLatex* lat3 = new TLatex(xv,yv,"Excluded");
  lat3->SetNDC(true);
  lat3->SetTextFont(52);
  lat3->SetTextSize(0.06);
  lat3->Draw("same");

  TGraph* noRegion = new TGraph(3);
  noRegion->SetPoint(0,TMath::Min(xmin,ymin),TMath::Min(xmin,ymin));
  noRegion->SetPoint(1,xmax,ymin);
  noRegion->SetPoint(2,TMath::Min(xmax,ymax),TMath::Min(xmax,ymax));
  noRegion->SetFillColor(16);

  TLatex* lat4 = new TLatex(0.7,0.25,"#tilde{g} NLSP");
  lat4->SetNDC(true);
  lat4->SetTextSize(0.04);

  if(bino.Contains("mNScan")){
    noRegion->Draw("same f");
    lat4->Draw("same");
  }

  can_excl02->RedrawAxis();


  if(print) {
    can_excl02->Print("",".gif");
    can_excl02->Print("",".pdf");
  }

  fout->cd();
  fout->Write();

  can_acc->Write();
  can_xs->Write();
  can_limit->Write();
  can_excl01->Write();
  can_excl02->Write();

  for(int i=0; i<nlimit; i++){
    for(int j=0; j<nxs-1; j++){
      curv[i][j]->Write();
      curvS[i][j]->Write();
    }
  }


}




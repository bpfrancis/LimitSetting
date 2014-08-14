#include "makeTemplate_new.h"

void makeTemplate_new() {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  TString hist_dir = "inputHists";

  TFile * f_xsec = new TFile("xsecdat/stop-bino_xsecs.root", "READ");
  TH2D * h_xsec = (TH2D*)f_xsec->Get("real_xsec");
  TH2D * h_xsec_errors = (TH2D*)f_xsec->Get("real_errors");

  TFile * f_ele_signal = new TFile(hist_dir + "/signalLimits_ele_bjj.root", "READ");
  TFile * f_muon_signal = new TFile(hist_dir + "/signalLimits_muon_bjj.root", "READ");

  TFile * fElectron = new TFile(hist_dir + "/limitInputs_ele_bjj.root", "READ");
  TFile * fMuon = new TFile(hist_dir + "/limitInputs_muon_bjj.root", "READ");

  TH1D * ele_data = (TH1D*)fElectron->Get("pfMET_gg_ele_bjj");
  TH1D * muon_data = (TH1D*)fMuon->Get("pfMET_gg_muon_bjj");

  vector<ObservedData> data;
  data.push_back(ObservedData("ele", ele_data));
  data.push_back(ObservedData("muon", muon_data));

  vector<TH1D*> central_ele,
    btagWeightUp_ele, btagWeightDown_ele,
    puWeightUp_ele, puWeightDown_ele,
    topPtUp_ele, topPtDown_ele,
    JECup_ele, JECdown_ele,
    leptonSFup_ele, leptonSFdown_ele,
    photonSFup_ele, photonSFdown_ele,
    scaleUp_ele, scaleDown_ele,
    pdfUp_ele, pdfDown_ele;

  GetBackgroundHistograms(fElectron, "pfMET", "ele_bjj",
			  central_ele,
			  btagWeightUp_ele, btagWeightDown_ele,
			  puWeightUp_ele, puWeightDown_ele,
			  topPtUp_ele, topPtDown_ele,
			  JECup_ele, JECdown_ele,
			  leptonSFup_ele, leptonSFdown_ele,
			  photonSFup_ele, photonSFdown_ele,
			  scaleUp_ele, scaleDown_ele,
			  pdfUp_ele, pdfDown_ele);

  vector<TH1D*> central_muon,
    btagWeightUp_muon, btagWeightDown_muon,
    puWeightUp_muon, puWeightDown_muon,
    topPtUp_muon, topPtDown_muon,
    JECup_muon, JECdown_muon,
    leptonSFup_muon, leptonSFdown_muon,
    photonSFup_muon, photonSFdown_muon,
    scaleUp_muon, scaleDown_muon,
    pdfUp_muon, pdfDown_muon;

  GetBackgroundHistograms(fMuon, "pfMET", "muon_bjj",
			  central_muon,
			  btagWeightUp_muon, btagWeightDown_muon,
			  puWeightUp_muon, puWeightDown_muon,
			  topPtUp_muon, topPtDown_muon,
			  JECup_muon, JECdown_muon,
			  leptonSFup_muon, leptonSFdown_muon,
			  photonSFup_muon, photonSFdown_muon,
			  scaleUp_muon, scaleDown_muon,
			  pdfUp_muon, pdfDown_muon);
  
  vector<BackgroundProcess> ttbar;
  ttbar.push_back(BackgroundProcess(0, "ttbar", "ele",
				    central_ele,
				    btagWeightUp_ele, btagWeightDown_ele,
				    puWeightUp_ele, puWeightDown_ele,
				    topPtUp_ele, topPtDown_ele,
				    JECup_ele, JECdown_ele,
				    leptonSFup_ele, leptonSFdown_ele,
				    photonSFup_ele, photonSFdown_ele,
				    scaleUp_ele, scaleDown_ele,
				    pdfUp_ele, pdfDown_ele));

  ttbar.push_back(BackgroundProcess(0, "ttbar", "muon",
				    central_muon,
				    btagWeightUp_muon, btagWeightDown_muon,
				    puWeightUp_muon, puWeightDown_muon,
				    topPtUp_muon, topPtDown_muon,
				    JECup_muon, JECdown_muon,
				    leptonSFup_muon, leptonSFdown_muon,
				    photonSFup_muon, photonSFdown_muon,
				    scaleUp_muon, scaleDown_muon,
				    pdfUp_muon, pdfDown_muon));
				  
  vector<BackgroundProcess> ttgamma;
  ttgamma.push_back(BackgroundProcess(6, "ttgamma", "ele",
				      central_ele,
				      btagWeightUp_ele, btagWeightDown_ele,
				      puWeightUp_ele, puWeightDown_ele,
				      topPtUp_ele, topPtDown_ele,
				      JECup_ele, JECdown_ele,
				      leptonSFup_ele, leptonSFdown_ele,
				      photonSFup_ele, photonSFdown_ele,
				      scaleUp_ele, scaleDown_ele,
				      pdfUp_ele, pdfDown_ele));
  
  ttgamma.push_back(BackgroundProcess(6, "ttgamma", "muon",
				      central_muon,
				      btagWeightUp_muon, btagWeightDown_muon,
				      puWeightUp_muon, puWeightDown_muon,
				      topPtUp_muon, topPtDown_muon,
				      JECup_muon, JECdown_muon,
				      leptonSFup_muon, leptonSFdown_muon,
				      photonSFup_muon, photonSFdown_muon,
				      scaleUp_muon, scaleDown_muon,
				      pdfUp_muon, pdfDown_muon));

  GridPoint grid;
  grid.data = data;
  grid.ttbar = ttbar;
  grid.ttgamma = ttgamma;

  int mst[29] = {110, 160, 185, 210, 235, 260, 285, 310, 335, 360, 385, 410, 460, 510, 560, 610, 660, 710, 810, 910, 1010, 1110, 1210, 1310, 1410, 1510, 1710, 2010, 5010};
  int mBino[31] = {25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 375, 425, 475, 525, 575, 625, 675, 725, 825, 925, 1025, 1125, 1225, 1325, 1425, 1525, 1725, 2025};
  
  Double_t xbins[31];
  xbins[0] = 0;
  xbins[1] = 55;
  for(int i = 1; i < 29; i++) xbins[i+1] = (mst[i] + mst[i-1])/2.;
  xbins[30] = 6510;

  Double_t ybins[33];
  ybins[0] = 0;
  ybins[1] = 12.5;
  for(int i = 1; i < 31; i++) ybins[i+1] = (mBino[i] + mBino[i-1])/2.;
  ybins[32] = 2175;

  TH2D * h_yield_ele = new TH2D("yield_ele", "yield_ele", 30, xbins, 32, ybins);
  h_yield_ele->GetXaxis()->SetTitle("Stop mass (GeV/c^{2})");
  h_yield_ele->GetYaxis()->SetTitle("Bino mass (GeV/c^{2})");
  h_yield_ele->GetXaxis()->SetRangeUser(222.5, 960);
  h_yield_ele->GetYaxis()->SetRangeUser(137.5, 775);
  h_yield_ele->GetZaxis()->SetLabelSize(0.02);
  
  TH2D * h_yield_muon = (TH2D*)h_yield_ele->Clone("yield_muon");

  TH2D * h_stat_ele = (TH2D*)h_yield_ele->Clone("stat_ele");
  TH2D * h_stat_muon = (TH2D*)h_yield_ele->Clone("stat_muon");

  TH2D * h_btag_ele = (TH2D*)h_yield_ele->Clone("btag_ele");
  TH2D * h_btag_muon = (TH2D*)h_yield_ele->Clone("btag_muon");

  TH2D * h_pileup_ele = (TH2D*)h_yield_ele->Clone("pileup_ele");
  TH2D * h_pileup_muon = (TH2D*)h_yield_ele->Clone("pileup_muon");

  TH2D * h_jec_ele = (TH2D*)h_yield_ele->Clone("jec_ele");
  TH2D * h_jec_muon = (TH2D*)h_yield_ele->Clone("jec_muon");

  TH2D * h_leptonSF_ele = (TH2D*)h_yield_ele->Clone("leptonSF_ele");
  TH2D * h_leptonSF_muon = (TH2D*)h_yield_ele->Clone("leptonSF_muon");

  TH2D * h_photonSF_ele = (TH2D*)h_yield_ele->Clone("photonSF_ele");
  TH2D * h_photonSF_muon = (TH2D*)h_yield_ele->Clone("photonSF_muon");

  for(int i = 0; i < 899; i++) {
    
    TH1D *hsig_ele,
      *hsig_ele_btagWeightUp, *hsig_ele_btagWeightDown,
      *hsig_ele_puWeightUp, *hsig_ele_puWeightDown,
      *hsig_ele_topPtUp, *hsig_ele_topPtDown,
      *hsig_ele_JECup, *hsig_ele_JECdown,
      *hsig_ele_leptonSFup, *hsig_ele_leptonSFdown,
      *hsig_ele_photonSFup, *hsig_ele_photonSFdown;

    bool foundPoint = GetSignalHistograms(f_ele_signal, "pfMET", "ele_bjj",
					  mst[int(i)/31], mBino[int(i)%31],
					  hsig_ele,
					  hsig_ele_btagWeightUp, hsig_ele_btagWeightDown,
					  hsig_ele_puWeightUp, hsig_ele_puWeightDown,
					  hsig_ele_topPtUp, hsig_ele_topPtDown,
					  hsig_ele_JECup, hsig_ele_JECdown,
					  hsig_ele_leptonSFup, hsig_ele_leptonSFdown,
					  hsig_ele_photonSFup, hsig_ele_photonSFdown);

    TH1D *hsig_muon,
      *hsig_muon_btagWeightUp, *hsig_muon_btagWeightDown,
      *hsig_muon_puWeightUp, *hsig_muon_puWeightDown,
      *hsig_muon_topPtUp, *hsig_muon_topPtDown,
      *hsig_muon_JECup, *hsig_muon_JECdown,
      *hsig_muon_leptonSFup, *hsig_muon_leptonSFdown,
      *hsig_muon_photonSFup, *hsig_muon_photonSFdown;

    foundPoint &= GetSignalHistograms(f_muon_signal, "pfMET", "muon_bjj",
				      mst[int(i)/31], mBino[int(i)%31],
				      hsig_muon,
				      hsig_muon_btagWeightUp, hsig_muon_btagWeightDown,
				      hsig_muon_puWeightUp, hsig_muon_puWeightDown,
				      hsig_muon_topPtUp, hsig_muon_topPtDown,
				      hsig_muon_JECup, hsig_muon_JECdown,
				      hsig_muon_leptonSFup, hsig_muon_leptonSFdown,
				      hsig_muon_photonSFup, hsig_muon_photonSFdown);
    
    if(!foundPoint) continue;

    vector<SignalYield> signal;

    double xsec = h_xsec->GetBinContent(h_xsec->FindBin(mst[int(i)/31], mBino[int(i)%31]));
    double xsecError = h_xsec_errors->GetBinContent(h_xsec_errors->FindBin(mst[int(i)/31], mBino[int(i)%31]));

    signal.push_back(SignalYield("ele",
				 hsig_ele,
				 hsig_ele_btagWeightUp, hsig_ele_btagWeightDown,
				 hsig_ele_puWeightUp, hsig_ele_puWeightDown,
				 hsig_ele_topPtUp, hsig_ele_topPtDown,
				 hsig_ele_JECup, hsig_ele_JECdown,
				 hsig_ele_leptonSFup, hsig_ele_leptonSFdown,
				 hsig_ele_photonSFup, hsig_ele_photonSFdown,
				 xsec, xsecError));

    signal[0].FillHistograms(mst[int(i)/31], mBino[int(i)%31],
			     h_yield_ele,
			     h_stat_ele,
			     h_btag_ele,
			     h_pileup_ele,
			     h_jec_ele,
			     h_leptonSF_ele,
			     h_photonSF_ele);

    signal.push_back(SignalYield("muon",
				 hsig_muon,
				 hsig_muon_btagWeightUp, hsig_muon_btagWeightDown,
				 hsig_muon_puWeightUp, hsig_muon_puWeightDown,
				 hsig_muon_topPtUp, hsig_muon_topPtDown,
				 hsig_muon_JECup, hsig_muon_JECdown,
				 hsig_muon_leptonSFup, hsig_muon_leptonSFdown,
				 hsig_muon_photonSFup, hsig_muon_photonSFdown,
				 xsec, xsecError));

    signal[1].FillHistograms(mst[int(i)/31], mBino[int(i)%31],
			     h_yield_muon,
			     h_stat_muon,
			     h_btag_muon,
			     h_pileup_muon,
			     h_jec_muon,
			     h_leptonSF_muon,
			     h_photonSF_muon);

    grid.signal = signal;
    grid.mStop = mst[int(i)/31];
    grid.mBino = mBino[int(i)%31];
    grid.xsec = xsec;
    grid.xsecError = xsecError;

    grid.Print();

  }

  TCanvas * can = new TCanvas("can", "can", 10, 10, 2000, 2000);

  h_yield_ele->Draw("colz");
  can->SaveAs("yield_ele.png");

  h_yield_muon->Draw("colz");
  can->SaveAs("yield_muon.png");

  h_yield_ele->Divide(h_yield_muon);
  h_yield_ele->Draw("colz");
  can->SaveAs("yield_ratio.png");

  h_stat_ele->Draw("colz");
  can->SaveAs("stat_ele.png");

  h_stat_muon->Draw("colz");
  can->SaveAs("stat_muon.png");

  h_stat_ele->Divide(h_stat_muon);
  h_stat_ele->Draw("colz");
  can->SaveAs("stat_ratio.png");

  h_btag_ele->Draw("colz");
  can->SaveAs("btag_ele.png");

  h_btag_muon->Draw("colz");
  can->SaveAs("btag_muon.png");

  h_btag_ele->Divide(h_btag_muon);
  h_btag_ele->Draw("colz");
  can->SaveAs("btag_ratio.png");

  h_pileup_ele->Draw("colz");
  can->SaveAs("pileup_ele.png");

  h_pileup_muon->Draw("colz");
  can->SaveAs("pileup_muon.png");

  h_pileup_ele->Divide(h_pileup_muon);
  h_pileup_ele->Draw("colz");
  can->SaveAs("pileup_ratio.png");

  h_jec_ele->Draw("colz");
  can->SaveAs("jec_ele.png");

  h_jec_muon->Draw("colz");
  can->SaveAs("jec_muon.png");

  h_jec_ele->Divide(h_jec_muon);
  h_jec_ele->Draw("colz");
  can->SaveAs("jec_ratio.png");

  h_leptonSF_ele->Draw("colz");
  can->SaveAs("leptonSF_ele.png");

  h_leptonSF_muon->Draw("colz");
  can->SaveAs("leptonSF_muon.png");

  h_leptonSF_ele->Divide(h_leptonSF_muon);
  h_leptonSF_ele->Draw("colz");
  can->SaveAs("leptonSF_ratio.png");

  h_photonSF_ele->Draw("colz");
  can->SaveAs("photonSF_ele.png");

  h_photonSF_muon->Draw("colz");
  can->SaveAs("photonSF_muon.png");

  h_photonSF_ele->Divide(h_photonSF_muon);
  h_photonSF_ele->Draw("colz");
  can->SaveAs("photonSF_ratio.png");

  delete can;

  f_xsec->Close();
  f_ele_signal->Close();
  f_muon_signal->Close();
  fElectron->Close();
  fMuon->Close();

}


#include "makeTemplate_new.h"

void makeTemplate_new() {

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
  ttgamma.push_back(BackgroundProcess(6, "ttgamma",
				      central_ele,
				      btagWeightUp_ele, btagWeightDown_ele,
				      puWeightUp_ele, puWeightDown_ele,
				      topPtUp_ele, topPtDown_ele,
				      JECup_ele, JECdown_ele,
				      leptonSFup_ele, leptonSFdown_ele,
				      photonSFup_ele, photonSFdown_ele,
				      scaleUp_ele, scaleDown_ele,
				      pdfUp_ele, pdfDown_ele));
  
  ttgamma.push_back(BackgroundProcess(6, "ttgamma",
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
  
  for(int i = 0; i < 899; i++) {
    
    TH1D *hsig_ele,
      *hsig_ele_btagWeightUp, *hsig_ele_btagWeightDown,
      *hsig_ele_puWeightUp, *hsig_ele_puWeightDown,
      *hsig_ele_topPtUp, *hsig_ele_topPtDown,
      *hsig_ele_JECup, *hsig_ele_JECdown,
      *hsig_ele_leptonSFup, *hsig_ele_leptonSFdown,
      *hsig_ele_photonSFup, *hsig_ele_photonSFdown;

    bool foundPoint = GetSignalHistograms(f_ele_signal, "pfMET", "ele_bjj",
					  mst[i], mBino[i],
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
				      mst[i], mBino[i],
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

    signal.push_back(SignalYield("muon",
				 hsig_muon,
				 hsig_muon_btagWeightUp, hsig_muon_btagWeightDown,
				 hsig_muon_puWeightUp, hsig_muon_puWeightDown,
				 hsig_muon_topPtUp, hsig_muon_topPtDown,
				 hsig_muon_JECup, hsig_muon_JECdown,
				 hsig_muon_leptonSFup, hsig_muon_leptonSFdown,
				 hsig_muon_photonSFup, hsig_muon_photonSFdown,
				 xsec, xsecError));

    grid.signal = signal;
    grid.mStop = mst[int(i)/31];
    grid.mBino[int(i)%31];
    grid.xsec = xsec;
    grid.xsecError = xsecError;

    grid.Print();

  }

  f_xsec->Close();
  f_ele_signal->Close();
  f_muon_signal->Close();
  fElectron->Close();
  fMuon->Close();

}


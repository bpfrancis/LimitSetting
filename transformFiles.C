#include <vector>

using namespace std;

TH1D * get_ttjets(TFile * file, TString channel, TString systematic) {

  TH1D * ttjets = (TH1D*)file->Get("pfMET_ttJetsHadronic_"+channel+systematic);
  ttjets->Add((TH1D*)file->Get("pfMET_ttJetsSemiLep_"+channel+systematic));
  ttjets->Add((TH1D*)file->Get("pfMET_ttJetsHadronic_"+channel+systematic));
  
  return ttjets;
}

TH1D * get_ttgamma(TFile * file, TString channel, TString systematic) {

  TH1D * ttgamma = (TH1D*)file->Get("pfMET_ttA_2to5_"+channel+systematic);
  ttgamma->Add((TH1D*)file->Get("pfMET_ttA_2to5_"+channel+systematic));
  ttgamma->Add((TH1D*)file->Get("pfMET_ttA_2to5_"+channel+systematic));
  
  return ttgamma;
}

void transformFiles() {

  const int nSystematics = 17;

  TString systematics[nSystematics] = {"",
				       "_btagWeightUp", "_btagWeightDown",
				       "_puWeightUp", "_puWeightDown",
				       "_topPtUp", "_topPtDown",
				       "_JECup", "_JECdown",
				       "_leptonSFup", "_leptonSFdown",
				       "_photonSFup", "_photonSFdown",
				       "_scaleUp", "_scaleDown",
				       "_pdfUp", "_pdfDown"};

  TString fixedNames[nSystematics] = {"",
				      "_btagWeightUp", "_btagWeightDown",
				      "_puWeightUp", "_puWeightDown",
				      "_topPtUp", "_topPtDown",
				      "_JECUp", "_JECDown",
				      "_leptonSFUp", "_leptonSFDown",
				      "_photonSFUp", "_photonSFDown",
				      "_scaleUp", "_scaleDown",
				      "_pdfUp", "_pdfDown"};

  TFile * fDataEle = new TFile("inputHists/limitInputs_ele_bjj.root", "READ");
 
  vector<TH1D*> h_ele_ttjets, h_ele_ttgamma;
  for(int i = 0; i < nSystematics; i++) {
    h_ele_ttjets.push_back(get_ttjets(fDataEle, "ele_bjj", systematics[i]));
    h_ele_ttgamma.push_back(get_ttgamma(fDataEle, "ele_bjj", systematics[i]));
  }
			   
  TFile * fDataMuon = new TFile("inputHists/limitInputs_muon_bjj.root", "READ");

  vector<TH1D*> h_muon_ttjets, h_muon_ttgamma;
  for(int i = 0; i < nSystematics; i++) {
    h_muon_ttjets.push_back(get_ttjets(fDataMuon, "muon_bjj", systematics[i]));
    h_muon_ttgamma.push_back(get_ttgamma(fDataMuon, "muon_bjj", systematics[i]));
  }

  TFile * fSignalEle = new TFile("inputHists/signalLimits_ele_bjj.root", "READ");
  TFile * fSignalMuon = new TFile("inputHists/signalLimits_muon_bjj.root", "READ");

  const int nSystematics_signal = 13;

  TString systematics_signal[nSystematics_signal] = {"",
				       "_btagWeightUp", "_btagWeightDown",
				       "_puWeightUp", "_puWeightDown",
				       "_topPtUp", "_topPtDown",
				       "_JECup", "_JECdown",
				       "_leptonSFup", "_leptonSFdown",
				       "_photonSFup", "_photonSFdown"};

  TString fixedNames_signal[nSystematics_signal] = {"",
						    "_btagWeightUp", "_btagWeightDown",
						    "_puWeightUp", "_puWeightDown",
						    "_topPtUp", "_topPtDown",
						    "_JECUp", "_JECDown",
						    "_leptonSFUp", "_leptonSFDown",
						    "_photonSFUp", "_photonSFDown"};

  Double_t mst[29] = {110, 160, 185, 210, 235, 260, 285, 310, 335, 360, 385, 410, 460, 510, 560, 610, 660, 710, 810, 910, 1010, 1110, 1210, 1310, 1410, 1510, 1710, 2010, 5010};
  Double_t mBino[31] = {25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 375, 425, 475, 525, 575, 625, 675, 725, 825, 925, 1025, 1125, 1225, 1325, 1425, 1525, 1725, 2025};

  vector<vector<TH1D*> > h_ele_signal, h_muon_signal;

   for(int i = 0; i < 899; i++) {

    index1 = mst[int(i)/31];
    index2 = mBino[int(i)%31];

    if(index1 < index2) continue;

    sprintf(code, "_mst_%d_m1_%d", index1, index2);
    TString code_t = code;

    TH1D * h_ele = (TH1D*)fSignalEle->Get("pfMET_gg_ele_bjj"+code);
    TH1D * h_muon = (TH1D*)fSignalMuon->Get("pfMET_gg_muon_bjj"+code);

    if(!h_ele || !h_muon) continue;

    vector<TH1D*> v_ele, v_muon;
    
    for(int i = 0; i < nSystematics_signal; i++) {
      v_ele.push_back((TH1D*)fSignalEle->Get("pfMET_gg_ele_bjj"+code+systematics_signal[i]));
      v_muon.push_back((TH1D*)fSignalMuon->Get("pfMET_gg_muon_bjj"+code+systematics_signal[i]));
    }

    h_ele_signal.push_back(v_ele);
    h_muon_signal.push_back(v_muon);

   }

  

  TFile	* fOut = new TFile("inputHists/stop-bino_shapes.root", "RECREATE");
  fOut->cd();

  for(unsigned int j = 0; j < h_ele_ttjets.size(); j++) {
    TString newName = "ttjets_ele"+fixedNames[j];
    h_ele_ttjets[j]->Write(newName.Data());
  }

  for(unsigned int j = 0; j < h_ele_ttgamma.size(); j++) {
    TString newName = "ttgamma_ele"+fixedNames[j];
    h_ele_ttgamma[j]->Write(newName.Data());
  }

  for(unsigned int j = 0; j < h_muon_ttjets.size(); j++) {
    TString newName = "ttjets_muon"+fixedNames[j];
    h_muon_ttjets[j]->Write(newName.Data());
  }

  for(unsigned int j = 0; j < h_muon_ttgamma.size(); j++) {
    TString newName = "ttgamma_muon"+fixedNames[j];
    h_muon_ttgamma[j]->Write(newName.Data());
  }

  fOut->Close();

  fDataEle->Close();
  fDataMuon->Close();
}

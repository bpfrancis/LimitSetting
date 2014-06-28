#include "durp.h"

void durp(TString scan="stop-bino") {

  TString output_dir = "hist";
  gSystem->mkdir(output_dir,true);

  int legendFillColor = 16;
  TString option2D = "COL Z";

  PlotMaker * pMaker = new PlotMaker(scan, legendFillColor, option2D);
  pMaker->InitStyle();
  pMaker->SetRange(222.5, 960, 137.5, 775);
  pMaker->SetAxisTitles("m_{Stop} [GeV]", "m_{Bino} [GeV]");

  TString datafile = "table/" + scan + ".table";

  ifstream fin;
  fin.open(datafile.Data());

  while(1){
    //mStop mBino xsec xsecError obsLimit expLimit exp_m1s exp_m2s exp_p1s exp_p2s

    int ms, mb;
    double xsec, xsecError, obsLimit, expLimit, exp_m1s, exp_m2s, exp_p1s, exp_p2s;
    fin >> ms >> mb >> xsec >> xsecError >> obsLimit >> expLimit >> exp_m1s >> exp_m2s >> exp_p1s >> exp_p2s;
    if(!fin.good()) break;

    pMaker->FillBin(ms, mb,
		    xsec, xsecError,
		    obsLimit,
		    expLimit, exp_m1s, exp_m2s, exp_p1s, exp_p2s);
		 
  }// while
  fin.close();

  pMaker->FillPotHoles();

  pMaker->DrawCrossSection();
  pMaker->DrawUpperLimit();

  pMaker->GetContours();
  pMaker->SmoothContours();

  pMaker->SetExclusion();

  pMaker->DrawExclusion();
  pMaker->DrawExclusionOnLimit();

  pMaker->Save(output_dir);
  
  delete pMaker;

}




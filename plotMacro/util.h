#include <TH2D.h>
#include <TGraph.h>
#include <TString.h>
#include <TMath.h>

const double epsilon = 1e-20;

// diagonal 0 : whole plane
// diagonal 1 : upper diagonal part of plane will be excluded 
// diagonal 2 : lower diagonal part of plane will be excluded 

void getMinMaxValues(TH2D *h, double& minVal, double& maxVal, int diagonal=0) {

  maxVal = 0;
  minVal = 9999999;

  int nbinsX = h->GetXaxis()->GetNbins();
  int nbinsY = h->GetYaxis()->GetNbins();

  for(int ix=1; ix <= nbinsX; ix++) {
    for(int iy=1; iy <= nbinsY - 2; iy++) {
      if(diagonal > 0) {
	double xx = h->GetXaxis()->GetBinCenter(ix);
	double yy = h->GetYaxis()->GetBinCenter(iy + 2);
	if(diagonal == 1 && xx < yy) continue;
	if(diagonal == 2 && xx > yy) continue;
      }
      double val = h->GetBinContent(ix,iy);
      if(TMath::IsNaN(val)) continue;
      if(val > maxVal) maxVal = val;
      if(val < minVal) minVal = val;
    }
  }
}


void fillPotHoles(TH2D *h, int diagonal=0) {

  // fill cells which have empty value

  int nbinsX = h->GetXaxis()->GetNbins();
  int nbinsY = h->GetYaxis()->GetNbins();

  for(int ix=1; ix <= nbinsX; ix++) {
    for(int iy=1; iy <= nbinsY - 2; iy++) {
      if(diagonal > 0) {
	double xx = h->GetXaxis()->GetBinCenter(ix);
	double yy = h->GetYaxis()->GetBinCenter(iy + 2);
	if(diagonal == 1 && xx < yy) continue;
	if(diagonal == 2 && xx > yy) continue;
      }
      double val = h->GetBinContent(ix,iy);
      if(TMath::IsNaN(val)) h->SetBinContent(ix,iy,0); // checking for NAN
    }
  }

  for(int ix=1; ix <= nbinsX; ix++) {
    for(int iy=1; iy <= nbinsY - 2; iy++) {
      if(diagonal > 0) {
	double xx = h->GetXaxis()->GetBinCenter(ix);
	double yy = h->GetYaxis()->GetBinCenter(iy + 2);
	if(diagonal == 1 && xx < yy) continue;
	if(diagonal == 2 && xx > yy) continue;
      }
      double val = h->GetBinContent(ix,iy);
      if(fabs(val) > epsilon) continue;
      int ncnt = 0;
      double sum = 0;
      double sumErr = 0;
      double up    = h->GetBinContent(ix,iy+1);
      if(fabs(up) > epsilon && iy < nbinsY){
	sum += up;
	sumErr += h->GetBinError(ix,iy+1)*h->GetBinError(ix,iy+1);
	ncnt++;
      }
      double down  = h->GetBinContent(ix,iy-1);
      if(fabs(down) > epsilon && iy > 1){
	sum += down;
	sumErr += h->GetBinError(ix,iy-1)*h->GetBinError(ix,iy-1);
	ncnt++;
      }
      double left  = h->GetBinContent(ix-1,iy);
      if(fabs(left) > epsilon && ix > 1){
	sum += left;
	sumErr += h->GetBinError(ix-1,iy)*h->GetBinError(ix-1,iy);
	ncnt++;
      }
      double right = h->GetBinContent(ix+1,iy);
      if(fabs(right) > epsilon && ix < nbinsX){
	sum += right;
	sumErr += h->GetBinError(ix+1,iy)*h->GetBinError(ix+1,iy);
	ncnt++;
      }
      if(ncnt > 0) {
	h->SetBinContent(ix,iy,sum/ncnt);
	h->SetBinError(ix,iy,TMath::Sqrt(sumErr)/ncnt);
      }
    } // for iy
  } // for ix

}



void fixBadCells(TH2D* h, int diagonal=0) {

  // fix bad cells which have wrong sign compared to 4 surrounding cells.
  // then assign average value from 4 cells

  int nbinsX = h->GetXaxis()->GetNbins();
  int nbinsY = h->GetYaxis()->GetNbins();

  for(int ix=1; ix <= nbinsX; ix++) {
    for(int iy=1; iy <= nbinsY - 2; iy++) {
      if(diagonal > 0) {
	double xx = h->GetXaxis()->GetBinCenter(ix);
	double yy = h->GetYaxis()->GetBinCenter(iy + 2);
	if(diagonal == 1 && xx < yy) continue;
	if(diagonal == 2 && xx > yy) continue;
      }
      double val = h->GetBinContent(ix,iy);
      if(val < epsilon) {
        int ncnt = 0;
        double up    = h->GetBinContent(ix,iy+1);
        if(up > epsilon) ncnt++;
        double down  = h->GetBinContent(ix,iy-1);
        if(down > epsilon) ncnt++;
        double left  = h->GetBinContent(ix-1,iy);
        if(left > epsilon) ncnt++;
        double right = h->GetBinContent(ix+1,iy);
        if(right > epsilon) ncnt++;
        if(ncnt == 4){
          val = (up+down+left+right)/ncnt;
          h->SetBinContent(ix,iy,val);
          up    = h->GetBinError(ix,iy+1);
          down  = h->GetBinError(ix,iy-1);
          left  = h->GetBinError(ix-1,iy);
          right = h->GetBinError(ix+1,iy);
          val = TMath::Sqrt(up*up + down*down + left*left + right*right)/ncnt;
          h->SetBinError(ix,iy,val);
        }
      }
      else {
        int ncnt = 0;
        double up    = h->GetBinContent(ix,iy+1);
        if(up < epsilon) ncnt++;
        double down  = h->GetBinContent(ix,iy-1);
        if(down < epsilon) ncnt++;
        double left  = h->GetBinContent(ix-1,iy);
        if(left < epsilon) ncnt++;
        double right = h->GetBinContent(ix+1,iy);
        if(right < epsilon) ncnt++;
        if(ncnt == 4){
          val = (up+down+left+right)/ncnt;
          h->SetBinContent(ix,iy,val);
          up    = h->GetBinError(ix,iy+1);
          down  = h->GetBinError(ix,iy-1);
          left  = h->GetBinError(ix-1,iy);
          right = h->GetBinError(ix+1,iy);
          val = TMath::Sqrt(up*up + down*down + left*left + right*right)/ncnt;
          h->SetBinError(ix,iy,val);
        }
      }
    } // for iy
  } // for ix
}



void RemovePoints(TGraph* g) {
  int N = g->GetN();
  int last = -1;
  for(int i=0; i<N; i++) {
    double x,y;
    g->GetPoint(i,x,y);
    //    if(abs(x-y) < 0.1) last = i;
    if(x > 600) {
      last = i+1;
      break;
    }
  }
  for(int i=N-1; i>=last; i--) g->RemovePoint(i);
}


void RemoveOffOrderedPoints(TGraph* g) {
  int i=g->GetN();
  while(i >= 1){
    double x,y;
    double xpre, ypre;
    g->GetPoint(i,x,y);
    g->GetPoint(i-1,xpre,ypre);
    if(x < xpre) g->RemovePoint(i);
    i--;
  }

}


void PrintPoints(TGraph* g) {
  int N = g->GetN();
  for(int i=0; i<N; i++) {
    double x,y;
    g->GetPoint(i,x,y);
    std::cout << "x, y = " << x << ", " << y << std::endl;
  }
}


void Print2DHist(TH2D* h) {

  int nbinsX = h->GetXaxis()->GetNbins();
  int nbinsY = h->GetYaxis()->GetNbins();

  for(int ix=1; ix <= nbinsX; ix++) {
    for(int iy=1; iy <= nbinsY; iy++) {
      double val = h->GetBinContent(ix,iy);
      std::cout << "(" << ix << "," << iy << ") = " << val << std::endl;
    }
  }

}


TGraph* makeBandGraph(TGraph* g1, TGraph* g2) {

  int n1 = g1->GetN();
  int n2 = g2->GetN();

  TGraph* graph = new TGraph(n1+n2);

  for(int i=0; i<n1; i++){
    double x,y;
    g1->GetPoint(i,x,y);
    graph->SetPoint(i,x,y);
  }
  for(int i=0; i<n2; i++){
    double x,y;
    g2->GetPoint(n2-1-i,x,y);
    graph->SetPoint(n1+i,x,y);
  }

  return graph;
}


TGraph* getContour(TH2D* h, TString name) {

  TGraph* graph = new TGraph(100);
  graph->SetName(name);
  int ip = 0;
  int nx = h->GetXaxis()->GetNbins();
  int ny = h->GetYaxis()->GetNbins();

  // for y>x
  int ix = -1;
  for(int j=ny;true; j--) {
    int k = -1;
    for(int i=2; i<nx-1; i++) {
      if(h->GetBinContent(i,j) < 0) {
	std::cout << "i,j,z : " << i << ", " << j << ", " << h->GetBinContent(i,j) << std::endl;
        k = i;
        break;
      }
    }// for i
    if(k<0) continue;
    double y = h->GetYaxis()->GetBinCenter(j);
    double x1 = h->GetXaxis()->GetBinCenter(k-1);
    double x2 = h->GetXaxis()->GetBinCenter(k);
    double z1 = h->GetBinContent(k-1,j);
    double z2 = h->GetBinContent(k,j);
    double x = x1 + (x2-x1)*fabs(z1)/fabs(z2-z1);
    std::cout << "y, x1, x2, z1, z2, x : " << y << ", " << x1 << ", " << x2 << ", " << x << ", " << z1 << ", " << z2 << std::endl;
    graph->SetPoint(ip++,x,y);

    if(h->GetYaxis()->GetBinCenter(j) < h->GetXaxis()->GetBinCenter(k)) {
      ix = k;
      break;
    }
  }// for j

  if(ix < 0) std::cout << "Something wrong...." << std::endl;

  // for y<x
  for(int i=ix; i<=nx; i++) {
    int k = -1;
    for(int j=2; j<ny-1; j++) {
      if(h->GetBinContent(i,j) < 0) {
        k = j;
        break;
      }
    }// for j
    if(k<0) continue;
    double x = h->GetXaxis()->GetBinCenter(i);
    double y1 = h->GetYaxis()->GetBinCenter(k-1);
    double y2 = h->GetYaxis()->GetBinCenter(k);
    double z1 = h->GetBinContent(i,k-1);
    double z2 = h->GetBinContent(i,k);
    double y = y1 + (y2-y1)*fabs(z1)/fabs(z2-z1);
    std::cout << "x, y1, y2, z1, z2, y : " << x << ", " << y1 << ", " << y2 << ", " << y << ", " << z1 << ", " << z2 << std::endl;
    graph->SetPoint(ip++,x,y);
  }// for i

  ip = graph->GetN()-1;
  while(1) {
    double x, y;
    graph->GetPoint(ip,x,y);
    if(x>1) break;
    else graph->RemovePoint(ip);
    ip--;
  }

  return graph;
}


void RemoveBadCells(TH2D* h, int diagonal=0) {

  // fix bad cells which have wrong sign compared to 4 surrounding cells.
  // then assign average value from 4 cells

  int nbinsX = h->GetXaxis()->GetNbins();
  int nbinsY = h->GetYaxis()->GetNbins();

  for(int ix=1; ix <= nbinsX; ix++) {
    for(int iy=1; iy <= nbinsY - 2; iy++) {
      double xx = h->GetXaxis()->GetBinCenter(ix);
      double yy = h->GetYaxis()->GetBinCenter(iy + 2);
      if(diagonal > 0) {
	if(diagonal == 1 && xx < yy) continue;
	if(diagonal == 2 && xx > yy) continue;
      }
      if(TMath::Abs(xx-yy) > 25) continue;
      double val = h->GetBinContent(ix,iy);
      double val_left = h->GetBinContent(ix-1,iy);
      double val_up = h->GetBinContent(ix,iy+1);
      if(val > 1 && val_left < 1 && val_up < 1) h->SetBinContent(ix,iy,(val_up+val_left)/2.0);
    }// for iy
  }// for ix

}

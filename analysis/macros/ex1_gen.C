/*
  ex1_gen.C
*/
#include <iostream>
#include <cstdio>
#include <cmath>
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"

Double_t f0(Double_t* x, Double_t* p) {
  Double_t y = 0.0;
  if (std::fabs(x[0])<1.0E-15) x[0] = 1.0E-15;
  Double_t a = p[0]/(x[0]*x[0]*x[0]*x[0]*x[0]);
  Double_t b = p[1]/x[0];
  Double_t c= std::exp(b)-1.0;
  y = a/c;
  return y;
}

Double_t f1(Double_t* x, Double_t* p) {
  Double_t y = 0.0;
  if (std::fabs(x[0])<1.0E-15) x[0] = 1.0E-15;
  Double_t a = p[0]/(x[0]*x[0]*x[0]*x[0]*x[0]);
  Double_t b = p[1]/x[0];
  y = a*std::exp(-b);
  return y;
}

void setStyle() {
  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.03);
  
}

void ex1_gen() {
  setStyle();
  std::string fn1 = "work/ex1.dat";
  Double_t pars[2] = { 1.0, 10.0 };
  Double_t xmin=0.0, xmax=20.0;
  Double_t ymin = 1.0E-10;
  Double_t ymax = 0.4;
  
  TF1* func = new TF1("f0", f0, xmin, xmax, 2);
  func->SetParameters(pars);
  Double_t s = func->Integral(xmin, xmax);
  func->SetParameter(0, 1.0/s);

  TCanvas* c = new TCanvas("c", "", 500, 500);
  TH1* hframe = c->DrawFrame(xmin, ymin, xmax, ymax);
  func->Draw("same");

  TF1* func1 = new TF1("f1", f1, xmin, xmax, 2);
  for (int i=0; i<2; ++i) {
    func1->SetParameter(i, func->GetParameter(i));
  }
  func1->SetLineColor(kGreen);
  func1->SetLineStyle(2);
  func1->Draw("same");
  
  unsigned int seed = 20220930;
  TRandom3 random(seed);
  //  Double_t esyst=0.0015;
  Double_t esyst=0.001;
  
  int n=20;
  int n1=0;
  Double_t dx = (xmax-xmin)/n;
  Double_t ax[20];
  Double_t ay[20];
  Double_t aex[20];
  Double_t aey[20];
  for (int i=0; i<n-3; ++i) {
    Double_t x = dx * (i+1);
    Double_t y0 = func->Eval(x);
    Double_t estat = y0 * 0.1;
    Double_t etot = std::sqrt(esyst*esyst + estat*estat);
    // std::cout << "x,y,e=" << x << ", " << y0 << ", " << etot
    // 	      << "[ " << estat << ", " << esyst << "]" << std::endl;
    Double_t y = y0 + random.Gaus(0.0, etot);
    std::printf("%5.3e %5.3e\n", x, y);
    ax[i] = x;
    ay[i] = y;
    aex[i] = 0.0;
    aey[i] = etot;
    n1 ++;
  }
  TGraph* g = new TGraphErrors(n1, ax, ay, aex, aey);
  g->SetMarkerColor(kBlack);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.8);
  g->Draw("p");
}


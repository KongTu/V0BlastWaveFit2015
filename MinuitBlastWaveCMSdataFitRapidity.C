#include <iostream>
#include <fstream>
#include <math.h>
#include "TMath.h"
#include <stdio.h>
#include <iomanip>
#include <vector>
#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "makeMultiPanelCanvas.C"

using namespace std;

vector<double> x,ex, y,ey,x1,ex1,y1,ey1;

double getChi2;
int getNdf;

const double R = 1.0;
const double dr = 1e-2; // FIXME

double integral(double beta_T, double T, double n, double pt, double mt)
{
  
  double s = 0.;
  for(double r = dr/2; r < R; r += dr)
  {
    double beta_r = (beta_T*(n+2)/2) * TMath::Power((r/R),n);
    double rho = TMath::ATanH(beta_r);

    s += r * dr * TMath::BesselK1((mt*cosh(rho))/T) * TMath::BesselI0((pt*sinh(rho))/T);
  }

  return s;
}

void function(int &npar, double *gin, double &f, double *par, int flag)
{
  double aka    = par[0];
  double aka1   = par[4];
  double T      = par[2];
  double beta_T = par[1];
  double n = par[3];

  double chi2 = 0.;

  for(int i = 0; i < int(x.size()); i++)
  {
    double m = 0.497;

    double pt = x[i];
    double mt = sqrt(m*m + pt*pt);

    double a = 0.;
    a = aka;
  
    double theo = a * mt * integral(beta_T, T, n, pt, mt);
    double q = (y[i] - theo)/ey[i];

    chi2 += q*q;
  }

  for(i = 0; i < int(x1.size()); i++){

  	double m1 = 1.115;

  	double pt1 = x1[i];
  	double mt1 = sqrt(m1*m1 + pt1*pt1);

  	double a1 = 0.;
  	a1 = aka1;

  	double theo1 = a1 * mt1 * integral(beta_T,T,n,pt1,mt1);
  	double q1 = (y1[i]-theo1)/ey1[i];

  	chi2 += q1*q1;
  }

  f = chi2;

  getChi2 = f;
  getNdf  = x.size() + x1.size() - 5 - 1;
}

double MyFunc( double *pt, double *p){

  double mass = 0.497;

  double temp = 0.;
    double mt = sqrt(pt[0]*pt[0]+mass*mass);

    temp = p[2] * mt * integral(p[0],p[1],p[3],pt[0],mt);

    return temp;
  
}

double MyFunc1( double *pt, double *p){

  double mass = 1.115;

  double temp = 0.;
    double mt = sqrt(pt[0]*pt[0]+mass*mass);

    temp = p[2] * mt * integral(p[0],p[1],p[3],pt[0],mt);

    return temp;
  
}

double beta_T(double beta_s,double n){

  return (beta_s*2)/(n+2);

}


void MinuitBlastWaveCMSdataFitRapidity(){

  gStyle->SetErrorX(0);


  TFile* file = new TFile("/Users/kongkong/2015Research/Spectra/pPb2013/rpyDependence/files/220plus_rpyDependence_noEffError_v3.root");
  
  TH1D* ksSpectra[5];
  TH1D* laSpectra[5];

  stringstream ksHistName;
  stringstream laHistName;

  for (int rpy = 0; rpy < 5; rpy++){

    ksHistName.str("");
    laHistName.str("");

    ksHistName << "ksSpectra_vtx_";
    ksHistName << rpy+1;

    laHistName << "laSpectra_vtx_";
    laHistName << rpy+1;

    ksSpectra[rpy] = (TH1D*)file->Get(ksHistName.str().c_str());
    laSpectra[rpy] = (TH1D*)file->Get(laHistName.str().c_str());
    

  }

  double ks_ptbins[26] = {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double ks_ptbincenter[25] = {0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};
  double la_ptbins[21] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double la_ptbincenter[20] = {0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};

  TCanvas* c1 = new TCanvas();
    c1->Divide(2,3,0,0);

  double beta_T[5],Tkin[5],aka[5],aka1[5],n[5];
  double ebeta_T[5],eTkin[5],eaka[5],eaka1[5],en[5];

  stringstream f1Name;
  stringstream f2Name;

  TF1* f1[5];
  TF1* f2[5];

  TH1D* f1_hist[8];
  TH1D* f2_hist[8];

  stringstream f1HistName;
  stringstream f2HistName;

  double ks_ptbins_histo[28] = {0.0,0.1,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
  double la_ptbins_histo[22] = {0.0,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};


  for(int mult = 0; mult < 8; mult++){

    f1HistName.str("");
    f2HistName.str("");

    f1HistName << "ks_fitFuncHist_";
    f1HistName << mult+1;

    f2HistName << "la_fitFuncHist_";
    f2HistName << mult+1;

    f1_hist[mult] = new TH1D(f1HistName.str().c_str(),f1HistName.str().c_str(),27,ks_ptbins_histo);
    f2_hist[mult] = new TH1D(f2HistName.str().c_str(),f2HistName.str().c_str(),21,la_ptbins_histo);

  }

  TLine* l1 = new TLine(0,1,9.0,1.0);
  l1->SetLineWidth(2);
  l1->SetLineColor(kRed);
  l1->SetLineStyle(2);

  TLatex* ratio[5];
  ratio[0] = new TLatex(4,0.17,"-2.4 < y^{}_{cm} < -1.5");
  ratio[1] = new TLatex(4,0.17,"-1.5 < y^{}_{cm} < -0.8");
  ratio[2] = new TLatex(4,0.17,"-0.8 < y^{}_{cm} < 0");
  ratio[3] = new TLatex(4,0.17,"0 < y^{}_{cm} < 0.8");
  ratio[4] = new TLatex(4,0.17,"0.8 < y^{}_{cm} < 1.5");

  TLegend *w1 = new TLegend(0.25,0.4,0.5,0.5);
  w1->SetLineColor(kWhite);
  w1->SetFillColor(0);

  w1->AddEntry(ksSpectra[4],"K^{0}_{s}");
  w1->AddEntry(laSpectra[4],"#Lambda/#bar{#Lambda}");
  
  TGraph* test[5];

  for(int rpy = 0; rpy < 5; rpy++){

      x.clear();
      ex.clear();
      y.clear();
      ey.clear();

      x1.clear();
      ex1.clear();
      y1.clear();
      ey1.clear();

      for (int pt = 2; pt < 10; pt++){

        x.push_back( ks_ptbincenter[pt] );
        ex.push_back(0.0);
        y.push_back( ksSpectra[rpy]->GetBinContent(pt+3) );
        ey.push_back( ksSpectra[rpy]->GetBinError(pt+3) );
        
      }

      for(pt = 0; pt < 14; pt++){

        x1.push_back( la_ptbincenter[pt] );
        ex1.push_back(0.0);
        y1.push_back( laSpectra[rpy]->GetBinContent(pt+2) );
        ey1.push_back( laSpectra[rpy]->GetBinError(pt+2) );
        
      }

  /**
   * simultaneous fit for K0s and Lambda:
   */

      TMinuit * gMinuit[5];
      gMinuit[rpy] = new TMinuit(5);

      //set the function to minimize with minuit;
      gMinuit[rpy]->SetFCN(function);

      double arglist[10];
      arglist[0] = -1;
      gMinuit[rpy]->mnexcm("SET ERR", arglist, 1, 0);

      gMinuit[rpy]->mnparm(0, "aka",    10,    0.1,  1,    100000, 0);
      gMinuit[rpy]->mnparm(1, "beta_T", 0.4,  0.01, 0.2,   0.70,   0);
      gMinuit[rpy]->mnparm(2, "Tkin",   0.16, 0.01,  0.10, 0.3,   0);
      gMinuit[rpy]->mnparm(3, "n",      1.0,   0.01, 0.1,  10.0,  0);
      gMinuit[rpy]->mnparm(4, "aka1",   10,    0.1,  1,    100000, 0);

      gMinuit[rpy]->mnexcm("MIGRAD",    arglist,  1,   0);
      gMinuit[rpy]->mnexcm("CALL FCN",  arglist,  1,   0);
      //gMinuit[rpy]->mnexcm("HESSE",     arglist,  1,   0);

      gMinuit[rpy]->GetParameter(0, aka[rpy],    eaka[rpy]);
      gMinuit[rpy]->GetParameter(1, beta_T[rpy], ebeta_T[rpy]);
      gMinuit[rpy]->GetParameter(2, Tkin[rpy],   eTkin[rpy]);
      gMinuit[rpy]->GetParameter(3, n[rpy],      en[rpy]);
      gMinuit[rpy]->GetParameter(4, aka1[rpy],   eaka1[rpy]);

      gMinuit[rpy]->SetErrorDef(1);
      test[rpy] = (TGraph*)gMinuit[rpy]->Contour(100,1,2);

  /**
   * define the fit function by using the parameters from fit:
   */

      double xmin = 0.3;
      double xmax = 1.6;

      double xmin1 = 0.6;
      double xmax1 = 3.0;

      f1Name.str("");
      f2Name.str("");

      f1Name << "f1_";
      f1Name << rpy+1;

      f2Name << "f2_";
      f2Name << rpy+1;
  
      f1[rpy] = new TF1(f1Name.str().c_str(),MyFunc,xmin,xmax,4);
      f1[rpy]->SetParameter(0,beta_T[rpy]);
      f1[rpy]->SetParameter(1,Tkin[rpy]);
      f1[rpy]->SetParameter(2,aka[rpy]);
      f1[rpy]->SetParameter(3,n[rpy]);

      f2[rpy] = new TF1(f2Name.str().c_str(),MyFunc1,xmin1,xmax1,4);
      f2[rpy]->SetParameter(0,beta_T[rpy]);
      f2[rpy]->SetParameter(1,Tkin[rpy]);
      f2[rpy]->SetParameter(2,aka1[rpy]);
      f2[rpy]->SetParameter(3,n[rpy]);

      c1->cd(rpy+1);
      gPad->SetLogy(1);

      ksSpectra[rpy]->SetMarkerSize(1.3);
      ksSpectra[rpy]->SetMarkerStyle(20);
      ksSpectra[rpy]->SetMarkerColor(kBlue);
      ksSpectra[rpy]->SetStats(kFALSE);
      ksSpectra[rpy]->SetTitle("");
      ksSpectra[rpy]->SetXTitle("P^{}_{T,V0}(GeV/c)");
      ksSpectra[rpy]->SetYTitle("1/N^{}_{ev}1/(2#PiP^{}_{T}d^{2}N/(dP^{}_{T}dy) [(GeV/c)^{-2}]");

      ksSpectra[rpy]->GetYaxis()->SetRangeUser(0.0000001,1000);

      laSpectra[rpy]->SetMarkerSize(1.3);
      laSpectra[rpy]->SetMarkerStyle(20);
      laSpectra[rpy]->SetMarkerColor(kBlack);

      ksSpectra[rpy]->Draw("P");
      laSpectra[rpy]->Draw("Psame");
      ratio[rpy]->Draw("same");

      f2[rpy]->Draw("same");
      f1[rpy]->Draw("same");

      for(int pt = 0; pt < 10; pt++){

        double mass = 0.497;

        double temp = 0.;
        double mt = sqrt(ks_ptbincenter[pt]*ks_ptbincenter[pt]+mass*mass);

        temp = aka[rpy] * mt * integral(beta_T[rpy],Tkin[rpy],n[rpy],ks_ptbincenter[pt],mt);

        f1_hist[rpy]->SetBinContent(pt+3, temp);
        f1_hist[rpy]->SetBinError(pt+3, 0.0);

      }

      for(int pt = 0; pt < 12; pt++){

        double mass = 1.115;

        double temp = 0.;
        double mt = sqrt(la_ptbincenter[pt]*la_ptbincenter[pt]+mass*mass);

        temp = aka1[rpy] * mt * integral(beta_T[rpy],Tkin[rpy],n[rpy],la_ptbincenter[pt],mt);

        f2_hist[rpy]->SetBinContent(pt+2, temp);
        f2_hist[rpy]->SetBinError(pt+2, 0.0);

      }


  }

  w1->Draw("same");

  TCanvas *p1 = new TCanvas("p1","p1",1,1,1300,400);

  p1->Range(0,0,1,1);
  TPad* pad1[10];

  pad1[0] = new TPad("pad10", "pad10",0.0,   0.2, 0.2,   1);
  pad1[1] = new TPad("pad11", "pad11",0.2,   0.2, 0.4,   1);
  pad1[2] = new TPad("pad12", "pad12",0.4,   0.2, 0.6,   1);
  pad1[3] = new TPad("pad13", "pad13",0.6,   0.2, 0.8,   1);
  pad1[4] = new TPad("pad13", "pad13",0.8,   0.2, 1.0,   1);

  pad1[5] = new TPad("pad18", "pad18",0.0,   0.0, 0.2,   0.2);
  pad1[6] = new TPad("pad19", "pad19",0.2,   0.0, 0.4,   0.2);
  pad1[7] = new TPad("pad110", "pad110",0.4, 0.0, 0.6,   0.2);
  pad1[8] = new TPad("pad111", "pad111",0.6, 0.0, 0.8,   0.2);
  pad1[9] = new TPad("pad111", "pad111",0.8, 0.0, 1.0,   0.2);

  for(int i = 0; i < 10; i++){

  pad1[i]->SetLeftMargin(0.0);
  pad1[i]->SetRightMargin(0);
  pad1[i]->SetTopMargin(0.0);
  pad1[i]->SetBottomMargin(0);

  //pad1[i]->SetLogy(1);
  pad1[i]->Draw();

  }

  pad1[0]->SetLeftMargin(0.165);
  pad1[5]->SetLeftMargin(0.165);

  pad1[4]->SetRightMargin(0.05);
  pad1[9]->SetRightMargin(0.05);

  pad1[0]->SetTopMargin(0.02);
  pad1[1]->SetTopMargin(0.02);
  pad1[2]->SetTopMargin(0.02);
  pad1[3]->SetTopMargin(0.02);
  pad1[4]->SetTopMargin(0.02);

  pad1[5]->SetBottomMargin(0.20);
  pad1[6]->SetBottomMargin(0.20);
  pad1[7]->SetBottomMargin(0.20);
  pad1[8]->SetBottomMargin(0.20);
  pad1[9]->SetBottomMargin(0.20);

  for(rpy = 0; rpy < 5; rpy++){

    pad1[rpy]->cd();
    pad1[rpy]->SetLogy(1);
    ksSpectra[rpy]->Draw();
    laSpectra[rpy]->Draw("same");
    f1[rpy]->Draw("same");
    f2[rpy]->Draw("same");
    ratio[rpy]->Draw("same");
  }

  w1->Draw("same");

  for(rpy = 0; rpy < 5; rpy++){

    pad1[rpy+5]->cd();
    f1_hist[rpy]->Divide( ksSpectra[rpy] );
    f1_hist[rpy]->SetMarkerStyle(20);
    f1_hist[rpy]->SetMarkerColor(kBlue);
    f1_hist[rpy]->SetMarkerSize(1.3);
    f1_hist[rpy]->SetLineColor(kBlue);
    f1_hist[rpy]->SetTitle(" ");
    f1_hist[rpy]->SetStats(kFALSE);
    f1_hist[rpy]->GetYaxis()->SetTitle("Fit/data");
    f1_hist[rpy]->GetXaxis()->SetTitle("p^{}_{T} (GeV/c)");
    f1_hist[rpy]->GetYaxis()->SetTitleSize(0.09);
    f1_hist[rpy]->GetYaxis()->SetTitleOffset(0.63);
    f1_hist[rpy]->GetYaxis()->SetLabelSize(0.09);
    f1_hist[rpy]->GetYaxis()->SetLabelOffset(0.025);
    f1_hist[rpy]->GetYaxis()->SetRangeUser(0.6,1.4);

    f1_hist[rpy]->GetXaxis()->SetTitleSize(0.09);
    f1_hist[rpy]->GetXaxis()->SetLabelSize(0.10);

    f2_hist[rpy]->Divide( laSpectra[rpy] );
    f2_hist[rpy]->SetMarkerStyle(20);
    f2_hist[rpy]->SetMarkerColor(kBlack);
    f2_hist[rpy]->SetMarkerSize(1.3);
    f2_hist[rpy]->SetLineColor(kBlack);
    f2_hist[rpy]->SetTitle(" ");
    f2_hist[rpy]->SetStats(kFALSE);
    f2_hist[rpy]->GetYaxis()->SetTitle("Fit/data");
    f2_hist[rpy]->GetXaxis()->SetTitle("p^{}_{T} (GeV/c)");
    f2_hist[rpy]->GetYaxis()->SetLabelSize(0.09);
    f2_hist[rpy]->GetYaxis()->SetLabelOffset(0.025);
    f2_hist[rpy]->GetYaxis()->SetRangeUser(0.6,1.4);
    
    f1_hist[rpy]->Draw("P");
    f2_hist[rpy]->Draw("Psame");
    l1->Draw("same");
  }

  TCanvas* s1 = new TCanvas();

  TH1D* hist = new TH1D("h1","h1",100,0.2,0.7);
  hist->GetYaxis()->SetRangeUser(0.14,0.23);
  hist->GetXaxis()->SetRangeUser(0.4,0.7);
  hist->SetYTitle("T^{}_{kin}(GeV)");
  hist->SetXTitle("<#beta^{}_{T}>");

  TGraph* g = new TGraph(5);

  g->SetMarkerStyle(20);

  for(rpy = 0; rpy < 5; rpy++){

    g->SetPoint(rpy,beta_T[rpy],Tkin[rpy]);

  }

  hist->Draw("P");
  g->Draw("Psame");

  double ellipse_x[5][101];
  double ellipse_y[5][101];

  for(rpy = 0; rpy < 5; rpy++){

    for(int i = 0; i < 100; i++){

        test[rpy]->GetPoint(i,ellipse_x[rpy][i],ellipse_y[rpy][i]);

    } 

      ellipse_x[rpy][100] = ellipse_x[rpy][0];
      ellipse_y[rpy][100] = ellipse_y[rpy][0];  
  }

  TGraph* g1[5];

  for(rpy = 0; rpy < 5; rpy++){

    g1[rpy] = new TGraph(101);
      
      for(i = 0; i < 101; i++){
          
          g1[rpy]->SetPoint(i,ellipse_x[rpy][i],ellipse_y[rpy][i]);
          g1[rpy]->SetLineColor(kRed);
          g1[rpy]->SetMarkerColor(kRed);

      }

      g1[rpy]->Draw("same");
  }





}
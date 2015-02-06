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

#include "fitting.h"

using namespace RooFit;

using namespace std;

vector<double> x,ex, y,ey,x1,ex1,y1,ey1;

double getChi2;
int getNdf;

const double R = 1.0;
const double dr = 1e-2; // FIXME

double integral(double beta_s, double T, double n, double pt, double mt)
{
  
  double s = 0.;
  for(double r = dr/2; r < R; r += dr)
  {
    double beta_r = beta_s * TMath::Power((r/R),n);
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
  double beta_s = par[1];
  double n = par[3];

  double chi2 = 0.;

  for(int i = 0; i < int(x.size()); i++)
  {
    double m = 0.497;

    double pt = x[i];
    double mt = sqrt(m*m + pt*pt);

    double a = 0.;
    a = aka;
  
    double theo = a * mt * integral(beta_s, T, n, pt, mt);
    double q = (y[i] - theo)/ey[i];

    chi2 += q*q;
  }

  for(i = 0; i < int(x1.size()); i++){

  	double m1 = 1.115;

  	double pt1 = x1[i];
  	double mt1 = sqrt(m1*m1 + pt1*pt1);

  	double a1 = 0.;
  	a1 = aka1;

  	double theo1 = a1 * mt1 * integral(beta_s,T,n,pt1,mt1);
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


void MinuitBlastWaveCMSdataFit(){

  gStyle->SetErrorX(0);

  TFile* file = new TFile("/Users/kongkong/2015Research/Spectra/pPb2013/rpyDependence/files/MidRpy0_0.8_v1.root");
  
  TH1D* ksSpectra[8];
  TH1D* laSpectra[8];

  stringstream ksHistName;
  stringstream laHistName;

  for (int mult = 0; mult < 8; mult++){

    ksHistName.str("");
    laHistName.str("");

    ksHistName << "ksSpectra_vtx_";
    ksHistName << mult+1;

    laHistName << "laSpectra_vtx_";
    laHistName << mult+1;

    ksSpectra[mult] = (TH1D*)file->Get(ksHistName.str().c_str());
    laSpectra[mult] = (TH1D*)file->Get(laHistName.str().c_str());
    //laSpectra[mult]->Scale(8);

  }

  double ks_ptbins[26] = {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double ks_ptbincenter[25] = {0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};
  double la_ptbins[21] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double la_ptbincenter[20] = {0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};

  TCanvas* c1 = new TCanvas();
    c1->Divide(4,2,0,0);

  double beta_s[8],Tkin[8],aka[8],aka1[8],n[8];
  double ebeta_s[8],eTkin[8],eaka[8],eaka1[8],en[8];

  stringstream f1Name;
  stringstream f2Name;

  TF1* f1[8];
  TF1* f2[8];

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

  TLatex* ratio[8];
  ratio[0] = new TLatex(2.5,45,"0 < N^{offline}_{trk} < 35");
  ratio[1] = new TLatex(2.5,45,"35 < N^{offline}_{trk} < 60");
  ratio[2] = new TLatex(2.5,45,"60 < N^{offline}_{trk} < 90");
  ratio[3] = new TLatex(2.5,45,"90 < N^{offline}_{trk} < 120");
  ratio[4] = new TLatex(2.5,45,"120 < N^{offline}_{trk} < 150");
  ratio[5] = new TLatex(2.5,45,"150 < N^{offline}_{trk} < 185");
  ratio[6] = new TLatex(2.5,45,"185 < N^{offline}_{trk} < 220");
  ratio[7] = new TLatex(2.5,45,"220 < N^{offline}_{trk} < #infty");
  //ratio[8] = new TLatex(2.5,45,"N^{offline}_{trk} > 260");
  //
  TLine* l1 = new TLine(0,1,9.0,1.0);
  l1->SetLineWidth(2);
  l1->SetLineColor(kRed);
  l1->SetLineStyle(2);

  TLegend *w1 = new TLegend(0.25,0.4,0.5,0.5);
  w1->SetLineColor(kWhite);
  w1->SetFillColor(0);

  w1->AddEntry(ksSpectra[7],"K^{0}_{s}");
  w1->AddEntry(laSpectra[7],"#Lambda/#bar{#Lambda}");

  TGraph* test[8];
  

  for(int mult = 0; mult < 8; mult++){

      x.clear();
      ex.clear();
      y.clear();
      ey.clear();

      x1.clear();
      ex1.clear();
      y1.clear();
      ey1.clear();

      //ks pt strats at 0.3-0.4
      //la pt starts at 0.6-0.8
      //all correspond to pt = 0;

      for (int pt = 0; pt < 14; pt++){

        x.push_back( ks_ptbincenter[pt] );
        ex.push_back(0.0);
        y.push_back( ksSpectra[mult]->GetBinContent(pt+3) );
        ey.push_back( ksSpectra[mult]->GetBinError(pt+3) );
        //ey.push_back( 0.1 );
      }

      for(pt = 0; pt < 14; pt++){

          x1.push_back( la_ptbincenter[pt] );
          ex1.push_back(0.0);
          y1.push_back( laSpectra[mult]->GetBinContent(pt+2) );
          ey1.push_back( laSpectra[mult]->GetBinError(pt+2) );
          //ey1.push_back( 0.1 );
      }

  /**
   * simultaneous fit for K0s and Lambda:
   */

      TMinuit * gMinuit[9];
      gMinuit[mult] = new TMinuit(5);

      //set the function to minimize with minuit;
      gMinuit[mult]->SetFCN(function);

      double arglist[10];
      arglist[4] = 0.001;
      gMinuit[mult]->mnexcm("SET ERR", arglist, 1, 0);

      gMinuit[mult]->mnparm(0, "aka",    10,   0.1,  1,    10000, 0);
      gMinuit[mult]->mnparm(1, "beta_s", 0.70, 0.01, 0.15,   1.0,   0);
      gMinuit[mult]->mnparm(2, "Tkin",   0.15, 0.01, 0.05,   1.0,   0);
      gMinuit[mult]->mnparm(3, "n",      1.0,  0.01, 0.1,    10.0,  0);
      gMinuit[mult]->mnparm(4, "aka1",   10,   0.1,  1,    10000, 0);

      gMinuit[mult]->mnexcm("MIGRAD",    arglist,  1,   0);
      gMinuit[mult]->mnexcm("CALL FCN",  arglist,  1,   0);
      //gMinuit[mult]->mnexcm("HESSE",     arglist,  1,   0);

      gMinuit[mult]->GetParameter(0, aka[mult],    eaka[mult]);
      gMinuit[mult]->GetParameter(1, beta_s[mult], ebeta_s[mult]);
      gMinuit[mult]->GetParameter(2, Tkin[mult],   eTkin[mult]);
      gMinuit[mult]->GetParameter(3, n[mult],      en[mult]);
      gMinuit[mult]->GetParameter(4, aka1[mult],   eaka1[mult]);

      gMinuit[mult]->SetErrorDef(4);
      test[mult] = (TGraph*)gMinuit[mult]->Contour(100,1,2);
      //test->SetFillColor(42);

  /**
   * define the fit function by using the parameters from fit:
   */

      double xmin = 0.45;
      double xmax = 2.3;

      double xmin1 = 1.0;
      double xmax1 = 4.0;

      f1Name.str("");
      f2Name.str("");

      f1Name << "f1_";
      f1Name << mult+1;

      f2Name << "f2_";
      f2Name << mult+1;
  
      f1[mult] = new TF1(f1Name.str().c_str(),MyFunc,xmin,xmax,4);
      f1[mult]->SetParameter(0,beta_s[mult]);
      f1[mult]->SetParameter(1,Tkin[mult]);
      f1[mult]->SetParameter(2,aka[mult]);
      f1[mult]->SetParameter(3,n[mult]);

      f2[mult] = new TF1(f2Name.str().c_str(),MyFunc1,xmin1,xmax1,4);
      f2[mult]->SetParameter(0,beta_s[mult]);
      f2[mult]->SetParameter(1,Tkin[mult]);
      f2[mult]->SetParameter(2,aka1[mult]);
      f2[mult]->SetParameter(3,n[mult]);

      c1->cd(mult+1);
      //gPad->SetLogy(1);

      ksSpectra[mult]->SetMarkerSize(1.3);
      ksSpectra[mult]->SetMarkerStyle(20);
      ksSpectra[mult]->SetMarkerColor(kBlue);
      ksSpectra[mult]->SetStats(kFALSE);
      ksSpectra[mult]->SetTitle("");
      ksSpectra[mult]->SetXTitle("P^{}_{T,V0}(GeV/c)");
      ksSpectra[mult]->SetYTitle("1/N^{}_{ev}1/(2#PiP^{}_{T}d^{2}N/(dP^{}_{T}dy) [(GeV/c)^{-2}]");

      ksSpectra[mult]->GetYaxis()->SetRangeUser(0.0000001,1000);

      laSpectra[mult]->SetMarkerSize(1.3);
      laSpectra[mult]->SetMarkerStyle(20);
      laSpectra[mult]->SetMarkerColor(kBlack);

      ksSpectra[mult]->Draw("P");
      laSpectra[mult]->Draw("Psame");
      ratio[mult]->Draw("same");

      f2[mult]->Draw("same");
      f1[mult]->Draw("same");

      for(int pt = 0; pt < 14; pt++){
      
        double mass = 0.497;

        double temp = 0.;
        double mt = sqrt(ks_ptbincenter[pt]*ks_ptbincenter[pt]+mass*mass);

        temp = aka[mult] * mt * integral(beta_s[mult],Tkin[mult],n[mult],ks_ptbincenter[pt],mt);

        f1_hist[mult]->SetBinContent(pt+3, temp);

      }

      for(int pt = 0; pt < 14; pt++){
      
        double mass = 1.115;

        double temp = 0.;
        double mt = sqrt(la_ptbincenter[pt]*la_ptbincenter[pt]+mass*mass);

        temp = aka1[mult] * mt * integral(beta_s[mult],Tkin[mult],n[mult],la_ptbincenter[pt],mt);

        f2_hist[mult]->SetBinContent(pt+2, temp);

      }

      
  }


  w1->Draw("same");

  TCanvas *p1 = new TCanvas("p1","p1",1,1,1200,400);

  p1->Range(0,0,1,1);
  TPad* pad1[8];

  pad1[0] = new TPad("pad10", "pad10",0.0,   0.2, 0.3,   1);
  pad1[1] = new TPad("pad11", "pad11",0.3,   0.2, 0.533, 1);
  pad1[2] = new TPad("pad12", "pad12",0.533, 0.2, 0.766, 1);
  pad1[3] = new TPad("pad13", "pad13",0.766, 0.2, 1.0,   1);

  pad1[4] = new TPad("pad18", "pad18",0.0,      0.0, 0.3,   0.2);
  pad1[5] = new TPad("pad19", "pad19",0.3,      0.0, 0.533, 0.2);
  pad1[6] = new TPad("pad110", "pad110",0.533,  0.0, 0.766, 0.2);
  pad1[7] = new TPad("pad111", "pad111",0.766,  0.0, 1.0,   0.2);

  for(int i = 0; i < 8; i++){

  pad1[i]->SetLeftMargin(0.0);
  pad1[i]->SetRightMargin(0);
  pad1[i]->SetTopMargin(0.0);
  pad1[i]->SetBottomMargin(0);

  //pad1[i]->SetLogy(1);
  pad1[i]->Draw();

  }

  pad1[0]->SetLeftMargin(0.165);
  pad1[4]->SetLeftMargin(0.165);

  pad1[3]->SetRightMargin(0.05);
  pad1[7]->SetRightMargin(0.05);

  pad1[0]->SetTopMargin(0.02);
  pad1[1]->SetTopMargin(0.02);
  pad1[2]->SetTopMargin(0.02);
  pad1[3]->SetTopMargin(0.02);

  pad1[4]->SetBottomMargin(0.20);
  pad1[5]->SetBottomMargin(0.20);
  pad1[6]->SetBottomMargin(0.20);
  pad1[7]->SetBottomMargin(0.20);

  for(mult = 0; mult < 4; mult++){

    pad1[mult]->cd();
    pad1[mult]->SetLogy(1);
    ksSpectra[mult]->Draw();
    laSpectra[mult]->Draw("same");
    f1[mult]->Draw("same");
    f2[mult]->Draw("same");
    ratio[mult]->Draw("same");
  }

  w1->Draw("same");

  for(mult = 0; mult < 4; mult++){

    pad1[mult+4]->cd();
    f1_hist[mult]->Divide( ksSpectra[mult] );
    f1_hist[mult]->SetMarkerStyle(20);
    f1_hist[mult]->SetMarkerColor(kBlue);
    f1_hist[mult]->SetMarkerSize(1.3);
    f1_hist[mult]->SetLineColor(kWhite);
    f1_hist[mult]->SetTitle(" ");
    f1_hist[mult]->SetStats(kFALSE);
    f1_hist[mult]->GetYaxis()->SetTitle("Fit/data");
    f1_hist[mult]->GetXaxis()->SetTitle("p^{}_{T} (GeV/c)");
    f1_hist[mult]->GetYaxis()->SetTitleSize(0.09);
    f1_hist[mult]->GetYaxis()->SetTitleOffset(0.63);
    f1_hist[mult]->GetYaxis()->SetLabelSize(0.09);
    f1_hist[mult]->GetYaxis()->SetLabelOffset(0.025);
    f1_hist[mult]->GetYaxis()->SetRangeUser(0.6,1.4);

    f1_hist[mult]->GetXaxis()->SetTitleSize(0.09);
    f1_hist[mult]->GetXaxis()->SetLabelSize(0.10);

    f2_hist[mult]->Divide( laSpectra[mult] );
    f2_hist[mult]->SetMarkerStyle(20);
    f2_hist[mult]->SetMarkerColor(kBlack);
    f2_hist[mult]->SetMarkerSize(1.3);
    f2_hist[mult]->SetLineColor(kWhite);
    f2_hist[mult]->SetTitle(" ");
    f2_hist[mult]->SetStats(kFALSE);
    f2_hist[mult]->GetYaxis()->SetTitle("Fit/data");
    f2_hist[mult]->GetXaxis()->SetTitle("p^{}_{T} (GeV/c)");
    f2_hist[mult]->GetYaxis()->SetLabelSize(0.09);
    f2_hist[mult]->GetYaxis()->SetLabelOffset(0.025);
    f2_hist[mult]->GetYaxis()->SetRangeUser(-10,40);
    
    f1_hist[mult]->Draw("P");
    f2_hist[mult]->Draw("Psame");
    l1->Draw("same");
  }



  TCanvas *p2 = new TCanvas("p2","p2",1,1,1200,400);

  p2->Range(0,0,1,1);
  TPad* pad2[8];

  pad2[0] = new TPad("pad20", "pad20",0.0,   0.2, 0.3,   1);
  pad2[1] = new TPad("pad21", "pad21",0.3,   0.2, 0.533, 1);
  pad2[2] = new TPad("pad22", "pad22",0.533, 0.2, 0.766, 1);
  pad2[3] = new TPad("pad23", "pad23",0.766, 0.2, 1.0,   1);

  pad2[4] = new TPad("pad28", "pad28",0.0,      0.0, 0.3,   0.2);
  pad2[5] = new TPad("pad29", "pad29",0.3,      0.0, 0.533, 0.2);
  pad2[6] = new TPad("pad210", "pad210",0.533,  0.0, 0.766, 0.2);
  pad2[7] = new TPad("pad211", "pad211",0.766,  0.0, 1.0,   0.2);

  for(int i = 0; i < 8; i++){

  pad2[i]->SetLeftMargin(0.0);
  pad2[i]->SetRightMargin(0);
  pad2[i]->SetTopMargin(0.0);
  pad2[i]->SetBottomMargin(0);
  pad2[i]->Draw();

  }

  pad2[0]->SetLeftMargin(0.165);
  pad2[4]->SetLeftMargin(0.165);

  pad2[3]->SetRightMargin(0.05);
  pad2[7]->SetRightMargin(0.05);

  pad2[0]->SetTopMargin(0.02);
  pad2[1]->SetTopMargin(0.02);
  pad2[2]->SetTopMargin(0.02);
  pad2[3]->SetTopMargin(0.02);

  pad2[4]->SetBottomMargin(0.20);
  pad2[5]->SetBottomMargin(0.20);
  pad2[6]->SetBottomMargin(0.20);
  pad2[7]->SetBottomMargin(0.20);

  for(mult = 0; mult < 4; mult++){

    pad2[mult]->cd();
    pad2[mult]->SetLogy(1);
    ksSpectra[mult+4]->Draw();
    laSpectra[mult+4]->Draw("same");
    f1[mult+4]->Draw("same");
    f2[mult+4]->Draw("same");
    ratio[mult+4]->Draw("same");
  }

  w1->Draw("same");

  for(mult = 4; mult < 8; mult++){

    pad2[mult]->cd();
    f1_hist[mult]->Divide( ksSpectra[mult] );
    f1_hist[mult]->SetMarkerStyle(20);
    f1_hist[mult]->SetMarkerColor(kBlue);
    f1_hist[mult]->SetMarkerSize(1.3);
    f1_hist[mult]->SetLineColor(kWhite);
    f1_hist[mult]->SetTitle(" ");
    f1_hist[mult]->SetStats(kFALSE);
    f1_hist[mult]->GetYaxis()->SetTitle("Fit/data");
    f1_hist[mult]->GetXaxis()->SetTitle("p^{}_{T} (GeV/c)");
    f1_hist[mult]->GetYaxis()->SetTitleSize(0.09);
    f1_hist[mult]->GetYaxis()->SetTitleOffset(0.63);
    f1_hist[mult]->GetYaxis()->SetLabelSize(0.09);
    f1_hist[mult]->GetYaxis()->SetLabelOffset(0.025);
    f1_hist[mult]->GetYaxis()->SetRangeUser(0.6,1.4);

    f1_hist[mult]->GetXaxis()->SetTitleSize(0.09);
    f1_hist[mult]->GetXaxis()->SetLabelSize(0.10);

    f2_hist[mult]->Divide( laSpectra[mult] );
    f2_hist[mult]->SetMarkerStyle(20);
    f2_hist[mult]->SetMarkerColor(kBlack);
    f2_hist[mult]->SetMarkerSize(1.3);
    f2_hist[mult]->SetLineColor(kWhite);
    f2_hist[mult]->SetTitle(" ");
    f2_hist[mult]->SetStats(kFALSE);
    f2_hist[mult]->GetYaxis()->SetTitle("Fit/data");
    f2_hist[mult]->GetXaxis()->SetTitle("p^{}_{T} (GeV/c)");
    f2_hist[mult]->GetYaxis()->SetLabelSize(0.09);
    f2_hist[mult]->GetYaxis()->SetLabelOffset(0.025);
    f2_hist[mult]->GetYaxis()->SetRangeUser(-10,40);
    
    f1_hist[mult]->Draw("P");
    f2_hist[mult]->Draw("Psame");
    l1->Draw("same");

  }

  TCanvas* s1 = new TCanvas();

  TH1D* hist = new TH1D("h1","h1",100,0.2,0.7);
  hist->GetYaxis()->SetRangeUser(0.02,0.25);
  hist->GetXaxis()->SetRangeUser(0.2,0.7);

  TGraph* g = new TGraph(8);

  g->SetMarkerStyle(20);

  for(mult = 0; mult < 8; mult++){

    double temp = 0.0;
    temp = beta_T(beta_s[mult],n[mult]);

    g->SetPoint(mult,temp,Tkin[mult]);

  }

  hist->Draw("P");
  g->Draw("Psame");

  double ellipse_x[8][100];
  double ellipse_y[8][100];

  for(mult = 0; mult < 8; mult++){

    double temp = beta_T(beta_s[mult],n[mult]);
    double shift = beta_s[mult]-temp;

    for(int i = 0; i < 100; i++){

        test[mult]->GetPoint(i,ellipse_x[mult][i],ellipse_y[mult][i]);
        ellipse_x[mult][i] = ellipse_x[mult][i] - shift;
    }   
  }

  TGraph* g1[8];

  for(mult = 0; mult < 8; mult++){

    g1[mult] = new TGraph(100);
      
      for(i = 0; i < 100; i++){
          
          g1[mult]->SetPoint(i,ellipse_x[mult][i],ellipse_y[mult][i]);
          g1[mult]->SetLineColor(kRed);
          g1[mult]->SetMarkerColor(kRed);

      }

      g1[mult]->Draw("same");
  }


}
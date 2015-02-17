#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTree.h>
#include <iostream>
#include <TLine.h>
#include <math.h>
#include <TF1.h>
#include <TCut.h>
#include <TLatex.h>
#include "selectionCut.h"

void formatHist(TH1* h, int col = 1, double norm = 1);
void plotSubtracted(char* filename,char *title="",char *mycut= "",bool useMC = false)
{
   TFile *f= new TFile(filename);
   
   TTree * TrackletTree= dynamic_cast<TTree *>(f->Get("TrackletTree12"));

   TH2F *h1 = new TH2F("h1","",12,-3,3,50,0,5);
   TH2F *h2 = new TH2F("h2","",12,-3,3,50,0,5);
   TH2F *h2p = new TH2F("h2p","",12,-3,3,50,0,5);
   TH2F *h2MC = new TH2F("h2MC","",12,-3,3,50,0,5);
   TH2F *h3 = new TH2F("h3","",12,-3,3,50,0,5);
 
   selectionCut cut1(0);
   
   TCut cut = cut1.Cut && TCut(mycut);
   TCanvas *c = new TCanvas ("c","",400,400);
   TrackletTree->Draw("abs(deta):eta1>>h1","abs(dphi)<1.5"&&cut);

   h3->SetXTitle("#eta");
   h3->SetYTitle("#Delta#eta");

   TrackletTree->Draw("abs(deta):eta1>>h2","abs(dphi)>1.5&&abs(dphi)<3"&&cut);   

   h1->Sumw2();
   h2->Sumw2();
   h3->Sumw2();

   for (int x=0;x<80;x++) {
      for (int y=0;y<200;y++) {
         h3->SetBinContent(x,y, h1->GetBinContent(x,y)-h2->GetBinContent(x,y));
         h3->SetBinError(x,y, sqrt(h1->GetBinError(x,y)*h1->GetBinError(x,y)+h2->GetBinError(x,y)*h2->GetBinError(x,y)));
      }
   }
   
   h3->Draw("col");
   TCanvas *c2 = new TCanvas ("c2","",400,400);
   c2->SetLogy();
   TH1D *hx1 = h1->ProjectionX();
   TH1D *hx2 = h2->ProjectionX();
   TH1D *hx3 = h3->ProjectionX();
   TH1D *hy1 = h1->ProjectionY();
   TH1D *hy2 = h2->ProjectionY();
   TH1D *hy3 = h3->ProjectionY();

   int norm = h1->GetEntries();
   if (norm==0) norm=1;
   formatHist(hy1,1,norm);
   formatHist(hy2,4,norm);
   formatHist(hy3,2,norm);
   hy1->SetXTitle(Form("|#Delta#eta| %s",title));
   hy1->SetYTitle("Arbitrary Normalization");
   hy1->SetAxisRange(10e-7,2,"Y");
   hy1->Draw();
   hy2->SetMarkerColor(4);
   hy2->SetMarkerStyle(24);
   hy2->SetLineColor(4);
   hy2->Draw("same");
   hy3->SetMarkerColor(2);
   hy3->SetLineColor(2);
   hy3->SetMarkerStyle(25);
   hy3->Draw("same");
   TLegend * leg1 = new TLegend(0.28,0.72,0.50,0.85);
   leg1->SetFillStyle(0);  
   leg1->SetFillColor(0); 
   leg1->SetBorderSize(0);
   leg1->SetTextSize(0.03);
   
   leg1->AddEntry(hy1,"All tracklets","p");
   leg1->AddEntry(hy2,"Combinatorial Background","p");
   leg1->AddEntry(hy3,"Tracklets with Background Subtracted","p");

   leg1->Draw();
   TLatex *   tex = new TLatex(0.8,0.48,"CMS Preliminary");
   tex->SetLineWidth(2);
   tex->Draw();

   c2->SaveAs("Subtracted.gif");
   c2->SaveAs("Subtracted.C");
   c2->SaveAs("Subtracted.eps");

}

void formatHist(TH1* h, int col, double norm){

  h->Scale(1/norm);
  h->SetLineColor(col);
  h->SetMarkerColor(col);
  h->SetMarkerSize(1);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();

}

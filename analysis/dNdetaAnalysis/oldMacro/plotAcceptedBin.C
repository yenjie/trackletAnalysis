#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TArrow.h>
#include "selectionCut.h"
#include <TPaletteAxis.h>

void plotAcceptedBin(int ch=12)
{
//   TFile *inputFile = new TFile("correction/correction-12-Official.root","");
//   TFile *inputFile = new TFile(Form("plot/root/result-%d-Run124022-Official-BSC.root",ch),"");
//   TFile *inputFile = new TFile(Form("plot/root/result-%d-Run124120-yetkin-2360GeV-1224-HF2.root",ch),"");
   TFile *inputFile = new TFile(Form("plot/root/result-%d-Run132440-ATLAS-7TeV-HF1.root",ch),"");
//   TFile *inputFile = new TFile(Form("correction.root"),"");
//   TFile *inputFile = new TFile("./sample-1213/TrackletTree-900GeV-ATLAS.root","");
   TFile *inputFile2 = new TFile("./newMorePerfectSample/TrackletTree-Run132440-hfcuts-final.root");
//   TFile *inputFile2 = new TFile("./sample-official/TrackletTree-900GeV-D6T-20091210-KKVertex-removeDead2-newMult.root");
   
   TTree *TrackletTree = (TTree*)inputFile2->FindObjectAny(Form("TrackletTree%d",ch));
   TH2F *h = (TH2F*) inputFile->FindObjectAny("hAcceptance");
   TH2F *h2 = new TH2F("h2","",100,-3,3,100,-20,20);
   TCanvas *c = new TCanvas("c","",400,400);
   c->SetLeftMargin(0.1691919);
   c->SetRightMargin(0.1691919);
   c->SetTopMargin(0.0483871);
   c->SetBottomMargin(0.1505376);

   for (int i=2;i<=h->GetNbinsX()-1;i++){
      for (int j=0;j<=h->GetNbinsY();j++) {
         if (h->GetBinContent(i,j)!=0) h->SetBinContent(i,j,10e10); 
	 else                          { 
	    h->SetBinContent(i,j,0);
	    h->SetBinError(i,j,0);
         
	 }
      }
   }
   h->SetXTitle(Form("Pseudorapidity #eta (Tracklet %d)",ch));
   h->SetYTitle("Primary Vertex V_{z} (cm)");
   h->SetFillStyle(0);
   h->SetFillColor(15);
   h->SetLineColor(2);
   
   selectionCut myCut(0,0,10000);
   TrackletTree->Draw("vz[1]:eta1>>h2","abs(deta)<0.1&&abs(dphi)<1"&&myCut.Cut);
    
   h->Draw("box");
   c->SaveAs("plot/Acceptance/FigAcceptedBins.C");
   c->SaveAs("plot/Acceptance/FigAcceptedBins.eps");
   c->SaveAs("plot/Acceptance/FigAcceptedBins.gif");
   
   h2->SetXTitle(Form("Pseudorapidity #eta (Tracklet %d)",ch));
   h2->SetYTitle("Primary Vertex V_{z} (cm)");
/*
   TPaletteAxis *palette = new TPaletteAxis(3.068702,-19.86577,3.366412,20.13423,h2);
palette->SetLabelColor(1);
palette->SetLabelFont(42);
palette->SetLabelOffset(0.005);
palette->SetLabelSize(0.045);
palette->SetTitleOffset(1);
palette->SetTitleSize(0.04);
   palette->SetFillColor(100);
   palette->SetFillStyle(1001);
   h2->GetListOfFunctions()->Add(palette,"br");
*/
   h2->Draw("colz");
   h->Draw("box same");
   TLatex *   tex = new TLatex(-2.33,8.23,"CMS Preliminary");
   tex->SetLineWidth(2);
   tex->SetTextColor(0);
//   tex->Draw();
   c->SaveAs(Form("plot/Acceptance/FigAcceptedBinsOverlap%d.C",ch));
   c->SaveAs(Form("plot/Acceptance/FigAcceptedBinsOverlap%d.eps",ch));
   c->SaveAs(Form("plot/Acceptance/FigAcceptedBinsOverlap%d.gif",ch));

}


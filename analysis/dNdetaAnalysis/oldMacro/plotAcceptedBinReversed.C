#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TArrow.h>

void plotAcceptedBinReversed()
{
   TFile *inputFile = new TFile("correction/correction-12-10TeV.root","");
   TFile *inputFile2 = new TFile("sample/TrackletTree-10TeV-TV-2.4M.root");
//   TFile *inputFile2 = new TFile("sample/TrackletTree-10TeV-2-1L1.root");
   TTree *TrackletTree = (TTree*)inputFile2->FindObjectAny("TrackletTree12");
   TH2F *h1 = (TH2F*) inputFile->FindObjectAny("hAcceptance");
   TH2F *h =  new TH2F("h","",10,-10,10,12,-3,3);
   for (int i=0;i<=h1->GetNbinsX();i++){
      for (int j=0;j<=h1->GetNbinsY();j++) {
         if (h1->GetBinContent(i,j)!=0) h->SetBinContent(j,i,10e10);
      }
   }
   TH2F *h2 = new TH2F("h2","",200,-10,10,200,-3,3);
   TCanvas *c = new TCanvas("c","",400,400);
   for (int i=0;i<=h->GetNbinsX();i++){
      for (int j=0;j<=h->GetNbinsY();j++) {
         if (h->GetBinContent(i,j)!=0) h->SetBinContent(i,j,10e10);
      }
   }
   h->SetYTitle("Pseudorapidity (#eta)");
   h->SetXTitle("Primary Vertex V_{z} (cm)");
   h->SetFillStyle(0);
   h->SetFillColor(15);
   h->SetLineColor(2);
   
   TrackletTree->Draw("eta1:vz[3]>>h2","abs(deta)<0.1");
    
   h->Draw("box");
   c->SaveAs("plot/Acceptance/FigAcceptedBins.C");
   c->SaveAs("plot/Acceptance/FigAcceptedBins.eps");
   c->SaveAs("plot/Acceptance/FigAcceptedBins.gif");
   
   h2->SetYTitle("Pseudorapidity (#eta)");
   h2->SetXTitle("Primary Vertex V_{z} (cm)");
   h2->Draw("col");
   h->Draw("box same");
   TLatex *   tex = new TLatex(-2.33,8.23,"CMS Preliminary");
   tex->SetLineWidth(2);
   tex->SetTextColor(0);
   tex->Draw();
   c->SaveAs("plot/Acceptance/FigAcceptedBinsOverlap.C");
   c->SaveAs("plot/Acceptance/FigAcceptedBinsOverlap.eps");
   c->SaveAs("plot/Acceptance/FigAcceptedBinsOverlap.gif");

}


#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include <TText.h>

#include "selectionCut.h"

void plotDelta(char* infile,int type = 1)
{
   cout <<infile<<endl;
   TFile *f = new TFile(infile);
   TTree * myTree = (TTree *)(f->FindObjectAny("TrackletTree12"));
   TNtuple * myNtuple = (TNtuple *)(f->FindObjectAny("ntmatched"));

   selectionCut myCut(0);
   TCut dataCut = myCut.Cut;
   
   TCanvas *c = new TCanvas("c","#Delta#eta",400,400);
   c->SetLogy();
   TH1F *h = new TH1F("h","",200,-5,5);
   if (type) myTree->Draw("deta>>h",dataCut); else myNtuple->Draw("deta>>h");
   h->Scale(1./h->GetEntries());
   if (type) h->SetXTitle("#Delta#eta"); else h->SetXTitle("#Delta#eta (ToyMC)");
   h->SetTitleOffset(1.4,"Y");
   h->SetYTitle("Arbitrary Normalization");
   h->Draw();
   TText *cms = new TText(-4.659,0.1187,"CMS Preliminary");
   cms->Draw();

   c->SaveAs("plot/delta/FigDeltaEta.eps");
   c->SaveAs("plot/delta/FigDeltaEta.gif");
   c->SaveAs("plot/delta/FigDeltaEta.C");
   

   TCanvas *c2 = new TCanvas("c2","#Delta#phi",400,400);
   c2->SetLogy();
   TH1F *h2 = new TH1F("h2","",157,-3.14,3.14);
   if (type) myTree->Draw("dphi>>h2"); else myNtuple->Draw("dphi>>h2");
   h2->Scale(1./h2->GetEntries());
   if (type) h2->SetXTitle("#Delta#phi"); else h2->SetXTitle("#Delta#phi (ToyMC)");
   h2->SetYTitle("Arbitrary Normalization");
   h2->SetAxisRange(2e-4,0.3,"Y");
   h2->Draw();
   TText *cms2 = new TText(-2.85,0.15,"CMS Preliminary");
   cms2->Draw();
   c2->SaveAs("plot/delta/FigDeltaPhi.eps");
   c2->SaveAs("plot/delta/FigDeltaPhi.gif");
   c2->SaveAs("plot/delta/FigDeltaPhi.C");

   TCanvas *c3 = new TCanvas("c3","Scatter",400,400);
   c3->SetLogz();
   TH2F *h3 = new TH2F("h3","",50,-3.14,3.14,50,-5,5);
   if (type) myTree->Draw("deta:dphi>>h3",""); else myNtuple->Draw("deta:dphi>>h3");
  // h3->Scale(1./h3->GetEntries());
   if (type) h3->SetXTitle("#Delta#phi"); else h3->SetXTitle("#Delta#phi (ToyMC)");
   if (type) h3->SetYTitle("#Delta#eta"); else h3->SetYTitle("#Delta#eta (ToyMC)");
   h3->SetTitleOffset(1.4,"Y");
   h3->Draw("col");
   TText *cms3 = new TText(-2.85,3.777,"CMS Preliminary");
   cms3->SetTextColor(0);
   cms3->Draw();

   c3->SaveAs("plot/delta/FigPlotScatterPlot.eps");
   c3->SaveAs("plot/delta/FigPlotScatterPlot.gif");
   c3->SaveAs("plot/delta/FigPlotScatterPlot.C");

   TCanvas *c4 = new TCanvas("c4","#Delta#phi",400,400);
   c4->SetLogy();
   TH1F *h4 = new TH1F("h4","",157,-3.14,3.14);
   if (type) myTree->Draw("dphi>>h4","abs(deta)<1&&nhit1>50"); else myNtuple->Draw("dphi>>h4");
   h4->Scale(1./h4->GetEntries());
   if (type) h4->SetXTitle("#Delta#phi"); else h4->SetXTitle("#Delta#phi (ToyMC)");
   h4->SetYTitle("Arbitrary Normalization");
   h4->Draw();

   TCanvas *c5 = new TCanvas("c5","nhad>10",400,400);
   c5->SetLogy();
   TH1F *h5 = new TH1F("h5","",157,-3.14,3.14);
   if (type) myTree->Draw("dphi>>h5","abs(deta)>1&&nhit1>50"); else myNtuple->Draw("dphi>>h5");
   h5->Scale(1./h5->GetEntries());
   if (type) h5->SetXTitle("#Delta#phi"); else h5->SetXTitle("#Delta#phi (ToyMC)");
   h5->SetYTitle("Arbitrary Normalization");
   h5->Draw();



}

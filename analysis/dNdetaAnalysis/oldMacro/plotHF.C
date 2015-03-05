#include <iostream>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>

void plotHF(char *infname,char *infname2)
{  
   TFile *inf = new TFile(infname);
   TTree *t1 = (TTree*)inf->FindObjectAny("TrackletTree12");
   TFile *inf2 = new TFile(infname2);
   TTree *t2 = (TTree*)inf2->FindObjectAny("TrackletTree12");
   TH1F *p1 = new TH1F("p1","",100,0,100);
   TH1F *p2 = new TH1F("p2","",100,0,100);
   TH1F *p1noCut = new TH1F("p1noCut","",100,0,100);
   TH1F *p2noCut = new TH1F("p2noCut","",100,0,100);
   t1->Draw("Sum$(abs(eta)<2.5)>>p1","nHFn>0&&nHFp>0&&evtType!=92&&evtType!=93");
   t2->Draw("Sum$(abs(eta)<2.5)>>p2","nHFn>0&&nHFp>0&&evtType!=92&&evtType!=93");
   t1->Draw("Sum$(abs(eta)<2.5)>>p1noCut","evtType!=92&&evtType!=93");
   t2->Draw("Sum$(abs(eta)<2.5)>>p2noCut","evtType!=92&&evtType!=93");
   p1->Sumw2();
   p2->Sumw2();
   p1noCut->Sumw2();
   p2noCut->Sumw2();
   p1->Divide(p1noCut);
   p2->Divide(p2noCut);
   p1->SetXTitle("N_{ch}^{Gen} |#eta|<2.5");
   p1->SetYTitle("HF Selection Efficiency");
   

   TCanvas *c = new TCanvas("c","",600,600);
   p1->Draw();
   p2->SetLineColor(2);
   p2->SetMarkerColor(2);
   p2->Draw("same");
   TLegend *l = new TLegend(0.48,0.23,0.93,0.43);
   l->SetBorderSize(0);
   l->SetFillColor(0);
   l->SetLineColor(0);
   l->AddEntry(p1,"PYTHIA ATLAS 7TeV","pl");
   l->AddEntry(p2,"PYTHIA D6T 7TeV","pl");
   l->Draw();
   TCanvas *c2 = new TCanvas("c2","",600,600);
   TH1F *pRatio = (TH1F*)p1->Clone();
   pRatio->Divide(p2);
   pRatio->SetYTitle("Efficiency Ratio (ATLAS/D6T)");
   pRatio->Draw();   
}

#include <TFile.h>
#include <TH1F.h>
#include <TText.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>

void plotTriggerCorrection()
{
   TFile  *a = new TFile("correction/TriggerCorrection-10TeV.root");
   TFile  *b = new TFile("correction/TriggerCorrection-10TeV.root");
   TH1F *h1 = (TH1F*) a->FindObjectAny("h");
   h1->SetName("h1");
   TH1F *h2 = (TH1F*) b->FindObjectAny("h");
   h2->SetName("h2");
   TCanvas *c = new TCanvas ("c","",400,400);   
   h1->SetAxisRange(0.8,0.9,"Y");
   h1->SetXTitle("#eta");
   h1->SetYTitle("#xi");
 
   h1->Draw();
   h1->SetLineColor(2);
   h1->SetMarkerColor(2);
   h2->SetMarkerStyle(25);
   h2->Draw("same");

   TLegend * leg1 = new TLegend(0.43,0.2,0.92,0.34);
   leg1->SetFillStyle(0);
   leg1->SetFillColor(0);
   leg1->SetBorderSize(0);
   leg1->AddEntry(h2,"10TeV","pl");
   leg1->AddEntry(h1,"900GeV","pl");
  
   leg1->Draw();
   TLatex *   tex = new TLatex(-2.81,0.891,"CMS Preliminary");
   tex->SetLineWidth(2);
   tex->Draw();

   c->SaveAs("plot/Trigger/TriggerCorrection.gif");
   c->SaveAs("plot/Trigger/TriggerCorrection.C");
   c->SaveAs("plot/Trigger/TriggerCorrection.eps");
}

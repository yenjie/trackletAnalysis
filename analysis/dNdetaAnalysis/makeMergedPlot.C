#define dndetaRange 40
#define _NETABIN 30

#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TLegend.h>

#include "GraphErrorsBand.h"

int makeMergedPlot(const char* name = "EPOS-5TeV", const char* title = "") {
   TFile* inf12 = new TFile(Form("correction/correction-12-%s.root", name));
   TH1F* h12o = (TH1F*)inf12->FindObjectAny("hMeasuredFinal");
   h12o->SetName("h12o");
   h12o->SetAxisRange(0, dndetaRange, "Y");

   TFile* inf13 = new TFile(Form("correction/correction-13-%s.root", name));
   TH1F* h13o = (TH1F*)inf13->FindObjectAny("hMeasuredFinal");
   h13o->SetName("h13o");

   TFile* inf23 = new TFile(Form("correction/correction-23-%s.root", name));
   TH1F* h23o = (TH1F*)inf23->FindObjectAny("hMeasuredFinal");
   h23o->SetName("h23o");

   TFile* outfile = new TFile(Form("merged/merged-%s.root", name), "recreate");
   TCanvas* c1 = new TCanvas("c1", "", 600, 600);

   TH1F* h12 = (TH1F*)h12o->Clone();
   h12->SetName("h12");
   TH1F* h13 = (TH1F*)h13o->Clone();
   h13->SetName("h13");
   TH1F* h23 = (TH1F*)h23o->Clone();
   h23->SetName("h23");

   h13->SetMarkerStyle(26);
   h13->SetMarkerColor(1);
   h13->SetLineColor(1);
   h23->SetMarkerStyle(25);
   h23->SetMarkerColor(4);
   h23->SetLineColor(4);
   h12->SetMarkerSize(1);
   h13->SetMarkerSize(1);
   h23->SetMarkerSize(1);

   h12->SetXTitle("#eta");
   h12->SetYTitle("dN/d#eta");
   h12->Draw("same");
   h13->Draw("same");
   h23->Draw("same");

   TLegend* leg = new TLegend(0.2, 0.18, 0.9, 0.55);
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillStyle(0);
   leg->AddEntry("hTruth", title, "");
   leg->AddEntry("h12", "Reconstructed (1st+2nd layers)", "pl");
   leg->AddEntry("h13", "Reconstructed (1st+3rd layers)", "pl");
   leg->AddEntry("h23", "Reconstructed (2nd+3rd layers)", "pl");
   leg->Draw();
   c1->SaveAs(Form("merged/merged-%s.pdf", name));
   // c1->SaveAs(Form("merged/merged-%s.C", name));

   TCanvas* c2 = new TCanvas("c2", "", 600, 600);
   TH1F* hAvg = (TH1F*)h12->Clone();
   hAvg->SetName("hAvg");

   for (int i=1; i<=_NETABIN; i++) {
      double avg = 0;
      double avgErr = 0;
      avg += h12->GetBinContent(i);
      avgErr += h12->GetBinError(i)*h12->GetBinError(i);
      avg += h13->GetBinContent(i);
      avgErr += h13->GetBinError(i)*h13->GetBinError(i);
      avg += h23->GetBinContent(i);
      avgErr += h23->GetBinError(i)*h23->GetBinError(i);
      avgErr = sqrt(avgErr);
      if (i>6 && i<=_NETABIN-6) {
         avg /= 3.0;
         avgErr /= 3.0;
      }

      hAvg->SetBinContent(i, avg);
      hAvg->SetBinError(i, avgErr);
   }
   hAvg->Draw("p");

   TLegend* leg2 = new TLegend(0.2, 0.18, 0.9, 0.55);
   leg2->SetBorderSize(0);
   leg2->SetTextFont(62);
   leg2->SetLineColor(1);
   leg2->SetLineStyle(1);
   leg2->SetLineWidth(1);
   leg2->SetFillStyle(0);

   leg2->AddEntry("hTruth", title, "");
   leg2->AddEntry(hAvg, "Reconstructed Tracklets", "pl");
   leg2->Draw();
   c2->SaveAs(Form("merged/avg-%s.pdf", name));
   // c2->SaveAs(Form("merged/avg-%s.C", name));

   outfile->Write();

   return 0;
}

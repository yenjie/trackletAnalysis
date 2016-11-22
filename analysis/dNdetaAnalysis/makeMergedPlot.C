#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TLegend.h>

#include "GraphErrorsBand.h"

int makeMergedPlot(const char* name = "EPOS-8TeV", const char* title = "") {
   TFile* inf12 = new TFile(Form("correction/correction-12-%s.root", name));
   TH1F* h12 = (TH1F*)((TH1F*)inf12->FindObjectAny("hMeasuredFinal"))->Clone("h12");

   TFile* inf13 = new TFile(Form("correction/correction-13-%s.root", name));
   TH1F* h13 = (TH1F*)((TH1F*)inf13->FindObjectAny("hMeasuredFinal"))->Clone("h13");

   TFile* inf23 = new TFile(Form("correction/correction-23-%s.root", name));
   TH1F* h23 = (TH1F*)((TH1F*)inf23->FindObjectAny("hMeasuredFinal"))->Clone("h23");

   TFile* infepos = new TFile("correction/correction-12-EPOS-8TeV.root");
   TH1F* hepos = (TH1F*)((TH1F*)infepos->FindObjectAny("hTruthWOSelection"))->Clone("hepos");
   hepos->SetAxisRange(0, 27, "Y");

   TFile* outfile = new TFile(Form("merged/merged-%s.root", name), "recreate");
   TCanvas* c1 = new TCanvas("c1", "", 600, 600);

   h12->SetXTitle("#eta");
   h12->SetYTitle("dN/d#eta");
   h12->SetMarkerSize(1);

   h13->SetMarkerSize(1);
   h13->SetMarkerStyle(26);
   h13->SetMarkerColor(1);
   h13->SetLineColor(1);

   h23->SetMarkerSize(1);
   h23->SetMarkerStyle(25);
   h23->SetMarkerColor(4);
   h23->SetLineColor(4);

   hepos->SetMarkerColor(kGreen+2);
   hepos->SetLineColor(kGreen+2);

   hepos->Draw("hist c ][");

   h12->Draw("same");
   h13->Draw("same");
   h23->Draw("same");

   TLegend* l1 = new TLegend(0.2, 0.2, 0.8, 0.4);
   l1->SetBorderSize(0);
   l1->SetTextFont(62);
   l1->SetLineColor(1);
   l1->SetLineStyle(1);
   l1->SetLineWidth(1);
   l1->SetFillStyle(0);
   l1->AddEntry("hepos", "EPOS LHC 8 TeV", "l");
   l1->AddEntry("h12", "Reconstructed (1st+2nd layers)", "pl");
   l1->AddEntry("h13", "Reconstructed (1st+3rd layers)", "pl");
   l1->AddEntry("h23", "Reconstructed (2nd+3rd layers)", "pl");
   l1->Draw();

   c1->Draw();
   c1->SaveAs(Form("merged/merged-%s.pdf", name));

   TCanvas* c2 = new TCanvas("c2", "", 600, 600);
   TH1F* hAvg = (TH1F*)h12->Clone();
   hAvg->SetName("hAvg");

   for (int i=1; i<=30; i++) {
      double avg = 0;
      double avgErr = 0;
      avg += h12->GetBinContent(i);
      avgErr += h12->GetBinError(i)*h12->GetBinError(i);
      avg += h13->GetBinContent(i);
      avgErr += h13->GetBinError(i)*h13->GetBinError(i);
      avg += h23->GetBinContent(i);
      avgErr += h23->GetBinError(i)*h23->GetBinError(i);
      avgErr = sqrt(avgErr);
      if (i>6 && i<=24) {
         avg /= 3.0;
         avgErr /= 3.0;
      }

      hAvg->SetBinContent(i, avg);
      hAvg->SetBinError(i, avgErr);
   }

   hepos->Draw("hist ][");
   hAvg->Draw("p same");

   TLegend* l2 = new TLegend(0.2, 0.2, 0.8, 0.4);
   l2->SetBorderSize(0);
   l2->SetTextFont(62);
   l2->SetLineColor(1);
   l2->SetLineStyle(1);
   l2->SetLineWidth(1);
   l2->SetFillStyle(0);

   l2->AddEntry("hepos", "EPOS LHC 8TeV", "l");
   l2->AddEntry("hAvg", "Reconstructed Tracklets", "pl");
   l2->Draw();

   c2->Draw();
   c2->SaveAs(Form("merged/avg-%s.pdf", name));

   outfile->Write();

   return 0;
}

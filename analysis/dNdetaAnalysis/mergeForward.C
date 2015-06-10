#define dndetaRange 9

#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>

#include "GraphErrorsBand.h"

int mergeForward(const char* name = "PYTHIA_Monash13", double syserr = 0) {
   TFile* f12 = new TFile(Form("correction/correction-12-%s.root", name));
   TH1F* h12 = (TH1F*)f12->FindObjectAny("hMeasuredFinal");
   TH1F* hMC = (TH1F*)f12->FindObjectAny("hTruthWOSelection");
   h12->SetName("h12");
   h12->SetAxisRange(0, dndetaRange, "Y");

   TFile* f13 = new TFile(Form("correction/correction-13-%s.root", name));
   TH1F* h13 = (TH1F*)f13->FindObjectAny("hMeasuredFinal");
   h13->SetName("h13");

   TFile* f23 = new TFile(Form("correction/correction-23-%s.root", name));
   TH1F* h23 = (TH1F*)f23->FindObjectAny("hMeasuredFinal");
   h23->SetName("h23");

   TFile* f14 = new TFile(Form("correction/correction-14-%s.root", name));
   TH1F* h14 = (TH1F*)f14->FindObjectAny("hMeasuredFinal");
   h14->SetName("h14");

   TFile* f15 = new TFile(Form("correction/correction-15-%s.root", name));
   TH1F* h15 = (TH1F*)f15->FindObjectAny("hMeasuredFinal");
   h15->SetName("h15");

   TFile* f45 = new TFile(Form("correction/correction-45-%s.root", name));
   TH1F* h45 = (TH1F*)f45->FindObjectAny("hMeasuredFinal");
   h45->SetName("h45");

   TFile* outfile = new TFile(Form("merged/merged-%s.root", name), "recreate");
   TCanvas* c1 = new TCanvas("c1", "", 600, 600);

   h12->SetXTitle("#eta");
   h12->SetYTitle("dN/d#eta");
   h12->SetMarkerSize(1);
   h13->SetMarkerStyle(26);
   h13->SetMarkerSize(1);
   h13->SetMarkerColor(1);
   h13->SetLineColor(1);
   h23->SetMarkerStyle(25);
   h23->SetMarkerSize(1);
   h23->SetMarkerColor(4);
   h23->SetLineColor(4);
   h14->SetMarkerStyle(28);
   h14->SetMarkerSize(1);
   h14->SetMarkerColor(8);
   h14->SetLineColor(8);
   h15->SetMarkerStyle(30);
   h15->SetMarkerSize(1);
   h15->SetMarkerColor(7);
   h15->SetLineColor(7);
   h45->SetMarkerStyle(32);
   h45->SetMarkerSize(1);
   h45->SetMarkerColor(9);
   h45->SetLineColor(9);

   hMC->SetLineColor(1);
   hMC->SetMarkerColor(1);
   hMC->Draw("same");
   h12->Draw("same");
   h13->Draw("same");
   h23->Draw("same");
   h14->Draw("same");
   h15->Draw("same");
   h45->Draw("same");

   TLegend* leg = new TLegend(0.2, 0.18, 1, 0.48);
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->AddEntry("hTruth", name, "");
   leg->AddEntry("h12", "Reconstructed (1st+2nd layers)", "pl");
   leg->AddEntry("h13", "Reconstructed (1st+3rd layers)", "pl");
   leg->AddEntry("h23", "Reconstructed (2nd+3rd layers)", "pl");
   leg->AddEntry("h14", "Reconstructed (1st+4th layers)", "pl");
   leg->AddEntry("h15", "Reconstructed (1st+5th layers)", "pl");
   leg->AddEntry("h45", "Reconstructed (4th+5th layers)", "pl");
   leg->Draw();
   c1->SaveAs(Form("merged/merged-%s.pdf", name));

   TCanvas* c2 = new TCanvas("c2", "", 600, 600);
   TH1F* hAvg = (TH1F*)h12->Clone();
   hAvg->SetName("hAvg");

   double avg, avgErr;
   for (int i=4; i<6; i++) {
      avg = 0;
      avgErr = 0;
      avg += h12->GetBinContent(i);
      avgErr += h12->GetBinError(i)*h12->GetBinError(i);
      avg += h14->GetBinContent(i);
      avgErr += h14->GetBinError(i)*h14->GetBinError(i);
      avg += h15->GetBinContent(i);
      avgErr += h15->GetBinError(i)*h15->GetBinError(i);
      avg += h45->GetBinContent(i);
      avgErr += h45->GetBinError(i)*h45->GetBinError(i);
      avgErr = sqrt(avgErr);

      hAvg->SetBinContent(i, avg/4);
      hAvg->SetBinError(i, avgErr/4);

      avg = 0;
      avgErr = 0;
      avg += h12->GetBinContent(31-i);
      avgErr += h12->GetBinError(31-i)*h12->GetBinError(31-i);
      avg += h14->GetBinContent(31-i);
      avgErr += h14->GetBinError(31-i)*h14->GetBinError(31-i);
      avg += h15->GetBinContent(31-i);
      avgErr += h15->GetBinError(31-i)*h15->GetBinError(31-i);
      avg += h45->GetBinContent(31-i);
      avgErr += h45->GetBinError(31-i)*h45->GetBinError(31-i);
      avgErr = sqrt(avgErr);

      hAvg->SetBinContent(31-i, avg/4);
      hAvg->SetBinError(31-i, avgErr/4);
   }
   for (int i=6; i<7; i++) {
      avg = 0;
      avgErr = 0;
      avg += h12->GetBinContent(i);
      avgErr += h12->GetBinError(i)*h12->GetBinError(i);
      avg += h13->GetBinContent(i);
      avgErr += h13->GetBinError(i)*h13->GetBinError(i);
      avg += h23->GetBinContent(i);
      avgErr += h23->GetBinError(i)*h23->GetBinError(i);
      avg += h14->GetBinContent(i);
      avgErr += h14->GetBinError(i)*h14->GetBinError(i);
      avg += h15->GetBinContent(i);
      avgErr += h15->GetBinError(i)*h15->GetBinError(i);
      avg += h45->GetBinContent(i);
      avgErr += h45->GetBinError(i)*h45->GetBinError(i);
      avgErr = sqrt(avgErr);

      hAvg->SetBinContent(i, avg/6);
      hAvg->SetBinError(i, avgErr/6);

      avg = 0;
      avgErr = 0;
      avg += h12->GetBinContent(31-i);
      avgErr += h12->GetBinError(31-i)*h12->GetBinError(31-i);
      avg += h13->GetBinContent(31-i);
      avgErr += h13->GetBinError(31-i)*h13->GetBinError(31-i);
      avg += h23->GetBinContent(31-i);
      avgErr += h23->GetBinError(31-i)*h23->GetBinError(31-i);
      avg += h14->GetBinContent(31-i);
      avgErr += h14->GetBinError(31-i)*h14->GetBinError(31-i);
      avg += h15->GetBinContent(31-i);
      avgErr += h15->GetBinError(31-i)*h15->GetBinError(31-i);
      avg += h45->GetBinContent(31-i);
      avgErr += h45->GetBinError(31-i)*h45->GetBinError(31-i);
      avgErr = sqrt(avgErr);

      hAvg->SetBinContent(31-i, avg/6);
      hAvg->SetBinError(31-i, avgErr/6);
   }
   for (int i=7; i<8; i++) {
      avg = 0;
      avgErr = 0;
      avg += h12->GetBinContent(i);
      avgErr += h12->GetBinError(i)*h12->GetBinError(i);
      avg += h13->GetBinContent(i);
      avgErr += h13->GetBinError(i)*h13->GetBinError(i);
      avg += h23->GetBinContent(i);
      avgErr += h23->GetBinError(i)*h23->GetBinError(i);
      avg += h14->GetBinContent(i);
      avgErr += h14->GetBinError(i)*h14->GetBinError(i);
      avgErr = sqrt(avgErr);

      hAvg->SetBinContent(i, avg/4);
      hAvg->SetBinError(i, avgErr/4);

      avg = 0;
      avgErr = 0;
      avg += h12->GetBinContent(31-i);
      avgErr += h12->GetBinError(31-i)*h12->GetBinError(31-i);
      avg += h13->GetBinContent(31-i);
      avgErr += h13->GetBinError(31-i)*h13->GetBinError(31-i);
      avg += h23->GetBinContent(31-i);
      avgErr += h23->GetBinError(31-i)*h23->GetBinError(31-i);
      avg += h14->GetBinContent(31-i);
      avgErr += h14->GetBinError(31-i)*h14->GetBinError(31-i);
      avgErr = sqrt(avgErr);

      hAvg->SetBinContent(31-i, avg/4);
      hAvg->SetBinError(31-i, avgErr/4);
   }
   for (int i=8; i<16; i++) {
      avg = 0;
      avgErr = 0;
      avg += h12->GetBinContent(i);
      avgErr += h12->GetBinError(i)*h12->GetBinError(i);
      avg += h13->GetBinContent(i);
      avgErr += h13->GetBinError(i)*h13->GetBinError(i);
      avg += h23->GetBinContent(i);
      avgErr += h23->GetBinError(i)*h23->GetBinError(i);
      avgErr = sqrt(avgErr);

      hAvg->SetBinContent(i, avg/3);
      hAvg->SetBinError(i, avgErr/3);

      avg = 0;
      avgErr = 0;
      avg += h12->GetBinContent(31-i);
      avgErr += h12->GetBinError(31-i)*h12->GetBinError(31-i);
      avg += h13->GetBinContent(31-i);
      avgErr += h13->GetBinError(31-i)*h13->GetBinError(31-i);
      avg += h23->GetBinContent(31-i);
      avgErr += h23->GetBinError(31-i)*h23->GetBinError(31-i);
      avgErr = sqrt(avgErr);

      hAvg->SetBinContent(31-i, avg/3);
      hAvg->SetBinError(31-i, avgErr/3);
   }

   hAvg->Draw("p");
   TGraph* gErrorBand = GetErrorBand(hAvg, syserr, syserr, 0.1);
   gErrorBand->Draw("F");
   hMC->Draw("same");
   hAvg->Draw("p same");

   TLegend* leg2 = new TLegend(0.2, 0.18, 0.9, 0.35);
   leg2->SetBorderSize(0);
   leg2->SetTextFont(62);
   leg2->SetLineColor(1);
   leg2->SetLineStyle(1);
   leg2->SetLineWidth(1);
   leg2->SetFillColor(0);
   leg2->SetFillStyle(0);

   leg2->AddEntry("hTruth", name, "");
   leg2->AddEntry(hAvg, "Reconstructed Tracklets", "pl");
   leg2->Draw();
   c2->SaveAs(Form("merged/avg-%s.pdf", name));

   TCanvas* c3 = new TCanvas("c3", "", 600, 600);
   TH1F* hSym = (TH1F*)h12->Clone();
   hSym->SetName("hSym");

   for (int i=4; i<6; i++) {
      avg = 0;
      avgErr = 0;
      avg += h12->GetBinContent(i);
      avgErr += h12->GetBinError(i)*h12->GetBinError(i);
      avg += h14->GetBinContent(i);
      avgErr += h14->GetBinError(i)*h14->GetBinError(i);
      avg += h15->GetBinContent(i);
      avgErr += h15->GetBinError(i)*h15->GetBinError(i);
      avg += h45->GetBinContent(i);
      avgErr += h45->GetBinError(i)*h45->GetBinError(i);
      avg += h12->GetBinContent(31-i);
      avgErr += h12->GetBinError(31-i)*h12->GetBinError(31-i);
      avg += h14->GetBinContent(31-i);
      avgErr += h14->GetBinError(31-i)*h14->GetBinError(31-i);
      avg += h15->GetBinContent(31-i);
      avgErr += h15->GetBinError(31-i)*h15->GetBinError(31-i);
      avg += h45->GetBinContent(31-i);
      avgErr += h45->GetBinError(31-i)*h45->GetBinError(31-i);
      avgErr = sqrt(avgErr);

      hAvg->SetBinContent(i, avg/8);
      hAvg->SetBinError(i, avgErr/8);
      hAvg->SetBinContent(31-i, avg/8);
      hAvg->SetBinError(31-i, avgErr/8);
   }
   for (int i=6; i<7; i++) {
      avg = 0;
      avgErr = 0;
      avg += h12->GetBinContent(i);
      avgErr += h12->GetBinError(i)*h12->GetBinError(i);
      avg += h13->GetBinContent(i);
      avgErr += h13->GetBinError(i)*h13->GetBinError(i);
      avg += h23->GetBinContent(i);
      avgErr += h23->GetBinError(i)*h23->GetBinError(i);
      avg += h14->GetBinContent(i);
      avgErr += h14->GetBinError(i)*h14->GetBinError(i);
      avg += h15->GetBinContent(i);
      avgErr += h15->GetBinError(i)*h15->GetBinError(i);
      avg += h45->GetBinContent(i);
      avgErr += h45->GetBinError(i)*h45->GetBinError(i);
      avg += h12->GetBinContent(31-i);
      avgErr += h12->GetBinError(31-i)*h12->GetBinError(31-i);
      avg += h13->GetBinContent(31-i);
      avgErr += h13->GetBinError(31-i)*h13->GetBinError(31-i);
      avg += h23->GetBinContent(31-i);
      avgErr += h23->GetBinError(31-i)*h23->GetBinError(31-i);
      avg += h14->GetBinContent(31-i);
      avgErr += h14->GetBinError(31-i)*h14->GetBinError(31-i);
      avg += h15->GetBinContent(31-i);
      avgErr += h15->GetBinError(31-i)*h15->GetBinError(31-i);
      avg += h45->GetBinContent(31-i);
      avgErr += h45->GetBinError(31-i)*h45->GetBinError(31-i);
      avgErr = sqrt(avgErr);

      hAvg->SetBinContent(i, avg/12);
      hAvg->SetBinError(i, avgErr/12);
      hAvg->SetBinContent(31-i, avg/12);
      hAvg->SetBinError(31-i, avgErr/12);
   }
   for (int i=7; i<8; i++) {
      avg = 0;
      avgErr = 0;
      avg += h12->GetBinContent(i);
      avgErr += h12->GetBinError(i)*h12->GetBinError(i);
      avg += h13->GetBinContent(i);
      avgErr += h13->GetBinError(i)*h13->GetBinError(i);
      avg += h23->GetBinContent(i);
      avgErr += h23->GetBinError(i)*h23->GetBinError(i);
      avg += h14->GetBinContent(i);
      avgErr += h14->GetBinError(i)*h14->GetBinError(i);
      avg += h12->GetBinContent(31-i);
      avgErr += h12->GetBinError(31-i)*h12->GetBinError(31-i);
      avg += h13->GetBinContent(31-i);
      avgErr += h13->GetBinError(31-i)*h13->GetBinError(31-i);
      avg += h23->GetBinContent(31-i);
      avgErr += h23->GetBinError(31-i)*h23->GetBinError(31-i);
      avg += h14->GetBinContent(31-i);
      avgErr += h14->GetBinError(31-i)*h14->GetBinError(31-i);
      avgErr = sqrt(avgErr);

      hAvg->SetBinContent(i, avg/8);
      hAvg->SetBinError(i, avgErr/8);
      hAvg->SetBinContent(31-i, avg/8);
      hAvg->SetBinError(31-i, avgErr/8);
   }
   for (int i=8; i<16; i++) {
      avg = 0;
      avgErr = 0;
      avg += h12->GetBinContent(i);
      avgErr += h12->GetBinError(i)*h12->GetBinError(i);
      avg += h13->GetBinContent(i);
      avgErr += h13->GetBinError(i)*h13->GetBinError(i);
      avg += h23->GetBinContent(i);
      avgErr += h23->GetBinError(i)*h23->GetBinError(i);
      avg += h12->GetBinContent(31-i);
      avgErr += h12->GetBinError(31-i)*h12->GetBinError(31-i);
      avg += h13->GetBinContent(31-i);
      avgErr += h13->GetBinError(31-i)*h13->GetBinError(31-i);
      avg += h23->GetBinContent(31-i);
      avgErr += h23->GetBinError(31-i)*h23->GetBinError(31-i);
      avgErr = sqrt(avgErr);

      hAvg->SetBinContent(i, avg/6);
      hAvg->SetBinError(i, avgErr/6);
      hAvg->SetBinContent(31-i, avg/6);
      hAvg->SetBinError(31-i, avgErr/6);
   }

   hSym->Draw("p");
   TGraph* gErrorBand2 = GetErrorBand(hSym, syserr, syserr, 0.1);
   gErrorBand2->Draw("F");
   hMC->Draw("same");
   hSym->SetLineColor(4);
   hSym->SetMarkerColor(4);
   hSym->Draw("p same");

   TLegend* leg3 = new TLegend(0.2, 0.18, 0.9, 0.35);
   leg3->SetBorderSize(0);
   leg3->SetTextFont(62);
   leg3->SetLineColor(1);
   leg3->SetLineStyle(1);
   leg3->SetLineWidth(1);
   leg3->SetFillColor(0);
   leg3->SetFillStyle(0);

   leg3->AddEntry("hTruth", name, "");
   leg3->AddEntry(hSym, "Reconstructed Tracklets", "pl");
   leg3->Draw();
   c3->SaveAs(Form("merged/avgsym-%s.pdf", name));

   outfile->Write();

   return 0;
}

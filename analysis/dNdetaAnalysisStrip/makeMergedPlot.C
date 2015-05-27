#define dndetaRange 9
#define _NETABIN 12

#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TNtuple.h>
#include "GraphErrorsBand.h"

void clearNBins(int n, TH1F* h) {
   for (int i=1; i<=n; i++) {
      h->SetBinContent(i, 0);
      h->SetBinContent(_NETABIN+1-i, 0);
      h->SetBinError(i, 0);
      h->SetBinError(_NETABIN+1-i, 0);
   }
}

int makeMergedPlot(const char* name = "PYTHIA_Monash13") {

   TFile *infEPOS = new TFile("correction/correction-12-EPOS.root");
   TH1F* hMCEPOS = (TH1F*)infEPOS->FindObjectAny("hTruthWOSelection");
   TFile *infPYTHIA = new TFile("correction/correction-12-PU28.root");
   TH1F* hMCPYTHIA = (TH1F*)infPYTHIA->FindObjectAny("hTruthWOSelection");
   
   TFile *inf12 = new TFile(Form("correction/correction-12-%s.root", name));
   TH1F *h12 = (TH1F*)inf12->FindObjectAny("hMeasuredFinal");
   TH1F* hMC = (TH1F*)inf12->FindObjectAny("hTruthWOSelection");
   h12->SetName("h12");
   h12->SetAxisRange(0, dndetaRange, "Y");

   TFile *inf13 = new TFile(Form("correction/correction-13-%s.root", name));
   TH1F *h13 = (TH1F*)inf13->FindObjectAny("hMeasuredFinal");
   h13->SetName("h13");

   TFile *inf14 = new TFile(Form("correction/correction-14-%s.root", name));
   TH1F *h14 = (TH1F*)inf14->FindObjectAny("hMeasuredFinal");
   h14->SetName("h14");

   TFile *inf23 = new TFile(Form("correction/correction-23-%s.root", name));
   TH1F *h23 = (TH1F*)inf23->FindObjectAny("hMeasuredFinal");
   h23->SetName("h23");

   TFile *inf24 = new TFile(Form("correction/correction-24-%s.root", name));
   TH1F *h24 = (TH1F*)inf24->FindObjectAny("hMeasuredFinal");
   h24->SetName("h24");

   TFile *inf34 = new TFile(Form("correction/correction-34-%s.root", name));
   TH1F *h34 = (TH1F*)inf34->FindObjectAny("hMeasuredFinal");
   h34->SetName("h34");

   TFile *outfile = new TFile(Form("merged/merged-%s.root", name), "recreate");
   TCanvas *c = new TCanvas("c", "", 600, 600);

   h13->SetMarkerStyle(26);
   h13->SetMarkerColor(1);
   h13->SetLineColor(1);
   h14->SetMarkerStyle(27);
   h14->SetMarkerColor(2);
   h14->SetLineColor(2);
   h23->SetMarkerStyle(25);
   h23->SetMarkerColor(4);
   h23->SetLineColor(4);
   h24->SetMarkerStyle(28);
   h24->SetMarkerColor(1);
   h24->SetLineColor(1);
   h34->SetMarkerStyle(30);
   h34->SetMarkerColor(4);
   h34->SetLineColor(4);
   h12->SetMarkerSize(1);
   h13->SetMarkerSize(1);
   h14->SetMarkerSize(1);
   h23->SetMarkerSize(1);
   h24->SetMarkerSize(1);
   h34->SetMarkerSize(1);

   clearNBins(1, h12);
   clearNBins(2, h13);
   clearNBins(2, h14);
   clearNBins(2, h23);
   clearNBins(2, h24);
   clearNBins(2, h34);

   h12->SetXTitle("#eta");
   h12->SetYTitle("dN/d#eta");
   hMC->SetLineColor(1);
   hMC->SetMarkerColor(1);
   hMC->Draw("hist same");
   h12->Draw("same");
   h13->Draw("same");
   h14->Draw("same");
   h23->Draw("same");
   h24->Draw("same");
   h34->Draw("same");

   TLegend *leg = new TLegend(0.2, 0.18, 1, 0.45);
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
   leg->AddEntry("h14", "Reconstructed (1st+4th layers)", "pl");
   leg->AddEntry("h23", "Reconstructed (2nd+3rd layers)", "pl");
   leg->AddEntry("h24", "Reconstructed (2nd+4th layers)", "pl");
   leg->AddEntry("h34", "Reconstructed (3rd+4th layers)", "pl");
   leg->Draw();
   c->SaveAs(Form("merged/merged-%s.pdf", name));

   TCanvas *c2 = new TCanvas("c2", "", 600, 600);
   TH1F *hAvg = (TH1F*)h12->Clone();
   hAvg->SetName("hAvg");

   for (int i=1; i<=_NETABIN; i++) {
      double avg = 0;
      double avgErr = 0;
      avg += h12->GetBinContent(i);
      avgErr += h12->GetBinError(i)*h12->GetBinError(i);
      avg += h13->GetBinContent(i);
      avgErr += h13->GetBinError(i)*h13->GetBinError(i);
      avg += h14->GetBinContent(i);
      avgErr += h14->GetBinError(i)*h14->GetBinError(i);
      avg += h23->GetBinContent(i);
      avgErr += h23->GetBinError(i)*h23->GetBinError(i);
      avg += h24->GetBinContent(i);
      avgErr += h24->GetBinError(i)*h24->GetBinError(i);
      avg += h34->GetBinContent(i);
      avgErr += h34->GetBinError(i)*h34->GetBinError(i);
      avgErr = sqrt(avgErr);
      if (i>2&&i<=_NETABIN-2) {
         avg /= 6.0;
         avgErr /= 6.0;
      }

      hAvg->SetBinContent(i, avg);
      hAvg->SetBinError(i, avgErr);
   }
   hMC->Draw("same hist");
   hAvg->Draw("p same");

   TLegend *leg3 = new TLegend(0.2, 0.18, 0.9, 0.45);
   leg3->SetBorderSize(0);
   leg3->SetTextFont(62);
   leg3->SetLineColor(1);
   leg3->SetLineStyle(1);
   leg3->SetLineWidth(1);
   leg3->SetFillColor(0);
   leg3->SetFillStyle(0);

   leg3->AddEntry("hTruth", name, "");
   leg3->AddEntry(hAvg, "Reconstructed Tracklets", "pl");
   leg3->Draw();
   c2->SaveAs(Form("merged/avg-%s.pdf", name));

   TCanvas *c3 = new TCanvas("c3", "", 600, 600);
   TH1F* hAvg2 = (TH1F*)h12->Clone();
   hAvg2->SetName("hAvg2");

   for (int i=1; i<=_NETABIN/2; i++) {
      double avg = 0;
      double avgErr = 0;
      avg += h12->GetBinContent(i);
      avgErr += h12->GetBinError(i)*h12->GetBinError(i);
      avg += h13->GetBinContent(i);
      avgErr += h13->GetBinError(i)*h13->GetBinError(i);
      avg += h14->GetBinContent(i);
      avgErr += h14->GetBinError(i)*h14->GetBinError(i);
      avg += h23->GetBinContent(i);
      avgErr += h23->GetBinError(i)*h23->GetBinError(i);
      avg += h24->GetBinContent(i);
      avgErr += h24->GetBinError(i)*h24->GetBinError(i);
      avg += h34->GetBinContent(i);
      avgErr += h34->GetBinError(i)*h34->GetBinError(i);
      avg += h12->GetBinContent(_NETABIN+1-i);
      avgErr += h12->GetBinError(_NETABIN+1-i)*h12->GetBinError(_NETABIN+1-i);
      avg += h13->GetBinContent(_NETABIN+1-i);
      avgErr += h13->GetBinError(_NETABIN+1-i)*h13->GetBinError(_NETABIN+1-i);
      avg += h14->GetBinContent(_NETABIN+1-i);
      avgErr += h14->GetBinError(_NETABIN+1-i)*h14->GetBinError(_NETABIN+1-i);
      avg += h23->GetBinContent(_NETABIN+1-i);
      avgErr += h23->GetBinError(_NETABIN+1-i)*h23->GetBinError(_NETABIN+1-i);
      avg += h24->GetBinContent(_NETABIN+1-i);
      avgErr += h24->GetBinError(_NETABIN+1-i)*h24->GetBinError(_NETABIN+1-i);
      avg += h34->GetBinContent(_NETABIN+1-i);
      avgErr += h34->GetBinError(_NETABIN+1-i)*h34->GetBinError(_NETABIN+1-i);
      avgErr = sqrt(avgErr);
      if (i>2) {
         avg /= 12.0;
         avgErr /= 12.0;
      } else {
         avg /= 2.0;
         avgErr /= 2.0;
      }

      hAvg2->SetBinContent(i, avg);
      hAvg2->SetBinError(i, avgErr);
      hAvg2->SetBinContent(_NETABIN+1-i, avg);
      hAvg2->SetBinError(_NETABIN+1-i, avgErr);
      cout <<fabs(hAvg2->GetBinCenter(i))<<" "<<avg<<endl;
   }
   TH1D *hh = (TH1D*)hAvg2->Clone("hh");
   for (int i=0;i<=hAvg2->GetNbinsX();i++) {
      double val = hAvg2->GetBinContent(i);
      hh->SetBinContent(i,val);
      hh->SetBinError(i,val*0.06);
     // TLine *b = new TLine(hAvg2->GetBinLowEdge(i),val*1.1,hAvg2->GetBinLowEdge(i+1),val*0.9);
      
      //if (val!=0) b->Draw("same");
   
   }
   hh->SetLineWidth(10);
   hh->SetLineColor(kGray);
   hh->Draw("same");
   hMC->Draw("hist same");
   
   hAvg2->SetLineColor(4);
   hAvg2->SetMarkerColor(4);
   hAvg2->Draw("p same");

   TLegend *leg2 = new TLegend(0.2, 0.18, 0.9, 0.45);
   leg2->SetBorderSize(0);
   leg2->SetTextFont(62);
   leg2->SetLineColor(1);
   leg2->SetLineStyle(1);
   leg2->SetLineWidth(1);
   leg2->SetFillColor(0);
   leg2->SetFillStyle(0);

   leg2->AddEntry("hTruth", name, "");
   leg2->AddEntry(hAvg2, "Reconstructed Tracklets", "pl");
   leg2->Draw();
   c3->SaveAs(Form("merged/avgsym-%s.pdf", name));

   hMCEPOS->SetLineColor(6);
   hMCPYTHIA->SetLineColor(4);
   //hMCEPOS->Draw("hist same");
   //hMCPYTHIA->Draw("hist same");
   outfile->Write();

   return 0;
}

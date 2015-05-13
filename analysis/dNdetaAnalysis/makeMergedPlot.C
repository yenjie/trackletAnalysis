#define dndetaRange 8
#define _NETABIN 30

#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TNtuple.h>
#include "UA5Plot.h"
// #include "ResultPlots.h"
#include "GraphErrorsBand.h"

void clearNBins(int n, TH1F* h) {
   for (int i=1; i<=n; i++) {
      h->SetBinContent(i, 0);
      h->SetBinContent(_NETABIN+1-i, 0);
      h->SetBinError(i, 0);
      h->SetBinError(_NETABIN+1-i, 0);
   }
}

int makeMergedPlot(const char* name = "PYTHIA8_Monash13", int zerobins = 5, int UA5 = 0) {
   TFile *inf12 = new TFile(Form("correction/correction-12-%s.root", name));
   TH1F *h12 = (TH1F*)inf12->FindObjectAny("hMeasuredFinal");
   TH1F* hMC = (TH1F*)inf12->FindObjectAny("hTruthWOSelection");
   h12->SetName("h12");
   h12->SetAxisRange(0, dndetaRange, "Y");

   TFile *inf13 = new TFile(Form("correction/correction-13-%s.root", name));
   TH1F *h13 = (TH1F*)inf13->FindObjectAny("hMeasuredFinal");
   h13->SetName("h13");

   TFile *inf23 = new TFile(Form("correction/correction-23-%s.root", name));
   TH1F *h23 = (TH1F*)inf23->FindObjectAny("hMeasuredFinal");
   h23->SetName("h23");

   TFile *outfile = new TFile(Form("merged/merged-%s.root", name), "recreate");
   TCanvas *c = new TCanvas("c", "", 600, 600);

   h13->SetMarkerStyle(26);
   h13->SetMarkerColor(1);
   h13->SetLineColor(1);
   h23->SetMarkerStyle(25);
   h23->SetMarkerColor(4);
   h23->SetLineColor(4);
   h12->SetMarkerSize(1.2);
   h13->SetMarkerSize(1.2);
   h23->SetMarkerSize(1.2);

   clearNBins(3, h12);
   clearNBins(zerobins, h13);
   clearNBins(zerobins, h23);

   TH1F *hUA5 = getUA5NSD();
   TH1F *hUA5Scaled = getUA5NSDScaled();
   // TH1F *hTracklet900GeV = tracklet900GeV();
   // TH1F *hTracklet900GeV = tracklet900GeVHF1();
   // TH1F *hTracklet2360GeVHF1 = tracklet2360GeVHF1();
   // TH1F* hPixelCounting = getPixelCounting2360GeV();
   // TH1F *hPixelCounting7TeV = pixelCounting7TeV();
   // TH1F* hBSC2360GeV = getBSC2360GeV();
   // TH1F* hBSC900GeV = getBSC900GeV();
   // TH1F* hHF900GeV = getHF900GeV();

   h12->SetXTitle("#eta");
   h12->SetYTitle("dN/d#eta");
   if (UA5) hUA5->Draw("p");
   if (UA5>=4) hUA5Scaled->Draw("p same");
   hMC->Draw("same");
   h12->Draw("same");
   h13->Draw("same");
   h23->Draw("same");

   TLegend *leg = new TLegend(0.2, 0.18, 1, 0.35);
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
   leg->Draw();
   c->SaveAs(Form("merged/merged-%s.pdf", name));

   TCanvas *c2 = new TCanvas("c2", "", 600, 600);
   TH1F *hAvg = (TH1F*)h12->Clone();
   hAvg->SetName("hAvg");

   hMC->Draw();
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
      if (i>zerobins&&i<=_NETABIN-zerobins) {
         avg /= 3.0;
         avgErr /= 3.0;
      }

      hAvg->SetBinContent(i, avg);
      hAvg->SetBinError(i, avgErr);
   }
   hAvg->Draw("same");

   if (UA5) hUA5->Draw("p same");
   // if (UA5>=2) hTracklet900GeV->Draw("p same");
   // if (UA5>=3) hTracklet2360GeVHF1->Draw("p same");
   if (UA5>=4) hUA5Scaled->Draw("p same");
   // if (UA5>=3) hPixelCounting->Draw("same");
   // pixelCounting7TeV()->Draw("same");

   TLegend *leg3 = new TLegend(0.2, 0.18, 0.9, 0.35);
   leg3->SetBorderSize(0);
   leg3->SetTextFont(62);
   leg3->SetLineColor(1);
   leg3->SetLineStyle(1);
   leg3->SetLineWidth(1);
   leg3->SetFillColor(0);
   leg3->SetFillStyle(0);

   leg3->AddEntry("hTruth", name, "");
   leg3->AddEntry(hAvg, "Reconstructed Tracklets", "pl");
   // if (UA5>=3) leg3->AddEntry(hPixelCounting, "2.36 TeV p+#bar{p} by Pixel counting(CMS)", "pl");
   // if (UA5>=2) leg3->AddEntry(hTracklet900GeV, "0.9 TeV p+p by Tracklet (CMS)", "pl");
   if (UA5) leg3->AddEntry(hUA5, "0.9 TeV p+#bar{p} (UA5)", "pl");
   // leg3->AddEntry(hPixelCounting7TeV, "7.0 TeV p+p by PixelCounting (CMS)", "pl");
   leg3->Draw();
   c2->SaveAs(Form("merged/avg-%s.pdf", name));

   // TCanvas *c3 = new TCanvas("c3", "", 600, 600);
   // TH1F* hAvg2 = (TH1F*) h12->Clone();
   // hAvg2->SetName("hAvg2");

   // for (int i=1; i<=_NETABIN/2; i++) {
   //    double avg = 0;
   //    avg += h12->GetBinContent(i);
   //    avg += h13->GetBinContent(i);
   //    avg += h23->GetBinContent(i);
   //    avg += h12->GetBinContent(_NETABIN+1-i);
   //    avg += h13->GetBinContent(_NETABIN+1-i);
   //    avg += h23->GetBinContent(_NETABIN+1-i);

   //    if (i>zerobins) avg/=6.0; else avg/=2.0;
   //    double avgErr = avg*3.8/100.;

   //    hAvg2->SetBinContent(i, avg);
   //    hAvg2->SetBinError(i, 0, avgErr);
   //    hAvg2->SetBinContent(_NETABIN+1-i, avg);
   //    hAvg2->SetBinError(_NETABIN+1-i, avgErr);

   //    cout << i << " " << avg << " +- " << avgErr << endl;
   // }
   // hAvg2->Draw();

   // if (UA5) hUA5->Draw("p same");
   // // if (UA5>=2) hTracklet900GeV->Draw("p same");
   // // if (UA5>=3) hTracklet2360GeVHF1->Draw("p same");
   // if (UA5>=4) hUA5Scaled->Draw("p same");

   // hAvg2->SetLineColor(1);
   // hAvg2->SetMarkerColor(1);
   // hAvg2->Draw("p same");

   // TLegend *leg2 = new TLegend(0.2, 0.18, 0.9, 0.35);
   // leg2->SetBorderSize(0);
   // leg2->SetTextFont(62);
   // leg2->SetLineColor(1);
   // leg2->SetLineStyle(1);
   // leg2->SetLineWidth(1);
   // leg2->SetFillColor(0);
   // leg2->SetFillStyle(0);

   // // hTracklet900GeV->Draw("p same");
   // // hBSC2360GeV->Draw("same");
   // // hBSC900GeV->Draw("same");
   // // hHF900GeV->SetMarkerStyle(20);
   // // hHF900GeV->SetMarkerColor(2);
   // // hHF900GeV->SetLineColor(2);
   // // hHF900GeV->Draw("same");

   // // leg2->AddEntry("hTruth", name, "");
   // leg2->AddEntry(hAvg2, "7.0 TeV p+p by Tracklet (CMS)", "pl");
   // // leg2->AddEntry(hTracklet2360GeVHF1, "2.36 TeV p+p by Tracklet (CMS)", "pl");
   // // leg2->AddEntry(hAvg2, "0.9 TeV p+p by Tracklet HF1(CMS)", "pl");
   // // leg2->AddEntry(hBSC900GeV, "0.9 TeV p+p by Tracklet BSC(CMS)", "pl");
   // // leg2->AddEntry(hHF900GeV, "0.9 TeV p+p by Tracklet (CMS)", "pl");
   // // leg2->AddEntry(hAvg2, "Run 124023 p+p by Tracklet (CMS)", "pl");
   // // leg2->AddEntry(hTracklet900GeV, "Run 123596 p+p by Tracklet (CMS)", "pl");
   // if (UA5) leg2->AddEntry(hUA5, "0.9 TeV p+#bar{p} (UA5)", "pl");
   // leg2->Draw();

   outfile->Write();

   return 0;
}

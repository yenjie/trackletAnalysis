#define dndetaRange 8
#define _NETABIN 30

#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TLegend.h>

#include "GraphErrorsBand.h"

void clearNBins(int n, TH1F* h) {
   for (int i=1; i<=n; i++) {
      h->SetBinContent(i, 0);
      h->SetBinContent(_NETABIN+1-i, 0);
      h->SetBinError(i, 0);
      h->SetBinError(_NETABIN+1-i, 0);
   }
}

int makeMergedPlot(const char* name = "PYTHIA_Monash13", const char* title = "") {
   TFile *infPYTHIA = new TFile("gen/pythia_CUETP8M1.root");
   TH1F* hPYTHIA = (TH1F*)infPYTHIA->FindObjectAny("hEta");
   hPYTHIA->SetName("hPYTHIA");
   TFile *infPYTHIAM = new TFile("gen/pythia_Monash.root");
   TH1F* hPYTHIAM = (TH1F*)infPYTHIAM->FindObjectAny("hEta");
   hPYTHIAM->SetName("hPYTHIAM");
   TFile *infEPOS = new TFile("gen/EPOS.root");
   TH1F* hEPOS = (TH1F*)infEPOS->FindObjectAny("hEta");
   hEPOS->SetName("hEPOS");
   TFile *infQGSJet = new TFile("gen/QGSJet.root");
   TH1F* hQGSJet = (TH1F*)infQGSJet->FindObjectAny("hEta");
   hQGSJet->SetName("hQGSJet");
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

   clearNBins(4, h12);
   clearNBins(6, h13);
   clearNBins(6, h23);

   h12->SetXTitle("#eta");
   h12->SetYTitle("dN/d#eta");
   hPYTHIA->SetLineColor(4);
   hPYTHIA->SetMarkerColor(4);
   hPYTHIAM->SetLineColor(kGreen+2);
   hPYTHIAM->SetMarkerColor(kGreen+2);
   hQGSJet->SetLineColor(2);
   hQGSJet->SetMarkerColor(2);
   hEPOS->SetLineColor(6);
   hEPOS->SetMarkerColor(6);
   //hPYTHIA->Draw("hist c same");
   //hPYTHIAM->Draw("hist c same");
   //hEPOS->Draw("hist c same");
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
   //leg->AddEntry("hPYTHIA", "PYTHIA8 CUETP8M1", "l");
   //leg->AddEntry("hPYTHIAM", "PYTHIA8 Monash", "l");
   //leg->AddEntry("hEPOS", "EPOS LHC", "l");
   //leg->AddEntry("hQGSJet", "QGSJet-II", "l");
   leg->AddEntry("h12", "Reconstructed (1st+2nd layers)", "pl");
   leg->AddEntry("h13", "Reconstructed (1st+3rd layers)", "pl");
   leg->AddEntry("h23", "Reconstructed (2nd+3rd layers)", "pl");
   leg->Draw();
   c1->SaveAs(Form("merged/merged-%s.pdf", name));
   c1->SaveAs(Form("merged/merged-%s.C", name));

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
      if (i>6&&i<=_NETABIN-6) {
         avg /= 3.0;
         avgErr /= 3.0;
      }

      hAvg->SetBinContent(i, avg);
      hAvg->SetBinError(i, avgErr);
   }
   hAvg->Draw("p");
   Double_t erreta[30] = {     0,      0,      0,  0.048,  0.048,
                           0.048,  0.048,  0.048,  0.048,  0.048,
                           0.048,  0.048,  0.048,  0.048,  0.048,
                           0.048,  0.048,  0.048,  0.048,  0.048,
                           0.048,  0.048,  0.048,  0.048,  0.048,
                           0.048,  0.048,      0,      0,      0};
   TGraph* gErrorBand = GetErrorBand(hAvg, erreta, erreta, 0.1);
   gErrorBand->Draw("F");
   hPYTHIA->Draw("hist c same");
   hPYTHIAM->Draw("hist c same");
   hQGSJet->Draw("hist c same");
   hEPOS->Draw("hist c same");
   hAvg->Draw("p same");

   TLegend* leg2 = new TLegend(0.2, 0.18, 0.9, 0.55);
   leg2->SetBorderSize(0);
   leg2->SetTextFont(62);
   leg2->SetLineColor(1);
   leg2->SetLineStyle(1);
   leg2->SetLineWidth(1);
   leg2->SetFillStyle(0);

   leg2->AddEntry("hTruth", title, "");
   leg2->AddEntry("hPYTHIA", "PYTHIA8 CUETP8M1", "l");
   leg2->AddEntry("hPYTHIAM", "PYTHIA8 Monash", "l");
   leg2->AddEntry("hEPOS", "EPOS LHC", "l");
   leg2->AddEntry("hQGSJet", "QGSJet-II", "l");
   leg2->AddEntry(hAvg, "Reconstructed Tracklets", "pl");
   leg2->Draw();
   c2->SaveAs(Form("merged/avg-%s.pdf", name));
   c2->SaveAs(Form("merged/avg-%s.C", name));

   TCanvas* c3 = new TCanvas("c3", "", 600, 600);
   TH1F* hSym = (TH1F*)h12->Clone();
   hSym->SetName("hSym");

   for (int i=1; i<=_NETABIN/2; i++) {
      double avg = 0;
      double avgErr = 0;
      avg += h12->GetBinContent(i);
      avgErr += h12->GetBinError(i)*h12->GetBinError(i);
      avg += h13->GetBinContent(i);
      avgErr += h13->GetBinError(i)*h13->GetBinError(i);
      avg += h23->GetBinContent(i);
      avgErr += h23->GetBinError(i)*h23->GetBinError(i);
      avg += h12->GetBinContent(_NETABIN+1-i);
      avgErr += h12->GetBinError(_NETABIN+1-i)*h12->GetBinError(_NETABIN+1-i);
      avg += h13->GetBinContent(_NETABIN+1-i);
      avgErr += h13->GetBinError(_NETABIN+1-i)*h13->GetBinError(_NETABIN+1-i);
      avg += h23->GetBinContent(_NETABIN+1-i);
      avgErr += h23->GetBinError(_NETABIN+1-i)*h23->GetBinError(_NETABIN+1-i);
      avgErr = sqrt(avgErr);
      if (i>6) {
         avg /= 6.0;
         avgErr /= 6.0;
      } else {
         avg /= 2.0;
         avgErr /= 2.0;
      }

      hSym->SetBinContent(i, avg);
      hSym->SetBinError(i, avgErr);
      hSym->SetBinContent(_NETABIN+1-i, avg);
      hSym->SetBinError(_NETABIN+1-i, avgErr);
   }
   hSym->Draw("p");
   TGraph* gErrorBand2 = GetErrorBand(hSym, erreta, erreta, 0.1);
   gErrorBand2->Draw("F");
   hPYTHIA->Draw("hist c same");
   hPYTHIAM->Draw("hist c same");
   hQGSJet->Draw("hist c same");
   hEPOS->Draw("hist c same");
   hSym->SetLineColor(4);
   hSym->SetMarkerColor(4);
   hSym->Draw("p same");

   TLegend* leg3 = new TLegend(0.2, 0.18, 0.9, 0.55);
   leg3->SetBorderSize(0);
   leg3->SetTextFont(62);
   leg3->SetLineColor(1);
   leg3->SetLineStyle(1);
   leg3->SetLineWidth(1);
   leg3->SetFillStyle(0);

   leg3->AddEntry("hTruth", title, "");
   leg3->AddEntry("hPYTHIA", "PYTHIA8 CUETP8M1", "l");
   leg3->AddEntry("hPYTHIAM", "PYTHIA8 Monash", "l");
   leg3->AddEntry("hEPOS", "EPOS LHC", "l");
   leg3->AddEntry("hQGSJet", "QGSJet-II", "l");
   leg3->AddEntry(hSym, "Reconstructed Tracklets", "pl");
   leg3->Draw();
   c3->SaveAs(Form("merged/avgsym-%s.pdf", name));
   c3->SaveAs(Form("merged/avgsym-%s.C", name));

   outfile->Write();

   TCanvas *cRatio = new TCanvas("cRatio","",600,600);
   TH1D *h12Ratio = (TH1D*)h12->Clone("h12Ratio");
   TH1D *h13Ratio = (TH1D*)h13->Clone("h13Ratio");
   TH1D *h23Ratio = (TH1D*)h23->Clone("h23Ratio");
   h12Ratio->Divide(hSym);
   h13Ratio->Divide(hSym);
   h23Ratio->Divide(hSym);

   TH2D *hRatioTmp = new TH2D("hRatioTmp",";#eta;dN/d#eta / <dN/d#eta>",100,-3,3,100,0.8,1.1);
   hRatioTmp->Draw();
   h12Ratio->Draw("same");
   h13Ratio->Draw("same");
   h23Ratio->Draw("same");

   leg->Draw();

   TFile *infAcc12 = new TFile("correction/acceptance-12.root");
   TFile *infAcc13 = new TFile("correction/acceptance-13.root");
   TFile *infAcc23 = new TFile("correction/acceptance-23.root");

   TH1D *hRatio12 = (TH1D*)infAcc12->Get("hRatio");
   hRatio12->SetName("hRatio12");
   TH1D *hRatio13 = (TH1D*)infAcc13->Get("hRatio");
   hRatio13->SetName("hRatio13");
   TH1D *hRatio23 = (TH1D*)infAcc23->Get("hRatio");
   hRatio23->SetName("hRatio23");

   TH1D *hDoubleRatio13=(TH1D*)hRatio12->Clone("hDoubleRatio13");
   TH1D *hDoubleRatio23=(TH1D*)hRatio12->Clone("hDoubleRatio23");

   hDoubleRatio13->Divide(hRatio13);
   hDoubleRatio23->Divide(hRatio23);
//   hDoubleRatio13->Draw("same");
//   hDoubleRatio23->Draw("same");
   
   return 0;
}

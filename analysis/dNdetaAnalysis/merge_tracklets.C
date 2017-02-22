#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TLine.h>
#include <TLegend.h>

#include "GraphErrorsBand.h"

int merge_tracklets(const char* name, const char* title, const char* refname) {
   TFile* inf12 = new TFile(Form("correction/correction-12-%s.root", name));
   TH1F* h12 = (TH1F*)((TH1F*)inf12->FindObjectAny("hMeasuredFinal"))->Clone("h12");

   TFile* inf13 = new TFile(Form("correction/correction-13-%s.root", name));
   TH1F* h13 = (TH1F*)((TH1F*)inf13->FindObjectAny("hMeasuredFinal"))->Clone("h13");

   TFile* inf23 = new TFile(Form("correction/correction-23-%s.root", name));
   TH1F* h23 = (TH1F*)((TH1F*)inf23->FindObjectAny("hMeasuredFinal"))->Clone("h23");

   TFile* infref = new TFile(Form("correction/correction-12-%s.root", refname));
   TH1F* href = (TH1F*)((TH1F*)infref->FindObjectAny("hTruthWOSelection"))->Clone("href");
   href->SetAxisRange(0, 30, "Y");

   TFile* outfile = new TFile(Form("rootfiles/merged-%s.root", name), "recreate");
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

   href->SetMarkerSize(0);
   href->SetMarkerColor(kGreen+2);
   href->SetLineColor(kGreen+2);

   href->Draw("c hist");

   h12->Draw("same");
   h13->Draw("same");
   h23->Draw("same");

   TLegend* l1 = new TLegend(0.2, 0.24, 0.8, 0.4);
   l1->SetBorderSize(0);
   l1->SetTextFont(43);
   l1->SetTextSize(18);
   l1->SetLineColor(1);
   l1->SetLineStyle(1);
   l1->SetLineWidth(1);
   l1->SetFillStyle(0);
   l1->AddEntry("href", title, "p");
   l1->AddEntry("h12", "Reconstructed (1st+2nd layers)", "pl");
   l1->AddEntry("h13", "Reconstructed (1st+3rd layers)", "pl");
   l1->AddEntry("h23", "Reconstructed (2nd+3rd layers)", "pl");
   l1->Draw();

   c1->Draw();
   c1->SaveAs(Form("figs/merged/merged-%s.png", name));

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
      if (i>5 && i<26) {
         avg /= 3.0;
         avgErr /= 3.0;
      }

      hAvg->SetBinContent(i, avg);
      hAvg->SetBinError(i, avgErr);
   }

   href->Draw("c hist");
   hAvg->Draw("p same");

   TLegend* l2 = new TLegend(0.2, 0.36, 0.8, 0.45);
   l2->SetBorderSize(0);
   l2->SetTextFont(43);
   l2->SetTextSize(18);
   l2->SetLineColor(1);
   l2->SetLineStyle(1);
   l2->SetLineWidth(1);
   l2->SetFillStyle(0);

   l2->AddEntry("href", title, "p");
   l2->AddEntry("hAvg", "Reconstructed Tracklets", "pl");
   l2->Draw();

   c2->Draw();
   c2->SaveAs(Form("figs/merged/avg-%s.png", name));

   TCanvas* c3 = new TCanvas("c3", "", 600, 600);

   TH1D* hratio12 = (TH1D*)h12->Clone("hratio12");
   TH1D* hratio13 = (TH1D*)h13->Clone("hratio13");
   TH1D* hratio23 = (TH1D*)h23->Clone("hratio23");

   hratio12->Divide(hAvg);
   hratio13->Divide(hAvg);
   hratio23->Divide(hAvg);

   hratio12->SetAxisRange(0.8, 1.2, "Y");
   hratio12->SetYTitle("Ratio");
   hratio12->Draw();
   hratio13->Draw("same");
   hratio23->Draw("same");

   TLine* lup = new TLine(-3, 1.03, 3, 1.03);
   lup->SetLineStyle(2);
   lup->Draw();
   TLine* ldown = new TLine(-3, 0.97, 3, 0.97);
   ldown->SetLineStyle(2);
   ldown->Draw();

   TLegend* l3 = new TLegend(0.2, 0.24, 0.8, 0.4);
   l3->SetBorderSize(0);
   l3->SetTextFont(43);
   l3->SetTextSize(18);
   l3->SetLineColor(1);
   l3->SetLineStyle(1);
   l3->SetLineWidth(1);
   l3->SetFillStyle(0);
   l3->AddEntry(href, title, "p");
   l3->AddEntry("hratio12", "Reconstructed (1st+2nd layers)", "pl");
   l3->AddEntry("hratio13", "Reconstructed (1st+3rd layers)", "pl");
   l3->AddEntry("hratio23", "Reconstructed (2nd+3rd layers)", "pl");
   l3->Draw();

   c3->SaveAs(Form("figs/merged/ratio-%s.png", name));

   h12->Write();
   h13->Write();
   h23->Write();
   href->Write();

   outfile->Write("", TObject::kOverwrite);
   outfile->Close();

   return 0;
}

int main(int argc, char* argv[]) {
   if (argc == 4)
      return merge_tracklets(argv[1], argv[2], argv[3]);
   else
      return 1;
}

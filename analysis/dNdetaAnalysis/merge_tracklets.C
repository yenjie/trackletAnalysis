#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLine.h>
#include <TLegend.h>

int merge_tracklets(const char* label, const char* ref) {
   TFile* f12 = new TFile(Form("correction/correction-12-%s.root", label));
   TH1F* h12 = (TH1F*)((TH1F*)f12->FindObjectAny("hMeasuredFinal"))->Clone("h12");

   TFile* f13 = new TFile(Form("correction/correction-13-%s.root", label));
   TH1F* h13 = (TH1F*)((TH1F*)f13->FindObjectAny("hMeasuredFinal"))->Clone("h13");

   TFile* f23 = new TFile(Form("correction/correction-23-%s.root", label));
   TH1F* h23 = (TH1F*)((TH1F*)f23->FindObjectAny("hMeasuredFinal"))->Clone("h23");

   TFile* fref = new TFile(Form("correction/correction-12-%s.root", ref));
   TH1F* href = (TH1F*)((TH1F*)fref->FindObjectAny("hTruthWOSelection"))->Clone("href");
   href->SetAxisRange(0, 30, "Y");

   TFile* fout = new TFile(Form("rootfiles/merged-%s.root", label), "recreate");
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
   l1->AddEntry("href", label, "p");
   l1->AddEntry("h12", "Reconstructed (1st+2nd layers)", "pl");
   l1->AddEntry("h13", "Reconstructed (1st+3rd layers)", "pl");
   l1->AddEntry("h23", "Reconstructed (2nd+3rd layers)", "pl");
   l1->Draw();

   c1->Draw();
   c1->SaveAs(Form("figs/merged/merged-%s.png", label));

   TCanvas* c2 = new TCanvas("c2", "", 600, 600);
   TH1F* havg = (TH1F*)h12->Clone("havg");

   for (int i=1; i<=30; i++) {
      double avg = 0;
      double avg_err = 0;
      avg += h12->GetBinContent(i);
      avg_err += h12->GetBinError(i) * h12->GetBinError(i);
      avg += h13->GetBinContent(i);
      avg_err += h13->GetBinError(i) * h13->GetBinError(i);
      avg += h23->GetBinContent(i);
      avg_err += h23->GetBinError(i) * h23->GetBinError(i);
      avg_err = sqrt(avg_err);
      if (i>5 && i<26) {
         avg /= 3.0;
         avg_err /= 3.0;
      }

      havg->SetBinContent(i, avg);
      havg->SetBinError(i, avg_err);
   }

   href->Draw("c hist");
   havg->Draw("p same");

   TLegend* l2 = new TLegend(0.2, 0.36, 0.8, 0.45);
   l2->SetBorderSize(0);
   l2->SetTextFont(43);
   l2->SetTextSize(18);
   l2->SetLineColor(1);
   l2->SetLineStyle(1);
   l2->SetLineWidth(1);
   l2->SetFillStyle(0);

   l2->AddEntry("href", label, "p");
   l2->AddEntry("havg", "Reconstructed Tracklets", "pl");
   l2->Draw();

   c2->Draw();
   c2->SaveAs(Form("figs/merged/merged-%s-avg.png", label));

   TCanvas* c3 = new TCanvas("c3", "", 600, 600);

   TH1D* hratio12 = (TH1D*)h12->Clone("hratio12");
   TH1D* hratio13 = (TH1D*)h13->Clone("hratio13");
   TH1D* hratio23 = (TH1D*)h23->Clone("hratio23");

   hratio12->Divide(havg);
   hratio13->Divide(havg);
   hratio23->Divide(havg);

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
   l3->AddEntry("href", label, "p");
   l3->AddEntry("hratio12", "Reconstructed (1st+2nd layers)", "pl");
   l3->AddEntry("hratio13", "Reconstructed (1st+3rd layers)", "pl");
   l3->AddEntry("hratio23", "Reconstructed (2nd+3rd layers)", "pl");
   l3->Draw();

   c3->SaveAs(Form("figs/merged/merged-%s-ratio.png", label));

   h12->Write();
   h13->Write();
   h23->Write();
   href->Write();

   fout->Write("", TObject::kOverwrite);
   fout->Close();

   return 0;
}

int main(int argc, char* argv[]) {
   if (argc == 3)
      return merge_tracklets(argv[1], argv[2]);
   else
      return 1;
}

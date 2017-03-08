#include "TFile.h"
#include "TTree.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"

int vertex_eff_and_res(const char* fname, int energy) {
   TFile* f = new TFile(fname, "READ");
   TTree* t = (TTree*)f->Get("TrackletTree12");

   TProfile *heff = new TProfile("heff", "", 40, 0, 40);
   t->Draw("vz[1]>-99:nhit1>>heff", "nhit1>0 && abs(vz[0])<15 && passHLT", "goff");

   heff->SetMarkerStyle(8);
   heff->SetLineColor(4);
   heff->SetMarkerColor(4);

   TH2D* h2d_vz_nhit1 = new TH2D("h2d_vz_nhit1", "", 40, 0, 40, 100, -2, 2);
   t->Draw("vz[1]-vz[0]:nhit1>>h2d_vz_nhit1", "nhit1>0 && abs(vz[0])<15 && passHLT", "goff colz");
   h2d_vz_nhit1->FitSlicesY();
   TH1D* hres = (TH1D*)gDirectory->Get("h2d_vz_nhit1_2");

   hres->SetMarkerStyle(8);
   hres->SetLineColor(4);
   hres->SetMarkerColor(4);

   TCanvas* c1 = new TCanvas("c1", "", 600, 600);
   heff->SetAxisRange(0.0, 1.2, "Y");
   heff->SetXTitle("nhit1");
   heff->SetYTitle("Efficiency");
   heff->SetStats(0);
   heff->Draw();

   TLegend* l1 = new TLegend(0.5, 0.6, 0.9, 0.72);
   l1->SetBorderSize(0);
   l1->SetFillStyle(0);
   l1->AddEntry(heff, Form("EPOS LHC %i TeV", energy), "pl");
   l1->Draw();

   c1->SaveAs(Form("figs/vertex-eff-epos-%itev.png", energy));

   TCanvas* c2 = new TCanvas("c2", "", 600, 600);
   c2->SetLogy();
   hres->SetAxisRange(0.002, 0.2, "Y");
   hres->SetTitle("");
   hres->SetXTitle("nhit1");
   hres->SetYTitle("Resolution (cm)");
   hres->SetStats(0);
   hres->Draw("same");

   TLegend* l2 = new TLegend(0.5, 0.6, 0.9, 0.72);
   l2->SetBorderSize(0);
   l2->SetFillStyle(0);
   l2->AddEntry(hres, Form("EPOS LHC %i TeV", energy), "pl");
   l2->Draw();

   c2->SaveAs(Form("figs/vertex-res-epos-%itev.png", energy));

   return 0;
}

int main(int argc, char* argv[]) {
   if (argc == 3)
      return vertex_eff_and_res(argv[1], atoi(argv[2]));
   else
      return 1;
}

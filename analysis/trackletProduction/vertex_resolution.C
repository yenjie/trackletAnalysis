#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>

int vertex_resolution(const char* fname, int energy) {
   TFile* f = new TFile(fname, "READ");
   TTree* t = (TTree*)f->Get("TrackletTree12");
   TH2D* h2d_vz_nhit1 = new TH2D("h2d_vz_nhit1", "", 40, 0, 40, 100, -2, 2);
   t->Draw("vz[1]-vz[0]:nhit1>>h2d_vz_nhit1", "weight * (nhit1>0 && abs(vz[1])<20 && abs(vz[0])<20 && passHLT)", "goff colz");
   h2d_vz_nhit1->FitSlicesY();
   TH1D* hres = (TH1D*)gDirectory->Get("h2d_vz_nhit1_2");

   hres->SetMarkerStyle(8);
   hres->SetLineColor(4);
   hres->SetMarkerColor(4);

   TCanvas* c1 = new TCanvas("c1", "", 600, 600);
   c1->SetLogy();
   hres->SetAxisRange(0.002, 0.2, "Y");
   hres->SetTitle("");
   hres->SetXTitle("nhit1");
   hres->SetYTitle("Resolution (cm)");
   hres->SetStats(0);
   hres->Draw("same");

   TLegend* l = new TLegend(0.5, 0.6, 0.9, 0.72);
   l->SetBorderSize(0);
   l->SetFillStyle(0);
   l->AddEntry(hres, Form("EPOS LHC %i TeV", energy), "pl");
   l->Draw();

   c1->SaveAs(Form("figs/vertex-res-epos-%itev.png", energy));

   return 0;
}

int main(int argc, char* argv[]) {
   if (argc == 3)
      return vertex_resolution(argv[1], atoi(argv[2]));
   else
      return 1;
}

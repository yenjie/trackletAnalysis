#include <TProfile.h>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>

int vertex_efficiency(const char* fname, int energy) {
   TFile* f = new TFile(fname);
   TTree* t = (TTree*)f->Get("TrackletTree12");
   TProfile *heff = new TProfile("heff", "", 40, 0, 40);
   t->Draw("vz[1]>-99:nhit1>>heff", "weight * (abs(vz[0])<20 && passHLT)", "goff");

   heff->SetMarkerStyle(8);
   heff->SetLineColor(4);
   heff->SetMarkerColor(4);

   TCanvas* c1 = new TCanvas("c1", "", 600, 600);
   heff->SetAxisRange(0.0, 1.2, "Y");
   heff->SetXTitle("nhit1");
   heff->SetYTitle("Efficiency");
   heff->SetStats(0);
   heff->Draw();

   TLegend* l = new TLegend(0.5, 0.6, 0.9, 0.72);
   l->SetBorderSize(0);
   l->SetFillStyle(0);
   l->AddEntry(heff, Form("EPOS LHC %i TeV", energy), "pl");
   l->Draw();

   c1->SaveAs(Form("figs/vertex-eff-epos-%itev.png", energy));

   return 0;
}

int main(int argc, char* argv[]) {
   if (argc == 3)
      return vertex_efficiency(argv[1], atoi(argv[2]));
   else
      return 1;
}

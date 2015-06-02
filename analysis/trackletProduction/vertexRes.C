#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>

void vertexRes(const char* fname) {
   TFile* f = new TFile(fname, "READ");
   TTree* t = (TTree*)f->Get("TrackletTree12");
   TH2D* hTracklet = new TH2D("hTracklet", "", 25, 0, 50, 100, -2, 2);
   t->Draw("vz[1]-vz[0]:nhit1>>hTracklet", "vz[1]>-40 && nhit1>0 && abs(vz[1])<15 && abs(vz[0])<15", "colz");
   hTracklet->FitSlicesY();
   TH1D* hRes = (TH1D*)gDirectory->Get("hTracklet_2");

   hRes->SetMarkerStyle(8);
   hRes->SetLineColor(4);
   hRes->SetMarkerColor(4);

   TCanvas* c1 = new TCanvas("c1", "", 600, 600);
   c1->SetLogy();
   hRes->SetAxisRange(0.002, 2, "Y");
   hRes->SetXTitle("N_{Hits}");
   hRes->SetTitle("Vertex Resolution");
   hRes->SetYTitle("Resolution (cm)");
   hRes->SetStats(0);
   hRes->Draw("same");

   // TLegend *l = new TLegend(0.5, 0.6, 0.9, 0.7);
   // l->AddEntry(hRes,"4-layer strip detectors","pl");
   // l->Draw();

   c1->SaveAs("vertex-res.png");
}

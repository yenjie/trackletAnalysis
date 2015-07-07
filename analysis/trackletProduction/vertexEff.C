#include <TProfile.h>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>

void vertexEff(const char* fname) {
   TFile* f = new TFile(fname);
   TTree* t = (TTree*)f->Get("TrackletTree12");
   TProfile *hTracklet = new TProfile("hTracklet", "", 20, 0, 40);
   t->Draw("vz[1]>-40:nhit1>>hTracklet", "abs(vz[0])<20", "goff");

   TProfile *hTracking = new TProfile("hTracking", "", 20, 0, 40);
   t->Draw("vz[2]>-40:nhit1>>hTracking", "abs(vz[0])<20", "goff");

   hTracklet->SetMarkerStyle(8);
   hTracklet->SetLineColor(4);
   hTracklet->SetMarkerColor(4);

   TCanvas* c1 = new TCanvas("c1", "", 600, 600);
   hTracklet->SetAxisRange(0.0, 1.1, "Y");
   hTracklet->SetTitle("Vertex reconstruction efficiency");
   hTracklet->SetXTitle("Layer 1 N_{Hits}");
   hTracklet->SetYTitle("Efficiency");
   hTracklet->SetStats(0);
   hTracklet->Draw();
   hTracking->SetLineColor(2);
   hTracking->SetMarkerColor(2);
   hTracking->Draw("same");
   TLegend* l = new TLegend(0.5, 0.6, 0.9, 0.72);
   l->SetBorderSize(0);
   l->SetFillStyle(0);
   l->AddEntry(hTracklet,"Tracklet Vertex", "pl");
   l->AddEntry(hTracking,"Line Tracking Vertex", "pl");
   l->Draw();

   c1->SaveAs("vertex-eff.png");
}

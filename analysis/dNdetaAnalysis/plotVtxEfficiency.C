#include <iostream>
#include <TProfile.h>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>

void plotVtxEfficiency(char *infile = "sample-7000GeV/TrackletTree-ATLAS-trackletVtx-1.root")
{
   TFile *inf = new TFile(infile);
   TTree *t = (TTree*)inf->FindObjectAny("TrackletTree12");
   TProfile *hTracklet = new TProfile("hTracklet","",20,0,20);
   TProfile *hCluster = new TProfile("hCluster","",20,0,20);
   TProfile *hTracking = new TProfile("hTracking","",20,0,20);
   t->Draw("vz[1]>-99:nhit1>>hTracklet","evtType!=92&&evtType!=93&&abs(vz[0])<10");
   t->Draw("vz[2]>-99:nhit1>>hTracking","evtType!=92&&evtType!=93&&abs(vz[0])<10");
   t->Draw("vz[4]>-999:nhit1>>hCluster","evtType!=92&&evtType!=93&&abs(vz[0])<10");
   
   hTracklet->SetLineColor(4);
   hCluster->SetLineColor(2);
   hTracklet->SetMarkerColor(4);
   hCluster->SetMarkerColor(2);
   hTracklet->SetMarkerSize(0.7);
   hTracking->SetMarkerSize(0.7);
   hCluster->SetMarkerSize(0.7);
   
   TCanvas *c = new TCanvas("c","",600,600);
//   c->SetLogy();
   hCluster->SetAxisRange(0.00,1.1,"Y");
   hCluster->SetXTitle("N_{Hits}^{Layer1} (B=0)");
   hCluster->SetTitle("");
   hCluster->SetYTitle("Vertex reconstruction efficiency");
   hCluster->Draw();
   TLine *l1 = new TLine(0,1,20,1);
   l1->Draw();
   hTracklet->Draw("same");
//   hTracking->Draw("same");
   hCluster->Draw("same c");
   hTracklet->Draw("same c");
//   hTracking->Draw("same c");
   
   TLegend *l = new TLegend(0.5,0.45,0.9,0.64);
   l->SetBorderSize(0);
   l->SetFillStyle(0);
   l->AddEntry(hCluster,"Cluster Vertex","pl");
   l->AddEntry(hTracklet,"Tracklet Vertex","pl");
//   l->AddEntry(hTracking,"Agglomerative Vertex","pl");
   l->Draw();
}

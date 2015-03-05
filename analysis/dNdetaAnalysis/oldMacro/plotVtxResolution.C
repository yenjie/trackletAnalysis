#include <iostream>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>

void plotVtxResolution(char *infile = "sample-7000GeV/TrackletTree-ATLAS-trackletVtx-1.root")
{
   TFile *inf = new TFile(infile);
   TTree *t = (TTree*)inf->FindObjectAny("TrackletTree12");
   TH2D *hTracklet = new TH2D("hTracklet","",150,0,150,1000,-2,2);
   TH2D *hCluster = new TH2D("hCluster","",150,0,150,1000,-10,10);
   TH2D *hTracking = new TH2D("hTracking","",150,0,150,1000,-2,2);
   t->Draw("vz[1]-vz[0]:nhit1>>hTracklet","vz[1]>-99&&nhit1>0&&abs(vz[1])<15","col");
   t->Draw("vz[2]-vz[0]:nhit1>>hTracking","vz[2]>-99&&nhit1>0&&abs(vz[2])<15","col");
   t->Draw("vz[4]-vz[0]:nhit1>>hCluster","vz[4]>-99&&nhit1>0&&abs(vz[4])<15","col");
   hTracklet->FitSlicesY();
   hTracking->FitSlicesY();
   hCluster->FitSlicesY();
   TH1D *hTrackletRes = (TH1D*)gDirectory->Get("hTracklet_2");
   TH1D *hTrackingRes = (TH1D*)gDirectory->Get("hTracking_2");
   TH1D *hClusterRes = (TH1D*)gDirectory->Get("hCluster_2");
   
   hTrackletRes->SetLineColor(4);
   hClusterRes->SetLineColor(2);
   hTrackletRes->SetMarkerColor(4);
   hClusterRes->SetMarkerColor(2);
   hTrackletRes->SetMarkerSize(0.5);
   hTrackingRes->SetMarkerSize(0.5);
   hClusterRes->SetMarkerSize(0.5);
   
   TCanvas *c = new TCanvas("c","",600,600);
   c->SetLogy();
   hClusterRes->SetAxisRange(0.002,2,"Y");
   hClusterRes->SetXTitle("N_{Hits} (Layer 1) B=0 sample");
   hClusterRes->SetTitle("");
   hClusterRes->SetYTitle("Resolution (cm)");
   hClusterRes->Draw();
   hTrackletRes->Draw("same");
   //hTrackingRes->Draw("same");
   
   TLegend *l = new TLegend(0.5,0.45,0.9,0.64);
   l->SetBorderSize(0);
   l->SetFillStyle(0);
   l->AddEntry(hClusterRes,"Cluster Vertex","pl");
   l->AddEntry(hTrackletRes,"Tracklet Vertex","pl");
   //l->AddEntry(hTrackingRes,"Agglomerative Vertex","pl");
   l->Draw();
}

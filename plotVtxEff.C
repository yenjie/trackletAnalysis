#include <TFile.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TText.h>

void plotVtxEff(){
   TFile *a = new TFile("sample/TrackletTree-10TeV-TV-2.4M.root");
   TNtuple *t10TeV = (TNtuple*) a->FindObjectAny("TrackletTree12");
   TFile *b = new TFile("sample/TrackletTree-900GeV-TV-3.2M.root");   
   TNtuple *t900GeV = (TNtuple*) b->FindObjectAny("TrackletTree12");

   // Vz distribution plot
   TCanvas *c1 = new TCanvas("c1","",400,400);
   TH1F *hVz900GeV = new TH1F("hVz900GeV","",100,-10,10);
   TH1F *hVz10TeV = new TH1F("hVz10TeV","",100,-10,10);
   
   int norm900GeV = t900GeV->GetEntries("vz[3]>-9999&&abs(vz[0])<20");
   int norm10TeV = t10TeV->GetEntries("vz[3]>-9999&&abs(vz[0])<20");
   t900GeV->Draw("vz[3]>>hVz900GeV","abs(vz[0])<20");
   t10TeV->Draw("vz[3]>>hVz10TeV","abs(vz[0])<20");

   hVz900GeV->Scale(1./norm900GeV);
   hVz10TeV->Scale(1./norm10TeV);
   
   hVz900GeV->SetLineColor(2);
   hVz10TeV->SetYTitle("Arbitrary Normalization");
   hVz10TeV->SetXTitle("Reconstructed V_{z} (cm)");
   hVz10TeV->SetTitleOffset(1.6,"Y");
   hVz10TeV->SetAxisRange(0,0.024,"Y");
   hVz10TeV->Draw();
   hVz900GeV->Draw("same");
   

   TLegend * leg = new TLegend(0.37,0.18,0.84,0.34);
   leg->SetFillStyle(0);  
   leg->SetFillColor(0); 
   leg->SetBorderSize(0);
   leg->AddEntry(hVz10TeV,"PYTHIA 10TeV","l");
   leg->AddEntry(hVz900GeV,"PYTHIA 900GeV","l");
   leg->Draw();   
   
   TText *cms = new TText(-9.5,0.022,"CMS Preliminary");
   cms->Draw();

   TText *labelA = new TText(-9,0.020,"(a)");
   labelA->Draw();

   c1->SaveAs("plot/VtxEff/VtxSpectrum.C");
   c1->SaveAs("plot/VtxEff/VtxSpectrum.gif");
   c1->SaveAs("plot/VtxEff/VtxSpectrum.eps");

   // Vertex reconstruction efficiency
   TProfile *hTV900GeV = new TProfile("hTV900GeV","",20,0,20);
   TProfile *hPV900GeV = new TProfile("hPV900GeV","",20,0,20);
   TProfile *hTV10TeV = new TProfile("hTV10TeV","",20,0,20);
   TProfile *hPV10TeV = new TProfile("hPV10TeV","",20,0,20);

   TCanvas *c2 = new TCanvas("c2","",400,400);
   t900GeV->Draw("vz[3]!=-9999:nhit1>>hTV900GeV");
   t900GeV->Draw("vz[1]!=-9999:nhit1>>hPV900GeV");
   t10TeV->Draw("vz[3]!=-9999:nhit1>>hTV10TeV");
   t10TeV->Draw("vz[1]!=-9999:nhit1>>hPV10TeV");
 
   hTV900GeV->SetMarkerColor(2);
   hTV900GeV->SetMarkerStyle(4);
   hTV900GeV->SetLineColor(2);
   hPV900GeV->SetMarkerColor(2);
   hPV900GeV->SetLineColor(2);

   hPV10TeV->SetMarkerStyle(21);
   hTV10TeV->SetMarkerStyle(25);
   
   hTV900GeV->SetYTitle("Vertex Reconstruction Efficiency");
   hTV900GeV->SetXTitle("Number of First Layer Hits");
   hTV900GeV->Draw("P");
   hPV900GeV->Draw("P same");
   hTV10TeV->Draw("P same");
   hPV10TeV->Draw("P same");
   
   TLegend * leg2 = new TLegend(0.4,0.27,0.99,0.48);
   leg2->SetFillStyle(0);  
   leg2->SetFillColor(0); 
   leg2->SetBorderSize(0);
   leg2->AddEntry(hPV10TeV,"Pixel Triplet Vertexing (10TeV)","pl");
   leg2->AddEntry(hTV10TeV,"Pixel Tracklet Vertexing (10TeV)","pl");
   leg2->AddEntry(hPV900GeV,"Pixel Triplet Vertexing (900GeV)","pl");
   leg2->AddEntry(hTV900GeV,"Pixel Tracklet Vertexing (900GeV)","pl");

   leg2->Draw();   
   
   TText *cms2 = new TText(8.87,0.48,"CMS Preliminary");
   cms2->Draw();

   TText *labelB = new TText(0.6,0.95,"(b)");
   labelB->Draw();
   
   c2->SaveAs("plot/VtxEff/VtxEff.C");
   c2->SaveAs("plot/VtxEff/VtxEff.gif");
   c2->SaveAs("plot/VtxEff/VtxEff.eps");
   
}

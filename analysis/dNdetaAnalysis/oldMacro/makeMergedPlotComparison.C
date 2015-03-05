#include <TH1F.h>
#include "makeMergedPlot.C"

void makeMergedPlotComparison(char *filename = "Run123596",bool putUA5=0)
{
   TH1F *hBit40 = makeMergedPlot(Form("%s-Official",filename));
   hBit40->SetName("hBit40");

   TH1F *hBit41 = makeMergedPlot(Form("%s-Official-41",filename));
   hBit41->SetName("hBit41");
   TH1F *hBit34 = makeMergedPlot(Form("%s-Official-34",filename));
   hBit34->SetName("hBit34");
   
   hBit41->SetMarkerColor(1);
   hBit41->SetLineColor(1);
   hBit34->SetMarkerColor(4);
   hBit34->SetLineColor(4);
   hBit34->SetMarkerSize(1.5);
   hBit40->SetMarkerSize(1.5);
   hBit41->SetMarkerSize(1.5);
   hBit34->Draw();
   hBit40->Draw("same");
   hBit41->Draw("same");
   
   TH1F *UA5 = getUA5NSD();
   if (putUA5) UA5->Draw("same");
   
   TLegend *leg = new TLegend(0.3,0.18,1,0.35,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("hTruth",Form("Data-%s",filename),"");
   entry=leg->AddEntry("hBit34","Reconstructed with Bit 34","pl");
   entry=leg->AddEntry("hBit40","Reconstructed with Bit 40","pl");
   entry=leg->AddEntry("hBit41","Reconstructed with Bit 41","pl");
   leg->Draw(); 
   
}

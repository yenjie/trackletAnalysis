   // ROOT Library
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TLine.h>
#include <TF1.h>
#include <TCut.h>
#include <TPad.h>
#include <TText.h>

#include "selectionCut.h"
void plotDataDrivenEff2(int TrackletType=12, double type=1)
{
   // data input
   TFile *f = new TFile("newMorePerfectSample/TrackletTree-Run132440-veryloosecuts.root");
   TTree * tData= dynamic_cast<TTree*>(f->Get(Form("TrackletTree%d",TrackletType)));
   
   // MC input
   TFile *fMC = new TFile("sample-7000GeV/TrackletTree-ATLAS-BPIX20.root");
   TTree * tMC = dynamic_cast<TTree*>(fMC->Get(Form("TrackletTree%d",TrackletType)));

   // MCPHOJET input
   TFile *fMCPHOJET = new TFile("sample-7000GeV/TrackletTree-ATLAS.root");
   TTree * tMCPHOJET = dynamic_cast<TTree*>(fMCPHOJET->Get(Form("TrackletTree%d",TrackletType)));

   //selectionCut myCut(1,0,10000);
   const int nTrackletBin =20;// myCut.nTrackletBin;
   double TrackletBins[nTrackletBin+1] = {-5,1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,70,300};
   
   TFile *foutput = new TFile("output.root","recreate");
   
   TCut dataCut = "l1TBit[4]!=0&&l1TBit[36]!=1&&l1TBit[37]!=1&&l1TBit[38]!=1&&l1TBit[39]!=1&&l1TBit[0]==1&&vtxQualCut==1&&nHFn>0&&nHFp>0";
   TCut MCCutNoVtxCut = "l1TBit[36]!=1&&l1TBit[37]!=1&&l1TBit[38]!=1&&l1TBit[39]!=1";
   TCut MCCut = "l1TBit[36]!=1&&l1TBit[37]!=1&&l1TBit[38]!=1&&l1TBit[39]!=1&&nHFn>0&&nHFp>0";
   TCut HFCut;
   
   if (type==1) {
      cout <<"HF cuts"<<endl;
      HFCut = "vz[1]>-99"; 
   } else {
      cout <<"BSC Cuts"<<endl;
      HFCut = "(l1TBit[40]==1||l1TBit[41]==1)";
   }
   
   TH1F *hData = new TH1F("hData","",nTrackletBin,TrackletBins);
   TH1F *hDataHFCut = new TH1F("hDataHFCut","",nTrackletBin,TrackletBins);

   TH1F *hMC = new TH1F("hMC","",nTrackletBin,TrackletBins);
   TH1F *hMCHFCut = new TH1F("hMCHFCut","",nTrackletBin,TrackletBins);
   TH1F *hMCPHOJET = new TH1F("hMCPHOJET","",nTrackletBin,TrackletBins);
   TH1F *hMCPHOJETHFCut = new TH1F("hMCPHOJETHFCut","",nTrackletBin,TrackletBins);
   
   TCanvas *c = new TCanvas("c","",600,600);
   tData->Draw("mult2>>hDataHFCut",dataCut&&HFCut);
   tData->Draw("mult2>>hData",dataCut);
   tMC->Draw("mult2>>hMCHFCut",MCCut&&HFCut);
   tMC->Draw("mult2>>hMC",MCCut);
   tMCPHOJET->Draw("mult2>>hMCPHOJETHFCut",MCCut&&HFCut);
   tMCPHOJET->Draw("mult2>>hMCPHOJET",MCCut);

   hDataHFCut->SetXTitle("N_{Tracklets}");
   hDataHFCut->SetYTitle("Efficiency");
   hDataHFCut->Sumw2();
   hData->Sumw2();
   hDataHFCut->Divide(hDataHFCut,hData,1,1,"b"); 
   hDataHFCut->Draw();  

   hMCHFCut->Sumw2();
   hMC->Sumw2();
   hMCHFCut->Divide(hMCHFCut,hMC,1,1,"b"); 
   hMCHFCut->SetLineColor(2);
   hMCHFCut->SetMarkerColor(2);
   hMCHFCut->SetMarkerStyle(4);
   hMCHFCut->Draw("same");  

   hMCPHOJETHFCut->Sumw2();
   hMCPHOJET->Sumw2();
   hMCPHOJETHFCut->Divide(hMCPHOJETHFCut,hMCPHOJET,1,1,"b"); 
   hMCPHOJETHFCut->SetLineColor(4);
   hMCPHOJETHFCut->SetMarkerColor(4);
   hMCPHOJETHFCut->SetMarkerStyle(4);
   hMCPHOJETHFCut->Draw("same");  
   
   TLine *l1 = new TLine(0,1,300,1);
   l1->Draw();
      
   TLegend *leg = new TLegend(0.43,0.22,0.92,0.45);
   leg->SetBorderSize(0);
   leg->SetFillColor(0);
   if (type==1) {
      leg->AddEntry("hDataHFCut","Data driven HF efficiency","pl");
      leg->AddEntry("hMCHFCut","PYTHIA ATLAS HF efficiency","pl");
      leg->AddEntry("hMCPHOJETHFCut","PHOJET HF efficiency","pl");
   } else {
      leg->AddEntry("hDataHFCut","Data driven BSC efficiency","pl");
      leg->AddEntry("hMCHFCut","PYTHIA ATLAS BSC efficiency","pl");
      leg->AddEntry("hMCPHOJETHFCut","PHOJET BSC efficiency","pl");
   }
   leg->Draw();

   TCanvas *c2 = new TCanvas("c2","",600,600);
   TH1F *hRatio = (TH1F*) hDataHFCut->Clone();
   hRatio->SetName("hRatio");
   hRatio->SetLineColor(2);
   hRatio->SetMarkerColor(2);
   hRatio->SetMarkerStyle(4);
   hRatio->Divide(hMCHFCut);
   hRatio->SetYTitle("Efficiency Ratio");
   hRatio->Draw();
   TH1F *hRatioPHOJET = (TH1F*) hDataHFCut->Clone();
   hRatioPHOJET->SetName("hRatioPHOJET");
   hRatioPHOJET->Divide(hMCPHOJETHFCut);
   hRatioPHOJET->SetLineColor(4);
   hRatioPHOJET->SetMarkerColor(4);
   hRatioPHOJET->SetMarkerStyle(4);
   hRatioPHOJET->Draw("same");

   TLegend *leg2 = new TLegend(0.43,0.72,0.92,0.85);
   leg2->SetBorderSize(0);
   leg2->SetFillColor(0);
   
   if (type==1) {
      leg2->AddEntry("hRatio","Data / PYTHIA ATLAS HF efficiency","pl");
      leg2->AddEntry("hRatioPHOJET","Data /PHOJET HF efficiency","pl");
   } else {
      leg2->AddEntry("hRatio","Data / PYTHIA ATLAS BSC efficiency","pl");
      leg2->AddEntry("hRatioPHOJET","Data /PHOJET BSC efficiency","pl");
   }
   l1->Draw();
   leg2->Draw();
   
}

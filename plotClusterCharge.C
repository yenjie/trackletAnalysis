#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>

void plotClusterCharge()
{
   TFile *inf1 = new TFile("newMorePerfectSample/PixelTree-13324-hfcuts.root");
   TTree *tData = (TTree*) inf1->FindObjectAny("PixelTree");
      
   TFile *inf2 = new TFile("sample-7000GeV/PixelTree-7TeV-ATLAS_B0.root");
   TTree *tMC = (TTree*) inf2->FindObjectAny("PixelTree");   
   
   TH1F *hData = new TH1F("hData","",20,0,4000);
   TH1F *hMC = new TH1F("hMC","",20,0,4000);
   TH1F *hData2 = new TH1F("hData2","",20,0,4000);
   TH1F *hMC2 = new TH1F("hMC2","",20,0,4000);

  
   tData->Draw("ch2/cosh(eta2)/1000.>>hData","nhits1<100&&nHFn>0&&nHFp>0&&L1A[124]==1");
   tMC->Draw("ch2/cosh(eta2)/1000.>>hMC","nhits1<100&&nHFn>0&&nHFp>0&&L1A[124]==1");
   tData->Draw("ch2/cosh(eta2)/1000.>>hData2","nhits1>100&&nHFn>0&&nHFp>0&&L1A[124]==1");
   tMC->Draw("ch2/cosh(eta2)/1000.>>hMC2","nhits1>100&&nHFn>0&&nHFp>0&&L1A[124]==1");
   
   int nData = hData->GetEntries();
   int nMC = hMC->GetEntries();
   int nData2 = hData2->GetEntries();
   int nMC2 = hMC2->GetEntries();
   
   hData->Sumw2();
   hMC->Sumw2();
   hData2->Sumw2();
   hMC2->Sumw2();
   
   hData->Scale(1./nData);
   hMC->Scale(1./nMC);
   hData2->Scale(1./nData2);
   hMC2->Scale(1./nMC2);
   
   hData->SetLineColor(2);
   hData->SetMarkerColor(2);
   
//   hData->Draw();
   hData->Divide(hMC);
   hData2->Divide(hMC2);
   hData->Divide(hData2);
   hData->Draw();
//   hData2->Draw("same");
//   hMC->Draw("same");
   
}



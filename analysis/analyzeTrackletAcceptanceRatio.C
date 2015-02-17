#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TF1.h>
#include "Math/Vector3D.h"
#include <TMath.h>
#include <TCut.h>
#include <TChain.h>

/*
double getRandomVz()
{
   int flag=0;
   while (flag==0) {
      double myVz = gRandom->Rndm()*-20+20;
      if (gRandom->Rndm()<TMath::Gaus(myVz,-2.702,3.627,1)) return myVz;
   }
}*/

void normalize(TH2F *hData,int nEtaBin,int nPhiBin)
{
   double entries = hData->GetEntries();
   entries/=(nEtaBin*nPhiBin*10);
   entries=0;
   for (int x = 0; x <= nEtaBin+1; x++)
   {
      for (int y = 0; y <= nPhiBin+1; y++)
      {
         if (hData->GetBinContent(x,y)>entries
	     &&x!=0&&x<=nEtaBin
	     &&y!=0&&y<=nPhiBin
	    ) 
	 {
	    hData->SetBinContent(x,y,1);
	 } else {
	    hData->SetBinContent(x,y,0);
	 }
      }
   }
}

void analyzeTrackletAcceptanceRatio(int ch=13)

{
//  infNameData = "newBestDataSample/TrackletTree-Run124120.root";
//  infNameData = "./sample-official2/Pythia_MinBias_D6T_900GeV_d20091210-BSCTuned_v1.root";
//   infNameMC = "./sample-official2/TrackletTree-900GeV-D6T-BSCTuned.root";
//   cout <<infNameData<<" "<<infNameMC<<endl;
/*
   TChain *tMC = new TChain(Form("TrackletTree%d",ch));
   tMC->AddFile("sample-1224/TrackletTree-2360GeV-D6T.root");
   tMC->AddFile("sample-1224/TrackletRndVtx-1.root");
   tMC->AddFile("sample-1224/TrackletRndVtx-2.root");
   tMC->AddFile("sample-1224/TrackletRndVtx-3.root");
   tMC->AddFile("sample-1224/TrackletRndVtx-4.root");
//   tData->SetName("tData");
   TChain *tData = new TChain(Form("TrackletTree%d",ch));
   tData->AddFile("newBestDataSample/TrackletTree-Run124022-partial-minBias.root");
//   tData->AddFile("newBestDataSample/Tracklet-Data-RandomVertex.root");
   tData->AddFile("newBestDataSample/Tracklet-RndVtx-1.root");
   tData->AddFile("newBestDataSample/Tracklet-RndVtx-2.root");
   tData->AddFile("newBestDataSample/Tracklet-RndVtx-3.root");
   tData->AddFile("newBestDataSample/Tracklet-RndVtx-4.root");
   tData->AddFile("newBestDataSample/Tracklet-RndVtx-5.root");
   tData->SetMarkerSize(0.1);
*/

   TChain *tMC = new TChain(Form("TrackletTree%d",ch));
   tMC->AddFile("sample-7000GeV/TrackletTree-D6T-randomVtx-1.root");
   tMC->AddFile("sample-7000GeV/TrackletTree-D6T-randomVtx-2.root");
   tMC->AddFile("sample-7000GeV/TrackletTree-D6T-randomVtx-3.root");
   tMC->AddFile("sample-7000GeV/TrackletTree-D6T-randomVtx-4.root");
   tMC->AddFile("sample-7000GeV/TrackletTree-D6T-randomVtx-5.root");
   tMC->AddFile("sample-7000GeV/TrackletTree-D6T-randomVtx-6.root");
   tMC->AddFile("sample-7000GeV/TrackletTree-D6T-randomVtx-7.root");
//   tData->SetName("tData");
   TChain *tData = new TChain(Form("TrackletTree%d",ch));
   tData->AddFile("newUltimateSample-Random/TrackletTree-Run132440-Random-1.root");
   tData->AddFile("newUltimateSample-Random/TrackletTree-Run132440-Random-2.root");
   tData->AddFile("newUltimateSample-Random/TrackletTree-Run132440-Random-3.root");
   tData->AddFile("newUltimateSample-Random/TrackletTree-Run132440-Random-4.root");
   tData->AddFile("newUltimateSample-Random/TrackletTree-Run132440-Random-5.root");
   tData->AddFile("newUltimateSample-Random/TrackletTree-Run132440-Random-6.root");
   

//   tMC->SetName("tMC");
   tMC->SetMarkerSize(0.1);
   
   int nEtaBin=300*4;
   int nPhiBin=316*4;
   int nVzBin=300*4;
   int VzRangeL=-20;
   int VzRangeH=20;
   
   TCut myCut = "abs(deta)<0.1&&abs(dphi)<1&&vz[1]>-99";
   TFile *outfile = new TFile(Form("TrackletAccStudy-%d.root",ch),"recreate");
   
   TH2F *hData = new TH2F("hData","",nEtaBin,-3,3,nVzBin,VzRangeL,VzRangeH);
   TH2F *hMC = new TH2F("hMC","",nEtaBin,-3,3,nVzBin,VzRangeL,VzRangeH);
   TH2F *hAccData = new TH2F("hAccData","",nEtaBin,-3,3,nVzBin,VzRangeL,VzRangeH);
   TH2F *hAccMC = new TH2F("hAccMC","",nEtaBin,-3,3,nVzBin,VzRangeL,VzRangeH);
   tData->Project("hData",Form("vz[1]:eta1"),myCut,"");
   tMC->Project("hMC",Form("vz[1]:eta1"),myCut,"");

   normalize(hData,nEtaBin,nVzBin);
   normalize(hMC,nEtaBin,nVzBin);
   

   
   TCanvas *c = new TCanvas("c","Data",600,600);      
   hData->SetLineColor(2);
   hData->Draw("box");
//   tData->Draw(Form("vz[1]:eta1"),myCut,"same");
   TCanvas *c2 = new TCanvas("c2","MC",600,600);      
   hMC->SetLineColor(2);
   hMC->Draw("box");
//   tMC->Draw(Form("vz[1]:eta1"),myCut,"same",200000);
   
   TCanvas *c3 = new TCanvas("c3","Data & MC",600,600);
   
   TH1F *hDataEta =  (TH1F*)hData->ProjectionX();
   TH1F *hMCEta =  (TH1F*)hMC->ProjectionX();
  

   hDataEta->SetLineColor(1);
   hDataEta->SetXTitle("#eta");
   hDataEta->Rebin(nEtaBin/12);
   hMCEta->SetLineColor(2);
   hMCEta->Rebin(nEtaBin/12);
   hDataEta->Draw("");
   hMCEta->Draw("same");
   
   TCanvas *c4 = new TCanvas("c4","Ratio",600,600);
   
   TH1F* hRatio = (TH1F*) hMCEta->Clone();
   hRatio->SetName("hRatio");
   hRatio->Divide(hDataEta);
   hRatio->Draw();

   TH2F *hDataAcc = (TH2F*)hData->Clone();
   hDataAcc->SetName("hDataAcc");
   hDataAcc->RebinX(nEtaBin/12.);
   hDataAcc->RebinY(nVzBin/20.);
   TH2F *hMCAcc = (TH2F*)hMC->Clone();
   hMCAcc->SetName("hMCAcc");
   hMCAcc->RebinX(nEtaBin/12.);
   hMCAcc->RebinY(nVzBin/20.);

   outfile->Write();

} 

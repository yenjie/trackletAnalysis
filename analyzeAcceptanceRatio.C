#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TF1.h>
#include "Math/Vector3D.h"
#include <TMath.h>


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
   for (int x = 1; x <= nEtaBin; x++)
   {
      for (int y = 1; y <= nPhiBin; y++)
      {
         if (hData->GetBinContent(x,y)>entries) hData->SetBinContent(x,y,1);
	 else 	 hData->SetBinContent(x,y,0);
      }
   }
}

void analyzeAcceptanceRatio(int ch=1, 
//char *infNameData = "DataSample/pixeltree_run123596.root", 
char *infNameData = "newBestDataSample/pixel-96-22.root", 
char *infNameMC = "sample-official/Pythia_MinBias_D6T_9000GeV_d20091210.root")

{
//  infNameData = "./sample-1213/Pythia_MinBias_ATLAS_900GeV_d20091208.root";
//  infNameData = "./sample-official2/Pythia_MinBias_D6T_900GeV_d20091210-BSCTuned_v1.root";

   cout <<infNameData<<" "<<infNameMC<<endl;
   TFile *infData = new TFile (infNameData);
   TTree *tData = (TTree*)infData->FindObjectAny("PixelTree");
   tData->SetName("tData");
   
   TFile *infMC = new TFile (infNameMC);
   TTree *tMC = (TTree*)infMC->FindObjectAny("PixelTree");
   tMC->SetName("tMC");
   
   int nEtaBin=300*4;
   int nPhiBin=316*2;
   int nVzBin=100;
   

   TFile *outfile = new TFile(Form("acceptanceStudy-%d.root",ch),"recreate");
   
   TH2F *hData = new TH2F("hData","",nEtaBin,-3,3,nPhiBin,-3.16,3.16);
   TH2F *hMC = new TH2F("hMC","",nEtaBin,-3,3,nPhiBin,-3.16,3.16);
   TH2F *hAccData = new TH2F("hAccData","",nEtaBin,-3,3,nVzBin,-20,20);
   TH2F *hAccMC = new TH2F("hAccMC","",nEtaBin,-3,3,nVzBin,-20,20);
   tData->Project("hData",Form("phi%d:eta%d",ch,ch,"","",10000000));
   tMC->Project("hMC",Form("phi%d:eta%d",ch,ch),"","",10000000);

   normalize(hData,nEtaBin,nPhiBin);
   normalize(hMC,nEtaBin,nPhiBin);
   

   
   TCanvas *c = new TCanvas("c","Data",600,600);      
   hData->SetLineColor(2);
   hData->Draw("box");
   tData->Draw(Form("phi%d:eta%d",ch,ch),"","same");
   TCanvas *c2 = new TCanvas("c2","MC",600,600);      
   hMC->SetLineColor(2);
   hMC->Draw("box");
   tMC->Draw(Form("phi%d:eta%d",ch,ch),"","same",200000);
   
   TCanvas *c3 = new TCanvas("c3","Data & MC",600,600);
   
   TH1F *hDataEta =  (TH1F*)hData->ProjectionX();
   TH1F *hMCEta =  (TH1F*)hMC->ProjectionX();
  

   hDataEta->SetLineColor(1);
   hDataEta->SetXTitle("#eta");
//   hDataEta->Rebin(nEtaBin/12);
   hMCEta->SetLineColor(2);
//   hMCEta->Rebin(nEtaBin/12);
   hDataEta->Draw("");
   hMCEta->Draw("same");
   
   TCanvas *c4 = new TCanvas("c4","Ratio",600,600);
   
   TH1F* hRatio = (TH1F*) hMCEta->Clone();
   hRatio->SetName("hRatio");
   hRatio->Divide(hDataEta);
   hRatio->Draw();
   
   int nToyMC = 50;
   
   double rL1 = 4.39;
   double rL2 = 7.27;
   double rL3 = 10.19;
   
   double r=rL1;
   if (ch==2) r=rL2;
   if (ch==3) r=rL3;
   
   double VxGen = 0.174562;
   double VyGen = 0.144887;

   double Vx = 0.205124;
   double Vy = 0.164012;

   TF1 *fVz = new TF1("f","gaus",-30,30);
   fVz->SetParameters(1,-2.702,3.627);
   
   double cache[1000][1000];
   for (int x=0;x<1000;x++){
      for (int y=0;y<1000;y++){
         cache[x][y]=0;
      }
   }
   
   for (int i=0; i<nToyMC; i++){
      cout <<i<<endl;
      double myVz = fVz->GetRandom();
      
      for (int xx=1;xx<=nEtaBin;xx++){
         for (int yy=1;yy<=nPhiBin;yy++) {
	  
            if (hData->GetBinContent(xx,yy)==0&&hMC->GetBinContent(xx,yy)==0) continue;
	    double eta = -3+6./nEtaBin*xx;
	    double phi = -3.16+6.32/nPhiBin*xx;
	    double x = r*cos(phi);
	    double y = r*sin(phi);
	    double z = r/tan(atan(exp(-eta))*2)-myVz;
	    ROOT::Math::XYZVector tmpVector(x-Vx,y-Vy,z);
	    double cEta = tmpVector.eta();
	    ROOT::Math::XYZVector tmpVector2(x-VxGen,y-VyGen,z);
	    double cEta2 = tmpVector2.eta();
	    
	    //double cPhi = tmpVector.phi();
	    //cout <<cEta<<" "<<myVz<<endl;
	    int bin = hAccData->FindBin(cEta,myVz);
	    if (hData->GetBinContent(xx,yy)!=0) hAccData->Fill(cEta,myVz);
	    if (hMC->GetBinContent(xx,yy)!=0) hAccMC->Fill(cEta2,myVz);
	 }
      }
   }
   
   TCanvas *c5 = new TCanvas("c5","hAccData",600,600);   
   hAccData->Draw("col");

   TCanvas *c6 = new TCanvas("c6","hAccMC",600,600);   
   hAccMC->Draw("col");

   TCanvas *c7 = new TCanvas("c7","AccRatio",600,600);
   TH1F *hProjData = (TH1F*)hAccData->ProjectionX();
   hProjData->SetName("hProjData");
   TH1F *hProjMC = (TH1F*)hAccMC->ProjectionX();
   hProjMC->SetName("hProjMC");
   hProjData->Rebin(600/12.);
   hProjMC->Rebin(600/12.);
   hProjData->Divide(hProjMC);
   hProjData->SetName("hProjData");
   hProjData->Draw();
   outfile->Write();
   //outfile->Close();
} 

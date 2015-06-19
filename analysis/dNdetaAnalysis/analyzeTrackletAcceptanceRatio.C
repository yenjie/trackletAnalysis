#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TCut.h>

void normalize(TH2F *hData, int nEtaBin, int nVzBin, bool reweight = 1) {
/*
   1  Constant     1.86973e-01   5.69091e-04   2.02006e-06   3.48340e-01
   2  Mean        -2.14145e+00   1.08618e-02   4.84445e-05  -5.12350e-04
   3  Sigma        4.30854e+00 
  */
   TH1D *hTmp = (TH1D*)hData->ProjectionY("hTmp");
   for (int x=1; x<=nEtaBin; x++) {
      for (int y=1; y<=nVzBin; y++) {
         double myVz = hTmp->GetBinCenter(y);
         double DataPdf = TMath::Gaus(myVz, -2.14145, 4.30854, 1);
         if (!reweight) DataPdf=1;
         if (hData->GetBinContent(x, y)>0 && x!=0 && x<=nEtaBin && y!=0 && y<=nVzBin) {
            hData->SetBinContent(x, y, DataPdf);
         } else {
            hData->SetBinContent(x, y, 0);
         }
      }
   }
}

void analyzeTrackletAcceptanceRatio(int TrackletType, const char* fnMC="/data/biran/trackletAnalysis/analysis/trackletProduction/TrackletTree-PYTHIA8-OFFICIAL-ACCEPTANCE.root", const char* fnData="/data/biran/trackletAnalysis/analysis/trackletProduction/TrackletTree-Run247324-PromptReco-ACCEPTANCE.root") {
   TFile* fMC = new TFile(fnMC, "READ");
   TTree* tMC = (TTree*)fMC->Get(Form("TrackletTree%i", TrackletType));
   TFile* fData = new TFile(fnData, "READ");
   TTree* tData = (TTree*)fData->Get(Form("TrackletTree%i", TrackletType));

   int nEtaBin = 1200;
   int nVzBin = 1100;
   int VzRangeL = -13;
   int VzRangeH = 9;

   TFile *outfile = new TFile(Form("acceptance-%d.root", TrackletType), "recreate");

   TCut myCut = "abs(deta)<0.04 && abs(dphi)<0.04 && vz[1]>-13 && vz[1]<9";
   TH2F *hData = new TH2F("hData", "", nEtaBin, -3, 3, nVzBin, VzRangeL, VzRangeH);
   TH2F *hMC = new TH2F("hMC", "", nEtaBin, -3, 3, nVzBin, VzRangeL, VzRangeH);
   TH2F *hAccData = new TH2F("hAccData", "", nEtaBin, -3, 3, nVzBin, VzRangeL, VzRangeH);
   TH2F *hAccMC = new TH2F("hAccMC", "", nEtaBin, -3, 3, nVzBin, VzRangeL, VzRangeH);
   tData->Project("hData", Form("vz[1]:eta1"), myCut);
   tMC->Project("hMC", Form("vz[1]:eta1"), myCut);

   normalize(hData, nEtaBin, nVzBin);
   normalize(hMC, nEtaBin, nVzBin);

   TCanvas *c1 = new TCanvas("c1", "Data", 600, 600);
   hData->SetLineColor(2);
   hData->Draw("box");

   TCanvas *c2 = new TCanvas("c2", "MC", 600, 600);
   hMC->SetLineColor(4);
   hMC->Draw("box");

   TCanvas *c3 = new TCanvas("c3", "Data & MC", 600, 600);
   TH1F *hDataEta = (TH1F*)hData->ProjectionX();
   TH1F *hMCEta = (TH1F*)hMC->ProjectionX();

   hDataEta->SetLineColor(1);
   hDataEta->SetXTitle("#eta");
   hDataEta->Rebin(nEtaBin/30);
   hMCEta->SetLineColor(2);
   hMCEta->Rebin(nEtaBin/30);
   hDataEta->Draw();
   hMCEta->Draw("same");

   TCanvas *c4 = new TCanvas("c4", "Ratio", 600, 600);
   TH1F* hRatio = (TH1F*) hMCEta->Clone();
   hRatio->SetName("hRatio");
   hRatio->Divide(hDataEta);
   hRatio->Draw();

   TH2F *hDataAcc = (TH2F*)hData->Clone();
   hDataAcc->SetName("hDataAcc");
   hDataAcc->RebinX(nEtaBin/30);
   hDataAcc->RebinY(nVzBin/11);

   TH2F *hMCAcc = (TH2F*)hMC->Clone();
   hMCAcc->SetName("hMCAcc");
   hMCAcc->RebinX(nEtaBin/30);
   hMCAcc->RebinY(nVzBin/11);

   outfile->Write();
}

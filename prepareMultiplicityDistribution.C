#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TF1.h"
#include "TTimeStamp.h"
#include <TCanvas.h>

void prepareMultiplicityDistribution()
{
    char *mult = "mult";
    TFile *inf = new TFile("newMorePerfectSample/TrackletTree-Run132440-hfcuts-final.root");
    TTree *tData = (TTree*) inf->FindObjectAny("TrackletTree12");
//    TFile *inf2 = new TFile("sample-7000GeV/TrackletTree-ATLAS.root");
    TFile *inf2 = new TFile("test.root");
    TTree *tMC = (TTree*) inf2->FindObjectAny("TrackletTree12");
    
    TH1F *hData = new TH1F("hData","",61,-5,300);
    TH1F *hMC = new TH1F("hMC","",61,-5,300);
    
    TCanvas *c = new TCanvas("c","",600,600);

    tData->Draw(Form("%s>>hData",mult),"vz[1]>-99&&nHFn>0&&nHFp>0");
    tMC->Draw(Form("%s>>hMC",mult),"vz[1]>-99&&nHFn>0&&nHFp>0");
    
    hMC->Sumw2();
    hData->Sumw2();
    hMC->Scale(1./hMC->GetEntries());
    hData->Scale(1./hData->GetEntries());

    hData->SetLineColor(2);
    hData->SetMarkerColor(2);
    
    hData->Draw();
    hMC->Draw("same");
    
    TCanvas *c2 = new TCanvas("c2","",600,600);
    TH1F *ratio = (TH1F*)hData->Clone();
    ratio->Divide(hMC);
    ratio->Draw();
}

// Plot final results

// Standard library
#include <math.h>
#include <iostream>
#include <fstream>

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

void analyzeBetaCorrection(int TrackletType = 12)
{
   // Definition of Vz, Eta, Hit bins
   const int nTrackletBin =14;// myCut.nTrackletBin;
   const int nEtaBin =12;// myCut.nEtaBin;
   const int nVzBin  =20;// myCut.nVzBin;
   int VzRangeL = -20;
   int VzRangeH = 20;
   double HitBins[nTrackletBin+1] = {-5,5,10,15,20,25,30,35,40,45,50,55,60,65,100};
   double EtaBins[nEtaBin+1];
   double VzBins[nVzBin+1];
   for (int i=0;i<=nEtaBin;i++) EtaBins[i] = (double)i*6.0/(double)nEtaBin-3.0;;
   for (int i=0;i<=nVzBin;i++) VzBins[i] = (double)i*(VzRangeH-VzRangeL)/(double)nVzBin+VzRangeL;

   
   TFile *myFile = new TFile(Form("correction/alphaBeta-%d.root",TrackletType));
   TNtuple *betaCorrection = (TNtuple*)myFile->FindObjectAny("correction");
   
   TFile *outFile = new TFile(Form("correction/alphaBetaCoeff-%d.root",TrackletType),"recreate");
   TNtuple * correction = new TNtuple("correction","","eta:nTracklet:vz:a:b");

   TH3F *hAlphaA =
         new TH3F("hAlphaA","",nEtaBin,EtaBins,nTrackletBin,HitBins,nVzBin,VzBins);
   TH3F *hAlphaB =
         new TH3F("hAlphaB","",nEtaBin,EtaBins,nTrackletBin,HitBins,nVzBin,VzBins);

   TH2F *hAlphaBeta = new TH2F("hAlphaBeta","",1000,0,0.5,1000,0,10);

   if (!betaCorrection) cout <<"Error reading betaCorrection!"<<endl;
   betaCorrection->SetName("betaCorrection");
   for (int x=1;x<=nEtaBin;x++) {
      for (int y=1;y<=nTrackletBin;y++) {
         for (int z=1;z<=nVzBin;z++) {
   	       double a=1,b=1;
	       cout <<Form("eta==%f&&nTracklet==%f&&vz==%f&&alpha!=0",(EtaBins[x]+EtaBins[x-1])/2,(HitBins[y]+HitBins[y-1])/2,(VzBins[z]+VzBins[z-1])/2)<<endl;
	       betaCorrection->Project("hAlphaBeta","alpha:beta",Form("eta==%f&&nTracklet==%f&&vz==%f&&alpha!=0",(EtaBins[x]+EtaBins[x-1])/2,(HitBins[y]+HitBins[y-1])/2,(VzBins[z]+VzBins[z-1])/2));
               hAlphaBeta->Fit("pol1","m E Q");
	       TF1 *fLine = hAlphaBeta->GetFunction("pol1");
               if (fLine) {
	          a=fLine->GetParameter(0);
	          b=fLine->GetParameter(1);
	       }
	       hAlphaA->SetBinContent(x,y,z,a);
	       hAlphaB->SetBinContent(x,y,z,b);
    	       correction->Fill((EtaBins[x]+EtaBins[x-1])/2,(HitBins[y]+HitBins[y-1])/2,(VzBins[z]+VzBins[z-1])/2,a,b);

         }
      }
   }
   outFile->Write();
   outFile->Close();   
}

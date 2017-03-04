#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TCanvas.h"

void normalize(TH2F* hdata, int nEtaBin, int nVzBin, bool reweight = 1) {
   TH1D* hvz = (TH1D*)hdata->ProjectionY("hvz");

   for (int x=1; x<=nEtaBin; x++) {
      for (int y=1; y<=nVzBin; y++) {
         double vz = hvz->GetBinCenter(y);
         // Run 285832
         double data_pdf = TMath::Gaus(vz, 1.00333, 4.65240, 1);
         if (!reweight) data_pdf = 1;
         if (hdata->GetBinContent(x, y) != 0) {
            hdata->SetBinContent(x, y, data_pdf);
         } else {
            hdata->SetBinContent(x, y, 0);
         }
      }
   }

   hvz->Delete();
}

int tracklet_acceptances(int TrackletType, const char* mc_file, const char* data_file, int nbins = 500) {
   TFile* fmc = new TFile(mc_file, "READ");
   TTree* tmc = (TTree*)fmc->Get(Form("TrackletTree%i", TrackletType));
   TFile* fdata = new TFile(data_file, "READ");
   TTree* tdata = (TTree*)fdata->Get(Form("TrackletTree%i", TrackletType));

   int nEtaBin = nbins * 12;
   int nVzBin = nbins * 14;
   int VzRangeL = -13;
   int VzRangeH = 15;

   printf("with %i eta bins, %i vz bins, %i < vz < %i\n", nEtaBin, nVzBin, VzRangeL, VzRangeH);

   TFile* fout = new TFile(Form("acceptance-%d.root", TrackletType), "recreate");
   TH2F* hdata = new TH2F("hdata", "", nEtaBin, -3, 3, nVzBin, VzRangeL, VzRangeH);
   TH2F* hmc = new TH2F("hmc", "", nEtaBin, -3, 3, nVzBin, VzRangeL, VzRangeH);

   TCut signal_cut = "abs(deta)<0.1 && abs(dphi)<1 && vz[1]>-13 && vz[1]<15";
   printf("projecting...\n");
   tdata->Project("hdata", Form("vz[1]:eta1"), signal_cut);
   tmc->Project("hmc", Form("vz[1]:eta1"), signal_cut);

   printf("calculating acceptances...\n");
   normalize(hdata, nEtaBin, nVzBin);
   normalize(hmc, nEtaBin, nVzBin);

   TH2F* haccep_data = (TH2F*)hdata->Clone("haccep_data");
   haccep_data->RebinX(nEtaBin/30);
   haccep_data->RebinY(nVzBin/14);

   TH2F* haccep_mc = (TH2F*)hmc->Clone("haccep_mc");
   haccep_mc->RebinX(nEtaBin/30);
   haccep_mc->RebinY(nVzBin/14);

   TCanvas* c1 = new TCanvas("c1", "", 600, 600);
   TH2F* hratio = (TH2F*)haccep_mc->Clone("hratio");
   hratio->Divide(haccep_data);
   hratio->Draw("colz");

   fout->Write("", TObject::kOverwrite);

   return 0;
}

int main(int argc, char* argv[]) {
   if (argc == 4)
      return tracklet_acceptances(atoi(argv[1]), argv[2], argv[3]);
   else if (argc == 5)
      return tracklet_acceptances(atoi(argv[1]), argv[2], argv[3], atoi(argv[4]));
   else
      return 1;
}

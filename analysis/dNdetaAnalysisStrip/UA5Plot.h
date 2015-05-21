#include <TH1F.h>

TH1F* getUA5() {
   TH1F* hEta_UA5 = new TH1F("hEta_UA5", "hEta_UA5", 600, -3., 3.);

   hEta_UA5->SetBinContent(hEta_UA5->FindBin(0.125), 3.1);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(0.375), 3.05);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(0.625), 3.1);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(0.875), 3.3);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(1.125), 3.3);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(1.375), 3.5);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(1.625), 3.45);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(1.875), 3.4);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(2.125), 3.45);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(2.375), 3.4);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(2.625), 3.3);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(2.875), 3.05);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(3.125), 2.9);

   hEta_UA5->SetBinContent(hEta_UA5->FindBin(-0.125), 3.1);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(-0.375), 3.05);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(-0.625), 3.1);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(-0.875), 3.3);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(-1.125), 3.3);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(-1.375), 3.5);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(-1.625), 3.45);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(-1.875), 3.4);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(-2.125), 3.45);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(-2.375), 3.4);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(-2.625), 3.3);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(-2.875), 3.05);
   hEta_UA5->SetBinContent(hEta_UA5->FindBin(-3.125), 2.9);
   hEta_UA5->SetMarkerColor(kGreen-2);
   hEta_UA5->SetMarkerSize(0.5);
   hEta_UA5->SetLineColor(2);
   hEta_UA5->SetMarkerStyle(4);

   hEta_UA5->SetLineColor(2);
   hEta_UA5->SetMarkerStyle(4);
   return hEta_UA5;
}

TH1F *getUA5NSD() {
   TH1F* hEta_UA5_NSD = new TH1F("hEta_UA5_NSD", ";#eta;dN/d#eta", 50, -3, 3);

   // positive eta
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(0.125), 3.48);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(0.375), 3.38);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(0.625), 3.52);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(0.875), 3.68);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.125), 3.71);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.375), 3.86);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.625), 3.76);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.875), 3.66);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(2.125), 3.72);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(2.375), 3.69);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(2.625), 3.56);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(2.875), 3.41);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(3.125), 3.15);

   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(0.125), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(0.375), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(0.625), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(0.875), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.125), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.375), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.625), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.875), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(2.125), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(2.375), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(2.625), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(2.875), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(3.125), 0.07);

   //negative eta
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(-0.125), 3.48);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(-0.375), 3.38);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(-0.625), 3.52);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(-0.875), 3.68);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(-1.125), 3.71);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(-1.375), 3.86);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(-1.625), 3.76);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(-1.875), 3.66);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(-2.125), 3.72);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(-2.375), 3.69);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(-2.625), 3.56);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(-2.875), 3.41);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(-3.125), 3.15);

   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(-0.125), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(-0.375), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(-0.625), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(-0.875), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(-1.125), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(-1.375), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(-1.625), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(-1.875), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(-2.125), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(-2.375), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(-2.625), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(-2.875), 0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(-3.125), 0.07);

   hEta_UA5_NSD->SetMarkerStyle(25);
   hEta_UA5_NSD->SetMarkerSize(1.1);
   hEta_UA5_NSD->SetMarkerColor(kGreen-2);
   hEta_UA5_NSD->SetMarkerSize(0.5);
   hEta_UA5_NSD->SetLineColor(kGreen-2);
   hEta_UA5_NSD->SetMarkerStyle(4);

   return hEta_UA5_NSD;
}

TH1F* getUA5NSDScaled() {
   TH1F* hEta_UA5_NSD = new TH1F("hEta_UA5_NSDScale", ";#eta;dN/d#eta", 50, -3, 3);
   // scaleFactor = 1.23625;
   // double scaleFactor = 4.23/3.423;
   double scaleFactor = 5.6/3.423;

   // positive eta
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*0.125), scaleFactor*3.48);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*0.375), scaleFactor*3.38);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*0.625), scaleFactor*3.52);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*0.875), scaleFactor*3.68);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*1.125), scaleFactor*3.71);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*1.375), scaleFactor*3.86);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*1.625), scaleFactor*3.76);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*1.875), scaleFactor*3.66);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*2.125), scaleFactor*3.72);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*2.375), scaleFactor*3.69);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*2.625), scaleFactor*3.56);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*2.875), scaleFactor*3.41);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*3.125), scaleFactor*3.15);

   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*0.125), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*0.375), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*0.625), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*0.875), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*1.125), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*1.375), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*1.625), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*1.875), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*2.125), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*2.375), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*2.625), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*2.875), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*3.125), scaleFactor*0.07);

   //negative eta
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*-0.125), scaleFactor*3.48);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*-0.375), scaleFactor*3.38);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*-0.625), scaleFactor*3.52);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*-0.875), scaleFactor*3.68);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*-1.125), scaleFactor*3.71);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*-1.375), scaleFactor*3.86);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*-1.625), scaleFactor*3.76);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*-1.875), scaleFactor*3.66);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*-2.125), scaleFactor*3.72);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*-2.375), scaleFactor*3.69);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*-2.625), scaleFactor*3.56);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*-2.875), scaleFactor*3.41);
   hEta_UA5_NSD->SetBinContent(hEta_UA5_NSD->FindBin(1.301*-3.125), scaleFactor*3.15);

   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*-0.125), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*-0.375), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*-0.625), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*-0.875), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*-1.125), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*-1.375), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*-1.625), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*-1.875), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*-2.125), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*-2.375), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*-2.625), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*-2.875), scaleFactor*0.07);
   hEta_UA5_NSD->SetBinError(hEta_UA5_NSD->FindBin(1.301*-3.125), scaleFactor*0.07);

   // hEta_UA5_NSD->SetMarkerStyle(25);
   // hEta_UA5_NSD->SetMarkerSize(1.1);
   hEta_UA5_NSD->SetMarkerColor(kBlue-2);
   hEta_UA5_NSD->SetMarkerSize(0.5);
   hEta_UA5_NSD->SetLineColor(kBlue-2);
   hEta_UA5_NSD->SetMarkerStyle(4);

   return hEta_UA5_NSD;
}

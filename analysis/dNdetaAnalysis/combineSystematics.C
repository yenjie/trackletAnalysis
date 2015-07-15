#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

#include "GraphErrorsBand.h"

int combineSystematics(const char* fname = "merged/merged-DATA-C-NEWEST-0BG.root") {
   TFile *infPYTHIA = new TFile("../dNdetaAnalysis/gen/pythia_CUETP8M1.root");
   TH1F* hPYTHIA = (TH1F*)infPYTHIA->FindObjectAny("hEta");
   hPYTHIA->SetName("hPYTHIA");
   TFile *infPYTHIAM = new TFile("../dNdetaAnalysis/gen/pythia_Monash.root");
   TH1F* hPYTHIAM = (TH1F*)infPYTHIAM->FindObjectAny("hEta");
   hPYTHIAM->SetName("hPYTHIAM");
   TFile *infEPOS = new TFile("../dNdetaAnalysis/gen/EPOS.root");
   TH1F* hEPOS = (TH1F*)infEPOS->FindObjectAny("hEta");
   hEPOS->SetName("hEPOS");
   TFile *infQGSJet = new TFile("../dNdetaAnalysis/gen/QGSJet.root");
   TH1F* hQGSJet = (TH1F*)infQGSJet->FindObjectAny("hEta");
   hQGSJet->SetName("hQGSJet");

   hPYTHIA->SetLineColor(4);
   hPYTHIA->SetMarkerColor(4);
   hPYTHIAM->SetLineColor(kGreen+2);
   hPYTHIAM->SetMarkerColor(kGreen+2);
   hQGSJet->SetLineColor(2);
   hQGSJet->SetMarkerColor(2);
   hEPOS->SetLineColor(6);
   hEPOS->SetMarkerColor(6);

   TFile* resultf = new TFile(fname, "read");

   const char* fsys[7] = {"systematics/systematics-DATA-C-WITH-BG.root",
                          "systematics/systematics-DATA-C-STD2-MULT.root",
                          "systematics/systematics-DATA-C-STD2-SB.root",
                          "systematics/systematics-DATA-C-STD2-SPLIT.root",
                          "systematics/systematics-DATA-DROP-C-STD2.root",
                          "systematics/systematics-DATA-SMEARED-C-STD2.root",
                          "systematics/systematics-PYTHIA8-C-EPOS.root"};
   const char* descrip[7] = {"Noise",
                             "Efficiency correction parametrization",
                             "Side-band region definition",
                             "Pixel splitting",
                             "Pixel hit reconstruction efficiency",
                             "Misalignment",
                             "EPOS vs PYTHIA8"};
   const char* namesys[7] = {"systematics-background.pdf",
                             "systematics-multiplicity.pdf",
                             "systematics-sideband.pdf",
                             "systematics-splitting.pdf",
                             "systematics-hitrecoeff.pdf",
                             "systematics-misalignment.pdf",
                             "systematics-eposvspythia.pdf"};

   TFile* fadd[7];
   for (int j=0; j<7; j++)
     fadd[j] = new TFile(fsys[j], "read");

   TFile* ftotal = new TFile("systematics/systematics-total.root", "recreate");
   TH1F* htotal = new TH1F("htotal", "", 30, -3, 3);

   TH1F* hadd[7];
   TF1* haddfit[7];
   TCanvas* canvii[7];

   for (int j=0; j<7; j++) {
      hadd[j] = (TH1F*)fadd[j]->Get("hMeasuredFinal")->Clone();
      hadd[j]->SetBinContent(4, 0);
      hadd[j]->SetBinError(4, 0);
      hadd[j]->SetBinContent(27, 0);
      hadd[j]->SetBinError(27, 0);
      hadd[j]->SetName(Form("hadd%i", j));
      for (int i=1; i<=30; i++) {
         if (hadd[j]->GetBinContent(i) < 1)
            hadd[j]->SetBinContent(i, 2.0-hadd[j]->GetBinContent(i));
      }

      hadd[j]->Fit("pol2", "LL", "", -2.1, 2.1);
      haddfit[j] = hadd[j]->GetFunction("pol2");
      if (j==3) {
         hadd[j]->Fit("pol0", "LL", "", -2, 2);
         haddfit[j] = hadd[j]->GetFunction("pol0");
      }
      haddfit[j]->Write(Form("haddfit%i", j));

      for (int i=1; i<=30; i++) {
         double base = htotal->GetBinContent(i);
         double add = haddfit[j]->Eval(i*0.2-3.1) - 1;
         htotal->SetBinContent(i, sqrt(base*base + add*add));
      }
   }

   // for (int j=0; j<7; j++) {
   //    canvii[j] = new TCanvas(Form("c%i", j), "", 600, 600);
   //    hadd[j]->Draw();
   //    canvii[j]->SaveAs(namesys[j]);
   // }

   TH1F* smooth[7];

   TCanvas* ctotal = new TCanvas("ctotal", "", 600, 600);
   TH1F* frame = ctotal->DrawFrame(-3, 0, 3, 10);
   frame->GetXaxis()->SetTitle("#eta");
   frame->GetXaxis()->CenterTitle(true);
   frame->GetXaxis()->SetLabelFont(42);
   frame->GetXaxis()->SetLabelSize(0.035);
   frame->GetXaxis()->SetTitleSize(0.035);
   frame->GetXaxis()->SetTitleFont(42);
   frame->GetYaxis()->SetTitle("Systematic uncertainty (%)");
   frame->GetYaxis()->CenterTitle(true);
   frame->GetYaxis()->SetLabelFont(42);
   frame->GetYaxis()->SetLabelSize(0.035);
   frame->GetYaxis()->SetTitleSize(0.035);
   frame->GetYaxis()->SetTitleOffset(1.25);
   frame->GetYaxis()->SetTitleFont(42);
   frame->SetTitle("Summary of systematic uncertainties");

   for (int j=0; j<7; j++) {
      smooth[j] = new TH1F(Form("smooth%i", j), "", 30, -3, 3);
      for (int i=1; i<=30; i++) {
         smooth[j]->SetBinContent(i, haddfit[j]->Eval(i*0.2-3.1) - 1);
         smooth[j]->SetBinContent(i, smooth[j]->GetBinContent(i) * 100);
      }
      smooth[j]->SetMarkerColor(j+2);
      smooth[j]->SetLineColor(j+2);

      smooth[j]->Draw("same lf2");
   }

   for (int i=1; i<=30; i++)
     htotal->SetBinContent(i, htotal->GetBinContent(i) * 100);
   htotal->SetLineColor(1);
   htotal->SetLineWidth(2);
   htotal->Draw("same l");

   TLegend* l1 = new TLegend(0.24, 0.5, 0.84, 0.8);
   l1->SetBorderSize(0);
   for (int j=0; j<7; j++)
      l1->AddEntry(smooth[j], descrip[j]);
   l1->AddEntry(htotal, "Total systematics");
   l1->Draw();

   ctotal->SaveAs("results/combined-systematics.png");
   ctotal->SaveAs("results/combined-systematics.pdf");
   ctotal->SaveAs("results/combined-systematics.C");

   Double_t erreta[30];
   for (int i=0; i<30; i++)
      erreta[i] = (htotal->GetBinContent(i+1) + htotal->GetBinContent(30-i)) / 200;

   TCanvas* c3 = new TCanvas("c3", "", 600, 600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);

   TH1F* hSym = (TH1F*)resultf->Get("hSym");
   TH1F* h12 = (TH1F*)resultf->Get("h12");
   TH1F* h13 = (TH1F*)resultf->Get("h13");
   TH1F* h23 = (TH1F*)resultf->Get("h23");

   hSym->Draw("p");
   TGraph* gErrorBand = GetErrorBand(hSym, erreta, erreta, 0.1);
   gErrorBand->Draw("F");

   // Double_t xAxis5[31] = {-3, -2.8, -2.6, -2.4, -2.2, -2, -1.8, -1.6, -1.4, -1.2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3};
   // TH1F *hSym10 = new TH1F("hSym10","",30, xAxis5);
   // hSym10->SetBinContent(5,6.086305);
   // hSym10->SetBinContent(6,6.250434);
   // hSym10->SetBinContent(7,6.158051);
   // hSym10->SetBinContent(8,6.133353);
   // hSym10->SetBinContent(9,6.096857);
   // hSym10->SetBinContent(10,6.014641);
   // hSym10->SetBinContent(11,5.905);
   // hSym10->SetBinContent(12,5.799025);
   // hSym10->SetBinContent(13,5.66402);
   // hSym10->SetBinContent(14,5.560308);
   // hSym10->SetBinContent(15,5.456539);
   // hSym10->SetBinContent(16,5.456539);
   // hSym10->SetBinContent(17,5.560308);
   // hSym10->SetBinContent(18,5.66402);
   // hSym10->SetBinContent(19,5.799025);
   // hSym10->SetBinContent(20,5.905);
   // hSym10->SetBinContent(21,6.014641);
   // hSym10->SetBinContent(22,6.096857);
   // hSym10->SetBinContent(23,6.133353);
   // hSym10->SetBinContent(24,6.158051);
   // hSym10->SetBinContent(25,6.250434);
   // hSym10->SetBinContent(26,6.086305);
   // hSym10->SetBinError(5,0.03495868);
   // hSym10->SetBinError(6,0.01760187);
   // hSym10->SetBinError(7,0.01119936);
   // hSym10->SetBinError(8,0.008740966);
   // hSym10->SetBinError(9,0.008136519);
   // hSym10->SetBinError(10,0.007994242);
   // hSym10->SetBinError(11,0.007870157);
   // hSym10->SetBinError(12,0.007756479);
   // hSym10->SetBinError(13,0.007611605);
   // hSym10->SetBinError(14,0.007501756);
   // hSym10->SetBinError(15,0.007376401);
   // hSym10->SetBinError(16,0.007376401);
   // hSym10->SetBinError(17,0.007501756);
   // hSym10->SetBinError(18,0.007611605);
   // hSym10->SetBinError(19,0.007756479);
   // hSym10->SetBinError(20,0.007870157);
   // hSym10->SetBinError(21,0.007994242);
   // hSym10->SetBinError(22,0.008136519);
   // hSym10->SetBinError(23,0.008740966);
   // hSym10->SetBinError(24,0.01119936);
   // hSym10->SetBinError(25,0.01760187);
   // hSym10->SetBinError(26,0.03495868);
   // hSym10->SetMinimum(0);
   // hSym10->SetMaximum(8);
   // hSym10->SetEntries(1323801);
   // hSym10->SetStats(0);
   // hSym10->SetLineColor(4);
   // hSym10->SetLineStyle(0);
   // hSym10->SetMarkerColor(4);
   // hSym10->SetMarkerStyle(20);
   // hSym10->GetXaxis()->SetTitle("#eta");
   // hSym10->GetXaxis()->CenterTitle(true);
   // hSym10->GetXaxis()->SetLabelFont(42);
   // hSym10->GetXaxis()->SetLabelSize(0.035);
   // hSym10->GetXaxis()->SetTitleSize(0.035);
   // hSym10->GetXaxis()->SetTitleFont(42);
   // hSym10->GetYaxis()->SetTitle("dN/d#eta");
   // hSym10->GetYaxis()->CenterTitle(true);
   // hSym10->GetYaxis()->SetLabelFont(42);
   // hSym10->GetYaxis()->SetLabelSize(0.035);
   // hSym10->GetYaxis()->SetTitleSize(0.035);
   // hSym10->GetYaxis()->SetTitleOffset(1.25);
   // hSym10->GetYaxis()->SetTitleFont(42);
   // hSym10->GetZaxis()->SetLabelFont(42);
   // hSym10->GetZaxis()->SetLabelSize(0.035);
   // hSym10->GetZaxis()->SetTitleSize(0.035);
   // hSym10->GetZaxis()->SetTitleFont(42);
   // hSym10->Draw("p");

   // TGraph* gErrorBand = GetErrorBand(hSym10, erreta, erreta, 0.1);
   // gErrorBand->Draw("F");

   // TH1F *PureCutCorrEta = new TH1F("PureCutCorrEta","PureCutCorrEta",30,-3,3);
   // PureCutCorrEta->SetBinContent(3,2.671707);
   // PureCutCorrEta->SetBinContent(4,5.506016);
   // PureCutCorrEta->SetBinContent(5,5.494589);
   // PureCutCorrEta->SetBinContent(6,5.808212);
   // PureCutCorrEta->SetBinContent(7,5.60169);
   // PureCutCorrEta->SetBinContent(8,5.742654);
   // PureCutCorrEta->SetBinContent(9,5.689229);
   // PureCutCorrEta->SetBinContent(10,5.506961);
   // PureCutCorrEta->SetBinContent(11,5.819585);
   // PureCutCorrEta->SetBinContent(12,5.422974);
   // PureCutCorrEta->SetBinContent(13,5.371057);
   // PureCutCorrEta->SetBinContent(14,5.343599);
   // PureCutCorrEta->SetBinContent(15,5.375769);
   // PureCutCorrEta->SetBinContent(16,4.743305);
   // PureCutCorrEta->SetBinContent(17,5.082912);
   // PureCutCorrEta->SetBinContent(18,5.209239);
   // PureCutCorrEta->SetBinContent(19,5.09033);
   // PureCutCorrEta->SetBinContent(20,5.296446);
   // PureCutCorrEta->SetBinContent(21,5.631082);
   // PureCutCorrEta->SetBinContent(22,5.546448);
   // PureCutCorrEta->SetBinContent(23,5.477513);
   // PureCutCorrEta->SetBinContent(24,5.607662);
   // PureCutCorrEta->SetBinContent(25,5.617325);
   // PureCutCorrEta->SetBinContent(26,5.600933);
   // PureCutCorrEta->SetBinContent(27,5.370979);
   // PureCutCorrEta->SetBinContent(28,2.6269);
   // PureCutCorrEta->SetEntries(136.2551);
   // PureCutCorrEta->Draw("same");

   hPYTHIA->Draw("hist c same");
   hPYTHIAM->Draw("hist c same");
   hQGSJet->Draw("hist c same");
   hEPOS->Draw("hist c same");
   // hSym10->Draw("p same");
   hSym->Draw("p same");

   TLegend* leg3 = new TLegend(0.3, 0.24, 0.84, 0.48);
   leg3->SetBorderSize(0);
   leg3->AddEntry("NULL", "Run 247324 PromptReco", "");
   leg3->AddEntry("hPYTHIA", "PYTHIA8 CUETP8M1", "l");
   leg3->AddEntry("hPYTHIAM", "PYTHIA8 Monash13", "l");
   leg3->AddEntry("hEPOS", "EPOS LHC", "l");
   leg3->AddEntry("hQGSJet", "QGSJet-II", "l");
   // leg3->AddEntry(hSym10, "Reconstructed Tracklets", "pl");
   leg3->AddEntry(hSym, "Reconstructed Tracklets", "pl");
   // leg3->AddEntry(PureCutCorrEta, "Pixel Cluster Analysis", "l");
   leg3->Draw();

   c3->SaveAs("results/sym-results.png");
   c3->SaveAs("results/sym-results.pdf");
   c3->SaveAs("results/sym-results.root");
   c3->SaveAs("results/sym-results.C");

   TCanvas *c1 = new TCanvas("c1", "", 600, 600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);

   Double_t xAxis1[31] = {-3, -2.8, -2.6, -2.4, -2.2, -2, -1.8, -1.6, -1.4, -1.2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3};
   TH1D *hTruthWOSelection = new TH1D("hTruthWOSelection","",30, xAxis1);
   hTruthWOSelection->SetBinContent(1,5.623362);
   hTruthWOSelection->SetBinContent(2,5.709071);
   hTruthWOSelection->SetBinContent(3,5.763876);
   hTruthWOSelection->SetBinContent(4,5.793779);
   hTruthWOSelection->SetBinContent(5,5.817132);
   hTruthWOSelection->SetBinContent(6,5.850085);
   hTruthWOSelection->SetBinContent(7,5.816032);
   hTruthWOSelection->SetBinContent(8,5.776728);
   hTruthWOSelection->SetBinContent(9,5.774527);
   hTruthWOSelection->SetBinContent(10,5.70172);
   hTruthWOSelection->SetBinContent(11,5.59891);
   hTruthWOSelection->SetBinContent(12,5.505151);
   hTruthWOSelection->SetBinContent(13,5.394889);
   hTruthWOSelection->SetBinContent(14,5.334683);
   hTruthWOSelection->SetBinContent(15,5.329433);
   hTruthWOSelection->SetBinContent(16,5.347535);
   hTruthWOSelection->SetBinContent(17,5.340134);
   hTruthWOSelection->SetBinContent(18,5.40414);
   hTruthWOSelection->SetBinContent(19,5.492799);
   hTruthWOSelection->SetBinContent(20,5.538404);
   hTruthWOSelection->SetBinContent(21,5.663316);
   hTruthWOSelection->SetBinContent(22,5.732423);
   hTruthWOSelection->SetBinContent(23,5.794329);
   hTruthWOSelection->SetBinContent(24,5.846335);
   hTruthWOSelection->SetBinContent(25,5.820932);
   hTruthWOSelection->SetBinContent(26,5.874587);
   hTruthWOSelection->SetBinContent(27,5.828433);
   hTruthWOSelection->SetBinContent(28,5.759976);
   hTruthWOSelection->SetBinContent(29,5.693469);
   hTruthWOSelection->SetBinContent(30,5.616712);
   hTruthWOSelection->SetBinError(1,0.01676891);
   hTruthWOSelection->SetBinError(2,0.01689622);
   hTruthWOSelection->SetBinError(3,0.01697712);
   hTruthWOSelection->SetBinError(4,0.0170211);
   hTruthWOSelection->SetBinError(5,0.01705537);
   hTruthWOSelection->SetBinError(6,0.01710361);
   hTruthWOSelection->SetBinError(7,0.01705376);
   hTruthWOSelection->SetBinError(8,0.01699604);
   hTruthWOSelection->SetBinError(9,0.0169928);
   hTruthWOSelection->SetBinError(10,0.01688533);
   hTruthWOSelection->SetBinError(11,0.01673241);
   hTruthWOSelection->SetBinError(12,0.01659172);
   hTruthWOSelection->SetBinError(13,0.01642472);
   hTruthWOSelection->SetBinError(14,0.01633282);
   hTruthWOSelection->SetBinError(15,0.01632478);
   hTruthWOSelection->SetBinError(16,0.01635248);
   hTruthWOSelection->SetBinError(17,0.01634116);
   hTruthWOSelection->SetBinError(18,0.0164388);
   hTruthWOSelection->SetBinError(19,0.01657309);
   hTruthWOSelection->SetBinError(20,0.01664175);
   hTruthWOSelection->SetBinError(21,0.01682837);
   hTruthWOSelection->SetBinError(22,0.01693074);
   hTruthWOSelection->SetBinError(23,0.01702191);
   hTruthWOSelection->SetBinError(24,0.01709813);
   hTruthWOSelection->SetBinError(25,0.01706094);
   hTruthWOSelection->SetBinError(26,0.01713939);
   hTruthWOSelection->SetBinError(27,0.01707193);
   hTruthWOSelection->SetBinError(28,0.01697138);
   hTruthWOSelection->SetBinError(29,0.01687311);
   hTruthWOSelection->SetBinError(30,0.01675899);
   hTruthWOSelection->SetMinimum(0);
   hTruthWOSelection->SetMaximum(9);
   hTruthWOSelection->SetEntries(3390519);
   hTruthWOSelection->SetStats(0);
   hTruthWOSelection->GetXaxis()->SetTitle("#eta");
   hTruthWOSelection->GetXaxis()->CenterTitle(true);
   hTruthWOSelection->GetXaxis()->SetLabelFont(42);
   hTruthWOSelection->GetXaxis()->SetLabelSize(0.035);
   hTruthWOSelection->GetXaxis()->SetTitleSize(0.035);
   hTruthWOSelection->GetXaxis()->SetTitleFont(42);
   hTruthWOSelection->GetYaxis()->SetTitle("dN/d#eta");
   hTruthWOSelection->GetYaxis()->CenterTitle(true);
   hTruthWOSelection->GetYaxis()->SetLabelFont(42);
   hTruthWOSelection->GetYaxis()->SetLabelSize(0.035);
   hTruthWOSelection->GetYaxis()->SetTitleSize(0.035);
   hTruthWOSelection->GetYaxis()->SetTitleOffset(1.25);
   hTruthWOSelection->GetYaxis()->SetTitleFont(42);
   hTruthWOSelection->GetZaxis()->SetLabelFont(42);
   hTruthWOSelection->GetZaxis()->SetLabelSize(0.035);
   hTruthWOSelection->GetZaxis()->SetTitleSize(0.035);
   hTruthWOSelection->GetZaxis()->SetTitleFont(42);
   hTruthWOSelection->Draw("");

   for (int i=0; i<30; i++)
      h12->SetBinError(i+1, h12->GetBinContent(i+1)*erreta[i]);
   for (int i=0; i<30; i++)
      h13->SetBinError(i+1, h13->GetBinContent(i+1)*erreta[i]);
   for (int i=0; i<30; i++)
      h23->SetBinError(i+1, h23->GetBinContent(i+1)*erreta[i]);
   h12->Draw("same");
   h13->Draw("same");
   h23->Draw("same");

   // Double_t xAxis2[31] = {-3, -2.8, -2.6, -2.4, -2.2, -2, -1.8, -1.6, -1.4, -1.2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3};
   // TH1D *h121 = new TH1D("h121","",30, xAxis2);
   // h121->SetBinContent(5,6.002668);
   // h121->SetBinContent(6,6.186875);
   // h121->SetBinContent(7,6.142438);
   // h121->SetBinContent(8,6.15766);
   // h121->SetBinContent(9,6.131638);
   // h121->SetBinContent(10,6.040274);
   // h121->SetBinContent(11,5.9062);
   // h121->SetBinContent(12,5.760449);
   // h121->SetBinContent(13,5.626106);
   // h121->SetBinContent(14,5.531895);
   // h121->SetBinContent(15,5.423155);
   // h121->SetBinContent(16,5.379206);
   // h121->SetBinContent(17,5.473795);
   // h121->SetBinContent(18,5.571939);
   // h121->SetBinContent(19,5.727229);
   // h121->SetBinContent(20,5.806399);
   // h121->SetBinContent(21,5.907567);
   // h121->SetBinContent(22,6.02488);
   // h121->SetBinContent(23,6.15553);
   // h121->SetBinContent(24,6.26597);
   // h121->SetBinContent(25,6.313993);
   // h121->SetBinContent(26,6.169942);
   // h121->SetBinError(5,0.06147559);
   // h121->SetBinError(6,0.02709557);
   // h121->SetBinError(7,0.02079968);
   // h121->SetBinError(8,0.02000931);
   // h121->SetBinError(9,0.01986279);
   // h121->SetBinError(10,0.01958364);
   // h121->SetBinError(11,0.0192278);
   // h121->SetBinError(12,0.01886721);
   // h121->SetBinError(13,0.01852957);
   // h121->SetBinError(14,0.01830285);
   // h121->SetBinError(15,0.01797491);
   // h121->SetBinError(16,0.01785103);
   // h121->SetBinError(17,0.01811963);
   // h121->SetBinError(18,0.01837166);
   // h121->SetBinError(19,0.01878757);
   // h121->SetBinError(20,0.01899623);
   // h121->SetBinError(21,0.01928077);
   // h121->SetBinError(22,0.01961992);
   // h121->SetBinError(23,0.02006417);
   // h121->SetBinError(24,0.02160401);
   // h121->SetBinError(25,0.02247519);
   // h121->SetBinError(26,0.03330451);
   // h121->SetMinimum(0);
   // h121->SetMaximum(8);
   // h121->SetEntries(1323771);
   // h121->SetStats(0);
   // h121->SetLineColor(2);
   // h121->SetLineStyle(0);
   // h121->SetMarkerColor(2);
   // h121->SetMarkerStyle(20);
   // h121->GetXaxis()->SetTitle("#eta");
   // h121->GetXaxis()->CenterTitle(true);
   // h121->GetXaxis()->SetLabelFont(42);
   // h121->GetXaxis()->SetLabelOffset(0.01);
   // h121->GetXaxis()->SetLabelSize(0.045);
   // h121->GetXaxis()->SetTitleSize(0.055);
   // h121->GetXaxis()->SetTitleFont(42);
   // h121->GetYaxis()->SetTitle("dN/d#eta");
   // h121->GetYaxis()->CenterTitle(true);
   // h121->GetYaxis()->SetLabelFont(42);
   // h121->GetYaxis()->SetLabelOffset(0.01);
   // h121->GetYaxis()->SetLabelSize(0.045);
   // h121->GetYaxis()->SetTitleSize(0.055);
   // h121->GetYaxis()->SetTitleOffset(1.25);
   // h121->GetYaxis()->SetTitleFont(42);
   // h121->GetZaxis()->SetLabelFont(42);
   // h121->GetZaxis()->SetLabelSize(0.045);
   // h121->GetZaxis()->SetTitleSize(0.035);
   // h121->GetZaxis()->SetTitleFont(42);
   // h121->Draw("same");

   // Double_t xAxis3[31] = {-3, -2.8, -2.6, -2.4, -2.2, -2, -1.8, -1.6, -1.4, -1.2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3};
   // TH1D *h132 = new TH1D("h132","",30, xAxis3);
   // h132->SetBinContent(7,6.072267);
   // h132->SetBinContent(8,6.085483);
   // h132->SetBinContent(9,6.042458);
   // h132->SetBinContent(10,5.988057);
   // h132->SetBinContent(11,5.851167);
   // h132->SetBinContent(12,5.75098);
   // h132->SetBinContent(13,5.652775);
   // h132->SetBinContent(14,5.567731);
   // h132->SetBinContent(15,5.479251);
   // h132->SetBinContent(16,5.406099);
   // h132->SetBinContent(17,5.492217);
   // h132->SetBinContent(18,5.61308);
   // h132->SetBinContent(19,5.750385);
   // h132->SetBinContent(20,5.853657);
   // h132->SetBinContent(21,5.930417);
   // h132->SetBinContent(22,5.974352);
   // h132->SetBinContent(23,5.979696);
   // h132->SetBinContent(24,5.918695);
   // h132->SetBinError(7,0.03511494);
   // h132->SetBinError(8,0.024653);
   // h132->SetBinError(9,0.0203168);
   // h132->SetBinError(10,0.01978024);
   // h132->SetBinError(11,0.01935581);
   // h132->SetBinError(12,0.0190935);
   // h132->SetBinError(13,0.01882554);
   // h132->SetBinError(14,0.0186145);
   // h132->SetBinError(15,0.01836858);
   // h132->SetBinError(16,0.01816108);
   // h132->SetBinError(17,0.01842205);
   // h132->SetBinError(18,0.01877707);
   // h132->SetBinError(19,0.01916837);
   // h132->SetBinError(20,0.0194949);
   // h132->SetBinError(21,0.01972415);
   // h132->SetBinError(22,0.01987194);
   // h132->SetBinError(23,0.02062681);
   // h132->SetBinError(24,0.02445253);
   // h132->SetEntries(1367837);
   // h132->SetStats(0);
   // h132->SetLineStyle(0);
   // h132->SetMarkerStyle(26);
   // h132->Draw("same");
   // Double_t xAxis4[31] = {-3, -2.8, -2.6, -2.4, -2.2, -2, -1.8, -1.6, -1.4, -1.2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3};

   // TH1D *h233 = new TH1D("h233","",30, xAxis4);
   // h233->SetBinContent(7,6.253702);
   // h233->SetBinContent(8,6.160283);
   // h233->SetBinContent(9,6.167112);
   // h233->SetBinContent(10,6.093606);
   // h233->SetBinContent(11,6.003658);
   // h233->SetBinContent(12,5.84794);
   // h233->SetBinContent(13,5.750849);
   // h233->SetBinContent(14,5.648552);
   // h233->SetBinContent(15,5.536242);
   // h233->SetBinContent(16,5.515283);
   // h233->SetBinContent(17,5.647658);
   // h233->SetBinContent(18,5.769369);
   // h233->SetBinContent(19,5.957166);
   // h233->SetBinContent(20,6.00892);
   // h233->SetBinContent(21,6.127924);
   // h233->SetBinContent(22,6.240699);
   // h233->SetBinContent(23,6.261465);
   // h233->SetBinContent(24,6.295232);
   // h233->SetBinError(7,0.03438706);
   // h233->SetBinError(8,0.02212884);
   // h233->SetBinError(9,0.01996381);
   // h233->SetBinError(10,0.01947797);
   // h233->SetBinError(11,0.01926018);
   // h233->SetBinError(12,0.01884819);
   // h233->SetBinError(13,0.01861958);
   // h233->SetBinError(14,0.01838648);
   // h233->SetBinError(15,0.0180428);
   // h233->SetBinError(16,0.01800771);
   // h233->SetBinError(17,0.01840378);
   // h233->SetBinError(18,0.01873991);
   // h233->SetBinError(19,0.01922706);
   // h233->SetBinError(20,0.01932872);
   // h233->SetBinError(21,0.01963993);
   // h233->SetBinError(22,0.01994026);
   // h233->SetBinError(23,0.02060818);
   // h233->SetBinError(24,0.02454597);
   // h233->SetEntries(1488622);
   // h233->SetStats(0);
   // h233->SetLineColor(4);
   // h233->SetLineStyle(0);
   // h233->SetMarkerColor(4);
   // h233->SetMarkerStyle(25);
   // h233->Draw("same");

   TLegend* leg1 = new TLegend(0.24, 0.3, 0.84, 0.45);
   leg1->SetBorderSize(0);
   leg1->AddEntry("NULL", "Run 247324 PromptReco", "");
   // leg1->AddEntry("h121", "Reconstructed (1st+2nd layers)", "pl");
   // leg1->AddEntry("h132", "Reconstructed (1st+3rd layers)", "pl");
   // leg1->AddEntry("h233", "Reconstructed (2nd+3rd layers)", "pl");
   leg1->AddEntry("h12", "Reconstructed (1st+2nd layers)", "pl");
   leg1->AddEntry("h13", "Reconstructed (1st+3rd layers)", "pl");
   leg1->AddEntry("h23", "Reconstructed (2nd+3rd layers)", "pl");
   leg1->Draw();

   c1->SaveAs("results/merged-results.png");
   c1->SaveAs("results/merged-results.pdf");
   c1->SaveAs("results/merged-results.root");
   c1->SaveAs("results/merged-results.C");

   ftotal->Write();

   return 0;
}

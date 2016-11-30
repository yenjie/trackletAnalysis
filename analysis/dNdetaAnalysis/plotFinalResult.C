// Plot final results
#define canvasSizeX 600
#define canvasSizeY 600
#define dndetaRange 30.0
#define SDFactor 1

// Standard library
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

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
#include "TChain.h"
#include <TLine.h>
#include <TF1.h>
#include <TCut.h>
#include <TPad.h>
#include <TText.h>

// For plotting
#include "GraphErrorsBand.h"

void formatHist(TH1* h, int col = 1, double norm = 1, double msize = 1);

// ============================================================================
// Main Routine
// ============================================================================
int plotFinalResult(int TrackletType,
                    std::string file_name,
                    const char* title,       // plot title
                    bool useCorrectionFile,  // use correction file
                    const char* corr_name,   // correction file name
                    int selection = 1,       // MC selection
                    int mult_selection = 0,  // multiplicity
                    bool doAccepCorr = 0,    // do acceptance correction
                    bool doTriggerCorr = 1,  // do trigger eff correction
                    bool verbose = 0,        // set verbose level
                    bool plotAlphaBeta = 0,  // make alpha/beta plots
                    bool useExtBetaCorr = 0, // do beta correction
                    bool useDR = 0)
{
   std::vector<std::string> file_list;
   if (file_name.substr(file_name.find_last_of(".") + 1) == "root") {
      file_list.push_back(file_name);
   } else {
      ifstream file_stream(file_name);
      if (!file_stream) return 1;
      std::string line;
      while (std::getline(file_stream, line))
         file_list.push_back(line);
   }

   // Input trackletTree
   TChain* TrackletTree = new TChain(Form("TrackletTree%d", TrackletType));
   for (std::size_t n=0; n<file_list.size(); ++n)
      TrackletTree->Add(file_list[n].c_str());

   bool isMC = false;
   if (TrackletTree->GetEntries("npart!=0")!=0) {
      isMC = true;
      cout << "This is a Monte Carlo study." << endl;
      doAccepCorr = 0;
   } else {
      cout << "This is a data analysis." << endl;
   }

   if (!doTriggerCorr)
      cout << "Trigger correction off!!!" << endl;

   // Choose multiplicity handle
   const char* multiplicity;
   if (mult_selection == 1) {
      multiplicity = "mult2";
      cout << "Use # of clusters as event multiplicity" << endl;
   } else if (mult_selection == 2) {
      multiplicity = "nTracklets";
      cout << "Use # of tracklets as event multiplicity" << endl;
   } else {
      multiplicity = "mult";
      cout << "Use # of tracklets after background subtraction as event multiplicity" << endl;
   }

   TH1::SetDefaultSumw2();

   // Read alpha, beta, geometry correction from file.
   TFile* fCorrection = 0;
   if (useCorrectionFile) {
      const char* correctionfname = Form("correction/correction-%d-%s.root", TrackletType, corr_name);
      fCorrection = new TFile(correctionfname);
      printf("Use correction file: %s\n", correctionfname);
   }

   TFile* fAcceptance = 0;
   if (useCorrectionFile && doAccepCorr) {
      const char* acceptancefname = Form("correction/acceptance-%d.root", TrackletType);
      fAcceptance = new TFile(acceptancefname);
      printf("Use acceptance file: %s\n", acceptancefname);
   }

   TH3F* hAlphaA = 0;
   TH3F* hAlphaB = 0;
   if (useExtBetaCorr) {
      TFile* myFile = new TFile(Form("correction/alphaBetaCoeff-%d.root", TrackletType));
      hAlphaA = (TH3F*)myFile->FindObjectAny("hAlphaA");
      hAlphaB = (TH3F*)myFile->FindObjectAny("hAlphaB");
   }

   int VzRangeL = -15;
   int VzRangeH = 15;

   // Definition of Vz, Eta, Hit bins
   const int nTrackletBin = 12;
   const int nEtaBin = 30;
   const int nVzBin = 15;

   double TrackletBins[nTrackletBin+1] = {-5, 2, 10, 15, 20, 25, 30, 36, 42, 50, 60, 72, 300};
   double EtaBins[nEtaBin+1];
   for (int i=0; i<=nEtaBin; i++)
      EtaBins[i] = (double)i*6.0/(double)nEtaBin-3.0;
   double VzBins[nVzBin+1];
   for (int i=0; i<=nVzBin; i++)
      VzBins[i] = (double)i*(VzRangeH-VzRangeL)/(double)nVzBin+VzRangeL;

   // Signal and Sideband regions =============================================
   double signalRegionCut = 1.0;   // dphi cut for signal region
   double sidebandRegionCut = 2.0; // dphi cut for side-band region
   double detaCut = 0.1;           // deta cut

   TCut signalRegion   = Form("abs(dphi)<%f&&abs(deta)<%f", signalRegionCut, detaCut);
   TCut sidebandRegion = Form("abs(dphi)>%f&&abs(dphi)<%f&&abs(deta)<%f", signalRegionCut, sidebandRegionCut, detaCut);
   if (useDR) {
      cout << "Use dR for analysis" << endl;
      signalRegion   = "dR<0.1";
      sidebandRegion = "dR>0.1&&dR<0.2";
   }

   TString vtxCut = "(vz[1]<15 && vz[1]>-15)";
   TCut MCSelection;
   TString offlineSelection;
   TCut evtSelection;

   switch (selection) {
      case 0:
         MCSelection = "(evtType!=102)";
         offlineSelection = "1";
         printf("---------- INELASTIC definition\n");
         break;
      case 1:
         MCSelection = "(evtType!=102&&evtType!=103&&evtType!=104)";
         offlineSelection = "1"; // effectively HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1
         printf("---------------- NSD definition\n");
         break;
   }
   if (!isMC) MCSelection = "1";
   evtSelection = TCut(vtxCut + "&&" + offlineSelection);

   // Output file =============================================================
   TFile* outf = new TFile(Form("correction-%i-%s.root", TrackletType, title), "recreate");

   TNtuple* betas = new TNtuple("betas", "", "eta:nTracklet:vz:beta:betaErr");
   TNtuple* alphas = new TNtuple("alphas", "", "eta:nTracklet:vz:alpha:alphaErr");
   TNtuple* correction = new TNtuple("correction", "", "eta:nTracklet:vz:alpha:alphaErr:beta:betaErr:obs:gen:err:betaMC");

   TH3F* hEverything = new TH3F("hEverything", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* hReproducedBackground = new TH3F("hReproducedBackground", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* hSubtracted = new TH3F("hSubtracted", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* hHadron = new TH3F("hHadron", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* hHadronAccepted = new TH3F("hHadronAccepted", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* hHadronWOSelection = new TH3F("hHadronWOSelection", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* hCorrected = new TH3F("hCorrected", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* hAcceptance2D = new TH3F("hAcceptance2D", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);

   TH2F* hEtaVzRatio = new TH2F("hEtaVzatio", "#eta:Vz", nEtaBin, EtaBins, nVzBin, VzBins);
   TH2F* hEtaHitRatio = new TH2F("hEtaHitatio", "#eta:N_{Hit}", nEtaBin, EtaBins, nTrackletBin, TrackletBins);
   TH2F* hAcceptance1D = new TH2F("hAcceptance1D", "", nEtaBin, EtaBins, nVzBin, VzBins);

   // Acceptance
   TH2F* hDataAcc = new TH2F("hDataAcc", "", nEtaBin, EtaBins, nVzBin, VzBins);
   TH2F* hMCAcc = new TH2F("hMCAcc", "", nEtaBin, EtaBins, nVzBin, VzBins);

   TH2F* hVzNTracklet = new TH2F("hVzNTracklet", "", nTrackletBin, TrackletBins, nVzBin, VzBins);

   TH1F* hCorrectedEtaBin = new TH1F("hCorrectedEtaBin", "Corrected", nEtaBin, -3, 3);
   TH1F* hDNDEtaVertexed = new TH1F("hDNDEtaVertexed", "", nEtaBin, -3, 3);
   TH1F* hDNDEtaNoVertexed = new TH1F("hDNDEtaNoVertexed", "", nEtaBin, -3, 3);
   TH1F* hTrigEff = new TH1F("hTrigEff", "", nTrackletBin, TrackletBins);
   TH1F* hTrigEffNoCut = new TH1F("hTrigEffNoCut", "", nTrackletBin, TrackletBins);
   TH1F* hSD = new TH1F("hSD", "", nTrackletBin, TrackletBins);
   TH1F* hSDFrac = new TH1F("hSDFrac", "", nTrackletBin, TrackletBins);
   TH1F* hVz = new TH1F("hVz", "", nVzBin, VzBins);
   TH1F* hnTracklet = new TH1F("hnTracklet", "", nTrackletBin, TrackletBins);

   TH1F* alphaPlots[nEtaBin][nVzBin];
   TH1F* betaPlots[nEtaBin][nVzBin];
   TH1F* betaMCPlots[nEtaBin][nVzBin];
   TH1F* alphaErrPlots[nEtaBin][nVzBin];
   TH1F* betaErrPlots[nEtaBin][nVzBin];

   // Prepare histograms ======================================================
   for (int i=0; i<nEtaBin; i++) {
      for (int j=0; j<nVzBin; j++) {
         alphaPlots[i][j] = new TH1F(Form("alpha%dVz%d", i, j), "", nTrackletBin, TrackletBins);
         alphaErrPlots[i][j] = new TH1F(Form("alphaErr%dVz%d", i, j), "", nTrackletBin, TrackletBins);
         betaPlots[i][j] = new TH1F(Form("beta%dVz%d", i, j), "", nTrackletBin, TrackletBins);
         betaErrPlots[i][j] = new TH1F(Form("betaErr%dVz%d", i, j), "", nTrackletBin, TrackletBins);
         betaMCPlots[i][j] = new TH1F(Form("betaMC%dVz%d", i, j), "", nTrackletBin, TrackletBins);
      }
   }

   // Fit functions of Beta and Alpha =========================================
   TF1* funAlpha[nEtaBin][nVzBin];
   TF1* funAlphaErr[nEtaBin][nVzBin];
   TF1* fAlpha[nEtaBin][nVzBin];
   TF1* fAlphaErr[nEtaBin][nVzBin];

   TH1F* hTriggerCorrection;
   TH1F* hEmptyEvtCorrection;

   TF1* fitEmptyEvt;
   TF1* fEmptyEvt;

   // Number of events
   int nevent = TrackletTree->Draw("vz[1]", evtSelection && MCSelection, "goff");
   int neventWOSelection = TrackletTree->Draw("vz[1]", MCSelection, "goff");
   cout << "Number of events: " << nevent << endl;
   if (nevent < 1) {
      cout << "No event survived. Abort." << endl;
      return 0;
   }

   // Acceptance calculation ==================================================
   TrackletTree->Project("hnTracklet", Form("%s", multiplicity), evtSelection);
   TrackletTree->Project("hVzNTracklet", Form("vz[1]:%s", multiplicity), evtSelection);
   TrackletTree->Project("hVz", "vz[1]", evtSelection);

   TCanvas* cVz = new TCanvas("cVz", "Vz distribution", canvasSizeX, canvasSizeY);
   hVz->Sumw2();
   hVz->Scale(1./hVz->GetEntries());
   hVz->Fit("gaus");
   hVz->SetXTitle("v_{z} (cm)");
   hVz->Draw();
   // cVz->SaveAs(Form("figs/vz/vz-%s-%d.pdf", title, TrackletType));

   // Define the acceptance region to avoid large correction factors
   // End point in z (cm)
   bool accepRegion[nEtaBin][nVzBin];
   memset(accepRegion, 0, sizeof(bool)*nEtaBin*nVzBin);

   double endpoint2 = 30.0; // 26.66 (old)
   double rho = 7.6; // Second layer rho
   double etaLimit = 2.3;
   if (TrackletType % 10 == 3) {
      rho = 10.5; // Third layer rho
      etaLimit = 1.9;
   }
   for (int i=0; i<nEtaBin; i++) {
      for (int j=0; j<nVzBin; j++) {
         double minEta = EtaBins[i];
         double maxEta = EtaBins[i+1];
         double maxEdge = VzBins[j+1]-rho/tan(atan(exp(maxEta-0.2))*2);
         double minEdge = VzBins[j]-rho/tan(atan(exp(minEta+0.2))*2);
         if (maxEdge>-endpoint2 && minEdge<endpoint2 && maxEta<etaLimit && minEta>-etaLimit)
            accepRegion[i][j] = 1;
      }
   }

   for (int i=0; i<nEtaBin; i++) {
      for (int j=0; j<nVzBin; j++) {
         if (accepRegion[i][j]) {
            if (verbose) cout << " Selected! " << endl;
            hAcceptance1D->SetBinContent(i+1, j+1, hVz->GetBinContent(j+1));
            hAcceptance1D->SetBinError(i+1, j+1, 0);
            for (int k=0; k<nTrackletBin; k++)
               hAcceptance2D->SetBinContent(i+1, k+1, j+1, hVzNTracklet->GetBinContent(k+1, j+1));
         } else {
            hAcceptance1D->SetBinContent(i+1, j+1, 0);
            for (int k=0; k<nTrackletBin; k++)
               hAcceptance2D->SetBinContent(i+1, k+1, j+1, 0);
         }
      }
   }

   TCanvas* cAcc = new TCanvas("cAcc", "Acceptance", canvasSizeX, canvasSizeY);
   hAcceptance1D->ProjectionX()->Draw();

   // Charged Hadrons =========================================================
   hHadron->SetXTitle("#eta");
   hHadron->SetYTitle("N_{hit}^{Layer1} |#eta|<3");

   TrackletTree->Project("hHadron", Form("vz[1]:%s:eta", multiplicity), "abs(eta)<3" && evtSelection);
   TrackletTree->Project("hHadronWOSelection", Form("vz[1]:%s:eta", multiplicity), "abs(eta)<3" && MCSelection);
   hHadron->Sumw2();
   hHadronWOSelection->Sumw2();

   TrackletTree->Project("hHadronAccepted", Form("vz[1]:%s:eta", multiplicity), "abs(eta)<3" && evtSelection);
   hHadronAccepted = (TH3F*)hHadron->Clone();
   hHadronAccepted->SetName("hHadronAccepted");

   // Prepare Tracklet Three-Dimensional Histogram ============================
   // Signal region && evtSelection
   TrackletTree->Project("hEverything", Form("vz[1]:%s:eta1", multiplicity), signalRegion && evtSelection);
   hEverything->Sumw2();

   // Side-band region && evtSelection
   TrackletTree->Project("hReproducedBackground", Form("vz[1]:%s:eta1", multiplicity), sidebandRegion && evtSelection);
   hReproducedBackground->Sumw2();

   // Read Acceptance =========================================================
   if (doAccepCorr) {
      hMCAcc = (TH2F*)fAcceptance->FindObjectAny("hMCAcc");
      hDataAcc = (TH2F*)fAcceptance->FindObjectAny("hDataAcc");
   }

   TCanvas* c1 = new TCanvas("c1", "scratch", 400, 400);

   // Beta calculation (background fraction) ==================================
   for (int x=1; x<=nEtaBin; x++) {
      for (int y=1; y<=nTrackletBin; y++) {
         for (int z=1; z<=nVzBin; z++) {
            double beta = 0;
            betaPlots[x-1][z-1]->SetBinContent(y, 0);
            betaPlots[x-1][z-1]->SetBinError(y, 0);
            betaErrPlots[x-1][z-1]->SetBinContent(y, 0);

            if (hAcceptance1D->GetBinContent(x, z) == 0) continue;
            if (hEverything->GetBinContent(x, y, z) != 0) {
               beta = hReproducedBackground->GetBinContent(x, y, z)/hEverything->GetBinContent(x, y, z);
               double e1 = hEverything->GetBinError(x, y, z)/hEverything->GetBinContent(x, y, z);
               double e2 = hReproducedBackground->GetBinError(x, y, z)/hReproducedBackground->GetBinContent(x, y, z);
               double betaErr = beta*sqrt(e2*e2);
               if (beta/betaErr>-10) {
                  betas->Fill((EtaBins[x]+EtaBins[x-1])/2, (TrackletBins[y]+TrackletBins[y-1])/2, (VzBins[z]+VzBins[z-1])/2, beta, betaErr);
                  betaPlots[x-1][z-1]->SetBinContent(y, beta);
                  betaPlots[x-1][z-1]->SetBinError(y, betaErr);
                  betaErrPlots[x-1][z-1]->SetBinContent(y, betaErr);
               }
            }
         }
      }
   }

   for (int j=0; j<nVzBin; j++) {
      for (int i=0; i<nEtaBin; i++) {
         // c[i]= new TCanvas (Form("c%d", i), "", canvasSizeX, canvasSizeY);
         // p1->cd(i+1);
         formatHist(betaPlots[i][j], 2, 1);
         double etaMin = i*6.0/nEtaBin-3;
         double etaMax = (i+1)*6.0/nEtaBin-3;
         betaPlots[i][j]->SetXTitle("N_{Hits}");
         betaPlots[i][j]->SetYTitle(Form("#beta %.1f < #eta < %.1f", etaMin, etaMax));
         betaPlots[i][j]->SetAxisRange(0, 1, "Y");
         betaPlots[i][j]->SetAxisRange(0, 100, "X");
         // betaPlots[i][j]->Draw("p");
      }
   }

   // alpha calculation (efficiency correction) ===============================
   if (!useCorrectionFile) {
      for (int x=1; x<=nEtaBin; x++) {
         for (int y=1; y<=nTrackletBin; y++) {
            for (int z=1; z<=nVzBin; z++) {
               alphaPlots[x-1][z-1]->SetBinContent(y, 0);
               alphaPlots[x-1][z-1]->SetBinError(y, 0);
               alphaErrPlots[x-1][z-1]->SetBinContent(y, 0);
               if (hAcceptance1D->GetBinContent(x, z) == 0) continue;

               if (hEverything->GetBinContent(x, y, z)!=0 && hHadron->GetBinContent(x, y, z)!=0) {
                  double val = hEverything->GetBinContent(x, y, z);
                  double beta = betaPlots[x-1][z-1]->GetBinContent(y);
                  double e1 = hEverything->GetBinError(x, y, z);
                  double e2 = hReproducedBackground->GetBinError(x, y, z);
                  double nsig = val*(1-beta);
                  double valErr = sqrt(e1*e1 + e2*e2);
                  double truth = hHadron->GetBinContent(x, y, z);
                  double truthErr = hHadron->GetBinError(x, y, z);
                  double alpha = truth/nsig;
                  if (verbose) cout << "alpha calc: " << x << " " << y << " " << z << " " << truth << " " << val << " " << (1-beta) << " " << endl;
                  double alphaErr = truth/nsig * sqrt(valErr/nsig*valErr/nsig + truthErr/truth*truthErr/truth * 0);
                  if (beta!=1 && alpha/alphaErr>-3 && alpha<3 && alpha>0) {
                     alphas->Fill((EtaBins[x]+EtaBins[x-1])/2, (TrackletBins[y]+TrackletBins[y-1])/2, (VzBins[z]+VzBins[z-1])/2, alpha, alphaErr);
                     alphaPlots[x-1][z-1]->SetBinContent(y, alpha);
                     alphaPlots[x-1][z-1]->SetBinError(y, alphaErr);
                     alphaErrPlots[x-1][z-1]->SetBinContent(y, alphaErr);
                  }
               }
            }
         }
      }
   }

   // Alpha correction calculation ============================================
   if (useCorrectionFile) {
      hTriggerCorrection = (TH1F*)fCorrection->FindObjectAny("hTriggerCorrection");
      hTriggerCorrection->SetName("hTriggerCorrection");
      // Use the alpha value obtained from the Correction file.
      for (int i=0; i<nEtaBin; i++) {
         for (int j=0; j<nVzBin; j++) {
            fAlpha[i][j] = (TF1*)fCorrection->FindObjectAny(Form("funAlpha%dVz%d", i, j));
            fAlphaErr[i][j] = (TF1*)fCorrection->FindObjectAny(Form("funAlphaErr%dVz%d", i, j));
            alphaPlots[i][j] = (TH1F*)fCorrection->FindObjectAny(Form("alpha%dVz%d", i, j));
            betaMCPlots[i][j] = (TH1F*)fCorrection->FindObjectAny(Form("beta%dVz%d", i, j));
         }
      }
   } else {
      for (int j=0; j<nVzBin; j++) {
         for (int i=0; i<nEtaBin; i++) {
            // c[i]= new TCanvas (Form("c%d", i), "", canvasSizeX, canvasSizeY);
            formatHist(alphaPlots[i][j], 2, 1);
            funAlpha[i][j] = new TF1(Form("funAlpha%dVz%d", i, j), "[1]/(x+[3]+0.5)+[2]/(x+0.5)/(x+0.5)+[0]", 0, 100);
            funAlphaErr[i][j] = new TF1(Form("funAlphaErr%dVz%d", i, j), "[0]+[1]/(x+0.5)+[2]*exp([3]*x)", 0, 100);
            double etaMin = i*0.2-3;
            double etaMax = i*0.2-3+0.2;
            alphaPlots[i][j]->Fit(Form("funAlpha%dVz%d", i, j), "M E Q", "", 0, 200);
            alphaPlots[i][j]->SetXTitle("N_{Hits}");
            alphaPlots[i][j]->SetYTitle(Form("#alpha %.1f < #eta < %.1f", etaMin, etaMax));
            alphaPlots[i][j]->SetAxisRange(0.9, 2.0, "Y");
            alphaPlots[i][j]->SetAxisRange(0, 100, "X");
            if (i<2 || i>9) alphaPlots[i][j]->SetAxisRange(4, 10, "y");
            alphaErrPlots[i][j]->Fit(Form("funAlphaErr%dVz%d", i, j), "M E Q", "", 0, 200);
            // alphaPlots[i][j]->Draw("p");
            fAlpha[i][j] = alphaPlots[i][j]->GetFunction(Form("funAlpha%dVz%d", i, j));
            fAlphaErr[i][j] = alphaErrPlots[i][j]->GetFunction(Form("funAlphaErr%dVz%d", i, j));
         }
      }

      // Vertex and event selection efficiency
      // TCanvas *cTrigEff = new TCanvas("cTrigEff", "TrigEff", canvasSizeX, canvasSizeY);
      TrackletTree->Project("hTrigEff", Form("%s", multiplicity), TCut(offlineSelection) && "vz[1]>-99" && MCSelection);
      hTrigEff->Sumw2();
      TrackletTree->Project("hTrigEffNoCut", Form("%s", multiplicity), MCSelection && "vz[1]>-99");
      hTrigEffNoCut->Sumw2();
      hTrigEff->Divide(hTrigEffNoCut);
      // hTrigEff->Draw();
      // cTrigEff->SaveAs(Form("figs/TrigEff-%s-%d.png", title, TrackletType));

      // Calculate SD'/IN'
      // TCanvas *cSDFrac = new TCanvas("cSDFrac", "SD Fraction After Cut", canvasSizeX, canvasSizeY);
      TrackletTree->Project("hSDFrac", Form("%s", multiplicity), evtSelection && !(MCSelection));
      hSDFrac->Sumw2();
      TrackletTree->Project("hSD", Form("%s", multiplicity), evtSelection);
      hSD->Sumw2();
      hSDFrac->Divide(hSD);
      // hSDFrac->Draw();
      // cSDFrac->SaveAs(Form("figs/SDFrac-%s-%d.png", title, TrackletType));

      // Calculate Vertexed dN/deta
      // TCanvas *cdNdEtaVertexed = new TCanvas("cdNdEtaVertexed", "dNdEta after vertexing", canvasSizeX, canvasSizeY);
      int nEvtAfterVtx = TrackletTree->Draw("nhit1", TCut("vz[1]>-99") && (MCSelection), "goff");
      TrackletTree->Project("hDNDEtaVertexed", "eta", TCut("vz[1]>-99") && (MCSelection));
      hDNDEtaVertexed->SetXTitle("#eta Truth after Vtx");
      hDNDEtaVertexed->SetYTitle("dN/#eta");
      hDNDEtaVertexed->Scale(2./nEvtAfterVtx);
      // hDNDEtaVertexed->Draw();

      // Calculate NoVertexed dN/deta
      // TCanvas *cdNdEtaNoVertexed = new TCanvas("cdNdEtaNoVertexed", "dNdEta after vertexing", canvasSizeX, canvasSizeY);
      int nEvtBeforeVtx = TrackletTree->Draw("nhit1", (MCSelection), "goff");
      TrackletTree->Project("hDNDEtaNoVertexed", "eta", (MCSelection));
      hDNDEtaNoVertexed->SetXTitle("#eta Truth before Vtx");
      hDNDEtaNoVertexed->SetYTitle("dN/#eta");
      hDNDEtaNoVertexed->Scale(2./nEvtBeforeVtx);
      // hDNDEtaNoVertexed->Draw();

      // Calculate Xi (Trigger correction)
      // TCanvas *cTriggerCorrection = new TCanvas("cTriggerCorrection", "Xi", canvasSizeX, canvasSizeY);
      int nEvtAfterEvtCut = TrackletTree->Draw("nhit1", evtSelection, "goff");
      TH1F* hdNdetaWithEvtCut = new TH1F("hdNdetaWithEvtCut", "", nEtaBin, -3, 3);
      TrackletTree->Project("hdNdetaWithEvtCut", "eta", evtSelection);
      hdNdetaWithEvtCut->Sumw2();
      hdNdetaWithEvtCut->Scale(1./nEvtAfterEvtCut);

      hTriggerCorrection = (TH1F*)hHadronWOSelection->Project3D("x");
      hTriggerCorrection->SetName("hTriggerCorrection");
      hTriggerCorrection->Sumw2();
      hTriggerCorrection->Scale(1./neventWOSelection);
      hTriggerCorrection->Divide(hdNdetaWithEvtCut);
      // hTriggerCorrection->Draw();
      // cTriggerCorrection->SaveAs(Form("figs/Xi-%s-%d.png", title, TrackletType));
   }

   // Make beta and alpha plot
   // if (plotAlphaBeta) {
   //    // Beta plot
   //    for (int z=1; z<=nVzBin; z++) {
   //       TLegend * l1 = new TLegend(0.63, 0.75, 0.93, 0.93);
   //       l1->SetFillStyle(0);
   //       l1->SetFillColor(0);
   //       l1->SetBorderSize(0);
   //       l1->SetTextSize(0.04);

   //       TLegend * l3 = new TLegend(0.26, 0.75, 0.56, 0.93);
   //       l3->SetFillStyle(0);
   //       l3->SetFillColor(0);
   //       l3->SetBorderSize(0);
   //       l3->SetTextSize(0.04);

   //       TCanvas *cc = new TCanvas(Form("cBetaPlot%d", z), "", canvasSizeX, canvasSizeY);
   //       for (int x=2; x<nEtaBin; x++) {
   //          int color = x-1;
   //          int mtype = 4;
   //          if (color>nEtaBin/2-1) {
   //             color = nEtaBin-1-color;
   //             mtype = 20;
   //          }
   //          formatHist(betaPlots[x-1][z-1], color);
   //          if (betaPlots[x-1][z-1]->GetEntries()) {
   //             betaPlots[x-1][z-1]->SetMarkerStyle(mtype);
   //             betaPlots[x-1][z-1]->SetAxisRange(0, 0.5, "Y");
   //             betaPlots[x-1][z-1]->Draw("L hist same");
   //             betaPlots[x-1][z-1]->Draw("e same");
   //             betaPlots[x-1][z-1]->SetYTitle("#beta");
   //             betaPlots[x-1][z-1]->SetXTitle(Form("N_{Hits} (%.0f #leq V_{z} < %.0f cm)", VzBins[z-1], VzBins[z]));
   //             if (x-1<nEtaBin/2)
   //                l3->AddEntry(betaPlots[x-1][z-1], Form("%.1f #leq #eta<%.1f", EtaBins[x-1], EtaBins[x]), "pl");
   //             else
   //                l1->AddEntry(betaPlots[x-1][z-1], Form("%.1f #leq #eta<%.1f", EtaBins[x-1], EtaBins[x]), "pl");
   //          }
   //       }
   //       l1->Draw();
   //       l3->Draw();
   //       cc->SaveAs(Form("figs/betaPlot/betaPlot-%s-%d-%d.png", title, z, TrackletType));
   //       cc->Close();
   //    }

   //    for (int z=1; z<=nVzBin; z++) {
   //       TLegend * l1 = new TLegend(0.63, 0.55, 0.93, 0.93);
   //       l1->SetFillStyle(0);
   //       l1->SetFillColor(0);
   //       l1->SetBorderSize(0);
   //       l1->SetTextSize(0.04);

   //       TLegend * l3 = new TLegend(0.26, 0.55, 0.56, 0.93);
   //       l3->SetFillStyle(0);
   //       l3->SetFillColor(0);
   //       l3->SetBorderSize(0);
   //       l3->SetTextSize(0.04);

   //       TCanvas *cc = new TCanvas(Form("cAlphaPlot%d", z), "", canvasSizeX, canvasSizeY);
   //       cc->SetLogy();
   //       for (int x=2; x<nEtaBin; x++) {
   //          int color = x-1;
   //          int mtype = 4;
   //          if (color>nEtaBin/2-1) {
   //             color = nEtaBin-1-color;
   //             mtype = 20;
   //          }
   //          formatHist(alphaPlots[x-1][z-1], color);
   //          if (alphaPlots[x-1][z-1]->GetEntries()) {
   //             alphaPlots[x-1][z-1]->SetMarkerStyle(mtype);
   //             alphaPlots[x-1][z-1]->SetYTitle("#alpha");
   //             alphaPlots[x-1][z-1]->SetXTitle(Form("N_{Hits} (%.0f #leq V_{z} < %.0f cm)", VzBins[z-1], VzBins[z]));
   //             alphaPlots[x-1][z-1]->SetStats(0);
   //             alphaPlots[x-1][z-1]->GetYaxis()->SetMoreLogLabels(1);
   //             alphaPlots[x-1][z-1]->SetAxisRange(0.9, 12, "Y");
   //             alphaPlots[x-1][z-1]->GetYaxis()->SetNoExponent(1);
   //             alphaPlots[x-1][z-1]->Draw("L hist same");
   //             alphaPlots[x-1][z-1]->Draw("hist error p same");
   //             if (x-1<nEtaBin/2)
   //                l3->AddEntry(alphaPlots[x-1][z-1], Form("%.1f #leq #eta<%.1f", EtaBins[x-1], EtaBins[x]), "pl");
   //             else
   //                l1->AddEntry(alphaPlots[x-1][z-1], Form("%.1f #leq #eta<%.1f", EtaBins[x-1], EtaBins[x]), "pl");
   //          }
   //       }
   //       l1->Draw();
   //       l3->Draw();
   //       cc->SaveAs(Form("figs/alphaPlot/alphaPlot-%s-%d-%d.png", title, z, TrackletType));
   //       cc->Close();
   //    }
   // }

   // Apply correction ========================================================
   for (int x=1; x<=nEtaBin; x++) {
      double totalN = 0;
      double totalNErr = 0;
      for (int y=1; y<=nTrackletBin; y++) {
         for (int z=1; z<=nVzBin; z++) {
            if (hAcceptance1D->GetBinContent(x, z) == 0) {
               hHadronAccepted->SetBinContent(x, y, z, 0);
               hHadronAccepted->SetBinError(x, y, z, 0);
               continue;
            }
            double val = hEverything->GetBinContent(x, y, z);
            double beta = betaPlots[x-1][z-1]->GetBinContent(y);
            double betaErr = betaPlots[x-1][z-1]->GetBinError(y);

            double alpha = alphaPlots[x-1][z-1]->GetBinContent(y);
            double alphaErr = alphaPlots[x-1][z-1]->GetBinError(y);

            double betaMC = betaMCPlots[x-1][z-1]->GetBinContent(y);

            if (doAccepCorr) {
               double accData = hDataAcc->GetBinContent(x, z);
               double accMC = hMCAcc->GetBinContent(x, z);
               if (accData==0 && accMC==0) {
                  cout << "Error in Acceptance Correction!!!! " << x << " " << z << endl;
               } else {
                  alpha = alpha * accMC / accData;
                  alphaErr = alphaErr * accMC / accData;
                  // cout << accMC/accData << endl;
               }
            }

            // Use extrapolated value if alpha is not available
            if (alpha==0 && fAlpha[x-1][z-1]!=0)
               alpha = fAlpha[x-1][z-1]->Eval(TrackletBins[y]);
            if (alpha == 0) {
               for (int k=0; k<y; k++) {
                  alpha = alphaPlots[x-1][z-1]->GetBinContent(y-k);
                  alphaErr = alphaPlots[x-1][z-1]->GetBinError(y-k);
                  if (alpha!=0) break;
               }
               // if (alpha==0) alpha = 1;
               // cout << (EtaBins[x]+EtaBins[x-1])/2 << " " << (TrackletBins[y]+TrackletBins[y-1])/2 << " " << (VzBins[z]+VzBins[z-1])/2 << "Used " << alpha << " " << endl;
               // cout << "Empty!!!" << endl;
            }

            if (alpha==0 || alpha>20) {
               hAcceptance1D->SetBinContent(x, z, 0);
               hHadronAccepted->SetBinContent(x, y, z, 0);
               hHadronAccepted->SetBinError(x, y, z, 0);
               continue;
            }

            double nCorrected = val*(1-beta)*alpha;
            hSubtracted->SetBinContent(x, y, z, nCorrected/alpha);
            hCorrected->SetBinContent(x, y, z, nCorrected);
            if (verbose) cout << "apply: " << x << " " << y << " " << z << " " << nCorrected << " " << alpha << " " << beta << " " << val << endl;
            double valErr = sqrt(alpha*(1-beta)*alpha*hEverything->GetBinContent(x, y, z)*(1-beta) +
                                 betaErr*betaErr*alpha*hEverything->GetBinContent(x, y, z)*alpha*hEverything->GetBinContent(x, y, z));
            hCorrected->SetBinError(x, y, z, valErr);
            correction->Fill((EtaBins[x]+EtaBins[x-1])/2, (TrackletBins[y]+TrackletBins[y-1])/2, (VzBins[z]+VzBins[z-1])/2, alpha, alphaErr, beta, betaErr, nCorrected, hHadronAccepted->GetBinContent(x, y, z), valErr, betaMC);
            totalN += nCorrected;
            totalNErr += valErr*valErr;
            if (verbose) cout << x << " " << y << " " << z << " " << nCorrected << " " << valErr << endl;

            hEtaVzRatio->Fill(x, z, nCorrected-hHadronAccepted->GetBinContent(x, y, z));
            hEtaHitRatio->Fill(x, y, nCorrected-hHadronAccepted->GetBinContent(x, y, z));
         }
      }
      hCorrectedEtaBin->SetBinContent(x, totalN);
      hCorrectedEtaBin->SetBinError(x, sqrt(totalNErr));
   }

   // Plot RawTracklet and Background Tracklet in nTracklet bin
   // TCanvas *cRawTrackletnTracklet = new TCanvas("cRawTrackletnTracklet", "Raw (nTracklet)", canvasSizeX, canvasSizeY);
   TH1F* hMCTruthnTracklet = (TH1F*)hHadronAccepted->Project3D("y");
   hMCTruthnTracklet->Scale(1./nevent);
   hMCTruthnTracklet->SetXTitle("N_{Hit1}|#eta|<1");
   // if (isMC) hMCTruthnTracklet->Draw("hist");

   TH1F* hRawTrackletnTracklet = (TH1F*)hEverything->Project3D("y");
   hRawTrackletnTracklet->Scale(1./nevent);
   // hRawTrackletnTracklet->Draw("same");

   TH1F* hBackgroundTrackletnTracklet = (TH1F*)hReproducedBackground->Project3D("y");
   hBackgroundTrackletnTracklet->SetLineColor(2);
   hBackgroundTrackletnTracklet->SetMarkerColor(2);
   hBackgroundTrackletnTracklet->Scale(1./nevent);
   // hBackgroundTrackletnTracklet->Draw("same");

   TH1F* hRawMinusBackgroundTrackletnTracklet = (TH1F*)hRawTrackletnTracklet->Clone();
   hRawMinusBackgroundTrackletnTracklet->SetLineColor(4);
   hRawMinusBackgroundTrackletnTracklet->SetMarkerColor(4);
   hRawMinusBackgroundTrackletnTracklet->Add(hBackgroundTrackletnTracklet, -1);
   // hRawMinusBackgroundTrackletnTracklet->Draw("same");
   // cRawTrackletnTracklet->Draw();

   // Plot RawTracklet and Background Tracklet in Vz bin
   // TCanvas *cRawTrackletVz = new TCanvas("cRawTrackletVz", "Raw (Vz)", canvasSizeX, canvasSizeY);
   TH1F* hMCTruthVz = (TH1F*)hHadronAccepted->Project3D("z");
   hMCTruthVz->Scale(1./nevent);
   hMCTruthVz->SetXTitle("V_{z}");
   // hMCTruthVz->Draw("hist");

   TH1F* hRawTrackletVz = (TH1F*)hEverything->Project3D("z");
   hRawTrackletVz->Scale(1./nevent);
   // hRawTrackletVz->Draw("same");

   TH1F* hBackgroundTrackletVz = (TH1F*)hReproducedBackground->Project3D("z");
   hBackgroundTrackletVz->SetLineColor(2);
   hBackgroundTrackletVz->SetMarkerColor(2);
   hBackgroundTrackletVz->Scale(1./nevent);
   // hBackgroundTrackletVz->Draw("same");

   TH1F* hRawMinusBackgroundTrackletVz = (TH1F*)hRawTrackletVz->Clone();
   hRawMinusBackgroundTrackletVz->SetLineColor(4);
   hRawMinusBackgroundTrackletVz->SetMarkerColor(4);
   hRawMinusBackgroundTrackletVz->Add(hBackgroundTrackletVz, -1);
   // hRawMinusBackgroundTrackletVz->Draw("same");
   // cRawTrackletVz->Draw();

   // Plot RawTracklet and Background Tracklet in eta bin
   // TCanvas *cRawTrackletEta = new TCanvas("cRawTrackletEta", "Raw (Eta)", canvasSizeX, canvasSizeY);
   TH1F* hMCTruthEta = (TH1F*)hHadronAccepted->Project3D("x");
   hMCTruthEta->Scale(2./nevent);
   hMCTruthEta->SetXTitle("#eta");
   // hMCTruthEta->Draw("hist");

   TH1F* hRawTrackletEta = (TH1F*)hEverything->Project3D("x");
   hRawTrackletEta->Scale(2./nevent);
   // hRawTrackletEta->Draw("same");

   TH1F* hBackgroundTrackletEta = (TH1F*)hReproducedBackground->Project3D("x");
   hBackgroundTrackletEta->SetLineColor(2);
   hBackgroundTrackletEta->SetMarkerColor(2);
   hBackgroundTrackletEta->Scale(2./nevent);
   // hBackgroundTrackletEta->Draw("same");

   TH1F* hRawMinusBackgroundTrackletEta = (TH1F*)hRawTrackletEta->Clone();
   hRawMinusBackgroundTrackletEta->SetLineColor(4);
   hRawMinusBackgroundTrackletEta->SetMarkerColor(4);
   hRawMinusBackgroundTrackletEta->Add(hBackgroundTrackletEta, -1);
   // hRawMinusBackgroundTrackletEta->Draw("same");
   // cRawTrackletEta->Draw();

   // Plot Before Acceptance correction
   // TCanvas *cDNdEtaC = new TCanvas("cDNdEtaC", "Before Acceptance correction", canvasSizeX, canvasSizeY);
   TH1F* hTruthAccepted = (TH1F*)hHadronAccepted->Project3D("x");
   hTruthAccepted->SetName("hTruthAccepted");

   TH1F* hTruth = (TH1F*)hTruthAccepted->Clone();
   hTruth->SetName("hTruth");
   hTruth->Sumw2();
   formatHist(hTruth, 1, nevent/nEtaBin*6);
   hTruth->Divide(hAcceptance1D->ProjectionX());

   TH1F* hTruthEvtCutCorrectedByXi = (TH1F*)hHadron->Project3D("x");
   hTruthEvtCutCorrectedByXi->SetName("hTruthEvtCutCorrectedByXi");
   hTruthEvtCutCorrectedByXi->Sumw2();
   formatHist(hTruthEvtCutCorrectedByXi, 1, nevent/nEtaBin*6);

   formatHist(hCorrectedEtaBin, 2, nevent/nEtaBin*6, 1);
   // hCorrectedEtaBin->Draw("e");

   hTruthAccepted->Sumw2();
   formatHist(hTruthAccepted, 1, nevent/nEtaBin*6);
   hTruthAccepted->SetAxisRange(0, dndetaRange, "y");
   hTruthAccepted->SetXTitle("#eta (Calculated Hand)");
   hTruthAccepted->SetYTitle("dN/d#eta");
   // hTruthAccepted->Draw("hist same");
   // cDNdEtaC->Draw();

   // Final dN/deta results ===================================================
   TH1F* hMeasured = (TH1F*)hCorrected->Project3D("x");
   hMeasured->SetName("hMeasured");
   hMeasured->Sumw2();
   formatHist(hMeasured, 4, nevent/nEtaBin*6, 1);
   hMeasured->Divide(hAcceptance1D->ProjectionX());

   TH1F* hUncorrected = (TH1F*)hEverything->Project3D("x");
   hUncorrected->SetName("hUncorrected");
   hUncorrected->Sumw2();
   formatHist(hUncorrected, 8, nevent/nEtaBin*6, 1);
   hUncorrected->Divide(hAcceptance1D->ProjectionX());

   TH1F* hBackgroundSubtracted = (TH1F*)hSubtracted->Project3D("x");
   hBackgroundSubtracted->SetName("hBackgroundSubtracted");
   hBackgroundSubtracted->Sumw2();
   formatHist(hBackgroundSubtracted, 9, nevent/nEtaBin*6, 1);
   hBackgroundSubtracted->Divide(hAcceptance1D->ProjectionX());

   TH1F* hTruthWOSelection = (TH1F*)hHadronWOSelection->Project3D("x");
   hTruthWOSelection->SetName("hTruthWOSelection");
   formatHist(hTruthWOSelection, 2, neventWOSelection/nEtaBin*6);

   TH1F* hTruthWOSelectionMult = (TH1F*)hHadronWOSelection->Project3D("y");
   hTruthWOSelectionMult->SetName("hTruthWOSelectionMult");
   formatHist(hTruthWOSelectionMult, 2, neventWOSelection);

   // Different calculation
   TH2F* hMeasuredEtanTracklet = new TH2F("hMeasuredEtanTracklet", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins);
   TH2F* hTruthEtanTracklet = new TH2F("hTruthEtanTracklet", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins);

   for (int x=1; x<=nEtaBin; x++) {
      for (int y=1; y<=nTrackletBin; y++) {
         double total = 0, totalErr = 0;
         double totalMC = 0, totalMCErr = 0;
         for (int z=1; z<=nVzBin; z++) {
            if (hAcceptance1D->GetBinContent(x, z)!=0) {
               total += hCorrected->GetBinContent(x, y, z);
               double err = hCorrected->GetBinError(x, y, z);
               totalErr = sqrt(totalErr*totalErr+err*err);
               totalMC += hHadronAccepted->GetBinContent(x, y, z);
               double errMC = hHadronAccepted->GetBinError(x, y, z);
               totalMCErr = sqrt(totalMCErr*totalMCErr + errMC*errMC);
               if (verbose) cout << x << " " << y << " " << z << " " << hCorrected->GetBinContent(x, y, z) << " " << hCorrected->GetBinError(x, y, z) << endl;
            }
         }
         hMeasuredEtanTracklet->SetBinContent(x, y, total);
         hMeasuredEtanTracklet->SetBinError(x, y, totalErr);
         hTruthEtanTracklet->SetBinContent(x, y, totalMC);
         hTruthEtanTracklet->SetBinError(x, y, totalMCErr);
      }
   }

   // Prepared dN/dhit/deta, apply vertexing correction
   double nEvt = 0;
   for (int y=1; y<=nTrackletBin; y++) {
      double SDFrac = hSDFrac->GetBinContent(y);
      double TrigEff = hTrigEff->GetBinContent(y);
      if (!doTriggerCorr) {
         SDFrac = 0;
         TrigEff = 1;
      }
      if (TrigEff != 0) nEvt += hnTracklet->GetBinContent(y)/TrigEff*(1-SDFrac*SDFactor);
      for (int x=1; x<=nEtaBin; x++) {
         for (int z=1; z<=nVzBin; z++) {
            if (hAcceptance1D->GetBinContent(x, z) != 0) {
               if (TrigEff != 0)
                  hAcceptance2D->SetBinContent(x, y, z, hAcceptance2D->GetBinContent(x, y, z)/TrigEff*(1-SDFrac*SDFactor));
            } else {
               hAcceptance2D->SetBinContent(x, y, z, 0);
            }
         }
      }
   }

   for (int x=1; x<=nEtaBin; x++) {
      for (int y=1; y<=nTrackletBin; y++) {
         double SDFrac = hSDFrac->GetBinContent(y);
         double TrigEff = hTrigEff->GetBinContent(y);
         if (!doTriggerCorr) {
            SDFrac = 0;
            TrigEff = 1;
         }
         if (TrigEff != 0) {
            hMeasuredEtanTracklet->SetBinContent(x, y, hMeasuredEtanTracklet->GetBinContent(x, y)/TrigEff*(1-SDFrac*SDFactor));
            hMeasuredEtanTracklet->SetBinError(x, y, hMeasuredEtanTracklet->GetBinError(x, y)/TrigEff*(1-SDFrac*SDFactor));
            hTruthEtanTracklet->SetBinContent(x, y, hTruthEtanTracklet->GetBinContent(x, y)/TrigEff*(1-SDFrac*SDFactor));
            hTruthEtanTracklet->SetBinError(x, y, hTruthEtanTracklet->GetBinError(x, y)/TrigEff*(1-SDFrac*SDFactor));
         }
      }
   }

   TH1F* hAcceptance2DEta = (TH1F*)hAcceptance2D->Project3D("x");
   hAcceptance2DEta->SetName("hAcceptance2DEta");

   TH1F* hAcceptance2DMult = (TH1F*)hAcceptance2D->Project3D("y");
   hAcceptance2DMult->SetName("hAcceptance2DMult");

   TH1F* hMeasuredTrigEffCorrected = (TH1F*)hMeasuredEtanTracklet->ProjectionX();
   hMeasuredTrigEffCorrected->SetName("hMeasuredTrigEffCorrected");
   formatHist(hMeasuredTrigEffCorrected, 2, 1./nEtaBin*6, 1);
   hMeasuredTrigEffCorrected->Divide(hAcceptance2DEta);
   hMeasuredTrigEffCorrected->SetMarkerStyle(4);

   TH1F* hMeasuredTrigEffCorrectedMult = (TH1F*)hMeasuredEtanTracklet->ProjectionY();
   hMeasuredTrigEffCorrectedMult->SetName("hMeasuredTrigEffCorrectedMult");
   formatHist(hMeasuredTrigEffCorrectedMult, 8, hAcceptance2DEta->Integral()/nEtaBin, 1);
   // hMeasuredTrigEffCorrectedMult->Divide(hAcceptance2DMult);

   TCanvas* cEmpty = new TCanvas("cEmpty", "Empty Correction", canvasSizeX, canvasSizeY);
   if (!useCorrectionFile) {
      hEmptyEvtCorrection = (TH1F*)hTruthWOSelection->Clone();
      hEmptyEvtCorrection->SetName("hEmptyEvtCorrection");
      hEmptyEvtCorrection->Divide(hMeasuredTrigEffCorrected);
      for (int x=1; x<=nEtaBin; x++)
         hEmptyEvtCorrection->SetBinError(x, 0);

      hEmptyEvtCorrection->Fit("pol2", "LL", "", -etaLimit+0.3, etaLimit-0.3);
      fEmptyEvt = hEmptyEvtCorrection->GetFunction("pol2");
      fEmptyEvt->SetName("fEmptyEvt");
   } else {
      hEmptyEvtCorrection = (TH1F*)fCorrection->FindObjectAny("hEmptyEvtCorrection");
      fEmptyEvt = (TF1*)fCorrection->FindObjectAny("fEmptyEvt");
   }
   hEmptyEvtCorrection->Draw();
   cEmpty->Draw();

   // Final result canvas
   TCanvas* cDNdEta = new TCanvas("cDNdEta", "Final result", canvasSizeX, canvasSizeY);

   hTruthWOSelection->SetAxisRange(0, dndetaRange, "y");
   hTruthWOSelection->SetXTitle("#eta");
   hTruthWOSelection->SetYTitle("dN/d#eta");
   hTruthWOSelection->SetStats(0);
   hTruthWOSelection->Draw("hist");

   TH1F* hMeasuredFinal = (TH1F*)hMeasuredTrigEffCorrected->Clone();
   hMeasuredFinal->SetName("hMeasuredFinal");
   hMeasuredFinal->SetMarkerStyle(20);
   if (doTriggerCorr) {
      for (int x=1; x<=nEtaBin; x++) {
         // double emptyCorrection = fEmptyEvt->Eval((EtaBins[x-1]+EtaBins[x])/2);
         double emptyCorrection = hEmptyEvtCorrection->GetBinContent(x);
         hMeasuredFinal->SetBinContent(x, hMeasuredFinal->GetBinContent(x)*emptyCorrection);
         hMeasuredFinal->SetBinError(x, hMeasuredFinal->GetBinError(x)*emptyCorrection);
      }
   }
   formatHist(hMeasuredFinal, 2, 1, 1.1);
   hMeasuredFinal->SetStats(0);
   hMeasuredFinal->Draw("e same");

   hTruthEvtCutCorrectedByXi->Multiply(hTriggerCorrection);
   hTruthEvtCutCorrectedByXi->SetAxisRange(0, dndetaRange, "y");
   hTruthEvtCutCorrectedByXi->SetXTitle("#eta");
   hTruthEvtCutCorrectedByXi->SetYTitle("dN/d#eta");
   // hTruthEvtCutCorrectedByXi->Draw("hist");

   // double systematicError13TeV[30];
   // TGraph* gErrorBand = GetErrorBand(hMeasuredFinal, systematicError13TeV, systematicError13TeV, 0.1);
   // gErrorBand->Draw("F");

   TLegend* l1 = new TLegend(0.3, 0.3, 0.8, 0.48);
   l1->SetFillStyle(0);
   l1->SetFillColor(0);
   l1->SetBorderSize(0);
   l1->SetTextFont(43);
   l1->SetTextSize(20);
   l1->AddEntry(hTruth, Form("%s", title));
   l1->AddEntry(hTruthWOSelection, "Truth", "l");
   l1->AddEntry(hMeasuredFinal, "Reconstructed", "pl");
   l1->Draw();

   // TText* text = new TText(-2.6, 5, "CMS Preliminary");
   // text->Draw();

   cDNdEta->Draw();
   cDNdEta->SaveAs(Form("figs/result/result-%s-%d.pdf", title, TrackletType));

   // Compare with Truth ======================================================
   TCanvas* cDNdEtaCompare = new TCanvas("cDNdEtaCompare", "Compare", canvasSizeX, canvasSizeY);

   hMeasuredFinal->SetAxisRange(0, dndetaRange, "y");
   hMeasuredFinal->Draw("e same");

   hTruthWOSelection->Draw("hist same");

   hUncorrected->SetMarkerStyle(25);
   hUncorrected->Draw("p same");

   hBackgroundSubtracted->SetMarkerStyle(26);
   hBackgroundSubtracted->Draw("p same");

   TH1F* hMeasuredNoCorrection = (TH1F*)hMeasured->Clone();
   hMeasuredNoCorrection->SetName("hMeasuredNoCorrection");
   hMeasuredNoCorrection->SetMarkerStyle(4);
   hMeasuredNoCorrection->SetMarkerColor(1);
   hMeasuredNoCorrection->SetLineColor(1);
   hMeasuredNoCorrection->Draw("e same");

   hMeasuredTrigEffCorrected->Draw("e same");

   TH1F* hTruthTrigEffCorrected = (TH1F*)hTruthEtanTracklet->ProjectionX();
   hTruthTrigEffCorrected->SetName("hTruthTrigEffCorrected");
   formatHist(hTruthTrigEffCorrected, 4, 1./nEtaBin*6, 1);
   hTruthTrigEffCorrected->Divide(hAcceptance2DEta);
   hTruthTrigEffCorrected->Draw("e same");

   TLegend* l2 = new TLegend(0.24, 0.16, 1, 0.45);
   l2->SetFillStyle(0);
   l2->SetFillColor(1001);
   l2->SetBorderSize(0);
   l2->SetTextFont(43);
   l2->SetTextSize(15);
   l2->AddEntry(hTruth, Form("%s", title));
   l2->AddEntry(hTruthWOSelection, "Truth", "l");
   l2->AddEntry(hTruthTrigEffCorrected, "Truth corrected for trigger eff", "pl");
   l2->AddEntry(hUncorrected, "Raw tracklets", "pl");
   l2->AddEntry(hBackgroundSubtracted, "Background subtracted", "pl");
   l2->AddEntry(hMeasuredNoCorrection, "Corrected for efficiency", "pl");
   l2->AddEntry(hMeasuredTrigEffCorrected, "Corrected for trigger eff", "pl");
   l2->AddEntry(hMeasuredFinal, "Final result", "pl");
   l2->Draw();

   cDNdEtaCompare->Draw();
   cDNdEtaCompare->SaveAs(Form("figs/compare/compare-%s-%i.png", title, TrackletType));

   if (isMC) {
      // Ratio between measured and truth =====================================
      TCanvas* cRatio = new TCanvas("cRatio", "Ratio", canvasSizeX, canvasSizeY);

      TH1F* hRatio = (TH1F*)hMeasuredFinal->Clone();
      hRatio->SetName("hRatio");
      hRatio->Divide(hTruthWOSelection);
      hRatio->SetXTitle("#eta");
      hRatio->SetYTitle("Ratio");
      hRatio->SetAxisRange(0.8, 1.2, "y");
      hRatio->SetStats(0);
      hRatio->Draw();

      // TH1F* hRatioStat = (TH1F*)hRatio->Clone();
      // hRatioStat->SetName("hRatioStat");
      // hRatioStat->SetMarkerSize(0);
      // hRatioStat->Draw("same");

      TLegend* l3 = new TLegend(0.18, 0.24, 0.9, 0.36);
      l3->SetFillStyle(0);
      l3->SetFillColor(0);
      l3->SetBorderSize(0);
      l3->AddEntry(hTruth, title);
      l3->AddEntry(hRatio, "Reconstructed / MC Truth", "pl");
      l3->Draw();

      TLine* line1 = new TLine(-3, 1, 3, 1);
      line1->Draw();

      // TText* text2 = new TText(-2.6, 1.165, "CMS Preliminary");
      // text2->Draw();

      cRatio->Draw();
      cRatio->SaveAs(Form("figs/ratio/ratio-%s-%d.pdf", title, TrackletType));

      // TCanvas *cDNdnTracklet = new TCanvas("cDNdnTracklet", "Measured vs mult", canvasSizeX, canvasSizeY);
      TH1F* hTruthHit = (TH1F*)hHadronAccepted->Project3D("y");
      hTruthHit->Sumw2();
      formatHist(hTruthHit, 1, nevent);
      hTruthHit->SetAxisRange(0, dndetaRange, "y");
      hTruthHit->SetXTitle("N_{Hit1} |#eta|<3");
      hTruthHit->SetYTitle("dN/dN_{Hit1}");
      // hTruthHit->Draw("hist");

      TH1F* hMeasuredHit = (TH1F*)hCorrected->Project3D("y");
      hMeasuredHit->Sumw2();
      formatHist(hMeasuredHit, 2, nevent, 1);
      // hMeasuredHit->Draw("e same");
      // cDNdnTracklet->Draw();

      // TCanvas *cRatioTracklet = new TCanvas("cRatioTracklet", "Ratio vs mult", canvasSizeX, canvasSizeY);
      TH1F* hRatioHit = (TH1F*)hMeasuredHit->Clone();
      hRatioHit->Divide(hTruthHit);
      hRatioHit->SetXTitle("#eta");
      hRatioHit->SetYTitle("Ratio");
      hRatioHit->SetAxisRange(0.7, 1.3, "y");
      // hRatioHit->Draw();

      // TLine* line2 = new TLine(0, 1, nTrackletBin, 1);
      // line2->Draw();
      // cRatioTracklet->Draw();
   }

   for (int i=0; i<nEtaBin; i++) {
      for (int j=0; j<nVzBin; j++) {
         if (fAlpha[i][j]) fAlpha[i][j]->Write();
         if (fAlphaErr[i][j]) fAlphaErr[i][j]->Write();
         if (verbose && !fAlpha[i][j]) cout << "no alpha! " << i << " " << j << endl;
         if (verbose && !fAlphaErr[i][j]) cout << "no alphaErr! " << i << " " << j << endl;
      }
   }
   fEmptyEvt->Write();

   cDNdEta->Write();
   outf->Write();

   return 0;
}

// ============================================================================
// Format Histogram
// ============================================================================
void formatHist(TH1* h, int col, double norm, double msize) {
   h->Scale(1/norm);
   h->SetLineColor(col);
   h->SetMarkerColor(col);
   h->SetMarkerSize(msize);
   h->GetYaxis()->SetTitleOffset(1.25);
   h->GetXaxis()->CenterTitle();
   h->GetYaxis()->CenterTitle();
   h->SetTitle("");
}

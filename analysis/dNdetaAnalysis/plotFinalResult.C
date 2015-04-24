// Plot final results
#define canvasSizeX 600
#define canvasSizeY 600
#define dndetaRange 8.0
#define SDFactor 1

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

// For plotting
#include "UA5Plot.h"
#include "GraphErrorsBand.h"

#include "selectionCut.h"

void formatHist(TH1* h, int col = 1, double norm = 1, double msize = 1.1);

// ===========================================================================
// Read Stored correction factors
// ===========================================================================

// Correction factor setting:
// 1 = 900 GeV
// 2 = 10 TeV
// 3 = 900 GeV (half of the sample)
// 4 = 10 TeV (half of the ample)

TFile* getAcceptanceFile(int TrackletType) {
   TFile* fAcceptance;
   const char* filename = Form("correction/Acceptance-%d.root", TrackletType);
   fAcceptance = new TFile(filename);
   cout << "Use Acceptance file : " << filename << endl;
   return fAcceptance;
}

TFile* getCorrectionFile(string correctionFileName, int TrackletType) {
   TFile* fCorrection;
   const char* filename = Form("correction/correction-%d-%s.root", TrackletType, correctionFileName.data());
   fCorrection = new TFile(filename);
   cout << "Use correction file : " << filename << endl;
   return fCorrection;
}

// ===========================================================================
// Main Routine
// ===========================================================================
int plotFinalResult(int TrackletType, const char* filename,
                    const char* myPlotTitle = "PYTHIA8_Monash13",               // Title of the plot
                    bool useCorrectionFile = 0,                                 // use Correction file
                    string correctionName = "PYTHIA8_Monash13",                 // Correction file name
                    Long64_t nentries = 1000000000,                             // Number of entries
                    Long64_t firstentry = 0,                                    // First Entry
                    int LumiL = 0,                                              // nLumi Lower cut
                    int LumiH = 1000,                                           // nLumi Higher cut
                    int verbose = 0,                                            // Set Verbose level
                    int makePlot = 0,                                           // make alpha plots
                    bool putUA5 = 1,                                            // overlap UA5 result
                    bool doAcceptanceCorrection = 0,                            // do acceptance correction
                    bool doBetaCorrection = 0,                                  // do beta correction
                    int doMult2 = 0,                                            // use # of cluster in the first layer
                    bool doTriggerCorrection = 1,                               // correction on trigger eff
                    int UseExternalSDEff = 0,
                    bool useDR = 0)
{
   // Output log file:
   FILE* logFile = fopen("correction.log", "w");
   TFile* f= new TFile(filename);

   // Input trackletTree
   TTree* TrackletTree = dynamic_cast<TTree*>(f->Get(Form("TrackletTree%d", TrackletType)));
   TrackletTree->SetName("TrackletTree");
   bool isMC = false;
   if (TrackletTree->GetEntries("npart!=0")!=0) {
      isMC = true;
      cout << "This is a Monte Carlo study." << endl;
      doAcceptanceCorrection = 0;
   } else {
      cout << "This is a data analysis." << endl;
   }

   if (!doTriggerCorrection)
      cout << "Trigger correction off!!!" << endl;

   // choose multiplcity handle
   const char* multiplicity;
   if (doMult2==1) {
      multiplicity = "mult2";   // use # of cluster
      cout << "Use # of clusters as event multiplicity" << endl;
   } else if (doMult2==2) {
      multiplicity = "nTracklets";    // use # of Background subtracted tracklets
      cout << "Use # of tracklets as event multiplicity" << endl;
   } else {
      multiplicity = "mult";    // use # of Background subtracted tracklets
      cout << "Use # of tracklets after background subtraction as event multiplicity" << endl;
   }

   // Read alpha, beta, geometry correction from file.
   TFile* fCorrection;
   TFile* fAcceptance;
   if (useCorrectionFile) fCorrection = getCorrectionFile(correctionName, TrackletType);
   if (useCorrectionFile && doAcceptanceCorrection) fAcceptance = getAcceptanceFile(TrackletType);
   TFile *fCorrectionExternal = new TFile(Form("correction/correction-%d-PHOJET-7TeV-HF1.root", TrackletType));

   TH3F *hAlphaA;
   TH3F *hAlphaB;
   if (doBetaCorrection) {
      TFile *myFile = new TFile(Form("correction/alphaBetaCoeff-%d.root", TrackletType));
      hAlphaA = (TH3F*)myFile->FindObjectAny("hAlphaA");
      hAlphaB = (TH3F*)myFile->FindObjectAny("hAlphaB");
   }

   // Definition of Vz, Eta, Hit bins
   selectionCut myCut(isMC, LumiL, LumiH);
   const int nTrackletBin = 20; // myCut.nTrackletBin;
   const int nEtaBin = 30; // myCut.nEtaBin;
   const int nVzBin = 20; // myCut.nVzBin;
   int VzRangeL = myCut.VzRangeL;
   int VzRangeH = myCut.VzRangeH;

   double TrackletBins[nTrackletBin+1] = {-5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 45, 60, 90, 300};
   double EtaBins[nEtaBin+1];
   double VzBins[nVzBin+1];

   for (int i=0; i<=nEtaBin; i++)
      EtaBins[i] = (double)i*6.0/(double)nEtaBin-3.0; ;
   for (int i=0; i<=nVzBin; i++)
      VzBins[i] = (double)i*(VzRangeH-VzRangeL)/(double)nVzBin+VzRangeL;

   // Signal and Sideband regions ==================================================================================================
   double signalRegionCut = 1.0;      //delta phi cut for signal region
   double sideBandRegionCut = 2.0;    //delta phi cut for sideband
   double detaCut = 0.1;              //delta eta cut

   TCut signalRegion                  = Form("abs(dphi)<%f&&abs(deta)<%f", signalRegionCut, detaCut);
   TCut sideBandRegionEtaSignalRegion = Form("abs(dphi)>%f&&abs(dphi)<%f&&abs(deta)<%f", signalRegionCut, sideBandRegionCut, detaCut);
   if (useDR) {
      cout << "Use dR for analysis" << endl;
      signalRegion                  = "dR<0.1";
      sideBandRegionEtaSignalRegion = "dR>0.1&&dR<0.2";
   }

   TCut evtSelection = myCut.Cut;   // cut on Z position
   // TCut NSDCut = "";
   // TCut NSDCut = "(evtType!=103&&evtType!=104)"; cout << "PYTHIA8 MC NSD definition" << endl;
   TCut NSDCut = "1"; cout << "PYTHIA8 MC INEL definition" << endl;
   // TCut NSDCut = "(evtType==103||evtType==104)"; cout << "PYTHIA8 MC SD definition" << endl;
   // TCut NSDCut = "evtType!=5&&evtType!=6"; cout << "PHOJET MC definition" << endl; // NSD
   if (!isMC) NSDCut = "";
   // evtSelection += NSDCut;

   // Output file =================================================================================================================
   TFile* outf = new TFile ("correction.root", "recreate");

   TNtuple* betas = new TNtuple("betas", "", "eta:nTracklet:vz:beta:betaErr");
   TNtuple* alphas = new TNtuple("alphas", "", "eta:nTracklet:vz:alpha:alphaErr");
   TNtuple* correction = new TNtuple("correction", "", "eta:nTracklet:vz:alpha:alphaErr:beta:betaErr:obs:gen:err:betaMC");

   TH3F* hEverything = new TH3F("hEverything", "Everything in the signal region with vz cut", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* hReproducedBackground = new TH3F("hReproducedBackground", "Reproduced Background", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* hSubtracted = new TH3F("hSubtracted", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* hHadron = new TH3F("hHadron", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* hHadronAccepted = new TH3F("hHadronAccepted", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* hHadronWOSelection = new TH3F("hHadronWOSelection", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* hCorrected = new TH3F("hCorrected", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* hnTracklet3D = new TH3F("hnTracklet3D", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);

   TH2F *hEtaVzRatio = new TH2F("hEtaVzatio", "#eta:Vz", nEtaBin, EtaBins, nVzBin, VzBins);
   TH2F *hEtaHitRatio = new TH2F("hEtaHitatio", "#eta:N_{Hit}", nEtaBin, EtaBins, nTrackletBin, TrackletBins);
   TH2F *hAcceptance = new TH2F("hAcceptance", "", nEtaBin, EtaBins, nVzBin, VzBins);

   // Acceptance
   TH2F *hDataAcc = new TH2F("hDataAcc", "", nEtaBin, EtaBins, nVzBin, VzBins);
   TH2F *hMCAcc = new TH2F("hMCAcc", "", nEtaBin, EtaBins, nVzBin, VzBins);

   TH2F *hVzNTracklet = new TH2F("hVzNTracklet", "", nTrackletBin, TrackletBins, nVzBin, VzBins);

   TH1F *hCorrectedEtaBin = new TH1F("hCorrectedEtaBin", "Corrected", nEtaBin, -3, 3);
   TH1F *hDNDEtaVertexed = new TH1F("hDNDEtaVertexed", "", nEtaBin, -3, 3);
   TH1F *hDNDEtaNoVertexed = new TH1F("hDNDEtaNoVertexed", "", nEtaBin, -3, 3);
   TH1F *hTrigEff = new TH1F("hTrigEff", "", nTrackletBin, TrackletBins);
   TH1F *hTrigEffNoCut = new TH1F("hTrigEffNoCut", "", nTrackletBin, TrackletBins);
   TH1F *hSD = new TH1F("hSD", "", nTrackletBin, TrackletBins);
   TH1F *hSDFrac = new TH1F("hSDFrac", "", nTrackletBin, TrackletBins);
   TH1F *hTruthRatio;
   TH1F *hVz = new TH1F("hVz", "", nVzBin, VzBins);
   TH1F *hnTracklet = new TH1F("hnTracklet", "", nTrackletBin, TrackletBins);

   TH1F *alphaPlots[nEtaBin][nVzBin];
   TH1F *betaPlots[nEtaBin][nVzBin];
   TH1F *betaMCPlots[nEtaBin][nVzBin];
   TH1F *alphaErrPlots[nEtaBin][nVzBin];
   TH1F *betaErrPlots[nEtaBin][nVzBin];

   TH1F *hTriggerCorrection;
   TH1F *hEmptyEvtCorrection;

   // Prepare histograms ==========================================================================================================
   for (int i=0; i<nEtaBin; i++) {
      for (int j=0; j<nVzBin; j++) {
         alphaPlots[i][j] = new TH1F(Form("alpha%dVz%d", i, j), "", nTrackletBin, TrackletBins);
         alphaErrPlots[i][j] = new TH1F(Form("alphaErr%dVz%d", i, j), "", nTrackletBin, TrackletBins);
         betaPlots[i][j] = new TH1F(Form("beta%dVz%d", i, j), "", nTrackletBin, TrackletBins);
         betaErrPlots[i][j] = new TH1F(Form("betaErr%dVz%d", i, j), "", nTrackletBin, TrackletBins);
         betaMCPlots[i][j] = new TH1F(Form("betaMC%dVz%d", i, j), "", nTrackletBin, TrackletBins);
      }
   }

   // Fit functions of Beta and Alpha =============================================================================================
   TF1 *funAlpha[nEtaBin][nVzBin];
   TF1 *funAlphaErr[nEtaBin][nVzBin];
   TF1 *fAlpha[nEtaBin][nVzBin];
   TF1 *fAlphaErr[nEtaBin][nVzBin];

   TCanvas *cVz = new TCanvas("cVz", "Vz distribution", canvasSizeX, canvasSizeY);

   // Number of events
   double nevent = TrackletTree->Draw("vz[1]", evtSelection&&NSDCut, "", nentries, firstentry);
   double neventWOSelection = TrackletTree->Draw("vz[1]", NSDCut, "", nentries, firstentry);
   if (nevent==0) {
      cout << "No event survived. Abort." << endl;
      return 0;
   }

   // Acceptance calculation ======================================================================================================
   TrackletTree->Project("hnTracklet", Form("%s", multiplicity), evtSelection, "", nentries, firstentry);
   TrackletTree->Project("hVzNTracklet", Form("vz[1]:%s", multiplicity), evtSelection, "", nentries, firstentry);
   TrackletTree->Project("hVz", "vz[1]", evtSelection, "", nentries, firstentry);

   hVz->Sumw2();
   hVz->Scale(1./hVz->GetEntries());
   hVz->Fit("gaus");
   hVz->SetXTitle("v_{z} (cm)");
   hVz->Draw();
   cVz->SaveAs(Form("plot/vz/plotVz-%s-%d.gif", myPlotTitle, TrackletType));
   cVz->SaveAs(Form("plot/vz/plotVz-%s-%d.eps", myPlotTitle, TrackletType));
   cVz->SaveAs(Form("plot/vz/plotVz-%s-%d.C", myPlotTitle, TrackletType));
   cout << "Number of events: " << nevent << endl;
   // Determine the acceptance region to avoid large correction factors
   // Second layer rho (cm)
   double rho2 = 7.6;
   // For Tracklet 1+3 or Tracklet 2+3
   if (TrackletType % 10 == 3) rho2 = 10.5;
   // End point in z (cm)
   // double endpoint2 = 26.66;
   double endpoint2 = 30.0;

   for (int i=0; i<nEtaBin; i++) {
      for (int j=0; j<nVzBin; j++) {
         double minEta = EtaBins[i];
         double maxEta = EtaBins[i+1];
         double maxEdge = VzBins[j+1]-rho2/tan(atan(exp(maxEta-0.1))*2);
         double minEdge = VzBins[j]-rho2/tan(atan(exp(minEta+0.1))*2);
         if (verbose) cout << minEta << " " << maxEta << " " << VzBins[j] << " " << maxEdge << " " << minEdge;

         // if (maxEdge>-endpoint2 && minEdge<endpoint2 && (maxEta)<=2.5 && (minEta)>=-2.5) {
         //    hAcceptance->SetBinContent(i+1, j+1, hVz->GetBinContent(j+1));
         //    hAcceptance->SetBinError(i+1, j+1, 0);
         //    for (int k=0; k<nTrackletBin; k++)
         //       hnTracklet3D->SetBinContent(i+1, k+1, j+1, hVzNTracklet->GetBinContent(k+1, j+1));
         //    if (verbose) cout << " Selected! " << endl;
         // } else {
         //    hAcceptance->SetBinContent(i+1, j+1, 0);
         //    for (int k=0; k<nTrackletBin; k++)
         //       hnTracklet3D->SetBinContent(i+1, k+1, j+1, 0);
         // }
         for (int k=0; k<nTrackletBin; k++) {
            if (maxEdge>-endpoint2 && minEdge<endpoint2 && (maxEta)<=2.5 && (minEta)>=-2.5) {
               // hAcceptance->SetBinContent(i+1, j+1, hAcceptance->GetBinContent(i+1, j+1)+hVzNTracklet->GetBinContent(k+1, j+1)/hVzNTracklet->GetEntries());
               hAcceptance->SetBinContent(i+1, j+1, hVz->GetBinContent(j+1));
               hnTracklet3D->SetBinContent(i+1, k+1, j+1, hVzNTracklet->GetBinContent(k+1, j+1));
               hAcceptance->SetBinError(i+1, j+1, 0);
               if (verbose) cout << " Selected!";
            } else {
               hnTracklet3D->SetBinContent(i+1, k+1, j+1, 0);
               hAcceptance->SetBinContent(i+1, j+1, 0);
            }
            if (verbose) cout << " " << endl;
         }
      }
   }

   TCanvas *cAcc = new TCanvas("cAcc", "Acc", canvasSizeX, canvasSizeY);
   hAcceptance->ProjectionX()->Draw();

   // Charged Hadrons =============================================================================================================
   hHadron->SetXTitle("#eta");
   hHadron->SetYTitle("N_{hit}^{Layer1} |#eta|<1");

   TrackletTree->Project("hHadron", Form("vz[1]:%s:eta", multiplicity), "abs(eta)<3"&&evtSelection, "", nentries, firstentry);
   TrackletTree->Project("hHadronWOSelection", Form("vz[1]:%s:eta", multiplicity), "abs(eta)<3"&&NSDCut, "", nentries, firstentry);
   hHadron->Sumw2();
   hHadronWOSelection->Sumw2();

   TrackletTree->Project("hHadronAccepted", Form("vz[1]:%s:eta", multiplicity), "abs(eta)<3"&&evtSelection, "", nentries, firstentry);
   hHadronAccepted = (TH3F*) hHadron->Clone();
   hHadronAccepted->SetName("hHadronAccepted");

   // Prepare Tracklet Three-Dimensional Histogram ================================================================================
   // Signal region && evtSelection //
   TrackletTree->Project("hEverything", Form("vz[1]:%s:eta1", multiplicity), signalRegion&&evtSelection, "", nentries, firstentry);
   hEverything->Sumw2();

   // deltaPhi sideband region (with eta signal region cut) && evtSelection
   TrackletTree->Project("hReproducedBackground", Form("vz[1]:%s:eta1", multiplicity), sideBandRegionEtaSignalRegion&&evtSelection, "", nentries, firstentry);
   hReproducedBackground->Sumw2();

   // Read Acceptance ==================================================================================================
   if (doAcceptanceCorrection) {
      hMCAcc = (TH2F*)fAcceptance->FindObjectAny("hMCAcc");
      hDataAcc = (TH2F*)fAcceptance->FindObjectAny("hDataAcc");
   }

   // Beta calculation (background fraction)================================================================================================
   for (int x=1; x<=nEtaBin; x++) {
      for (int y=1; y<=nTrackletBin; y++) {
         for (int z=1; z<=nVzBin; z++) {
            double beta = 0;
            betaPlots[x-1][z-1]->SetBinContent(y, 0);
            betaPlots[x-1][z-1]->SetBinError(y, 0);
            betaErrPlots[x-1][z-1]->SetBinContent(y, 0);

            if (hAcceptance->GetBinContent(x, z)==0) continue;
            if (hEverything->GetBinContent(x, y, z)!=0) {
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
               // beta = hEverything->GetBinContent(x, y, z)-hReproducedBackground->GetBinContent(x, y, z);
               // betaErr = sqrt(e1*e1+e2*e2);
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

   // alpha calculation (efficiency correction) ===============================================================================
   if (!useCorrectionFile) {
      for (int x=1; x<=nEtaBin; x++) {
         for (int y=1; y<=nTrackletBin; y++) {
            for (int z=1; z<=nVzBin; z++) {
               // clear
               alphaPlots[x-1][z-1]->SetBinContent(y, 0);
               alphaPlots[x-1][z-1]->SetBinError(y, 0);
               alphaErrPlots[x-1][z-1]->SetBinContent(y, 0);
               if (hAcceptance->GetBinContent(x, z)==0) continue;

               if (hEverything->GetBinContent(x, y, z)!=0 && hHadron->GetBinContent(x, y, z)!=0) {
                  double val = hEverything->GetBinContent(x, y, z);
                  double beta = betaPlots[x-1][z-1]->GetBinContent(y);
                  double e1 = hEverything->GetBinError(x, y, z);
                  double e2 = hReproducedBackground->GetBinError(x, y, z);
                  double nsig = val*(1-beta);
                  double valErr = sqrt(e1*e1+e2*e2);
                  double truth = hHadron->GetBinContent(x, y, z);
                  double truthErr = hHadron->GetBinError(x, y, z);
                  double alpha = truth/nsig;
                  if (verbose) cout << "alpha calc: " << x << " " << y << " " << z << " " << truth << " " << val << " " << (1-beta) << " " << endl;
                  double alphaErr = truth/nsig*sqrt(valErr/nsig*valErr/nsig+truthErr/truth*truthErr/truth*0);
                  if (beta!=1 && alpha/alphaErr>-3 && alpha<500 && alpha>0) {
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

   // Alpha correction calculation =======================================================================
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

      hTruthRatio = (TH1F*)fCorrection->FindObjectAny("hTruthRatio");
      if (UseExternalSDEff == 0) {
         hTrigEff = (TH1F*)fCorrection->FindObjectAny("hTrigEff");
         hSDFrac = (TH1F*)fCorrection->FindObjectAny("hSDFrac");
         hEmptyEvtCorrection = (TH1F*)fCorrection->FindObjectAny("hEmptyEvtCorrection");
      }
      else {
         hTrigEff = (TH1F*)fCorrectionExternal->FindObjectAny("hTrigEff");
         hSDFrac = (TH1F*)fCorrectionExternal->FindObjectAny("hSDFrac");
         hEmptyEvtCorrection = (TH1F*)fCorrectionExternal->FindObjectAny("hEmptyEvtCorrection");
      }

      TCanvas *cTrigEff = new TCanvas("cTrigEff", "TrigEff", canvasSizeX, canvasSizeY);
      hTrigEff->Draw();

      TCanvas *cSDFrac = new TCanvas("cSDFrac", "SD Fraction After Cut", canvasSizeX, canvasSizeY);
      hSDFrac->Draw();
   } else {
      for (int j=0; j<nVzBin; j++) {
         for (int i=0; i<nEtaBin; i++) {
            // c[i]= new TCanvas (Form("c%d", i), "", canvasSizeX, canvasSizeY);
            // p2->cd(i+1);
            formatHist(alphaPlots[i][j], 2, 1);
            funAlpha[i][j] = new TF1(Form("funAlpha%dVz%d", i, j), "[1]/(x+[3]+0.5)+[2]/(x+0.5)/(x+0.5)+[0]", 0, 100);
            funAlphaErr[i][j] = new TF1(Form("funAlphaErr%dVz%d", i, j), "[0]+[1]/(x+0.5)+[2]*exp([3]*x)", 0, 100);
            double etaMin = i*0.5-3;
            double etaMax = i*0.5-3+0.5;
            alphaPlots[i][j]->Fit(Form("funAlpha%dVz%d", i, j), "m E Q", "", 0, 200);
            alphaPlots[i][j]->SetXTitle("N_{Hits}");
            alphaPlots[i][j]->SetYTitle(Form("#alpha %.1f < #eta < %.1f", etaMin, etaMax));
            alphaPlots[i][j]->SetAxisRange(0.9, 2.0, "Y");
            alphaPlots[i][j]->SetAxisRange(0, 100, "X");
            if (i<2 || i>9) alphaPlots[i][j]->SetAxisRange(4, 10, "y");
            alphaErrPlots[i][j]->Fit(Form("funAlphaErr%dVz%d", i, j), "m E Q", "", 0, 200);
            // alphaPlots[i][j]->Draw("p");
            fAlpha[i][j] = alphaPlots[i][j]->GetFunction(Form("funAlpha%dVz%d", i, j));
            fAlphaErr[i][j] = alphaErrPlots[i][j]->GetFunction(Form("funAlphaErr%dVz%d", i, j));
         }
      }

      // Vertex and event selection efficiency
      TCanvas *cTrigEff = new TCanvas("cTrigEff", "TrigEff", canvasSizeX, canvasSizeY);
      TrackletTree->Project("hTrigEff", Form("%s", multiplicity), TCut(myCut.evtSelection)&&"vz[1]>-99"&&NSDCut, "", nentries, firstentry);
      hTrigEff->Sumw2();
      TrackletTree->Project("hTrigEffNoCut", Form("%s", multiplicity), NSDCut&&"vz[1]>-99", "", nentries, firstentry);
      hTrigEffNoCut->Sumw2();
      hTrigEff->Divide(hTrigEffNoCut);
      hTrigEff->Draw();
      cTrigEff->SaveAs(Form("plot/TrigEff-%s-%d.gif", myPlotTitle, TrackletType));

      // Calculate SD'/IN'
      TCanvas *cSDFrac = new TCanvas("cSDFrac", "SD Fraction After Cut", canvasSizeX, canvasSizeY);
      TrackletTree->Project("hSDFrac", Form("%s", multiplicity), evtSelection&&!(NSDCut), "", nentries, firstentry);
      hSDFrac->Sumw2();
      TrackletTree->Project("hSD", Form("%s", multiplicity), evtSelection, "", nentries, firstentry);
      hSD->Sumw2();
      hSDFrac->Divide(hSD);
      hSDFrac->Draw();
      cSDFrac->SaveAs(Form("plot/SDFrac-%s-%d.gif", myPlotTitle, TrackletType));

      // Calculate Vertexed dN/deta
      TCanvas *cdNdEtaVertexed = new TCanvas("cdNdEtaVertexed", "dNdEta after vertexing", canvasSizeX, canvasSizeY);
      TrackletTree->Project("hDNDEtaVertexed", "eta", TCut("vz[1]>-99")&&(NSDCut), "", nentries, firstentry);
      int nEvtAfterVtx = TrackletTree->Draw("nhit1", TCut("vz[1]>-99")&&(NSDCut), "", nentries, firstentry);
      hDNDEtaVertexed->SetXTitle("#eta Truth after Vtx");
      hDNDEtaVertexed->SetYTitle("dN/#eta");
      hDNDEtaVertexed->Scale(2./nEvtAfterVtx);
      hDNDEtaVertexed->Draw();

      // Calculate NoVertexed dN/deta
      TCanvas *cdNdEtaNoVertexed = new TCanvas("cdNdEtaNoVertexed", "dNdEta after vertexing", canvasSizeX, canvasSizeY);
      int nEvtBeforeVtx = TrackletTree->Draw("nhit1", (NSDCut), "", nentries, firstentry);
      TrackletTree->Project("hDNDEtaNoVertexed", "eta", (NSDCut), "", nentries, firstentry);
      hDNDEtaNoVertexed->SetXTitle("#eta Truth before Vtx");
      hDNDEtaNoVertexed->SetYTitle("dN/#eta");
      hDNDEtaNoVertexed->Scale(2./nEvtBeforeVtx);
      hDNDEtaNoVertexed->Draw();

      // Calculate Xi
      TCanvas *cTriggerCorrection = new TCanvas("cTriggerCorrection", "Xi", canvasSizeX, canvasSizeY);
      int nEvtAfterEvtCut = TrackletTree->Draw("nhit1", evtSelection, "", nentries, firstentry);
      TH1F *hdNdetaWithEvtCut = new TH1F("hdNdetaWithEvtCut", "", nEtaBin, -3, 3);
      TrackletTree->Project("hdNdetaWithEvtCut", "eta", evtSelection, "", nentries, firstentry);
      hdNdetaWithEvtCut->Sumw2();
      hdNdetaWithEvtCut->Scale(1./nEvtAfterEvtCut);
      hTriggerCorrection = (TH1F*)hHadronWOSelection->Project3D("x");
      hTriggerCorrection->SetName("hTriggerCorrection");
      hTriggerCorrection->Sumw2();
      hTriggerCorrection->Scale(1./neventWOSelection);
      hTriggerCorrection->Divide(hdNdetaWithEvtCut);
      hTriggerCorrection->Draw();
      cTriggerCorrection->SaveAs(Form("plot/Xi-%s-%d.gif", myPlotTitle, TrackletType));

      // Calculate Empty Event Correction
      /*
      hEmptyEvtCorrection = new TH1F("hEmptyEvtCorrection", "", 1, 0, 2);
      TH1F *hNonEmptyEvt = new TH1F("hNonEmptyEvt", "", 1, 0, 2);
      TrackletTree->Project("hEmptyEvtCorrection", "1", NSDCut, "", nentries, firstentry);
      TrackletTree->Project("hNonEmptyEvt", "1", NSDCut&&Form("%s>0", multiplicity), "", nentries, firstentry);
      hEmptyEvtCorrection->Sumw2();
      hNonEmptyEvt->Sumw2();
      hEmptyEvtCorrection->Divide(hNonEmptyEvt);
      cout << "EmptyEvent Ratio : " << hEmptyEvtCorrection->GetBinContent(1) << endl;
      */
   }

   // Make beta and alpha plot
   if (makePlot) {
      // Beta plot
      for (int z=1; z<=nVzBin; z++) {
         TLegend * leg1 = new TLegend(0.63, 0.75, 0.93, 0.93);
         leg1->SetFillStyle(0);
         leg1->SetFillColor(0);
         leg1->SetBorderSize(0);
         leg1->SetTextSize(0.04);

         TLegend * leg2 = new TLegend(0.26, 0.75, 0.56, 0.93);
         leg2->SetFillStyle(0);
         leg2->SetFillColor(0);
         leg2->SetBorderSize(0);
         leg2->SetTextSize(0.04);

         TCanvas *cc = new TCanvas(Form("cBetaPlot%d", z), "", canvasSizeX, canvasSizeY);
         for (int x=2; x<nEtaBin; x++) {
            int color = x-1;
            int mtype = 4;
            if (color>nEtaBin/2-1) {
               color = nEtaBin-1-color;
               mtype = 20;
            }
            formatHist(betaPlots[x-1][z-1], color);
            if (betaPlots[x-1][z-1]->GetEntries()) {
               betaPlots[x-1][z-1]->SetMarkerStyle(mtype);
               betaPlots[x-1][z-1]->SetAxisRange(0, 0.5, "Y");
               betaPlots[x-1][z-1]->Draw("L hist same");
               betaPlots[x-1][z-1]->Draw("e same");
               betaPlots[x-1][z-1]->SetYTitle("#beta");
               betaPlots[x-1][z-1]->SetXTitle(Form("N_{Hits} (%.0f #leq V_{z} < %.0f cm)", VzBins[z-1], VzBins[z]));
               if (x-1<nEtaBin/2)
                  leg2->AddEntry(betaPlots[x-1][z-1], Form("%.1f #leq #eta<%.1f", EtaBins[x-1], EtaBins[x]), "pl");
               else
                  leg1->AddEntry(betaPlots[x-1][z-1], Form("%.1f #leq #eta<%.1f", EtaBins[x-1], EtaBins[x]), "pl");
            }
         }
         leg1->Draw();
         leg2->Draw();
         cc->SaveAs(Form("plot/betaPlot/betaPlot-%s-%d-%d.gif", myPlotTitle, z, TrackletType));
         cc->SaveAs(Form("plot/betaPlot/betaPlot-%s-%d-%d.C", myPlotTitle, z, TrackletType));
         cc->SaveAs(Form("plot/betaPlot/betaPlot-%s-%d-%d.eps", myPlotTitle, z, TrackletType));
         cc->Close();
      }

      for (int z=1; z<=nVzBin; z++) {
         TLegend * leg1 = new TLegend(0.63, 0.75, 0.93, 0.93);
         leg1->SetFillStyle(0);
         leg1->SetFillColor(0);
         leg1->SetBorderSize(0);
         leg1->SetTextSize(0.04);

         TLegend * leg2 = new TLegend(0.26, 0.75, 0.56, 0.93);
         leg2->SetFillStyle(0);
         leg2->SetFillColor(0);
         leg2->SetBorderSize(0);
         leg2->SetTextSize(0.04);

         TCanvas *cc = new TCanvas(Form("cAlphaPlot%d", z), "", canvasSizeX, canvasSizeY);
         cc->SetLogy();
         for (int x=2; x<nEtaBin; x++) {
            int color = x-1;
            int mtype = 4;
            if (color>nEtaBin/2-1) {
               color = nEtaBin-1-color;
               mtype = 20;
            }
            formatHist(alphaPlots[x-1][z-1], color);
            if (alphaPlots[x-1][z-1]->GetEntries()) {
               alphaPlots[x-1][z-1]->SetMarkerStyle(mtype);
               alphaPlots[x-1][z-1]->SetYTitle("#alpha");
               alphaPlots[x-1][z-1]->SetXTitle(Form("N_{Hits} (%.0f #leq V_{z} < %.0f cm)", VzBins[z-1], VzBins[z]));
               alphaPlots[x-1][z-1]->SetStats(0);
               alphaPlots[x-1][z-1]->GetYaxis()->SetMoreLogLabels(1);
               alphaPlots[x-1][z-1]->SetAxisRange(0.9, 12, "Y");
               alphaPlots[x-1][z-1]->GetYaxis()->SetNoExponent(1);
               alphaPlots[x-1][z-1]->Draw("L hist same");
               alphaPlots[x-1][z-1]->Draw("hist error p same");
               if (x-1<nEtaBin/2)
                  leg2->AddEntry(alphaPlots[x-1][z-1], Form("%.1f #leq #eta<%.1f", EtaBins[x-1], EtaBins[x]), "pl");
               else
                  leg1->AddEntry(alphaPlots[x-1][z-1], Form("%.1f #leq #eta<%.1f", EtaBins[x-1], EtaBins[x]), "pl");
            }
         }
         leg1->Draw();
         leg2->Draw();
         cc->SaveAs(Form("plot/alphaPlot/alphaPlot-%s-%d-%d.gif", myPlotTitle, z, TrackletType));
         cc->SaveAs(Form("plot/alphaPlot/alphaPlot-%s-%d-%d.C", myPlotTitle, z, TrackletType));
         cc->SaveAs(Form("plot/alphaPlot/alphaPlot-%s-%d-%d.eps", myPlotTitle, z, TrackletType));
         cc->Close();
      }
   }

   // Apply correction ===============================================================================
   ofstream alphaCode;
   alphaCode.open("alpha.C");

   alphaCode << "  const int nEtaBin = 30; " << endl;
   alphaCode << "  const int nTrackletBin = 14; " << endl;
   alphaCode << "  const int nVzBin = 10; " << endl;
   alphaCode << endl;
   alphaCode << "  double TrackletBins[nTrackletBin+1] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100, 200, 700}; " << endl;
   alphaCode << endl;
   alphaCode << "  double EtaBins[nEtaBin+1]; " << endl;
   alphaCode << "  for (int i=0; i<=nEtaBin; i++)" << endl;
   alphaCode << "    EtaBins[i] = (double)i*6.0/(double)nEtaBin-3.0; " << endl;
   alphaCode << "  double VzBins[nVzBin+1]; " << endl;
   alphaCode << "  for (int i=0; i<=nVzBin; i++)" << endl;
   alphaCode << "    VzBins[i] = (double)i*20.0/(double)nVzBin-10.0; " << endl;
   alphaCode << endl;
   alphaCode << "  double alpha[nEtaBin][nTrackletBin][nVzBin]; " << endl;
   alphaCode << endl;
   alphaCode << "  AlphaTracklets12_ = new TH3F(\"hAlphaTracklets12\", \"\", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins); " << endl;
   alphaCode << "  AlphaTracklets12_->SetDirectory(0); " << endl;
   alphaCode << endl;

   // int alpha0flag = 0;

   for (int x=1; x<=nEtaBin; x++) {
      double totalN = 0;
      double totalNErr = 0;
      for (int y=1; y<=nTrackletBin; y++) {
         for (int z=1; z<=nVzBin; z++) {
            if (hAcceptance->GetBinContent(x, z)==0) {
               hHadronAccepted->SetBinContent(x, y, z, 0);
               hHadronAccepted->SetBinError(x, y, z, 0);
               continue;
            }
            double val = hEverything->GetBinContent(x, y, z);
            double beta = betaPlots[x-1][z-1]->GetBinContent(y);
            double betaErr = betaPlots[x-1][z-1]->GetBinError(y);

            int y1 = y;
            // if (y1>10) y1 = 10;
            if (y1>19) y1 = 19;
            double alpha = alphaPlots[x-1][z-1]->GetBinContent(y1);
            double alphaErr = alphaPlots[x-1][z-1]->GetBinError(y1);

            double betaMC = betaMCPlots[x-1][z-1]->GetBinContent(y);
            if (doBetaCorrection) {
               double a = hAlphaA->GetBinContent(x, y, z);
               double b = hAlphaB->GetBinContent(x, y, z);
               if (fabs(b)>20) b = 1;
               double A = a + b * beta;
               // double A = 0.9479+2.655*beta;
               double B = a + b * betaMC;
               // double B = 0.974407+1.77653*betaMC;
               cout << A << " " << B << " " << endl;
            }
            if (doAcceptanceCorrection) {
               double accData = hDataAcc->GetBinContent(x, z);
               double accMC = hMCAcc->GetBinContent(x, z);
               if (accData==0 && accMC==0) {
                  cout << "Error in Acceptance Correction!!!! " << x << " " << z << endl;
               } else {
                  alpha = alpha / accData * accMC;
                  alphaErr = alphaErr / accData * accMC;
                  // cout << accMC/accData << endl;
               }
            }
            alphaCode << "  AlphaTracklets12_->SetBinContent(" << x << ", " << y << ", " << z << ", " << alpha << "); " << endl;

            // Use extrapolated value if alpha is not available
            if (alpha==0 && fAlpha[x-1][z-1]!=0) alpha = fAlpha[x-1][z-1]->Eval(TrackletBins[y]);
            if (alpha==0) {
               for (int k=0; k<y; k++) {
                  alpha = alphaPlots[x-1][z-1]->GetBinContent(y-k);
                  alphaErr = alphaPlots[x-1][z-1]->GetBinError(y-k);
                  if (alpha!=0) break;
               }
               if (alpha==0) alpha = 1;
               // cout << (EtaBins[x]+EtaBins[x-1])/2 << " " << (TrackletBins[y]+TrackletBins[y-1])/2 << " " << (VzBins[z]+VzBins[z-1])/2 << "Used " << alpha << " " << endl;
               // cout << "Empty!!!" << endl;
            }

            double nCorrected = val*(1-beta)*alpha;
            hSubtracted->SetBinContent(x, y, z, nCorrected/alpha);
            hCorrected->SetBinContent(x, y, z, nCorrected);
            if (verbose) cout << "apply: " << x << " " << y << " " << z << " " << nCorrected << " " << alpha << " " << beta << endl;
            double valErr = sqrt(alpha*(1-beta)*alpha*hEverything->GetBinContent(x, y, z)*(1-beta) +
                            betaErr*betaErr*alpha*hEverything->GetBinContent(x, y, z)*alpha*hEverything->GetBinContent(x, y, z));
            hCorrected->SetBinError(x, y, z, valErr);
            correction->Fill((EtaBins[x]+EtaBins[x-1])/2, (TrackletBins[y]+TrackletBins[y-1])/2, (VzBins[z]+VzBins[z-1])/2, alpha, alphaErr, beta, betaErr, nCorrected, hHadronAccepted->GetBinContent(x, y, z), valErr, betaMC);
            totalN += nCorrected;
            totalNErr += valErr*valErr;
            if (verbose) cout << x << " " << y << " " << z << " " << nCorrected << " " << valErr << endl;

            hEtaVzRatio->Fill(x, z, nCorrected-hHadronAccepted->GetBinContent(x, y, z));
            hEtaHitRatio->Fill(x, y, nCorrected-hHadronAccepted->GetBinContent(x, y, z));

            if (fabs(valErr)/nCorrected>-1) {
               fprintf(logFile, " beta (%3.1f, %3.1f, %3.1f) = %.2f +- %.2f", EtaBins[x], TrackletBins[y], VzBins[z], beta, betaErr);
               fprintf(logFile, " alpha = %.2f +- %.2f ", alpha, alphaErr);
               fprintf(logFile, " Value = %.1f +- %.1f Truth = %.1f \n", nCorrected, valErr, (double)hHadronAccepted->GetBinContent(x, y, z));
            }
         }
      }
      hCorrectedEtaBin->SetBinContent(x, totalN);
      hCorrectedEtaBin->SetBinError(x, sqrt(totalNErr));
   }
   alphaCode.close();

   // Plot RawTracklet and Background Tracklet in nTracklet bin
   TCanvas *cRawTrackletnTracklet = new TCanvas("cRawTrackletnTracklet", "", canvasSizeX, canvasSizeY);
   TH1F *hMCTruthnTracklet = (TH1F*)hHadronAccepted->Project3D("y");
   hMCTruthnTracklet->Scale(1./nevent);
   hMCTruthnTracklet->SetXTitle("N_{Hit1}|#eta|<1");
   if (isMC) hMCTruthnTracklet->Draw("hist");

   TH1F *hRawTrackletnTracklet = (TH1F*)hEverything->Project3D("y");
   hRawTrackletnTracklet->Scale(1./nevent);
   hRawTrackletnTracklet->Draw("same");

   TH1F *hBackgroundTrackletnTracklet = (TH1F*)hReproducedBackground->Project3D("y");
   hBackgroundTrackletnTracklet->SetLineColor(2);
   hBackgroundTrackletnTracklet->SetMarkerColor(2);
   hBackgroundTrackletnTracklet->Scale(1./nevent);
   hBackgroundTrackletnTracklet->Draw("same");

   TH1F *hRawMinusBackgroundTrackletnTracklet = (TH1F*)hRawTrackletnTracklet->Clone();
   hRawMinusBackgroundTrackletnTracklet->SetLineColor(4);
   hRawMinusBackgroundTrackletnTracklet->SetMarkerColor(4);
   hRawMinusBackgroundTrackletnTracklet->Add(hBackgroundTrackletnTracklet, -1);
   hRawMinusBackgroundTrackletnTracklet->Draw("same");
   cRawTrackletnTracklet->Update();

   // Plot RawTracklet and Background Tracklet in Vz bin
   TCanvas *cRawTrackletVz = new TCanvas("cRawTrackletVz", "", canvasSizeX, canvasSizeY);
   TH1F *hMCTruthVz = (TH1F*)hHadronAccepted->Project3D("z");
   hMCTruthVz->Scale(1./nevent);
   hMCTruthVz->SetXTitle("V_{z}");
   hMCTruthVz->Draw("hist");

   TH1F *hRawTrackletVz = (TH1F*)hEverything->Project3D("z");
   hRawTrackletVz->Scale(1./nevent);
   hRawTrackletVz->Draw("same");

   TH1F *hBackgroundTrackletVz = (TH1F*)hReproducedBackground->Project3D("z");
   hBackgroundTrackletVz->SetLineColor(2);
   hBackgroundTrackletVz->SetMarkerColor(2);
   hBackgroundTrackletVz->Scale(1./nevent);
   hBackgroundTrackletVz->Draw("same");

   TH1F *hRawMinusBackgroundTrackletVz = (TH1F*)hRawTrackletVz->Clone();
   hRawMinusBackgroundTrackletVz->SetLineColor(4);
   hRawMinusBackgroundTrackletVz->SetMarkerColor(4);
   hRawMinusBackgroundTrackletVz->Add(hBackgroundTrackletVz, -1);
   hRawMinusBackgroundTrackletVz->Draw("same");
   cRawTrackletVz->Update();

   // Plot RawTracklet and Background Tracklet in eta bin
   TCanvas *cRawTrackletEta = new TCanvas("cRawTrackletEta", "", canvasSizeX, canvasSizeY);
   TH1F *hMCTruthEta = (TH1F*)hHadronAccepted->Project3D("x");
   hMCTruthEta->Scale(2./nevent);
   hMCTruthEta->SetXTitle("#eta");
   hMCTruthEta->Draw("hist");

   TH1F *hRawTrackletEta = (TH1F*)hEverything->Project3D("x");
   hRawTrackletEta->Scale(2./nevent);
   hRawTrackletEta->Draw("same");

   TH1F *hBackgroundTrackletEta = (TH1F*)hReproducedBackground->Project3D("x");
   hBackgroundTrackletEta->SetLineColor(2);
   hBackgroundTrackletEta->SetMarkerColor(2);
   hBackgroundTrackletEta->Scale(2./nevent);
   hBackgroundTrackletEta->Draw("same");

   TH1F *hRawMinusBackgroundTrackletEta = (TH1F*)hRawTrackletEta->Clone();
   hRawMinusBackgroundTrackletEta->SetLineColor(4);
   hRawMinusBackgroundTrackletEta->SetMarkerColor(4);
   hRawMinusBackgroundTrackletEta->Add(hBackgroundTrackletEta, -1);
   hRawMinusBackgroundTrackletEta->Draw("same");
   cRawTrackletEta->Update();

   // Plot Before Acceptance correction
   TCanvas *cDNdEtaC = new TCanvas("cDNdEtaC", "Before Acceptance correction", canvasSizeX, canvasSizeY);
   TH1F *hTruthAccepted = (TH1F*)hHadronAccepted->Project3D("x");
   hTruthAccepted->SetName("hTruthAccepted");

   TH1F *hTruth = (TH1F*)hTruthAccepted->Clone();
   hTruth->SetName("hTruth");
   hTruth->Sumw2();
   formatHist(hTruth, 1, nevent/nEtaBin*6);
   hTruth->Divide(hAcceptance->ProjectionX());  // calibrate the acceptance

   TH1F *hTruthEvtCutCorrectedByXi = (TH1F*)hHadron->Project3D("x");
   hTruthEvtCutCorrectedByXi->SetName("hTruthEvtCutCorrectedByXi");
   hTruthEvtCutCorrectedByXi->Sumw2();
   formatHist(hTruthEvtCutCorrectedByXi, 1, nevent/nEtaBin*6);

   formatHist(hCorrectedEtaBin, 2, nevent/nEtaBin*6, 1.1);
   hCorrectedEtaBin->Draw("e");

   hTruthAccepted->Sumw2();
   formatHist(hTruthAccepted, 1, nevent/nEtaBin*6);
   hTruthAccepted->SetAxisRange(0, dndetaRange, "y");
   hTruthAccepted->SetXTitle("#eta (Calculated Hand)");
   hTruthAccepted->SetYTitle("dN/d#eta");
   hTruthAccepted->Draw("hist same");
   cDNdEtaC->Update();

   // Final dN/deta results =======================================================================================================
   TH1F *hMeasured = (TH1F*)hCorrected->Project3D("x");
   hMeasured->Sumw2();
   formatHist(hMeasured, 4, nevent/nEtaBin*6, 1.1);
   hMeasured->Divide(hAcceptance->ProjectionX());  // calibrate the acceptance

   TH1F *hTruthWOSelection = (TH1F*)hHadronWOSelection->Project3D("x");
   hTruthWOSelection->SetName("hTruthWOSelection");
   formatHist(hTruthWOSelection, 2, neventWOSelection/nEtaBin*6);

   // Different calculation
   TH2F *hMeasuredEtanTracklet = new TH2F("hMeasuredEtanTracklet", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins);
   TH2F *hTruthEtanTracklet = new TH2F("hTruthEtanTracklet", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins);
   for (int x=1; x<=nEtaBin; x++) {
      for (int y=1; y<=nTrackletBin; y++) {
         double total = 0, totalErr = 0;
         double totalMC = 0, totalMCErr = 0;
         for (int z=1; z<=nVzBin; z++) {
            if (hAcceptance->GetBinContent(x, z)!=0) {
               total += hCorrected->GetBinContent(x, y, z);
               double err = hCorrected->GetBinError(x, y, z);
               totalErr = sqrt(totalErr*totalErr+err*err);
               totalMC += hHadronAccepted->GetBinContent(x, y, z);
               double errMC = hHadronAccepted->GetBinError(x, y, z);
               totalMCErr = sqrt(totalMCErr*totalMCErr+errMC*errMC);
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
      if (!doTriggerCorrection) {
         SDFrac = 0;
         TrigEff = 1;
      }
      if (TrigEff!=0) nEvt += hnTracklet->GetBinContent(y)/TrigEff*(1-SDFrac*SDFactor);
      for (int x=1; x<=nEtaBin; x++) {
         for (int z=1; z<=nVzBin; z++) {
            if (hAcceptance->GetBinContent(x, z)!=0) {
               if (TrigEff!=0) {
                  double tmp = hnTracklet3D->GetBinContent(x, y, z);
                  hnTracklet3D->SetBinContent(x, y, z, tmp/TrigEff*(1-SDFrac*SDFactor));
               }
            } else {
               hnTracklet3D->SetBinContent(x, y, z, 0);
            }
         }
      }
   }

   for (int x=1; x<=nEtaBin; x++) {
      for (int y=1; y<=nTrackletBin; y++) {
         double SDFrac = hSDFrac->GetBinContent(y);
         double TrigEff = hTrigEff->GetBinContent(y);
         if (!doTriggerCorrection) {
            SDFrac = 0;
            TrigEff = 1;
         }
         if (TrigEff!=0) {
            hMeasuredEtanTracklet->SetBinContent(x, y, hMeasuredEtanTracklet->GetBinContent(x, y)/TrigEff*(1-SDFrac*SDFactor));
            hMeasuredEtanTracklet->SetBinError(x, y, hMeasuredEtanTracklet->GetBinError(x, y)/TrigEff*(1-SDFrac*SDFactor));
            hTruthEtanTracklet->SetBinContent(x, y, hTruthEtanTracklet->GetBinContent(x, y)/TrigEff*(1-SDFrac*SDFactor));
            hTruthEtanTracklet->SetBinError(x, y, hTruthEtanTracklet->GetBinError(x, y)/TrigEff*(1-SDFrac*SDFactor));
         }
      }
   }

   TH1F *hAcceptance2 = (TH1F*)hnTracklet3D->Project3D("x");
   hAcceptance2->SetName("hAcceptance2");

   TH1F *hMeasuredTrigEffCorrected = (TH1F*)hMeasuredEtanTracklet->ProjectionX();
   hMeasuredTrigEffCorrected->SetName("hMeasuredTrigEffCorrected");
   formatHist(hMeasuredTrigEffCorrected, 2, 1./nEtaBin*6, 1.1);
   hMeasuredTrigEffCorrected->Divide(hAcceptance2);  // calibrate the acceptance
   hMeasuredTrigEffCorrected->SetMarkerStyle(4);

   TCanvas *cEmpty = new TCanvas("c", "", canvasSizeX, canvasSizeY);
   if (!useCorrectionFile) {
      hEmptyEvtCorrection = (TH1F*)hTruthWOSelection->Clone();
      hEmptyEvtCorrection->SetName("hEmptyEvtCorrection");
      hEmptyEvtCorrection->Divide(hMeasuredTrigEffCorrected);
      for (int x=1; x<=nEtaBin; x++)
         hEmptyEvtCorrection->SetBinError(x, 0);
   } else {
      hEmptyEvtCorrection = (TH1F*)fCorrection->FindObjectAny("hEmptyEvtCorrection");
   }
   hEmptyEvtCorrection->Draw();

   TCanvas *cDNdEta = new TCanvas("cDNdEta", "Final result", canvasSizeX, canvasSizeY);
   hTruthWOSelection->SetAxisRange(0, dndetaRange, "y");
   hTruthWOSelection->SetXTitle("#eta");
   hTruthWOSelection->SetYTitle("dN/d#eta");
   hTruthWOSelection->Draw("hist");

   TH1F *hMeasuredFinal = (TH1F*) hMeasuredTrigEffCorrected->Clone();
   hMeasuredFinal->SetName("hMeasuredFinal");
   hMeasuredFinal->SetMarkerStyle(20);
   for (int x=1; x<nEtaBin; x++) {
      double emptyCorrection = (hEmptyEvtCorrection->GetBinContent(nEtaBin/2)+hEmptyEvtCorrection->GetBinContent(nEtaBin/2+1))/2;
      if (!doTriggerCorrection) emptyCorrection = 1;
      double val = hMeasuredFinal->GetBinContent(x)*emptyCorrection;
      double valerr = hMeasuredFinal->GetBinError(x)*emptyCorrection;
      hMeasuredFinal->SetBinContent(x, val);
      hMeasuredFinal->SetBinError(x, valerr);
   }
   formatHist(hMeasuredFinal, 2, 1, 1.1);
   hMeasuredFinal->Draw("e same");

   TH1F *hMeasuredStat = (TH1F*)hMeasuredFinal->Clone();
   hMeasuredStat->SetName("hMeasuredStat");
   hMeasuredStat->SetMarkerSize(0);
   hMeasuredStat->Draw("e same");

   hTruthEvtCutCorrectedByXi->Multiply(hTriggerCorrection);
   hTruthEvtCutCorrectedByXi->SetAxisRange(0, dndetaRange, "y");
   hTruthEvtCutCorrectedByXi->SetXTitle("#eta");
   hTruthEvtCutCorrectedByXi->SetYTitle("dN/d#eta");
   // hTruthEvtCutCorrectedByXi->Draw("hist");

   // double systematicError10TeV[12] = {0.132, 0.132, 0.093, 0.085, 0.084, 0.083, 0.083, 0.084, 0.085, 0.093, 0.132, 0.132};
   // TGraph *gErrorBand = GetErrorBand((TH1F*)hMeasuredFinal, systematicError10TeV, systematicError10TeV, 0.25);
   // gErrorBand->Draw("F");

   TH1F* hUA5(0);
   if (putUA5) {
      hUA5 = getUA5NSD();
      hUA5->Draw("p same");
   }

   TH1F* hPYTHIA(0);
   // if (putPYTHIA) {
   //    hPYTHIA = getPythiaD6T();
   //    hPYTHIA->Draw("hist same");
   // }

   TLegend* leg1 = new TLegend(0.30, 0.18, 1, 0.3);
   leg1->SetFillStyle(0);
   leg1->SetFillColor(0);
   leg1->SetBorderSize(0);
   // leg1->SetTextSize(0.045);
   leg1->AddEntry(hTruth, Form("%s", myPlotTitle), "");
   leg1->AddEntry(hMeasuredFinal, "Reconstructed", "pl");
   // if (putPYTHIA) leg1->AddEntry(hPYTHIA, "PYTHIA 900GeV D6T", "l");
   leg1->Draw();
   // TText* tex = new TText(-2.6, 5, "CMS Preliminary");
   // tex->Draw();
   cDNdEta->SaveAs(Form("plot/result/result-%s-%d.gif", myPlotTitle, TrackletType));
   cDNdEta->SaveAs(Form("plot/result/result-%s-%d.eps", myPlotTitle, TrackletType));
   cDNdEta->SaveAs(Form("plot/result/result-%s-%d.C", myPlotTitle, TrackletType));

   // Compare with Truth ===========================================================================================
   TCanvas *cDNdEtaCompare = new TCanvas("cDNdEtaCompare", "Compare", canvasSizeX, canvasSizeY);
   hTruthWOSelection->Draw("hist");

   TH1F *hMeasuredNoCorrection = (TH1F*) hMeasured->Clone();
   hMeasuredNoCorrection->SetName("hMeasuredNoCorrection");
   hMeasuredNoCorrection->SetMarkerStyle(4);
   hMeasuredNoCorrection->SetMarkerColor(1);
   hMeasuredNoCorrection->SetLineColor(1);
   hMeasuredNoCorrection->Draw("e same");

   // hTruthEvtCutCorrectedByXi->Draw("hist same");

   hMeasured->Multiply(hTriggerCorrection);
   // hMeasured->Draw("e same");

   hMeasuredTrigEffCorrected->Draw("e same");

   TH1F *hTruthTrigEffCorrected = (TH1F*)hTruthEtanTracklet->ProjectionX();
   hTruthTrigEffCorrected->SetName("hTruthTrigEffCorrected");
   formatHist(hTruthTrigEffCorrected, 4, 1./nEtaBin*6, 1.1);
   hTruthTrigEffCorrected->Divide(hAcceptance2);  // calibrate the acceptance
   hTruthTrigEffCorrected->Draw("e same");

   // hTruth->Multiply(hTriggerCorrection);
   hTruth->SetLineColor(6);
   hTruth->Draw("hist same");

   hDNDEtaVertexed->Draw("hist same");

   hMeasuredFinal->Draw("e same");

   if (putUA5) hUA5->Draw("p same");

   // if (putPYTHIA) hPYTHIA->Draw("hist same");

   TLegend * leg11 = new TLegend(0.20, 0.18, 1, 0.5);
   leg11->SetFillStyle(0);
   leg11->SetFillColor(0);
   leg11->SetBorderSize(0);
   leg11->AddEntry(hTruth, Form("%s", myPlotTitle), "");
   leg11->AddEntry(hTruthWOSelection, "MC Truth W/O selection", "l");
   leg11->AddEntry(hTruth, "MC Truth hTruth", "l");
   leg11->AddEntry(hMeasuredNoCorrection, "No trig correction", "pl");
   // leg11->AddEntry(hTruthEvtCutCorrectedByXi, "MC Truth accepted corrected by xi", "l");
   // leg11->AddEntry(hMeasured, "Xi correction", "pl");
   leg11->AddEntry(hDNDEtaVertexed, "MC Truth with vtx", "l");
   leg11->AddEntry(hTruthTrigEffCorrected, "MC Truth TrigEff corrected", "pl");
   leg11->AddEntry(hMeasuredTrigEffCorrected, "TrigEff corrected", "pl");
   leg11->AddEntry(hMeasuredFinal, "TrigEff+Empty correction", "pl");
   leg11->Draw();
   cDNdEtaCompare->Update();

   if (isMC) {
      // Ratio between measured and truth ===========================================================================================
      TCanvas *cRatio = new TCanvas("cRatio", "Ratio", canvasSizeX, canvasSizeY);
      TH1F *hRatio = (TH1F*)hMeasuredFinal->Clone();
      hRatio->SetName("hRatio");
      hRatio->Divide(hTruthWOSelection);
      hRatio->SetXTitle("#eta");
      hRatio->SetYTitle("Ratio");
      hRatio->SetAxisRange(0.8, 1.2, "y");
      hRatio->Draw();

      TH1F* hRatioStat = (TH1F*)hRatio->Clone();
      hRatioStat->SetName("hRatioStat");
      hRatioStat->SetMarkerSize(0);
      hRatioStat->Draw("same");

      TLegend* leg2 = new TLegend(0.18, 0.18, 1, 0.3);
      leg2->SetFillStyle(0);
      leg2->SetFillColor(0);
      leg2->SetBorderSize(0);
      //leg2->SetTextSize(0.045);
      leg2->AddEntry(hTruth, myPlotTitle, "");
      leg2->AddEntry(hRatio, "Reconstructed / MC Truth", "pl");
      leg2->Draw();

      TLine* l1 = new TLine(-3, 1, 3, 1);
      l1->Draw();
      // TText* tex2 = new TText(-2.6, 1.165, "CMS Preliminary");
      // tex2->Draw();
      cRatio->SaveAs(Form("plot/ratio/ratio-%s-%d.gif", myPlotTitle, TrackletType));
      cRatio->SaveAs(Form("plot/ratio/ratio-%s-%d.eps", myPlotTitle, TrackletType));
      cRatio->SaveAs(Form("plot/ratio/ratio-%s-%d.C", myPlotTitle, TrackletType));

      TCanvas *cDNdnTracklet = new TCanvas("cDNdnTracklet", "Measured vs mult", canvasSizeX, canvasSizeY);
      TH1F *hTruthHit = (TH1F*)hHadronAccepted->Project3D("y");
      hTruthHit->Sumw2();
      formatHist(hTruthHit, 1, nevent);
      hTruthHit->SetAxisRange(0, dndetaRange, "y");
      hTruthHit->SetXTitle("N_{Hit1} |#eta|<1");
      hTruthHit->SetYTitle("dN/dN_{Hit1}");
      hTruthHit->Draw("hist");

      TH1F *hMeasuredHit = (TH1F*)hCorrected->Project3D("y");
      hMeasuredHit->Sumw2();
      formatHist(hMeasuredHit, 2, nevent, 1.1);
      hMeasuredHit->Draw("e same");
      cDNdnTracklet->Update();

      TCanvas *cRatioTracklet = new TCanvas("cRatioTracklet", "Ratio vs mult", canvasSizeX, canvasSizeY);
      TH1F *hRatioHit = (TH1F*)hMeasuredHit->Clone();
      hRatioHit->Divide(hTruthHit);
      hRatioHit->SetXTitle("#eta");
      hRatioHit->SetYTitle("Ratio");
      hRatioHit->SetAxisRange(0.7, 1.3, "y");
      hRatioHit->Draw();

      TLine *l2 = new TLine(0, 1, nTrackletBin, 1);
      l2->Draw();
      cRatioTracklet->Update();
   }

   // hCorrected->Scale(1./nevent);
   // hHadron->Scale(1./nevent);

   for (int i=0; i<nEtaBin; i++) {
      for (int j=0; j<nVzBin; j++) {
         if (fAlpha[i][j]) fAlpha[i][j]->Write();
         // Don't save beta anymore.
         // fBeta[i]->Write();
         // fBetaErr[i]->Write();
         if (fAlphaErr[i][j]) fAlphaErr[i][j]->Write();
         if (verbose && !fAlpha[i][j]) cout << "no alpha! " << i << " " << j << endl;
         if (verbose && !fAlphaErr[i][j]) cout << "no alphaErr! " << i << " " << j << endl;
      }
   }

   cDNdEta->Write();
   outf->Write();
   fclose(logFile);

   // if (alpha0flag==1) cout << "Some bins are lost!!!" << endl;
   return 0;
}

//===========================================================================
// Format Histogram
//===========================================================================
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

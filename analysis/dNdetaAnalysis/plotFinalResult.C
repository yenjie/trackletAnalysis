#define _CANVAS_W 600
#define _CANVAS_H 600
#define _PLOT_RANGE 30.0
#define _SD_FACTOR 1

#define _ENERGY 8

// standard library
#include <fstream>
#include <string>

// ROOT library
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TCut.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TText.h"

// For plotting
#include "error_bands.h"

// External histograms
#include "distributions.h"

void formatHist(TH1* h, int col = 1, double norm = 1, double msize = 1);

// ============================================================================
// Main Routine
// ============================================================================
int plotFinalResult(int TrackletType,
                    std::string file_name,
                    const char* title,           // plot title
                    bool apply_correction,       // apply external corrections
                    const char* corr_name,       // correction file name
                    bool apply_accep_corr = 1,   // apply acceptance correction
                    bool verbose = 0,            // enable debug output
                    bool apply_ext_accep = 1,    // use pre-defined acceptance
                    int selection = 1,           // dn/deta selection in MC
                    int mult_selection = 0,      // multiplicity handle
                    bool use_wide_dphi = 0,      // use wider dphi region
                    bool apply_trigger_corr = 1) // do trigger eff correction
{
   std::vector<std::string> file_list;
   if (file_name.substr(file_name.find_last_of(".") + 1) == "root") {
      file_list.push_back(file_name);
   } else {
      std::ifstream file_stream(file_name);
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
      printf("$ Monte Carlo analysis\n");
      apply_accep_corr = 0;
   } else {
      printf("$ data analysis\n");
   }
   printf(" # at %iTeV\n", _ENERGY);

   if (!apply_trigger_corr)
      printf(" # no trigger correction applied\n");

   // Choose multiplicity handle
   const char* multiplicity;
   if (mult_selection == 1) {
      multiplicity = "nhit1_cut";
      printf("$ event multiplicity handle: number of clusters with cluster size cut\n");
   } else if (mult_selection == 2) {
      multiplicity = "nTracklets";
      printf("$ event multiplicity handle: number of tracklets\n");
   } else {
      multiplicity = "mult";
      printf("$ event multiplicity handle: number of background-subtracted tracklets\n");
   }

   TH1::SetDefaultSumw2();

   // Read alpha, beta, geometry correction from file.
   TFile* fCorrection = 0;
   if (apply_correction) {
      const char* correctionfname = Form("correction/correction-%d-%s.root", TrackletType, corr_name);
      fCorrection = new TFile(correctionfname);
      printf("$ correction file: %s\n", correctionfname);
   }

   TFile* fAcceptance = 0;
   if (apply_correction && apply_accep_corr) {
      const char* acceptancefname = Form("correction/acceptance-%d.root", TrackletType);
      fAcceptance = new TFile(acceptancefname);
      printf("$ acceptance file: %s\n", acceptancefname);
   }

   // Definition of Vz, Eta, Hit bins
   const int nTrackletBin = 20;
   double TrackletBins[nTrackletBin+1] = {-10, 8, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 66, 72, 80, 90, 100, 110, 130, 150, 300};

   const int nEtaBin = 30;
   double EtaBins[nEtaBin+1];
   for (int i=0; i<=nEtaBin; i++)
      EtaBins[i] = i*6.0/nEtaBin-3.0;

   const int nVzBin = 14;
   double VzBins[nVzBin+1];
   double VzRangeH = 15;
   double VzRangeL = -13;
   for (int i=0; i<=nVzBin; i++)
      VzBins[i] = i*(VzRangeH-VzRangeL)/nVzBin+VzRangeL;

   // Signal and Sideband regions =============================================
   double signal_region = 1.0;   // dphi cut for signal region
   double sideband_region = 2.0; // dphi cut for sideband region
   double deta_cut = 0.1;        // deta cut

   if (use_wide_dphi) {
      signal_region = 1.5;
      sideband_region = 3.0;
   }

   TCut signal_region_cut   = Form("abs(dphi)<%f && abs(deta)<%f", signal_region, deta_cut);
   TCut sideband_region_cut = Form("abs(dphi)>%f && abs(dphi)<%f && abs(deta)<%f", signal_region, sideband_region, deta_cut);

   TCut vertex_selection = "(vz[1]<15 && vz[1]>-13)";

   TCut offline_selection;
   TCut gen_selection;
   switch (selection) {
      case 0:
         offline_selection = "(1)";
         gen_selection = "(evtType!=102)";
         printf("$ INELASTIC definition\n");
         break;
      case 1:
         // HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1
         offline_selection = "(passHLT)";
         gen_selection = "(evtType!=102&&evtType!=103&&evtType!=104)";
         printf("$ NSD definition\n");
         break;
   }
   if (!isMC) gen_selection = "1";

   TCut event_selection = vertex_selection && offline_selection;
   TCut gen_event_selection = event_selection && gen_selection;

   // Output file =============================================================
   TFile* outf = new TFile(Form("correction-%i-%s.root", TrackletType, title), "recreate");

   TNtuple* betas = new TNtuple("betas", "", "eta:ntl:vz:beta:betaErr");
   TNtuple* alphas = new TNtuple("alphas", "", "eta:ntl:vz:alpha:alphaErr");
   TH3F* hbeta = new TH3F("hbeta", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* halpha = new TH3F("halpha", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* halpha_applied = new TH3F("halpha_applied", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);

   TH3F* hEverything = new TH3F("hEverything", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* hReproducedBackground = new TH3F("hReproducedBackground", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* hSubtracted = new TH3F("hSubtracted", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* hCorrected = new TH3F("hCorrected", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);

   TH3F* hHadron = new TH3F("hHadron", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* hHadronAccepted = new TH3F("hHadronAccepted", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);
   TH3F* hHadronWOSelection = new TH3F("hHadronWOSelection", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);

   TH2F* hAcceptance1D = 0;
   if (apply_ext_accep)
      hAcceptance1D = get_acceptance(TrackletType, _ENERGY);
   else if (apply_correction)
      hAcceptance1D = (TH2F*)fCorrection->FindObjectAny("hAcceptance1D");
   else
      hAcceptance1D = new TH2F("hAcceptance1D", "", nEtaBin, EtaBins, nVzBin, VzBins);
   TH3F* hAcceptance2D = new TH3F("hAcceptance2D", "", nEtaBin, EtaBins, nTrackletBin, TrackletBins, nVzBin, VzBins);

   TH1F* hVz = new TH1F("hVz", "", nVzBin, VzBins);
   TH1F* hnTracklet = new TH1F("hnTracklet", "", nTrackletBin, TrackletBins);
   TH2F* hVzNTracklet = new TH2F("hVzNTracklet", "", nVzBin, VzBins, nTrackletBin, TrackletBins);

   TH1F* hCorrectedEtaBin = new TH1F("hCorrectedEtaBin", "Corrected", nEtaBin, EtaBins);
   TH1F* hDNDEtaVertexed = new TH1F("hDNDEtaVertexed", "", nEtaBin, EtaBins);
   TH1F* hDNDEtaNoVertexed = new TH1F("hDNDEtaNoVertexed", "", nEtaBin, EtaBins);
   TH1F* hdNdetaWithEvtCut = new TH1F("hdNdetaWithEvtCut", "", nEtaBin, EtaBins);
   TH1F* hTrigEff = 0;
   if (apply_correction)
      hTrigEff = (TH1F*)fCorrection->Get("hTrigEff");
   else
      hTrigEff = new TH1F("hTrigEff", "", nTrackletBin, TrackletBins);
   TH1F* hTrigEffNoCut = new TH1F("hTrigEffNoCut", "", nTrackletBin, TrackletBins);
   TH1F* hSD = new TH1F("hSD", "", nTrackletBin, TrackletBins);
   TH1F* hSDFrac = new TH1F("hSDFrac", "", nTrackletBin, TrackletBins);

   TH1F* alphaPlots[nEtaBin][nVzBin];
   TH1F* betaPlots[nEtaBin][nVzBin];
   TH1F* alphaErrPlots[nEtaBin][nVzBin];
   TH1F* betaErrPlots[nEtaBin][nVzBin];

   // Prepare histograms ======================================================
   for (int i=0; i<nEtaBin; i++) {
      for (int j=0; j<nVzBin; j++) {
         alphaPlots[i][j] = new TH1F(Form("alpha%dVz%d", i, j), "", nTrackletBin, TrackletBins);
         alphaErrPlots[i][j] = new TH1F(Form("alphaErr%dVz%d", i, j), "", nTrackletBin, TrackletBins);
         betaPlots[i][j] = new TH1F(Form("beta%dVz%d", i, j), "", nTrackletBin, TrackletBins);
         betaErrPlots[i][j] = new TH1F(Form("betaErr%dVz%d", i, j), "", nTrackletBin, TrackletBins);
      }
   }

   // Fit functions of Beta and Alpha =========================================
   TF1* funAlpha[nEtaBin][nVzBin];
   TF1* funAlphaErr[nEtaBin][nVzBin];
   TF1* fAlpha[nEtaBin][nVzBin];
   TF1* fAlphaErr[nEtaBin][nVzBin];

   TH1F* hTriggerCorrection;
   TH1F* hEmptyEvtCorrection;

   // Number of events
   TH1F* hnevent = new TH1F("hnevent", "", nVzBin, VzBins);
   int nevententries = TrackletTree->Draw("vz[1]>>hnevent", "weight" * gen_event_selection, "goff");
   float nevent = hnevent->Integral(0, hnevent->GetNbinsX() + 1);

   TH1F* hneventWOSelection = new TH1F("hneventWOSelection", "", nVzBin, VzBins);
   TrackletTree->Draw("vz[1]>>hneventWOSelection", "weight" * gen_selection, "goff");
   float neventWOSelection = hneventWOSelection->Integral(0, hneventWOSelection->GetNbinsX() + 1);

   printf("$ weighted events: %f, entries: %i\n", nevent, nevententries);
   if (nevententries < 1) {
      printf("  ! no events selected - stopping\n");
      return 0;
   }

   // Acceptance calculation ==================================================
   TrackletTree->Project("hnTracklet", Form("%s", multiplicity), "weight" * event_selection);
   TrackletTree->Project("hVzNTracklet", Form("%s:vz[1]", multiplicity), "weight" * event_selection);
   TrackletTree->Project("hVz", "vz[1]", "weight" * event_selection);

   TCanvas* cVz = new TCanvas("cVz", "Vz distribution", _CANVAS_W, _CANVAS_H);
   hVz->Sumw2();
   hVz->Scale(1./hVz->GetEntries());
   hVz->Fit("gaus");
   hVz->SetXTitle("v_{z} (cm)");
   hVz->SetStats(0);
   hVz->Draw();
   cVz->SaveAs(Form("figs/vz/vz-%s-%d.png", title, TrackletType));

   // Define the acceptance region to avoid large correction factors
   for (int i=1; i<=nEtaBin; i++) {
      for (int j=1; j<=nVzBin; j++) {
         if ((!apply_ext_accep && !apply_correction) || hAcceptance1D->GetBinContent(i, j) != 0) {
            hAcceptance1D->SetBinContent(i, j, hVz->GetBinContent(j));
            hAcceptance1D->SetBinError(i, j, 0);
            for (int k=1; k<=nTrackletBin; k++)
               hAcceptance2D->SetBinContent(i, k, j, hVzNTracklet->GetBinContent(j, k));
         }
      }
   }

   // Charged Hadrons =========================================================
   hHadron->SetXTitle("#eta");
   hHadron->SetYTitle("N_{hit}^{Layer1} |#eta|<3");

   TrackletTree->Project("hHadron", Form("vz[1]:%s:eta", multiplicity), "weight" * (event_selection && "abs(eta)<3"));
   TrackletTree->Project("hHadronWOSelection", Form("vz[1]:%s:eta", multiplicity), "weight" * (gen_selection && "abs(eta)<3"));
   hHadron->Sumw2();
   hHadronWOSelection->Sumw2();

   TrackletTree->Project("hHadronAccepted", Form("vz[1]:%s:eta", multiplicity), "weight" * (event_selection && "abs(eta)<3"));
   hHadronAccepted = (TH3F*)hHadron->Clone();
   hHadronAccepted->SetName("hHadronAccepted");

   // Prepare Tracklet Three-Dimensional Histogram ============================
   TrackletTree->Project("hEverything", Form("vz[1]:%s:eta1", multiplicity), "weight" * (signal_region_cut && event_selection));
   hEverything->Sumw2();

   TrackletTree->Project("hReproducedBackground", Form("vz[1]:%s:eta1", multiplicity), "weight" * (sideband_region_cut && event_selection));
   hReproducedBackground->Sumw2();

   // Read Acceptance =========================================================
   TH2F* haccep_mc = 0;
   TH2F* haccep_data = 0;
   if (apply_accep_corr) {
      haccep_mc = (TH2F*)fAcceptance->FindObjectAny("haccep_mc");
      haccep_data = (TH2F*)fAcceptance->FindObjectAny("haccep_data");
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

            if (hAcceptance1D->GetBinContent(x, z) == 0)
               continue;
            if (hEverything->GetBinContent(x, y, z) != 0) {
               beta = hReproducedBackground->GetBinContent(x, y, z)/hEverything->GetBinContent(x, y, z);
               double e1 = hEverything->GetBinError(x, y, z)/hEverything->GetBinContent(x, y, z);
               double e2 = hReproducedBackground->GetBinError(x, y, z)/hReproducedBackground->GetBinContent(x, y, z);
               double betaErr = beta * sqrt(e1*e1 + e2*e2);
               if (beta/betaErr > -10) {
                  betas->Fill(x, y, z, beta, betaErr);
                  hbeta->SetBinContent(x, y, z, beta);
                  hbeta->SetBinError(x, y, z, betaErr);
                  betaPlots[x-1][z-1]->SetBinContent(y, beta);
                  betaPlots[x-1][z-1]->SetBinError(y, betaErr);
                  betaErrPlots[x-1][z-1]->SetBinContent(y, betaErr);
               } else {
                  if (verbose) printf("   |  ! warning: beta not used: %f, eta: %i, ntl: %i, vz: %i, betaErr: %f\n", beta, x, y, z, betaErr);
               }
            }
         }
      }
   }

   for (int i=0; i<nEtaBin; i++) {
      for (int j=0; j<nVzBin; j++) {
         // c[i]= new TCanvas (Form("c%d", i), "", _CANVAS_W, _CANVAS_H);
         // p1->cd(i+1);
         betaPlots[i][j]->SetXTitle("N_{Hits}");
         betaPlots[i][j]->SetYTitle(Form("#beta %.1f < #eta < %.1f", EtaBins[i], EtaBins[i+1]));
         betaPlots[i][j]->SetAxisRange(0, 1, "Y");
         betaPlots[i][j]->SetAxisRange(0, 200, "X");
         formatHist(betaPlots[i][j], 2, 1);
         // betaPlots[i][j]->Draw("p");
      }
   }

   // alpha calculation (efficiency correction) ===============================
   if (!apply_correction) {
      for (int x=1; x<=nEtaBin; x++) {
         for (int y=1; y<=nTrackletBin; y++) {
            for (int z=1; z<=nVzBin; z++) {
               alphaPlots[x-1][z-1]->SetBinContent(y, 0);
               alphaPlots[x-1][z-1]->SetBinError(y, 0);
               alphaErrPlots[x-1][z-1]->SetBinContent(y, 0);

               if (hAcceptance1D->GetBinContent(x, z) == 0)
                  continue;
               if (hEverything->GetBinContent(x, y, z) != 0 && hHadron->GetBinContent(x, y, z) != 0) {
                  double val = hEverything->GetBinContent(x, y, z);
                  double beta = betaPlots[x-1][z-1]->GetBinContent(y);
                  double e1 = hEverything->GetBinError(x, y, z);
                  double e2 = hReproducedBackground->GetBinError(x, y, z);
                  double nsig = val*(1-beta);
                  double valErr = sqrt(e1*e1 + e2*e2);
                  double truth = hHadron->GetBinContent(x, y, z);
                  double truthErr = hHadron->GetBinError(x, y, z);
                  double alpha = truth/nsig;
                  double alphaErr = truth/nsig * sqrt(valErr/nsig*valErr/nsig + truthErr/truth*truthErr/truth);
                  if (verbose) printf("   | calculation - eta: %i, ntl: %i, vz: %i, val: %f, beta: %f, nsig: %f, truth: %f, alpha: %f, alphaErr: %f\n", x, y, z, val, beta, nsig, truth, alpha, alphaErr);

                  if (alpha > 10 || alpha < 0 || (beta > 0.5 && alpha > 5)) {
                     if (TrackletBins[y] < 101) {
                        hAcceptance1D->SetBinContent(x, z, 0);
                        continue;
                     }
                     printf("   |  ! possibly bad acceptance: eta: %i, ntl: %i, vz: %i, val: %f, alpha: %f, alphaErr: %f, beta: %f\n", x, y, z, val, alpha, alphaErr, beta);
                  } else if (alpha > 0 && ((beta != 1 && alpha/alphaErr > 5 && alpha < 5) || (alpha < 2))) {
                     alphas->Fill(x, y, z, alpha, alphaErr);
                     halpha->SetBinContent(x, y, z, alpha);
                     halpha->SetBinError(x, y, z, alphaErr);
                     alphaPlots[x-1][z-1]->SetBinContent(y, alpha);
                     alphaPlots[x-1][z-1]->SetBinError(y, alphaErr);
                     alphaErrPlots[x-1][z-1]->SetBinContent(y, alphaErr);
                  } else {
                     printf("   |  ! warning: calculated alpha not used: %f, eta: %i, ntl: %i, vz: %i, beta: %f, alphaErr: %f\n", alpha, x, y, z, beta, alphaErr);
                  }
               } else if (TrackletBins[y] < 101) {
                  hAcceptance1D->SetBinContent(x, z, 0);
               }
            }
         }
      }
   }

   // Alpha correction calculation ============================================
   if (apply_correction) {
      hTriggerCorrection = (TH1F*)fCorrection->FindObjectAny("hTriggerCorrection");
      hTriggerCorrection->SetName("hTriggerCorrection");

      for (int i=0; i<nEtaBin; i++) {
         for (int j=0; j<nVzBin; j++) {
            alphaPlots[i][j] = (TH1F*)fCorrection->FindObjectAny(Form("alpha%dVz%d", i, j));

            fAlpha[i][j] = (TF1*)fCorrection->FindObjectAny(Form("funAlpha%dVz%d", i, j));
            fAlphaErr[i][j] = (TF1*)fCorrection->FindObjectAny(Form("funAlphaErr%dVz%d", i, j));
         }
      }
   } else {
      for (int j=0; j<nVzBin; j++) {
         for (int i=0; i<nEtaBin; i++) {
            alphaPlots[i][j]->SetXTitle("N_{Hits}");
            alphaPlots[i][j]->SetYTitle(Form("#alpha %.1f < #eta < %.1f", EtaBins[i], EtaBins[i+1]));
            alphaPlots[i][j]->SetAxisRange(0.5, 2.5, "Y");
            alphaPlots[i][j]->SetAxisRange(0, 200, "X");
            if (i < 7 || i > 22) alphaPlots[i][j]->SetAxisRange(2.5, 12.5, "y");
            formatHist(alphaPlots[i][j], 2, 1);

            funAlpha[i][j] = new TF1(Form("funAlpha%dVz%d", i, j), "[1]/(x+[3]+0.5)+[2]/(x+0.5)/(x+0.5)+[0]", 0, 150);
            funAlphaErr[i][j] = new TF1(Form("funAlphaErr%dVz%d", i, j), "[0]+[1]/(x+0.5)+[2]*exp([3]*x)", 0, 150);
            alphaPlots[i][j]->Fit(Form("funAlpha%dVz%d", i, j), "M E Q", "", 0, 150);
            alphaErrPlots[i][j]->Fit(Form("funAlphaErr%dVz%d", i, j), "M E Q", "", 0, 150);

            fAlpha[i][j] = alphaPlots[i][j]->GetFunction(Form("funAlpha%dVz%d", i, j));
            fAlphaErr[i][j] = alphaErrPlots[i][j]->GetFunction(Form("funAlphaErr%dVz%d", i, j));
         }
      }

      // Vertex and event selection efficiency
      // TCanvas *cTrigEff = new TCanvas("cTrigEff", "TrigEff", _CANVAS_W, _CANVAS_H);
      TrackletTree->Project("hTrigEff", Form("%s", multiplicity), "weight" * (gen_selection && offline_selection && "vz[1]>-99"));
      hTrigEff->Sumw2();
      TrackletTree->Project("hTrigEffNoCut", Form("%s", multiplicity), "weight" * (gen_selection && "vz[1]>-99"));
      hTrigEffNoCut->Sumw2();
      hTrigEff->Divide(hTrigEffNoCut);
      // hTrigEff->Draw();
      // cTrigEff->SaveAs(Form("figs/TrigEff-%s-%d.png", title, TrackletType));

      // Calculate SD'/IN'
      TCanvas *cSDFrac = new TCanvas("cSDFrac", "SD Fraction After Cut", _CANVAS_W, _CANVAS_H);
      TrackletTree->Project("hSDFrac", Form("%s", multiplicity), "weight" * (event_selection && !gen_selection));
      hSDFrac->Sumw2();
      TrackletTree->Project("hSD", Form("%s", multiplicity), "weight" * event_selection);
      hSD->Sumw2();
      hSDFrac->Divide(hSD);

      hSDFrac->SetStats(0);
      hSDFrac->SetAxisRange(-0.01, 0.2, "Y");
      hSDFrac->Draw();
      cSDFrac->SaveAs(Form("figs/corrs/sdfrac-%s-%d.png", title, TrackletType));

      // Calculate Vertexed dN/deta
      // TCanvas *cdNdEtaVertexed = new TCanvas("cdNdEtaVertexed", "dNdEta after vertexing", _CANVAS_W, _CANVAS_H);
      TH1F* hneventWvertexSelection = new TH1F("hneventWvertexSelection", "", 1, 0, 1000);
      TrackletTree->Draw("nhit1>>hneventWvertexSelection", "weight" * (gen_selection && "vz[1]>-99"), "goff");
      float neventWvertexSelection = hneventWvertexSelection->Integral(0, hneventWvertexSelection->GetNbinsX() + 1);
      TrackletTree->Project("hDNDEtaVertexed", "eta", "weight" * (gen_selection && "vz[1]>-99"));
      hDNDEtaVertexed->SetXTitle("#eta Truth after Vtx");
      hDNDEtaVertexed->SetYTitle("dN/#eta");
      hDNDEtaVertexed->Scale(nEtaBin/6*neventWvertexSelection);
      // hDNDEtaVertexed->Draw();

      // Calculate NoVertexed dN/deta
      // TCanvas *cdNdEtaNoVertexed = new TCanvas("cdNdEtaNoVertexed", "dNdEta after vertexing", _CANVAS_W, _CANVAS_H);
      TH1F* hneventWOvertexSelection = new TH1F("hneventWOvertexSelection", "", 1, 0, 1000);
      TrackletTree->Draw("nhit1>>hneventWOvertexSelection", "weight" * gen_selection, "goff");
      float neventWOvertexSelection = hneventWOvertexSelection->Integral(0, hneventWOvertexSelection->GetNbinsX() + 1);
      TrackletTree->Project("hDNDEtaNoVertexed", "eta", "weight" * gen_selection);
      hDNDEtaNoVertexed->SetXTitle("#eta Truth before Vtx");
      hDNDEtaNoVertexed->SetYTitle("dN/#eta");
      hDNDEtaNoVertexed->Scale(nEtaBin/6*neventWOvertexSelection);
      // hDNDEtaNoVertexed->Draw();

      // Calculate Xi (Trigger correction)
      TCanvas *cTriggerCorrection = new TCanvas("cTriggerCorrection", "Xi", _CANVAS_W, _CANVAS_H);
      TH1F* hneventWEventSelection = new TH1F("hneventWEventSelection", "", 1, 0, 1000);
      TrackletTree->Draw("nhit1>>hneventWEventSelection", "weight" * event_selection, "goff");
      float neventWEventSelection = hneventWEventSelection->Integral(0, hneventWEventSelection->GetNbinsX() + 1);
      TrackletTree->Project("hdNdetaWithEvtCut", "eta", "weight" * event_selection);
      hdNdetaWithEvtCut->Sumw2();
      hdNdetaWithEvtCut->Scale(1./neventWEventSelection);

      hTriggerCorrection = (TH1F*)hHadronWOSelection->Project3D("x");
      hTriggerCorrection->SetName("hTriggerCorrection");
      hTriggerCorrection->Sumw2();
      hTriggerCorrection->Scale(1./neventWOSelection);
      hTriggerCorrection->Divide(hdNdetaWithEvtCut);

      hTriggerCorrection->SetStats(0);
      hTriggerCorrection->SetAxisRange(0.9, 1, "Y");
      hTriggerCorrection->Draw();
      cTriggerCorrection->SaveAs(Form("figs/corrs/xi-%s-%d.png", title, TrackletType));
   }

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

            if (apply_accep_corr) {
               double accep_data = haccep_data->GetBinContent(x, z);
               double accep_mc = haccep_mc->GetBinContent(x, z);
               if (accep_data == 0 || accep_mc == 0) {
                  if (hAcceptance1D->GetBinContent(x, z) != 0)
                     printf("  ! acceptance correction error: eta: %i, vz: %i\n", x, z);
               } else {
                  alpha = alpha * accep_mc / accep_data;
                  alphaErr = alphaErr * accep_mc / accep_data;
               }
            }

            if (verbose) printf("   | application - eta: %i, ntl: %i, vz: %i, val: %f, beta: %f, nsig: %f, alpha: %f, alphaErr: %f\n", x, y, z, val, beta, val*(1-beta), alpha, alphaErr);
            // Use extrapolated value if alpha is not available
            if (alpha == 0 && fAlpha[x-1][z-1] != 0) {
               alpha = fAlpha[x-1][z-1]->Eval(TrackletBins[y]);
               if (verbose) printf("   |  ! extrapolated alpha: %f\n", alpha);
            }
            if (alpha == 0) {
               for (int k=0; k<y; k++) {
                  alpha = alphaPlots[x-1][z-1]->GetBinContent(y-k);
                  alphaErr = alphaPlots[x-1][z-1]->GetBinError(y-k);
                  if (verbose) printf("   |  ! use alpha from another tracklet bin: %f, bin: %i\n", alpha, y-k);
                  if (alpha != 0) break;
               }
            }

            if (alpha == 0 || alpha > 10) {
               alpha = 1;
               if (verbose) printf("   |  ! alpha outside bounds: %f\n", alpha);
            }

            halpha_applied->SetBinContent(x, y, z, alpha);
            halpha_applied->SetBinError(x, y, z, alphaErr);

            double nCorrected = val*(1-beta)*alpha;
            double valErr = sqrt(alpha*alpha*(1-beta)*(1-beta)*val + betaErr*betaErr*alpha*alpha*val*val);

            hSubtracted->SetBinContent(x, y, z, nCorrected/alpha);
            hCorrected->SetBinContent(x, y, z, nCorrected);
            hCorrected->SetBinError(x, y, z, valErr);

            totalN += nCorrected;
            totalNErr += valErr*valErr;
         }
      }
      hCorrectedEtaBin->SetBinContent(x, totalN);
      hCorrectedEtaBin->SetBinError(x, sqrt(totalNErr));
   }

   TCanvas* cAcc = new TCanvas("cAcc", "Acceptance", _CANVAS_W, _CANVAS_H);
   hAcceptance1D->SetStats(0);
   hAcceptance1D->Draw();
   cAcc->SaveAs(Form("figs/accep/accep-%s-%i.png", title, TrackletType));

   TH1F* hAcceptance0D = (TH1F*)hAcceptance1D->ProjectionX();
   hAcceptance0D->SetName("hAcceptance0D");
   hAcceptance0D->Scale(1./hAcceptance0D->GetMaximum());

   // Plot RawTracklet and Background Tracklet in nTracklet bin
   // TCanvas *cRawTrackletnTracklet = new TCanvas("cRawTrackletnTracklet", "Raw (nTracklet)", _CANVAS_W, _CANVAS_H);
   TH1F* hMCTruthnTracklet = (TH1F*)hHadronAccepted->Project3D("y");
   hMCTruthnTracklet->Scale(1./nevent);
   hMCTruthnTracklet->SetXTitle("mult |#eta|<3");
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
   // TCanvas *cRawTrackletVz = new TCanvas("cRawTrackletVz", "Raw (Vz)", _CANVAS_W, _CANVAS_H);
   TH1F* hMCTruthVz = (TH1F*)hHadronAccepted->Project3D("z");
   hMCTruthVz->Scale(1./nevent);
   hMCTruthVz->SetXTitle("v_{z}");
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
   // TCanvas *cRawTrackletEta = new TCanvas("cRawTrackletEta", "Raw (Eta)", _CANVAS_W, _CANVAS_H);
   TH1F* hMCTruthEta = (TH1F*)hHadronAccepted->Project3D("x");
   hMCTruthEta->Scale(nEtaBin/6*nevent);
   hMCTruthEta->SetXTitle("#eta");
   // hMCTruthEta->Draw("hist");

   TH1F* hRawTrackletEta = (TH1F*)hEverything->Project3D("x");
   hRawTrackletEta->Scale(nEtaBin/6*nevent);
   // hRawTrackletEta->Draw("same");

   TH1F* hBackgroundTrackletEta = (TH1F*)hReproducedBackground->Project3D("x");
   hBackgroundTrackletEta->SetLineColor(2);
   hBackgroundTrackletEta->SetMarkerColor(2);
   hBackgroundTrackletEta->Scale(nEtaBin/6*nevent);
   // hBackgroundTrackletEta->Draw("same");

   TH1F* hRawMinusBackgroundTrackletEta = (TH1F*)hRawTrackletEta->Clone();
   hRawMinusBackgroundTrackletEta->SetLineColor(4);
   hRawMinusBackgroundTrackletEta->SetMarkerColor(4);
   hRawMinusBackgroundTrackletEta->Add(hBackgroundTrackletEta, -1);
   // hRawMinusBackgroundTrackletEta->Draw("same");
   // cRawTrackletEta->Draw();

   // Plot Before Acceptance correction
   // TCanvas *cDNdEtaC = new TCanvas("cDNdEtaC", "Before Acceptance correction", _CANVAS_W, _CANVAS_H);
   TH1F* hTruthAccepted = (TH1F*)hHadronAccepted->Project3D("x");
   hTruthAccepted->SetName("hTruthAccepted");

   TH1F* hTruth = (TH1F*)hTruthAccepted->Clone();
   hTruth->SetName("hTruth");
   hTruth->Sumw2();
   formatHist(hTruth, 1, nevent/nEtaBin*6);
   hTruth->Divide(hAcceptance0D);

   TH1F* hTruthEvtCutCorrectedByXi = (TH1F*)hHadron->Project3D("x");
   hTruthEvtCutCorrectedByXi->SetName("hTruthEvtCutCorrectedByXi");
   hTruthEvtCutCorrectedByXi->Sumw2();
   formatHist(hTruthEvtCutCorrectedByXi, 1, nevent/nEtaBin*6);

   formatHist(hCorrectedEtaBin, 2, nevent/nEtaBin*6, 1);
   // hCorrectedEtaBin->Draw("e");

   hTruthAccepted->Sumw2();
   formatHist(hTruthAccepted, 1, nevent/nEtaBin*6);
   hTruthAccepted->SetAxisRange(0, _PLOT_RANGE, "Y");
   hTruthAccepted->SetXTitle("#eta (Calculated Hand)");
   hTruthAccepted->SetYTitle("dN/d#eta");
   // hTruthAccepted->Draw("hist same");
   // cDNdEtaC->Draw();

   // Final dN/deta results ===================================================
   TH1F* hMeasured = (TH1F*)hCorrected->Project3D("x");
   hMeasured->SetName("hMeasured");
   hMeasured->Sumw2();
   formatHist(hMeasured, 4, nevent/nEtaBin*6, 1);
   hMeasured->Divide(hAcceptance0D);

   TH1F* hUncorrected = (TH1F*)hEverything->Project3D("x");
   hUncorrected->SetName("hUncorrected");
   hUncorrected->Sumw2();
   formatHist(hUncorrected, 8, nevent/nEtaBin*6, 1);
   hUncorrected->Divide(hAcceptance0D);

   TH1F* hBackgroundSubtracted = (TH1F*)hSubtracted->Project3D("x");
   hBackgroundSubtracted->SetName("hBackgroundSubtracted");
   hBackgroundSubtracted->Sumw2();
   formatHist(hBackgroundSubtracted, 9, nevent/nEtaBin*6, 1);
   hBackgroundSubtracted->Divide(hAcceptance0D);

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
            if (hAcceptance1D->GetBinContent(x, z) != 0) {
               total += hCorrected->GetBinContent(x, y, z);
               double err = hCorrected->GetBinError(x, y, z);
               totalErr = sqrt(totalErr*totalErr + err*err);

               totalMC += hHadronAccepted->GetBinContent(x, y, z);
               double errMC = hHadronAccepted->GetBinError(x, y, z);
               totalMCErr = sqrt(totalMCErr*totalMCErr + errMC*errMC);
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
      if (!apply_trigger_corr) {
         SDFrac = 0;
         TrigEff = 1;
      }
      if (TrigEff != 0)
         nEvt += hnTracklet->GetBinContent(y)/TrigEff*(1-SDFrac*_SD_FACTOR);
      for (int x=1; x<=nEtaBin; x++) {
         for (int z=1; z<=nVzBin; z++) {
            if (hAcceptance1D->GetBinContent(x, z) != 0) {
               if (TrigEff != 0)
                  hAcceptance2D->SetBinContent(x, y, z, hAcceptance2D->GetBinContent(x, y, z)/TrigEff*(1-SDFrac*_SD_FACTOR));
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
         if (!apply_trigger_corr) {
            SDFrac = 0;
            TrigEff = 1;
         }
         if (TrigEff != 0) {
            hMeasuredEtanTracklet->SetBinContent(x, y, hMeasuredEtanTracklet->GetBinContent(x, y)/TrigEff*(1-SDFrac*_SD_FACTOR));
            hMeasuredEtanTracklet->SetBinError(x, y, hMeasuredEtanTracklet->GetBinError(x, y)/TrigEff*(1-SDFrac*_SD_FACTOR));
            hTruthEtanTracklet->SetBinContent(x, y, hTruthEtanTracklet->GetBinContent(x, y)/TrigEff*(1-SDFrac*_SD_FACTOR));
            hTruthEtanTracklet->SetBinError(x, y, hTruthEtanTracklet->GetBinError(x, y)/TrigEff*(1-SDFrac*_SD_FACTOR));
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
   hMeasuredTrigEffCorrectedMult->Divide(hAcceptance2DMult);

   TCanvas* cEmpty = new TCanvas("cEmpty", "Empty Correction", _CANVAS_W, _CANVAS_H);
   if (!apply_correction) {
      hEmptyEvtCorrection = (TH1F*)hTruthWOSelection->Clone("hEmptyEvtCorrection");
      hEmptyEvtCorrection->Divide(hMeasuredTrigEffCorrected);
   } else {
      hEmptyEvtCorrection = (TH1F*)fCorrection->Get("hEmptyEvtCorrection");
   }

   hEmptyEvtCorrection->SetStats(0);
   hEmptyEvtCorrection->SetAxisRange(0.8, 1.2, "Y");
   hEmptyEvtCorrection->Draw();
   cEmpty->Draw();
   cEmpty->SaveAs(Form("figs/corrs/empty-event-%s-%d.png", title, TrackletType));

   // Final result canvas
   TCanvas* cDNdEta = new TCanvas("cDNdEta", "Final result", _CANVAS_W, _CANVAS_H);

   hTruthWOSelection->SetAxisRange(0, _PLOT_RANGE, "Y");
   hTruthWOSelection->SetXTitle("#eta");
   hTruthWOSelection->SetYTitle("dN/d#eta");
   hTruthWOSelection->SetStats(0);
   hTruthWOSelection->Draw("hist");

   TH1F* hMeasuredFinal = (TH1F*)hMeasuredTrigEffCorrected->Clone();
   hMeasuredFinal->SetName("hMeasuredFinal");
   hMeasuredFinal->SetMarkerStyle(20);
   if (apply_trigger_corr) {
      hMeasuredFinal->Multiply(hEmptyEvtCorrection);
   }
   formatHist(hMeasuredFinal, 2, 1, 1.1);
   hMeasuredFinal->SetStats(0);
   hMeasuredFinal->Draw("e same");

   hTruthEvtCutCorrectedByXi->Multiply(hTriggerCorrection);
   hTruthEvtCutCorrectedByXi->SetAxisRange(0, _PLOT_RANGE, "Y");
   hTruthEvtCutCorrectedByXi->SetXTitle("#eta");
   hTruthEvtCutCorrectedByXi->SetYTitle("dN/d#eta");
   // hTruthEvtCutCorrectedByXi->Draw("hist");

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

   cDNdEta->Draw();
   cDNdEta->SaveAs(Form("figs/result/result-%s-%d.png", title, TrackletType));

   // Compare with Truth ======================================================
   TCanvas* cDNdEtaCompare = new TCanvas("cDNdEtaCompare", "Compare", _CANVAS_W, _CANVAS_H);

   hMeasuredFinal->SetAxisRange(0, _PLOT_RANGE, "Y");
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

   for (int i=0; i<nEtaBin; i++) {
      for (int j=0; j<nVzBin; j++) {
         if (fAlpha[i][j]) fAlpha[i][j]->Write();
         if (fAlphaErr[i][j]) fAlphaErr[i][j]->Write();
      }
   }

   outf->Write("", TObject::kOverwrite);

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

int main(int argc, char* argv[]) {
   if (argc == 6)
      return plotFinalResult(atoi(argv[1]), argv[2], argv[3], atoi(argv[4]), argv[5]);
   else if (argc == 7)
      return plotFinalResult(atoi(argv[1]), argv[2], argv[3], atoi(argv[4]), argv[5], atoi(argv[6]));
   else if (argc == 8)
      return plotFinalResult(atoi(argv[1]), argv[2], argv[3], atoi(argv[4]), argv[5], atoi(argv[6]), atoi(argv[7]));
   else if (argc == 9)
      return plotFinalResult(atoi(argv[1]), argv[2], argv[3], atoi(argv[4]), argv[5], atoi(argv[6]), atoi(argv[7]), atoi(argv[8]));
   else if (argc == 10)
      return plotFinalResult(atoi(argv[1]), argv[2], argv[3], atoi(argv[4]), argv[5], atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atoi(argv[9]));
   else if (argc == 11)
      return plotFinalResult(atoi(argv[1]), argv[2], argv[3], atoi(argv[4]), argv[5], atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atoi(argv[9]), atoi(argv[10]));
   else if (argc == 12)
      return plotFinalResult(atoi(argv[1]), argv[2], argv[3], atoi(argv[4]), argv[5], atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atoi(argv[9]), atoi(argv[10]), atoi(argv[11]));
   else if (argc == 13)
      return plotFinalResult(atoi(argv[1]), argv[2], argv[3], atoi(argv[4]), argv[5], atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atoi(argv[9]), atoi(argv[10]), atoi(argv[11]), atoi(argv[12]));
   else
      return -1;
}

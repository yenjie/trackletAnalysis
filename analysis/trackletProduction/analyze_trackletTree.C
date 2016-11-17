#define _DZ 0.12
#define _DPHI 0.05

#define _PI 3.14159265358979

#include <deque>

#include <TFile.h>
#include <TNtuple.h>
#include <TTimeStamp.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TGraph.h>
#include <TFitResult.h>

#include "Tracklet.h"
#include "ClusterVertexCompat.h"
#include "pdfs.h"

using namespace std;

typedef struct Vertex {
   double vz;
   int nz;
   double vzmean;
   double sigma2;
   double chi2;
   double par0;
} Vertex;

double dphi(double phi1, double phi2);

bool sortphi(RecoHit h1, RecoHit h2);
bool sorteta(RecoHit h1, RecoHit h2);

bool sortvz(Vertex v1, Vertex v2);
bool sortnz(Vertex v1, Vertex v2);
bool sortsigma2(Vertex v1, Vertex v2);
bool sortchi2(Vertex v1, Vertex v2);

Double_t csfit(Double_t *x, Double_t *par) {
   return par[0]*cosh(x[0]);
}

int analyze_trackletTree(const char* infile = "PixelTree.root", // Input Pixel Tree file
                         const char* outfile = "output.root",   // Ouptut Tracklet Tree
                         long startEntry = 0,                   // Starting Entry number in the Pixel Tree
                         long endEntry = 1000000000,            // Ending Entry number in the Pixel Tree
                         double pileUp = 0,                     // Artifically overlap event to mimic pile-up, nBX=208: 0.0542, nBX!=208: 0.005
                         bool reWeight = 0,                     // Reweight vertex distribution to match data
                         bool useRandomVertex = 0,              // Use random vertex (instead of the reco one)
                         bool useForwardPixels = 0,             // Use forward pixel detector for vertexing
                         double BGfrac = 0,                     // Additional background level
                         bool addL1Bck = 0,                     // Add random background to first pixel layer
                         bool addL2Bck = 0,                     // Add random background to second pixel layer
                         bool addL3Bck = 0,                     // Add random background to third pixel layer
                         double splitProb = 0,                  // Splitting probability of the pixel hit
                         bool doTwoHitCorrelation = 0,          // Create a pixel counting tree instead of tracklet tree
                         bool smearPixels = 0,                  // Smear pixel hits
                         double dropProb = 0,                   // Emulate efficiency loss
                         bool cutOnClusterSize = 0,             // Cut on clusterSize to reduce background
                         bool reweightMultiplicity = 0,         // Reweight multiplicity distribution
                         bool putBeamHalo = false,              // Adding beam Halo
                         double beamHaloRatio = 0.0,
                         const char* beamHaloFile = "BeamHalo.root",
                         bool useKKVertex = 0,                  // Use vertex from other recoVtx collection
                         double smearVertex = 0)                // Add additional smearing to vertex position
{
   // Set Random Seed =========================================================
   TTimeStamp myTime;
   gRandom->SetSeed(myTime.GetNanoSec());
   cout << "Randomize " << gRandom->Rndm() << endl;

   // Input file ==============================================================
   TFile* inf = TFile::Open(infile);
   TTree* t = dynamic_cast<TTree*>(inf->Get("ana/PixelTree"));
   // TFile* beamHaloInf;
   // TTree* beamHaloTree;
   // if (putBeamHalo) {
   //    cout << "Add Beam Halo Background!!!!" << endl;
   //    cout << "Ratio = "<< beamHaloRatio << endl;
   //    cout << "File = "<< beamHaloFile << endl;
   //    beamHaloInf = new TFile(beamHaloFile);
   //    beamHaloTree = dynamic_cast<TTree*>(beamHaloInf->FindObjectAny("PixelTree"));
   // }

   // Output file =============================================================
   TFile* outf = new TFile(outfile, "recreate");
   TTree* trackletTree12 = new TTree("TrackletTree12", "Tree of Reconstructed Tracklets");
   TTree* trackletTree13 = new TTree("TrackletTree13", "Tree of Reconstructed Tracklets");
   TTree* trackletTree23 = new TTree("TrackletTree23", "Tree of Reconstructed Tracklets");
   // TTree* trackletTree14 = new TTree("TrackletTree14", "Tree of Reconstructed Tracklets");
   // TTree* trackletTree15 = new TTree("TrackletTree15", "Tree of Reconstructed Tracklets");
   // TTree* trackletTree45 = new TTree("TrackletTree45", "Tree of Reconstructed Tracklets");

   TTree* hltTree = ((TTree*)inf->Get("hltanalysis/HltTree"))->CloneTree();

   bool isMC = 0;
   double vzShift = 0;

   // Selection on Hits and events ============================================
   SelectionCriteria cuts;
   cuts.drCut   = 0.4;   // to remove double hit
   cuts.dPhiCut = 0.04;  // to remove double hit
   cuts.dEtaCut = 0.2;   // to remove double hit
   cuts.vzCut   = 10;    // vertex cut

   // Settings ================================================================
   cuts.verbose_          = false;
   cuts.useDeltaPhi_      = false;
   cuts.useDeltaRho_      = false;
   cuts.checkSecondLayer_ = true;

   // Tracklet Tree data format ===============================================
   TrackletData tdata12;
   TrackletData tdata13;
   TrackletData tdata23;
   // TrackletData tdata14;
   // TrackletData tdata15;
   // TrackletData tdata45;

   setTrackletTreeBranch(trackletTree12, tdata12);
   setTrackletTreeBranch(trackletTree13, tdata13);
   setTrackletTreeBranch(trackletTree23, tdata23);
   // setTrackletTreeBranch(trackletTree14, tdata14);
   // setTrackletTreeBranch(trackletTree15, tdata15);
   // setTrackletTreeBranch(trackletTree45, tdata45);

   if (t->GetEntries("nRun<10")!=0) {
      isMC = true;
      cout << "This is a Monte Carlo study." << endl;
      vzShift = -0.36;
      cout << "vzShift = " << vzShift << endl;
   } else {
      cout << "This is a data analysis." << endl;
      pileUp = 0;
      reWeight = 0;
      BGfrac = 0;
      addL1Bck = 0;
      addL2Bck = 0;
      addL3Bck = 0;
      smearVertex = 0;
   }

   TH3D* hLayer1Hit = new TH3D("hLayer1Hit", "", 75, 0, 15, 60, -3, 3, 64, -3.2, 3.2);
   TH3D* hLayer2Hit = new TH3D("hLayer2Hit", "", 75, 0, 15, 60, -3, 3, 64, -3.2, 3.2);
   TH3D* hLayer3Hit = new TH3D("hLayer3Hit", "", 75, 0, 15, 60, -3, 3, 64, -3.2, 3.2);

   // Prepare hit spectra for random hit
   if (addL1Bck) {
      cout << "Projecting...1" << endl;
      t->Project("hLayer1Hit", "phi1:eta1:r1");
   }
   if (addL2Bck) {
      cout << "Projecting...2" << endl;
      t->Project("hLayer2Hit", "phi2:eta2:r2");
   }
   if (addL3Bck) {
      cout << "Projecting...3" << endl;
      t->Project("hLayer3Hit", "phi3:eta3:r3");
   }
   cout << "Projecting...done" << endl;

   // Parameters for the tree =================================================
   Parameters par;
   getPixelTreeBranch(t, par);

   // int nBeamHalo = 0;
   // Parameters beamHaloPar;
   // if (putBeamHalo) getPixelTreeBranch(beamHaloTree, beamHaloPar);
   cout << "Number of Events: " << t->GetEntries() << endl;

   int nPileUp = 1;
   if (pileUp!=0)
      printf("@@@@@ Do pileup mixing with mean: %.4f\n", pileUp);
   if (addL1Bck)
      printf("$$$$$$$$$$$$$$$$$ ADDING NOISE $$$$$$$$$\n");
   if (reWeight)
      printf("######### REWEIGHT VZ ##################\n");
   if (useRandomVertex)
      printf("-------------------- Using Random Vertex\n");
   else if (useKKVertex)
      printf("------------- Using Reconstructed Vertex\n");
   else
      printf("------------------ Using Tracklet Vertex\n");
   printf("................................................................\n");

   TF1* csfitf = new TF1("csfit", csfit, -4, 4, 1);
   csfitf->SetParameters(1,8882, 0);
   csfitf->SetParLimits(0, 1.8882, 1.8882);

   // Main loop ===============================================================
   for (int i=startEntry; i<t->GetEntries()&&i<endEntry; i=i+nPileUp) {
      t->GetEntry(i);
      if (i % 1000 == 0) {
         cout << "Run " << par.nRun << " Event " << i << " "
              << trackletTree12->GetEntries() << " "
              << trackletTree13->GetEntries() << " "
              << trackletTree23->GetEntries() << " "
              // << trackletTree14->GetEntries() << " "
              // << trackletTree15->GetEntries() << " "
              // << trackletTree45->GetEntries()
              // << " Add Beam Halo: " << nBeamHalo << " " << nBeamHalo/(double)i
              << endl;
      }

      if (!isMC) {
         if (par.nLumi<89)
            continue;
      }

      // Reweight MC vertex distribution to match data
      if (reWeight) {
         double myVz = par.vz[1];
         if (myVz < -90) {
            // TF1* f = new TF1("f", "gaus", -30, 30);
            TF1* f = new TF1("f", "1", -20, 20);
            // f->SetParameters(1, -0.394086, 5.32670);
            myVz = f->GetRandom();
            delete f;
            myVz=par.vz[0];
         }

         // 13 TeV Run 246908
         // double DataPdf = TMath::Gaus(myVz, -2.03961-vzShift, 4.2783, 1);
         // 13 TeV Run 247324
         double DataPdf = TMath::Gaus(myVz, -2.14145, 4.30854, 1);
         // double MCPdf = TMath::Gaus(myVz, -0.394086, 5.32670, 1);
         double MCPdf = TMath::Gaus(myVz, -0.356413, 5.33904, 1);

         double Ratio = DataPdf / MCPdf;
         double x = gRandom->Rndm()*3;
         if (x > Ratio) continue;
      }

      // Beam Halo ============================================================
      // if (gRandom->Rndm()<beamHaloRatio && putBeamHalo) {
      //    nBeamHalo++;
      //    bool selectFlag = false;
      //    while (!selectFlag) {
      //       int nEntry = beamHaloTree->GetEntries();
      //       beamHaloTree->GetEntry(nEntry*gRandom->Rndm());
      //       if (par.hltBit[67]==1)
      //          selectFlag = true;
      //    }
      // }

      // Fill reco vertex information
      tdata12.nv = par.nv+1;
      tdata13.nv = par.nv+1;
      tdata23.nv = par.nv+1;
      // tdata14.nv = par.nv+1;
      // tdata15.nv = par.nv+1;
      // tdata45.nv = par.nv+1;
      for (int j=1; j<par.nv; j++) {
         tdata12.vz[j+1] = par.vz[j];
         tdata13.vz[j+1] = par.vz[j];
         tdata23.vz[j+1] = par.vz[j];
         tdata12.vx[j+1] = par.vx[j];
         tdata13.vx[j+1] = par.vx[j];
         tdata23.vx[j+1] = par.vx[j];
         tdata12.vy[j+1] = par.vy[j];
         tdata13.vy[j+1] = par.vy[j];
         tdata23.vy[j+1] = par.vy[j];
         // tdata14.vz[j+1] = par.vz[j];
         // tdata15.vz[j+1] = par.vz[j];
         // tdata45.vz[j+1] = par.vz[j];
         // tdata14.vx[j+1] = par.vx[j];
         // tdata15.vx[j+1] = par.vx[j];
         // tdata45.vx[j+1] = par.vx[j];
         // tdata14.vy[j+1] = par.vy[j];
         // tdata15.vy[j+1] = par.vy[j];
         // tdata45.vy[j+1] = par.vy[j];
      }
      // Fill MC vertex
      tdata12.vx[0] = par.vx[0];
      tdata13.vx[0] = par.vx[0];
      tdata23.vx[0] = par.vx[0];
      tdata12.vy[0] = par.vy[0];
      tdata13.vy[0] = par.vy[0];
      tdata23.vy[0] = par.vy[0];
      tdata12.vz[0] = par.vz[0];
      tdata13.vz[0] = par.vz[0];
      tdata23.vz[0] = par.vz[0];
      // tdata14.vx[0] = par.vx[0];
      // tdata15.vx[0] = par.vx[0];
      // tdata45.vx[0] = par.vx[0];
      // tdata14.vy[0] = par.vy[0];
      // tdata15.vy[0] = par.vy[0];
      // tdata45.vy[0] = par.vy[0];
      // tdata14.vz[0] = par.vz[0];
      // tdata15.vz[0] = par.vz[0];
      // tdata45.vz[0] = par.vz[0];

      // Add background hits
      int bckHits1 = par.nhits1 * BGfrac;
      int bckHits2 = par.nhits2 * BGfrac;
      int bckHits3 = par.nhits3 * BGfrac;

      if (addL1Bck!=0) {
         for (int i=par.nhits1; i<par.nhits1+bckHits1; i++) {
            double eta, phi, r;
            hLayer1Hit->GetRandom3(r, eta, phi);
            par.eta1[i] = eta;
            par.phi1[i] = phi;
            par.r1[i] = r;
         }
         par.nhits1 += bckHits1;
      }

      if (addL2Bck!=0) {
         for (int i=par.nhits2; i<par.nhits2+bckHits2; i++) {
            double eta, phi, r;
            hLayer2Hit->GetRandom3(r, eta, phi);
            par.eta2[i] = eta;
            par.phi2[i] = phi;
            par.r2[i] = r;
         }
         par.nhits2 += bckHits2;
      }

      if (addL3Bck!=0) {
         for (int i=par.nhits3; i<par.nhits3+bckHits3; i++) {
            double eta, phi, r;
            hLayer3Hit->GetRandom3(r, eta, phi);
            par.eta3[i] = eta;
            par.phi3[i] = phi;
            par.r3[i] = r;
         }
         par.nhits3 += bckHits3;
      }

      Parameters pupar[12];
      pupar[0] = par;
      double vzPileUp[12];
      vzPileUp[0] = par.vz[0];
      int recoPU = 0;

      vector<RecoHit> layer1raw, layer2raw, layer3raw, layer4raw, layer5raw;
      prepareHits(layer1raw, par, cuts, 1, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      prepareHits(layer2raw, par, cuts, 2, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      // prepareHits(layer3raw, par, cuts, 3, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      // prepareHits(layer4raw, par, cuts, 4, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      // prepareHits(layer5raw, par, cuts, 5, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);

      vector<RecoHit> allhits, fallhits;
      prepareHits(allhits, par, cuts, 1, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      prepareHits(allhits, par, cuts, 2, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      // prepareHits(allhits, par, cuts, 3, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);

      if (pileUp!=0) {
         nPileUp = 0;
         do {
            nPileUp = gRandom->Poisson(pileUp);
         } while (nPileUp == 0);

         for (int p=1; p<nPileUp; p++) {
            if (i + p == t->GetEntries()) {
               nPileUp = p;
               break;
            }
            t->GetEntry(i+p);
            vzPileUp[p] = par.vz[0];

            bckHits1 = par.nhits1 * BGfrac;
            bckHits2 = par.nhits2 * BGfrac;
            bckHits3 = par.nhits3 * BGfrac;

            if (addL1Bck!=0) {
               for (int i=par.nhits1; i<par.nhits1+bckHits1; i++) {
                  double eta, phi, r;
                  hLayer1Hit->GetRandom3(r, eta, phi);
                  par.eta1[i] = eta;
                  par.phi1[i] = phi;
                  par.r1[i] = r;
               }
               par.nhits1 += bckHits1;
            }

            if (addL2Bck!=0) {
               for (int i=par.nhits2; i<par.nhits2+bckHits2; i++) {
                  double eta, phi, r;
                  hLayer2Hit->GetRandom3(r, eta, phi);
                  par.eta2[i] = eta;
                  par.phi2[i] = phi;
                  par.r2[i] = r;
               }
               par.nhits2 += bckHits2;
            }

            if (addL3Bck!=0) {
               for (int i=par.nhits3; i<par.nhits3+bckHits3; i++) {
                  double eta, phi, r;
                  hLayer3Hit->GetRandom3(r, eta, phi);
                  par.eta3[i] = eta;
                  par.phi3[i] = phi;
                  par.r3[i] = r;
               }
               par.nhits3 += bckHits3;
            }
            pupar[p] = par;

            prepareHits(layer1raw, par, cuts, 1, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
            prepareHits(layer2raw, par, cuts, 2, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
            // prepareHits(layer3raw, par, cuts, 3, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
            // prepareHits(layer4raw, par, cuts, 4, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
            // prepareHits(layer5raw, par, cuts, 5, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);

            prepareHits(allhits, par, cuts, 1, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
            prepareHits(allhits, par, cuts, 2, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
            // prepareHits(allhits, par, cuts, 3, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
         }
         t->GetEntry(i);
         par = pupar[0];
      }

      // Vertex reconstruction ================================================
      double trackletVertex = -99;

      std::sort(allhits.begin(), allhits.end(), sortphi);
      for (unsigned int m=0; m<allhits.size() && allhits[m].phi<_DPHI-_PI; m++) {
         allhits.push_back(allhits[m]);
         allhits.back().phi += 2*_PI;
      }
      RecoHit fake(0, 4, 0, 0, 0, 1);
      allhits.push_back(fake);

      std::vector<Vertex> vertices;
      std::deque<RecoHit> trackcands[3];
      for (unsigned int a=0; a<allhits.size(); a++) {
         for (int q=0; q<2; q++) {
            while (!trackcands[q].empty() && allhits[a].phi>trackcands[q].front().phi+_DPHI) {
               for (std::deque<RecoHit>::iterator it = trackcands[1-q].begin(); it != trackcands[1-q].end(); it++) {
                  Vertex vertex;
                  double r1 = trackcands[q].front().r;
                  double z1 = r1/tan(2*atan(exp(-trackcands[q].front().eta)));
                  double r2 = (*it).r;
                  double z2 = r2/tan(2*atan(exp(-(*it).eta)));
                  vertex.vz = z1-(z2-z1)/(r2-r1)*r1;
                  if (fabs(vertex.vz)<200.0)
                     vertices.push_back(vertex);
               }
               trackcands[q].pop_front();
            }
         }
         trackcands[allhits[a].layer-1].push_back(allhits[a]);
      }

      // for (unsigned int a=0; a<allhits.size(); a++) {
      //    for (int q=0; q<3; q++) {
      //       while (!trackcands[q].empty() && allhits[a].phi>trackcands[q].front().phi+_DPHI) {
      //          for (int w=0; w<3; w++) {
      //             if (q == w) continue;
      //             for (std::deque<RecoHit>::iterator it = trackcands[w].begin(); it != trackcands[w].end(); it++) {
      //                Vertex vertex;
      //                double r1 = trackcands[q].front().r;
      //                double z1 = r1/tan(2*atan(exp(-trackcands[q].front().eta)));
      //                double r2 = (*it).r;
      //                double z2 = r2/tan(2*atan(exp(-(*it).eta)));
      //                vertex.vz = z1-(z2-z1)/(r2-r1)*r1;
      //                if (fabs(vertex.vz)<20.0)
      //                   vertices.push_back(vertex);
      //             }
      //          }
      //          trackcands[q].pop_front();
      //       }
      //    }
      //    trackcands[allhits[a].layer-1].push_back(allhits[a]);
      // }

      // if (useForwardPixels) {
      //    prepareHits(fallhits, par, cuts, 1, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      //    prepareHits(fallhits, par, cuts, 4, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      //    prepareHits(fallhits, par, cuts, 5, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);

      //    if (nPileUp>0) {
      //       for (int p=1; p<nPileUp; p++) {
      //          t->GetEntry(i+p);
      //          prepareHits(fallhits, par, cuts, 1, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      //          prepareHits(fallhits, par, cuts, 4, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      //          prepareHits(fallhits, par, cuts, 5, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      //       }
      //       t->GetEntry(i);
      //    }

      //    std::sort(fallhits.begin(), fallhits.end(), sortphi);
      //    for (unsigned int m=0; m<fallhits.size() && fallhits[m].phi<_DPHI-_PI; m++) {
      //       fallhits.push_back(fallhits[m]);
      //       fallhits.back().phi += 2*_PI;
      //    }
      //    fallhits.push_back(fake);

      //    std::deque<RecoHit> ftrackcands[5];
      //    for (unsigned int a=0; a<fallhits.size(); a++) {
      //       for (int q=0; q<5; q++) {
      //          while (!ftrackcands[q].empty() && fallhits[a].phi>ftrackcands[q].front().phi+_DPHI) {
      //             for (int w=0; w<5; w++) {
      //                if (q == w) continue;
      //                if (q && w) continue;
      //                for (std::deque<RecoHit>::iterator it = ftrackcands[w].begin(); it != ftrackcands[w].end(); it++) {
      //                   Vertex vertex;
      //                   double r1 = ftrackcands[q].front().r;
      //                   double z1 = r1/tan(2*atan(exp(-ftrackcands[q].front().eta)));
      //                   double r2 = (*it).r;
      //                   double z2 = r2/tan(2*atan(exp(-(*it).eta)));
      //                   vertex.vz = z1-(z2-z1)/(r2-r1)*r1;
      //                   if (fabs(vertex.vz)<20.0)
      //                      vertices.push_back(vertex);
      //                }
      //             }
      //             ftrackcands[q].pop_front();
      //          }
      //       }
      //       ftrackcands[fallhits[a].layer-1].push_back(fallhits[a]);
      //    }
      // }

      // Vertex clustering ====================================================
      if (vertices.size()) {
         std::sort(vertices.begin(), vertices.end(), sortvz);
         for (unsigned int z=0; z<vertices.size(); z++) {
            vertices[z].nz = 0;
            vertices[z].vzmean = 0;
            unsigned int y = 0;
            for (; y<vertices.size() && vertices[y].vz-vertices[z].vz<_DZ; y++) {
               if (fabs(vertices[y].vz-vertices[z].vz)<_DZ) {
                  vertices[z].nz++;
                  vertices[z].vzmean += vertices[y].vz;
               }
            }
            vertices[z].vzmean /= vertices[z].nz;
            vertices[z].sigma2 = 0;
            for (--y; y<vertices.size() && vertices[z].vz-vertices[y].vz<_DZ; y--) {
               if (fabs(vertices[y].vz-vertices[z].vz)<_DZ)
                  vertices[z].sigma2 += (vertices[y].vz-vertices[z].vzmean)*(vertices[y].vz-vertices[z].vzmean);
            }
            vertices[z].sigma2 /= vertices[z].nz;
         }

         std::stable_sort(vertices.begin(), vertices.end(), sortnz);
         if (vertices[0].nz == 1) {
            std::sort(vertices.begin(), vertices.end(), sortvz);
            for (unsigned int z=0; z<vertices.size(); z++) {
               vertices[z].nz = 0;
               vertices[z].vzmean = 0;
               unsigned int y = z;
               for (; y<vertices.size() && vertices[y].vz-vertices[z].vz<2*_DZ; y++) {
                  vertices[z].nz++;
                  vertices[z].vzmean += vertices[y].vz;
               }
               vertices[z].vzmean /= vertices[z].nz;
               vertices[z].sigma2 = 0;
               for (--y; y>=z && y<vertices.size(); y--) {
                  vertices[z].sigma2 += (vertices[y].vz-vertices[z].vzmean)*(vertices[y].vz-vertices[z].vzmean);
               }
               vertices[z].sigma2 /= vertices[z].nz;
            }
            std::stable_sort(vertices.begin(), vertices.end(), sortnz);
         }

         std::vector<Vertex> candidates;
         candidates.push_back(vertices[0]);
         for (unsigned int c=1; c<vertices.size(); c++) {
            if (vertices[c].nz + 1 < vertices[0].nz)
               break;
            if (vertices[c].vzmean != candidates.back().vzmean)
               candidates.push_back(vertices[c]);
         }

         std::vector<Vertex> pileupcands;
         pileupcands.push_back(vertices[0]);
         for (unsigned int p=1; p<vertices.size(); p++) {
            if (vertices[p].nz < 4)
               break;
            if (vertices[p].nz < pileupcands[0].nz*0.5)
               break;
            std::vector<Vertex>::iterator it = pileupcands.begin();
            for (; it!=pileupcands.end() && abs((*it).vzmean-vertices[p].vzmean)>0.06; it++);
            if (it!=pileupcands.end())
               continue;
            pileupcands.push_back(vertices[p]);
         }
         recoPU = pileupcands.size();

         // Use chi2 and sigma2 ===============================================
         for (unsigned int d=0; d<candidates.size(); d++) {
            TGraph* csveta_g = new TGraph(layer1raw.size()+layer2raw.size());

            int g = 0;
            for (unsigned int c1=0; c1<layer1raw.size(); c1++) {
               double r1 = layer1raw[c1].r;
               double z1 = r1/tan(2*atan(exp(-layer1raw[c1].eta)));
               double neweta = (z1 - candidates[d].vzmean > 0) ? -log(tan(atan(r1/(z1-candidates[d].vzmean))/2)) : log(tan(atan(r1/(candidates[d].vzmean-z1))/2));
               if (layer1raw[c1].cs<24 && layer1raw[c1].cs>=0.9*cosh(neweta))
                  csveta_g->SetPoint(g++, neweta, layer1raw[c1].cs);
            }
            for (unsigned int c2=0; c2<layer2raw.size(); c2++) {
               double r2 = layer2raw[c2].r;
               double z2 = r2/tan(2*atan(exp(-layer2raw[c2].eta)));
               double neweta = (z2 - candidates[d].vzmean > 0) ? -log(tan(atan(r2/(z2-candidates[d].vzmean))/2)) : log(tan(atan(r2/(candidates[d].vzmean-z2))/2));
               if (layer2raw[c2].cs<24)
                  csveta_g->SetPoint(g++, neweta, layer2raw[c2].cs);
            }
            csveta_g->Set(g);
            if (g) {
               TFitResultPtr csveta_g_r = csveta_g->Fit("csfit", "WBQS");
               candidates[d].chi2 = csveta_g_r->Chi2();
               candidates[d].par0 = csveta_g_r->Parameter(0);
            } else {
               candidates[d].chi2 = 9999.9;
               candidates[d].par0 = 9999.9;
            }
         }

         std::vector<Vertex>::iterator eqnz = candidates.begin();
         for (; eqnz!=candidates.end() && (*eqnz).nz==candidates[0].nz; eqnz++);
         std::sort(candidates.begin(), eqnz, sortchi2);
         std::sort(eqnz, candidates.end(), sortchi2);

         std::vector<Vertex>::iterator eqchi2 = candidates.begin();
         for (; eqchi2!=eqnz && (*eqchi2).chi2/candidates[0].chi2<1.12; eqchi2++);
         std::sort(candidates.begin(), eqchi2, sortsigma2);

         std::vector<Vertex>::iterator nzminusone = eqnz;
         for (; nzminusone!=candidates.end() && (*nzminusone).chi2<candidates[0].chi2/2; nzminusone++);
         if (eqnz != nzminusone) {
            std::stable_sort(eqnz, nzminusone, sortsigma2);
            candidates.insert(candidates.begin(), *eqnz);
         }

         // // Use sigma2 only ===================================================
         // std::vector<Vertex>::iterator eqnz = candidates.begin();
         // for (; eqnz!=candidates.end() && (*eqnz).nz==candidates[0].nz; eqnz++);
         // std::sort(candidates.begin(), eqnz, sortsigma2);

         trackletVertex = candidates[0].vzmean;
      }

      // For particle gun =====================================================
      //    if (i%1000 == 0) cout << "!!! USE GEN VERTEX (FOR PARTICLE GUN) " << endl;
      //    trackletVertex = par.vz[0];

      // double smear = 0;
      // if (smearVertex!=0) {
      //    if (i==1) cout << "Vertex smeared!" << endl;
      //    while (smear==0) {
      //       double x = gRandom->Rndm()*2-1;
      //       if (gRandom->Rndm()<TMath::Gaus(x, 0, smearVertex, 1))
      //          smear = x;
      //    }
      //    trackletVertex += smear;
      // }

      if (useRandomVertex) {
         tdata12.vz[1] = gRandom->Rndm()*22-13-vzShift;
         tdata13.vz[1] = tdata12.vz[1];
         tdata23.vz[1] = tdata12.vz[1];
         // tdata14.vz[1] = tdata12.vz[1];
         // tdata15.vz[1] = tdata12.vz[1];
         // tdata45.vz[1] = tdata12.vz[1];
      } else if (useKKVertex) {
         tdata12.vx[1] = par.vx[1];
         tdata12.vy[1] = par.vy[1];
         tdata12.vz[1] = par.vz[1];
         tdata13.vx[1] = par.vx[1];
         tdata13.vy[1] = par.vy[1];
         tdata13.vz[1] = par.vz[1];
         tdata23.vx[1] = par.vx[1];
         tdata23.vy[1] = par.vy[1];
         tdata23.vz[1] = par.vz[1];
         // tdata14.vx[1] = par.vx[1];
         // tdata14.vy[1] = par.vy[1];
         // tdata14.vz[1] = par.vz[1];
         // tdata15.vx[1] = par.vx[1];
         // tdata15.vy[1] = par.vy[1];
         // tdata15.vz[1] = par.vz[1];
         // tdata45.vx[1] = par.vx[1];
         // tdata45.vy[1] = par.vy[1];
         // tdata45.vz[1] = par.vz[1];
      } else {
         tdata12.vz[1] = trackletVertex;
         tdata13.vz[1] = trackletVertex;
         tdata23.vz[1] = trackletVertex;
         // tdata14.vz[1] = trackletVertex;
         // tdata15.vz[1] = trackletVertex;
         // tdata45.vz[1] = trackletVertex;
      }

      // Process hits with Vz constraint:
      vector<RecoHit> layer1Cut;
      prepareHits(layer1Cut, par, cuts, 1, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, 1, par.nRun, par.nLumi, 0);

      std::vector<RecoHit> layer1, layer2, layer3, layer4, layer5;
      prepareHits(layer1, par, cuts, 1, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      prepareHits(layer2, par, cuts, 2, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      prepareHits(layer3, par, cuts, 3, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      // prepareHits(layer4, par, cuts, 4, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      // prepareHits(layer5, par, cuts, 5, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);

      // prepareHits(layer2, par, cuts, 24, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      // prepareHits(layer2, par, cuts, 25, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      // prepareHits(layer3, par, cuts, 34, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      // prepareHits(layer3, par, cuts, 35, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);

      if (nPileUp>0) {
         for (int p=1; p<nPileUp; p++) {
            t->GetEntry(i+p);
            par = pupar[p];
            prepareHits(layer1, par, cuts, 1, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
            prepareHits(layer2, par, cuts, 2, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
            prepareHits(layer3, par, cuts, 3, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
            // prepareHits(layer4, par, cuts, 4, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
            // prepareHits(layer5, par, cuts, 5, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);

            // prepareHits(layer2, par, cuts, 24, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
            // prepareHits(layer2, par, cuts, 25, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
            // prepareHits(layer3, par, cuts, 34, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
            // prepareHits(layer3, par, cuts, 35, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
         }
         t->GetEntry(i);
         par = pupar[0];
      }

      std::vector<Tracklet> recoTracklets12, recoTracklets13, recoTracklets23, recoTracklets14, recoTracklets15, recoTracklets45;

      std::vector<RecoHit> combinedhits;
      if (doTwoHitCorrelation) {
         // Hit correlations in a single layer ================================
         combinedhits.reserve(layer3.size());
         combinedhits.insert(combinedhits.end(), layer1.begin(), layer1.end());
         std::sort(combinedhits.begin(), combinedhits.end(), sorteta);
         recoTracklets12 = recoTracklets(combinedhits, 1, 1);

         combinedhits.clear();
         combinedhits.insert(combinedhits.end(), layer2.begin(), layer2.end());
         std::sort(combinedhits.begin(), combinedhits.end(), sorteta);
         recoTracklets13 = recoTracklets(combinedhits, 2, 2);

         combinedhits.clear();
         combinedhits.insert(combinedhits.end(), layer3.begin(), layer3.end());
         std::sort(combinedhits.begin(), combinedhits.end(), sorteta);
         recoTracklets23 = recoTracklets(combinedhits, 3, 3);
      } else {
         // Tracklet Reconstruction ===========================================
         combinedhits.reserve(layer2.size() + layer3.size());
         combinedhits.insert(combinedhits.end(), layer1.begin(), layer1.end());
         combinedhits.insert(combinedhits.end(), layer2.begin(), layer2.end());
         std::sort(combinedhits.begin(), combinedhits.end(), sorteta);
         recoTracklets12 = recoTracklets(combinedhits, 1, 2);

         combinedhits.clear();
         combinedhits.insert(combinedhits.end(), layer1.begin(), layer1.end());
         combinedhits.insert(combinedhits.end(), layer3.begin(), layer3.end());
         std::sort(combinedhits.begin(), combinedhits.end(), sorteta);
         recoTracklets13 = recoTracklets(combinedhits, 1, 3);

         combinedhits.clear();
         combinedhits.insert(combinedhits.end(), layer2.begin(), layer2.end());
         combinedhits.insert(combinedhits.end(), layer3.begin(), layer3.end());
         std::sort(combinedhits.begin(), combinedhits.end(), sorteta);
         recoTracklets23 = recoTracklets(combinedhits, 2, 3);

         // combinedhits.clear();
         // combinedhits.insert(combinedhits.end(), layer1.begin(), layer1.end());
         // combinedhits.insert(combinedhits.end(), layer4.begin(), layer4.end());
         // std::sort(combinedhits.begin(), combinedhits.end(), sorteta);
         // recoTracklets14 = recoTracklets(combinedhits, 1, 4);

         // combinedhits.clear();
         // combinedhits.insert(combinedhits.end(), layer1.begin(), layer1.end());
         // combinedhits.insert(combinedhits.end(), layer5.begin(), layer5.end());
         // std::sort(combinedhits.begin(), combinedhits.end(), sorteta);
         // recoTracklets15 = recoTracklets(combinedhits, 1, 5);

         // combinedhits.clear();
         // combinedhits.insert(combinedhits.end(), layer4.begin(), layer4.end());
         // combinedhits.insert(combinedhits.end(), layer5.begin(), layer5.end());
         // std::sort(combinedhits.begin(), combinedhits.end(), sorteta);
         // recoTracklets45 = recoTracklets(combinedhits, 4, 5);
      }

      // Cluster Vertex Compatibility =========================================
      double clusVtxQual1 = 0;
      double clusVtxQual2 = 0;
      double clusVtxQual3 = 0;

      // double clusVtxQual1 = getClusVtxCompat(layer1, 1);
      // double clusVtxQual2 = getClusVtxCompat(layer2, 2);
      // double clusVtxQual3 = getClusVtxCompat(layer3, 3);

      // Move the Vertex back
      // if (smearVertex!=0) {
      //    tdata12.vz[1] = trackletVertex - smear;
      //    tdata13.vz[1] = trackletVertex - smear;
      //    tdata23.vz[1] = trackletVertex - smear;
      // }

      // Fill Ntuple
      tdata12.nTracklet  = recoTracklets12.size();
      tdata12.nhit1      = layer1.size();
      tdata12.nhit2      = layer2.size();
      tdata12.nRun       = par.nRun;
      tdata12.nEv        = par.nEv;
      tdata12.nLumi      = par.nLumi;
      tdata12.nBX        = par.nBX;
      tdata12.nHFn       = par.nHFp;
      tdata12.nHFp       = par.nHFn;
      tdata12.nHits      = layer1.size() + layer2.size() + layer3.size();
      tdata12.ntrks      = par.ntrks;
      tdata12.ntrksCut   = par.ntrksCut;
      tdata12.nPU        = nPileUp;
      tdata12.recoPU     = recoPU;
      tdata12.clusVtxQual1 = clusVtxQual1;
      tdata12.clusVtxQual2 = clusVtxQual2;
      tdata12.clusVtxQual3 = clusVtxQual3;

      for (int j=0; j<nPileUp; j++)
         tdata12.vzPU[j] = vzPileUp[j];

      int ntracklet12s = 0;
      int ntracklet12b = 0;
      for (int j=0; j<(int)tdata12.nTracklet; j++) {
         tdata12.eta1[j] = recoTracklets12[j].eta1();
         tdata12.eta2[j] = recoTracklets12[j].eta2();
         tdata12.r1[j]   = recoTracklets12[j].r1();
         tdata12.r2[j]   = recoTracklets12[j].r2();
         tdata12.cs1[j]  = recoTracklets12[j].cs1();
         tdata12.cs2[j]  = recoTracklets12[j].cs2();
         tdata12.ch1[j]  = recoTracklets12[j].ch1();
         tdata12.ch2[j]  = recoTracklets12[j].ch2();
         tdata12.phi1[j] = recoTracklets12[j].phi1();
         tdata12.phi2[j] = recoTracklets12[j].phi2();
         tdata12.deta[j] = recoTracklets12[j].deta();
         tdata12.dphi[j] = recoTracklets12[j].dphi();
         if (fabs(tdata12.deta[j])<0.1) {
            if (fabs(tdata12.dphi[j])<1.0)
               ntracklet12s++;
            if (fabs(tdata12.dphi[j])>1.0 && fabs(tdata12.dphi[j])<2.0)
               ntracklet12b++;
         }
      }
      tdata12.mult = ntracklet12s - ntracklet12b;
      tdata12.mult2 = layer1Cut.size();
      tdata12.npart = 0;

      // bool reWeightMultDropFlag = 0;
      // if (reweightMultiplicity) {
      //    TH1F* hRatio = getMultRatio();
      //    reWeightMultDropFlag = 0;
      //    double myVz = tdata12.vz[1];
      //    double Ratio = hRatio->GetBinContent(hRatio->FindBin(tdata12.mult));
      //    // cout << Ratio << endl;
      //    double x = gRandom->Rndm()*1.6;
      //    if (x > Ratio) reWeightMultDropFlag = 1;
      // }
      // if (reWeightMultDropFlag) continue;

      double pro1 = 0;
      double pro2 = 0;
      for (int j=0; j<12; j++) tdata12.nhad[j] = 0;
      for (int j=0; j<par.npart; j++) {
         if (fabs(par.pdg[j])==2212) {
            double momentum = par.pt[j] * cosh(par.eta[j]);
            if (momentum>pro2) {
               if (momentum>pro1) {
                  pro2 = pro1;
                  pro1 = momentum;
               } else {
                  pro2 = momentum;
               }
            }
         }
         if (fabs(par.eta[j])>3 || par.chg[j]==0 || abs(par.pdg[j])==11 || abs(par.pdg[j])==13) continue;
         tdata12.eta[tdata12.npart] = par.eta[j];
         tdata12.phi[tdata12.npart] = par.phi[j];
         tdata12.chg[tdata12.npart] = par.chg[j];
         tdata12.pdg[tdata12.npart] = par.pdg[j];
         if (abs(par.pdg[j]-11)<1) cout << "Oh no???" << endl;
         tdata12.pt[tdata12.npart] = par.pt[j];
         tdata12.npart++;
         int bin = (int)((par.eta[j]+3)*2);
         int pdg = (int)abs(par.pdg[j]);
         if (pdg==211 || pdg==321 || pdg==2212 || pdg==3122) tdata12.nhad[bin]++;
      }

      tdata12.evtType = par.evtType;
      tdata12.pro2 = pro2;
      for (int j=0; j<par.nv; j++)
         tdata12.vz[j] += vzShift;

      trackletTree12->Fill();

// IMPORTANT: Variable names MUST NOT contain q or w.
#define fillTrackletTree(q, w) {                               \
   tdata##q##w.nTracklet  = recoTracklets##q##w.size();        \
   tdata##q##w.nhit1      = layer##q.size();                   \
   tdata##q##w.nhit2      = layer##w.size();                   \
   tdata##q##w.nRun       = par.nRun;                          \
   tdata##q##w.nEv        = par.nEv;                           \
   tdata##q##w.nLumi      = par.nLumi;                         \
   tdata##q##w.nBX        = par.nBX;                           \
   tdata##q##w.nHFn       = par.nHFp;                          \
   tdata##q##w.nHFp       = par.nHFn;                          \
   tdata##q##w.nHits      = tdata12.nHits;                     \
   tdata##q##w.ntrks      = par.ntrks;                         \
   tdata##q##w.ntrksCut   = par.ntrksCut;                      \
   tdata##q##w.nPU        = nPileUp;                           \
   tdata##q##w.recoPU     = recoPU;                            \
   tdata##q##w.clusVtxQual1 = clusVtxQual1;                    \
   tdata##q##w.clusVtxQual2 = clusVtxQual2;                    \
   tdata##q##w.clusVtxQual3 = clusVtxQual3;                    \
                                                               \
   for (int j=0; j<nPileUp; j++)                               \
      tdata##q##w.vzPU[j] = vzPileUp[j];                       \
                                                               \
   int ntracklet##q##w##s = 0;                                 \
   int ntracklet##q##w##b = 0;                                 \
   for (int j=0; j<(int)tdata##q##w.nTracklet; j++) {          \
      tdata##q##w.eta1[j] = recoTracklets##q##w[j].eta1();     \
      tdata##q##w.eta2[j] = recoTracklets##q##w[j].eta2();     \
      tdata##q##w.r1[j]   = recoTracklets##q##w[j].r1();       \
      tdata##q##w.r2[j]   = recoTracklets##q##w[j].r2();       \
      tdata##q##w.cs1[j]  = recoTracklets##q##w[j].cs1();      \
      tdata##q##w.cs2[j]  = recoTracklets##q##w[j].cs2();      \
      tdata##q##w.ch1[j]  = recoTracklets##q##w[j].ch1();      \
      tdata##q##w.ch2[j]  = recoTracklets##q##w[j].ch2();      \
      tdata##q##w.phi1[j] = recoTracklets##q##w[j].phi1();     \
      tdata##q##w.phi2[j] = recoTracklets##q##w[j].phi2();     \
      tdata##q##w.deta[j] = recoTracklets##q##w[j].deta();     \
      tdata##q##w.dphi[j] = recoTracklets##q##w[j].dphi();     \
      if (fabs(tdata##q##w.deta[j])<0.1) {                     \
         if (fabs(tdata##q##w.dphi[j])<1.0)                    \
            ntracklet##q##w##s++;                              \
         if (fabs(tdata##q##w.dphi[j])>1.0 &&                  \
             fabs(tdata##q##w.dphi[j])<2.0)                    \
            ntracklet##q##w##b++;                              \
      }                                                        \
   }                                                           \
   tdata##q##w.mult = ntracklet##q##w##s - ntracklet##q##w##b; \
   tdata##q##w.mult2 = layer1Cut.size();                       \
   tdata##q##w.npart = 0;                                      \
   for (int j=0; j<12; j++) tdata##q##w.nhad[j] = 0;           \
   for (int j=0; j<par.npart; j++) {                           \
      if (fabs(par.eta[j])>3 || par.chg[j]==0 ||               \
          fabs(par.pdg[j])==11 || abs(par.pdg[j])==13)         \
         continue;                                             \
      tdata##q##w.eta[tdata##q##w.npart] = par.eta[j];         \
      tdata##q##w.phi[tdata##q##w.npart] = par.phi[j];         \
      tdata##q##w.chg[tdata##q##w.npart] = par.chg[j];         \
      tdata##q##w.pdg[tdata##q##w.npart] = par.pdg[j];         \
      tdata##q##w.pt[tdata##q##w.npart] = par.pt[j];           \
      tdata##q##w.npart++;                                     \
      int bin = (int)((par.eta[j]+3)*2);                       \
      int pdg = (int)abs(par.pdg[j]);                          \
      if (pdg==211 || pdg==321 || pdg==2212 || pdg==3122)      \
         tdata##q##w.nhad[bin]++;                              \
   }                                                           \
                                                               \
   tdata##q##w.evtType = par.evtType;                          \
   tdata##q##w.pro2 = pro2;                                    \
   for (int j=0; j<par.nv; j++)                                \
      tdata##q##w.vz[j] += vzShift;                            \
                                                               \
   trackletTree##q##w->Fill();                                 \
}
      fillTrackletTree(1, 3);
      fillTrackletTree(2, 3);
      /*
      fillTrackletTree(1, 4);
      fillTrackletTree(1, 5);
      fillTrackletTree(4, 5);
      */
   }

   inf->Close();

   hltTree->Write("", TObject::kOverwrite);
   outf->Write("", TObject::kOverwrite);
   outf->Close();

   return 0;
}

double dphi(double phi1, double phi2) {
   double pi = 3.14159265358979;
   double dphi = fabs(phi1 - phi2);

   if (dphi < pi)
      return dphi;
   else
      return 2 * pi - dphi;
}

bool sortphi(RecoHit h1, RecoHit h2) {
   return h1.phi < h2.phi;
}

bool sorteta(RecoHit h1, RecoHit h2) {
   return h1.eta > h2.eta;
}

bool sortvz(Vertex v1, Vertex v2) {
   return v1.vz < v2.vz;
}

bool sortnz(Vertex v1, Vertex v2) {
   return v1.nz > v2.nz;
}

bool sortsigma2(Vertex v1, Vertex v2) {
   return v1.sigma2 < v2.sigma2;
}

bool sortchi2(Vertex v1, Vertex v2) {
   return v1.chi2 < v2.chi2;
}

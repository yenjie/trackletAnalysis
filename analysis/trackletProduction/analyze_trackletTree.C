#define _DZ 0.12

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
                         int addL1Bck = 0,                      // Add random background to first pixel layer
                         int addL2Bck = 0,                      // Add random background to second pixel layer
                         int addL3Bck = 0,                      // Add random background to third pixel layer
                         bool smearPixels = 0,                  // Smear pixel hits
                         double smearVertex = 0,                // Add additional smearing to vertex position
                         bool reWeight = 0,                     // Reweight to Run 123596 vtx distribution
                         bool reweightMultiplicity = 0,         // Reweight the multiplicity distribution
                         bool cutOnClusterSize = 0,             // Cut on clusterSize to reduce background
                         int makeVzCut = 0,                     // Cut on Vz
                         double splitProb = 0,                  // Splitting probability of the pixel hit
                         double dropProb = 0,                   // Emulate efficiency loss
                         double pileUp = 0,                     // Artifically overlap event to mimic pile-up
                         bool putBeamHalo = false,              // Adding beam Halo
                         double beamHaloRatio = 0.0,
                         const char* beamHaloFile = "BeamHalo.root",
                         bool useKKVertex = 0,                  // Use vertex from other recoVtx collection
                         bool useRandomVertex = 0,              // Use random vertex (instead of the reco one)
                         bool mimicPixelCounting = 0,           // Create a pixel counting tree instead of tracklet tree
                         bool checkDuplicateEvent = 0)          // Check if we have duplicates in the sample (slow)
{
   // Set Random Seed =========================================================
   TTimeStamp myTime;
   gRandom->SetSeed(myTime.GetNanoSec());
   cout << "Randomize " << gRandom->Rndm() << endl;

   // Input file ==============================================================
   TFile* inf = new TFile(infile);
   TTree* t = dynamic_cast<TTree*>(inf->Get("ana/PixelTree"));
   TFile* beamHaloInf;
   TTree* beamHaloTree;
   if (putBeamHalo) {
      cout << "Add Beam Halo Background!!!!" << endl;
      cout << "Ratio = "<< beamHaloRatio << endl;
      cout << "File = "<< beamHaloFile << endl;
      beamHaloInf = new TFile(beamHaloFile);
      beamHaloTree = dynamic_cast<TTree*>(beamHaloInf->FindObjectAny("PixelTree"));
   }

   // Output file =============================================================
   TFile* outf = new TFile(outfile, "recreate");
   TTree* trackletTree12 = new TTree("TrackletTree12", "Tree of Reconstructed Tracklets");
   TTree* trackletTree13 = new TTree("TrackletTree13", "Tree of Reconstructed Tracklets");
   TTree* trackletTree23 = new TTree("TrackletTree23", "Tree of Reconstructed Tracklets");
   TTree* trackletTree14 = new TTree("TrackletTree14", "Tree of Reconstructed Tracklets");
   TTree* trackletTree15 = new TTree("TrackletTree15", "Tree of Reconstructed Tracklets");
   TTree* trackletTree45 = new TTree("TrackletTree45", "Tree of Reconstructed Tracklets");
   // TTree* hltTree = ((TTree*)inf->Get("hltanalysis/HltTree"))->CloneTree();
   // TTree* outTree = t->CloneTree();

   int vertexHitRegion = 500000;
   int nbins = 100;
   bool isMC = 0;
   double vzShift = 0;

   // Selection on Hits and events ============================================
   SelectionCriteria cuts;
   cuts.drCut   = 0.4;      // to remove double hit
   cuts.dPhiCut = 0.04;     // to remove double hit
   cuts.dEtaCut = 0.2;      // to remove double hit
   cuts.vzCut   = 10;       // vertex cut

   // Settings ================================================================
   cuts.verbose_          = false;
   cuts.useDeltaPhi_      = false;
   cuts.useDeltaRho_      = false;
   cuts.checkSecondLayer_ = true;

   // Tracklet Tree data format ===============================================
   TrackletData tdata12;
   TrackletData tdata13;
   TrackletData tdata23;
   TrackletData tdata14;
   TrackletData tdata15;
   TrackletData tdata45;

   setTrackletTreeBranch(trackletTree12, tdata12);
   setTrackletTreeBranch(trackletTree13, tdata13);
   setTrackletTreeBranch(trackletTree23, tdata23);
   setTrackletTreeBranch(trackletTree14, tdata14);
   setTrackletTreeBranch(trackletTree15, tdata15);
   setTrackletTreeBranch(trackletTree45, tdata45);

   // TH3D* hLayer1Hit = new TH3D("hLayer1Hit", "", 75, 0, 15, 60, -3, 3, 64, -3.2, 3.2);
   // TH3D* hLayer2Hit = new TH3D("hLayer2Hit", "", 75, 0, 15, 60, -3, 3, 64, -3.2, 3.2);
   // TH3D* hLayer3Hit = new TH3D("hLayer3Hit", "", 75, 0, 15, 60, -3, 3, 64, -3.2, 3.2);

   // Prepare hit spectra for random hit
   // cout << "Projecting...1" << endl;
   // if (addL1Bck!=0) t->Project("hLayer1Hit", "phi1:eta1:r1");
   // cout << "Projecting...2" << endl;
   // if (addL2Bck!=0) t->Project("hLayer2Hit", "phi2:eta2:r2");
   // cout << "Projecting...3" << endl;
   // if (addL3Bck!=0) t->Project("hLayer3Hit", "phi3:eta3:r3");
   // cout << "Projecting...done" << endl;

   if (t->FindBranch("npart")!=0) {
      isMC = true;
      cout << "This is a Monte Carlo study." << endl;
      vzShift = -0.4847;
      cout << "vzShift = " << vzShift << endl;
   } else {
      cout << "This is a data analysis." << endl;
      smearVertex = 0;
   }

   // Event record
   vector<int> events[500];
   // Parameters for the tree =================================================
   Parameters par;
   getPixelTreeBranch(t, par);
   // Parameters beamHaloPar;
   // if (putBeamHalo) getPixelTreeBranch(beamHaloTree, beamHaloPar);
   cout << "Number of Events: " << t->GetEntries() << endl;

   int nBeamHalo = 0;
   int nPileUp = 0;
   if (pileUp != 0)
      cout << "Do pileup! With probability of " << pileUp << endl;

   TF1* csfitf = new TF1("csfit", csfit, -4, 4, 1);
   csfitf->SetParameters(1,8882, 0);
   csfitf->SetParLimits(0, 1.8882, 1.8882);

   // Main loop ===============================================================
   for (int i=startEntry; i<t->GetEntries()&&i<endEntry; i=i+1+nPileUp) {
      t->GetEntry(i);
      if (i % 1000 == 0) {
         cout << "Run " << par.nRun << " Event " << i << " "
              << trackletTree12->GetEntries() << " "
              << trackletTree13->GetEntries() << " "
              << trackletTree23->GetEntries() << " "
              << trackletTree14->GetEntries() << " "
              << trackletTree15->GetEntries() << " "
              << trackletTree45->GetEntries()
              << " Add Beam Halo: " << nBeamHalo << " " << nBeamHalo/(double)i
              << endl;
         if (reWeight) cout << "Reweighted!!!!!!!" << endl;
      }

      // bool flagDuplicateEvent = 0;
      // if (checkDuplicateEvent) {
      //    for (unsigned int j=0; j<events[par.nLumi].size(); j++) {
      //       if (par.nEv==events[par.nLumi][j]) {
      //          flagDuplicateEvent = 1;
      //          continue;
      //       }
      //    }
      //    if (!flagDuplicateEvent) events[par.nLumi].push_back(par.nEv);
      // }
      // if (flagDuplicateEvent) continue;

      // bool reWeightDropFlag = 0;
      // // Reweight MC vertex distribution to be the same as data
      // if (reWeight) {
      //    reWeightDropFlag = 0;
      //    double myVz = par.vz[1];
      //    if (myVz<-90) {
      //       TF1 *f = new TF1("f", "gaus", -30, 30);
      //       f->SetParameters(1, -0.6536, 4.438);
      //       myVz = f->GetRandom();
      //       delete f;
      //    }

      //    // for early data 900 GeV
      //    // double MCPdf = TMath::Gaus(myVz,-2.709,4.551,1);
      //    // double DataPdf = TMath::Gaus(myVz,-2.702,3.627,1);

      //    // for early data 7000 GeV Run 132440
      //    double MCPdf = TMath::Gaus(myVz, -0.6536, 4.438, 1);
      //    double DataPdf = TMath::Gaus(myVz, 0.3533-vzShift, 2.161, 1);

      //    // double DataPdf = TMath::Gaus(myVz,-0.4623,2.731,1);

      //    double Ratio = DataPdf / MCPdf;
      //    double x = gRandom->Rndm()*2.5;

      //    if (x>Ratio) reWeightDropFlag = 1;
      // }
      // if (reWeightDropFlag) continue;

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

      // Selection on Events

      // Fill reco vertex information
      tdata12.nv = par.nv+1;
      tdata13.nv = par.nv+1;
      tdata23.nv = par.nv+1;
      tdata14.nv = par.nv+1;
      tdata15.nv = par.nv+1;
      tdata45.nv = par.nv+1;
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
         tdata14.vz[j+1] = par.vz[j];
         tdata15.vz[j+1] = par.vz[j];
         tdata45.vz[j+1] = par.vz[j];
         tdata14.vx[j+1] = par.vx[j];
         tdata15.vx[j+1] = par.vx[j];
         tdata45.vx[j+1] = par.vx[j];
         tdata14.vy[j+1] = par.vy[j];
         tdata15.vy[j+1] = par.vy[j];
         tdata45.vy[j+1] = par.vy[j];
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
      tdata14.vx[0] = par.vx[0];
      tdata15.vx[0] = par.vx[0];
      tdata45.vx[0] = par.vx[0];
      tdata14.vy[0] = par.vy[0];
      tdata15.vy[0] = par.vy[0];
      tdata45.vy[0] = par.vy[0];
      tdata14.vz[0] = par.vz[0];
      tdata15.vz[0] = par.vz[0];
      tdata45.vz[0] = par.vz[0];

      // Add background hits
      // int bckHits = 0;
      // if (addL1Bck!=0 || addL2Bck!=0 || addL3Bck!=0) {
      //    for (int i=0; i<par.nhits1; i++)
      //       if (gRandom->Rndm()<0.01)
      //          bckHits++;
      // }

      // if (addL1Bck!=0) {
      //    for (int i=par.nhits1; i<par.nhits1+bckHits; i++) {
      //       double eta, phi, r;
      //       hLayer1Hit->GetRandom3(r, eta, phi);
      //       par.eta1[i] = eta;
      //       par.phi1[i] = phi;
      //       par.r1[i] = r;
      //    }
      //    par.nhits1 += bckHits;
      // }

      // if (addL2Bck!=0) {
      //    for (int i=par.nhits2; i<par.nhits2+bckHits; i++) {
      //       double eta, phi, r;
      //       hLayer2Hit->GetRandom3(r, eta, phi);
      //       par.eta2[i] = eta;
      //       par.phi2[i] = phi;
      //       par.r2[i] = r;
      //    }
      //    par.nhits2 += bckHits;
      // }

      // if (addL3Bck!=0) {
      //    for (int i=par.nhits3; i<par.nhits3+bckHits; i++) {
      //       double eta, phi, r;
      //       hLayer3Hit->GetRandom3(r, eta, phi);
      //       par.eta3[i] = eta;
      //       par.phi3[i] = phi;
      //       par.r3[i] = r;
      //    }
      //    par.nhits3 += bckHits;
      // }

      // Add trackletVertex
      // if (tdata12.nv==2) tdata12.nv = 3;
      // if (tdata13.nv==2) tdata13.nv = 3;
      // if (tdata23.nv==2) tdata23.nv = 3;

      vector<RecoHit> layer1;
      prepareHits(layer1, par, cuts, 1, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      vector<RecoHit> layer2;
      prepareHits(layer2, par, cuts, 2, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      vector<RecoHit> layer3;
      prepareHits(layer3, par, cuts, 3, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      vector<RecoHit> layer4;
      prepareHits(layer4, par, cuts, 4, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      vector<RecoHit> layer5;
      prepareHits(layer5, par, cuts, 5, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);

      int recoPU = 0;
      double trackletVertex = -99;

      vector<RecoHit> allhits;
      prepareHits(allhits, par, cuts, 1, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      prepareHits(allhits, par, cuts, 2, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);

      if (pileUp!=0) {
         nPileUp = gRandom->Poisson(pileUp);
         for (int p=1; p<=nPileUp; p++) {
            t->GetEntry(i+p);
            prepareHits(allhits, par, cuts, 1, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
            prepareHits(allhits, par, cuts, 2, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
         }
         t->GetEntry(i);
      }

      std::sort(allhits.begin(), allhits.end(), sortphi);
      for (unsigned int m=0; m<allhits.size() && allhits[m].phi<0.08-_PI; m++) {
         allhits.push_back(allhits[m]);
         allhits.back().phi += 2*_PI;
      }
      RecoHit fake(0, 4, 0, 0, 1);
      allhits.push_back(fake);

      std::vector<Vertex> vertices;
      std::deque<RecoHit> trackcands[2];
      for (unsigned int a=0; a<allhits.size(); a++) {
         for (int q=0; q<2; q++) {
            while (!trackcands[q].empty() && allhits[a].phi>trackcands[q].front().phi+0.08) {
               for (std::deque<RecoHit>::iterator it = trackcands[1-q].begin(); it != trackcands[1-q].end(); it++) {
                  Vertex vertex;
                  double r1 = trackcands[q].front().r;
                  double z1 = r1/tan(2*atan(exp(-trackcands[q].front().eta)));
                  double r2 = (*it).r;
                  double z2 = r2/tan(2*atan(exp(-(*it).eta)));
                  vertex.vz = z1-(z2-z1)/(r2-r1)*r1;
                  if (fabs(vertex.vz)<20.0)
                     vertices.push_back(vertex);
               }
               trackcands[q].pop_front();
            }
         }
         trackcands[allhits[a].layer-1].push_back(allhits[a]);
      }

      // Using forward pixel hits =============================================
      // vector<RecoHit> fallhits;
      // prepareHits(fallhits, par, cuts, 1, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      // prepareHits(fallhits, par, cuts, 4, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      // prepareHits(fallhits, par, cuts, 5, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);

      // if (pileUp!=0) {
      //    nPileUp = gRandom->Poisson(pileUp);
      //    for (int p=1; p<=nPileUp; p++) {
      //       t->GetEntry(i+p);
      //       prepareHits(fallhits, par, cuts, 1, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      //       prepareHits(fallhits, par, cuts, 4, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      //       prepareHits(fallhits, par, cuts, 5, 0, 0, 0, splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      //    }
      //    t->GetEntry(i);
      // }

      // std::sort(fallhits.begin(), fallhits.end(), sortphi);
      // for (unsigned int m=0; m<fallhits.size() && fallhits[m].phi<0.08-_PI; m++) {
      //    fallhits.push_back(fallhits[m]);
      //    fallhits.back().phi += 2*_PI;
      // }
      // fallhits.push_back(fake);

      // std::deque<RecoHit> ftrackcands[3];
      // for (unsigned int a=0; a<fallhits.size(); a++) {
      //    for (int q=0; q<3; q++) {
      //       while (!ftrackcands[q].empty() && fallhits[a].phi>ftrackcands[q].front().phi+0.08) {
      //          for (int w=0; w<3; w++) {
      //             if (w == q || (q & w)) continue;
      //             for (std::deque<RecoHit>::iterator it = ftrackcands[w].begin(); it != ftrackcands[w].end(); it++) {
      //                Vertex vertex;
      //                double r1 = ftrackcands[q].front().r;
      //                double z1 = r1/tan(2*atan(exp(-ftrackcands[q].front().eta)));
      //                double r2 = (*it).r;
      //                double z2 = r2/tan(2*atan(exp(-(*it).eta)));
      //                vertex.vz = z1-(z2-z1)/(r2-r1)*r1;
      //                if (fabs(vertex.vz)<20.0)
      //                   vertices.push_back(vertex);
      //             }
      //          }
      //          ftrackcands[q].pop_front();
      //       }
      //    }
      //    int l = fallhits[a].layer - 1;
      //    if (l == 3) l = 1;
      //    if (l == 4) l = 2;
      //    ftrackcands[l].push_back(fallhits[a]);
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
            if (vertices[c].vzmean!=candidates.back().vzmean)
               candidates.push_back(vertices[c]);
         }

         std::vector<Vertex> pileupcands;
         pileupcands.push_back(vertices[0]);
         for (unsigned int p=1; p<vertices.size(); p++) {
            if (vertices[p].nz < 10)
               break;
            if (vertices[p].nz < pileupcands[0].nz*0.75)
               break;
            std::vector<Vertex>::iterator it = pileupcands.begin();
            for (; it!=pileupcands.end() && abs((*it).vzmean-vertices[p].vzmean)>2.0; it++);
            if (it!=pileupcands.end())
               continue;
            pileupcands.push_back(vertices[p]);
         }
         recoPU = pileupcands.size();

         // Use chi2 and sigma2 ===============================================
         // for (unsigned int d=0; d<candidates.size(); d++) {
         //    TGraph* csveta_g = new TGraph(layer1.size()+layer2.size());

         //    int g = 0;
         //    for (unsigned int c1=0; c1<layer1.size(); c1++) {
         //       double r1 = layer1[c1].r;
         //       double z1 = r1/tan(2*atan(exp(-layer1[c1].eta)));
         //       double neweta = (z1 - candidates[d].vzmean > 0) ? -log(tan(atan(r1/(z1-candidates[d].vzmean))/2)) : log(tan(atan(r1/(candidates[d].vzmean-z1))/2));
         //       if (layer1[c1].cs<24 && layer1[c1].cs>=0.9*cosh(neweta))
         //          csveta_g->SetPoint(g++, neweta, layer1[c1].cs);
         //    }
         //    for (unsigned int c2=0; c2<layer2.size(); c2++) {
         //       double r2 = layer2[c2].r;
         //       double z2 = r2/tan(2*atan(exp(-layer2[c2].eta)));
         //       double neweta = (z2 - candidates[d].vzmean > 0) ? -log(tan(atan(r2/(z2-candidates[d].vzmean))/2)) : log(tan(atan(r2/(candidates[d].vzmean-z2))/2));
         //       if (layer2[c2].cs<24)
         //          csveta_g->SetPoint(g++, neweta, layer2[c2].cs);
         //    }
         //    csveta_g->Set(g);
         //    if (g) {
         //       TFitResultPtr csveta_g_r = csveta_g->Fit("csfit", "WBQS");
         //       candidates[d].chi2 = csveta_g_r->Chi2();
         //       candidates[d].par0 = csveta_g_r->Parameter(0);
         //    } else {
         //       candidates[d].chi2 = 9999.9;
         //       candidates[d].par0 = 9999.9;
         //    }
         // }

         // std::vector<Vertex>::iterator eqnz = candidates.begin();
         // for (; eqnz!=candidates.end() && (*eqnz).nz==candidates[0].nz; eqnz++);
         // std::sort(candidates.begin(), eqnz, sortchi2);
         // std::sort(eqnz, candidates.end(), sortchi2);

         // std::vector<Vertex>::iterator eqchi2 = candidates.begin();
         // for (; eqchi2!=eqnz && (*eqchi2).chi2/candidates[0].chi2<1.12; eqchi2++);
         // std::sort(candidates.begin(), eqchi2, sortsigma2);

         // std::vector<Vertex>::iterator nzminusone = eqnz;
         // for (; nzminusone!=candidates.end() && (*nzminusone).chi2<candidates[0].chi2/2; nzminusone++);
         // if (eqnz != nzminusone) {
         //    std::stable_sort(eqnz, nzminusone, sortsigma2);
         //    candidates.insert(candidates.begin(), *eqnz);
         // }

         // Use sigma2 only ===================================================
         std::vector<Vertex>::iterator eqnz = candidates.begin();
         for (; eqnz!=candidates.end() && (*eqnz).nz==candidates[0].nz; eqnz++);
         std::sort(candidates.begin(), eqnz, sortsigma2);

         trackletVertex = candidates[0].vzmean;
      }

      // For particle gun ==========
      //    if (i%1000 == 0) cout << "!!! USE GEN VERTEX (FOR PARTICLE GUN) " << endl;
      //    trackletVertex = par.vz[0];
      // ===========================
      // vz[1] is always the selected algorithm

      //if (useKKVertex) trackletVertex = par.vz[1];
      // if (useKKVertex) trackletVertex = 0;

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

      tdata12.vz[1] = trackletVertex;
      tdata13.vz[1] = trackletVertex;
      tdata23.vz[1] = trackletVertex;
      tdata14.vz[1] = trackletVertex;
      tdata15.vz[1] = trackletVertex;
      tdata45.vz[1] = trackletVertex;

      if (useKKVertex) {
         if (i==0) cout << "Use Reconstructed Vertex " << endl;
         if (i==0) cout << par.vz[1] << endl;
         tdata12.vx[1] = par.vx[1];
         tdata12.vy[1] = par.vy[1];
         tdata13.vx[1] = par.vx[1];
         tdata13.vy[1] = par.vy[1];
         tdata23.vx[1] = par.vx[1];
         tdata23.vy[1] = par.vy[1];
         tdata14.vx[1] = par.vx[1];
         tdata14.vy[1] = par.vy[1];
         tdata15.vx[1] = par.vx[1];
         tdata15.vy[1] = par.vy[1];
         tdata45.vx[1] = par.vx[1];
         tdata45.vy[1] = par.vy[1];
      } else {
         if (i==0) cout << "Use Tracklet Vertex " << endl;
      }

      // if (useRandomVertex) {
      //    if (i==1) cout << "Random Vertex!!!" << endl;
      //    tdata12.vz[1] = gRandom->Rndm()*40-20;
      //    tdata13.vz[1] = tdata12.vz[1];
      //    tdata23.vz[1] = tdata12.vz[1];
      // }
      // // Use trackletVertex
      // if (fabs(tdata12.vz[1])>cuts.vzCut && makeVzCut) continue;

      vector<RecoHit> layer1Cut;
      prepareHits(layer1Cut, par, cuts, 1, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, 1, par.nRun, par.nLumi, 0);

      // Process hits with Vz constraint:
      // vector<RecoHit> layer1;
      // prepareHits(layer1, par, cuts, 1, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      // vector<RecoHit> layer2;
      // prepareHits(layer2, par, cuts, 2, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      // vector<RecoHit> layer3;
      // prepareHits(layer3, par, cuts, 3, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);

      // if (nPileUp>0) {
      //    for (int p=1; p<=nPileUp; p++) {
      //       t->GetEntry(i+p);
      //       prepareHits(layer1, par, cuts, 1, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize);
      //       prepareHits(layer2, par, cuts, 2, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize);
      //       prepareHits(layer3, par, cuts, 3, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize);
      //    }
      //    t->GetEntry(i);
      // }

      std::vector<RecoHit> combinedhits;
      prepareHits(combinedhits, par, cuts, 1, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      prepareHits(combinedhits, par, cuts, 2, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      std::sort(combinedhits.begin(), combinedhits.end(), sorteta);
      std::vector<Tracklet> recoTracklets12;
      recoTracklets12 = recoTracklets(combinedhits, 12);

      combinedhits.clear();
      prepareHits(combinedhits, par, cuts, 1, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      prepareHits(combinedhits, par, cuts, 3, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      std::sort(combinedhits.begin(), combinedhits.end(), sorteta);
      std::vector<Tracklet> recoTracklets13;
      recoTracklets13 = recoTracklets(combinedhits, 13);

      combinedhits.clear();
      prepareHits(combinedhits, par, cuts, 2, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      prepareHits(combinedhits, par, cuts, 3, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      std::sort(combinedhits.begin(), combinedhits.end(), sorteta);
      std::vector<Tracklet> recoTracklets23;
      recoTracklets23 = recoTracklets(combinedhits, 23);

      combinedhits.clear();
      prepareHits(combinedhits, par, cuts, 1, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      prepareHits(combinedhits, par, cuts, 4, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      std::sort(combinedhits.begin(), combinedhits.end(), sorteta);
      std::vector<Tracklet> recoTracklets14;
      recoTracklets14 = recoTracklets(combinedhits, 14);

      combinedhits.clear();
      prepareHits(combinedhits, par, cuts, 1, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      prepareHits(combinedhits, par, cuts, 5, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      std::sort(combinedhits.begin(), combinedhits.end(), sorteta);
      std::vector<Tracklet> recoTracklets15;
      recoTracklets15 = recoTracklets(combinedhits, 15);

      combinedhits.clear();
      prepareHits(combinedhits, par, cuts, 4, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      prepareHits(combinedhits, par, cuts, 5, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], splitProb, dropProb, cutOnClusterSize, par.nRun, par.nLumi, smearPixels);
      std::sort(combinedhits.begin(), combinedhits.end(), sorteta);
      std::vector<Tracklet> recoTracklets45;
      recoTracklets45 = recoTracklets(combinedhits, 45);

      // Form Tracklets
      // vector<Tracklet> protoTracklets12;
      // vector<Tracklet> protoTracklets13;
      // vector<Tracklet> protoTracklets23;
      // vector<Tracklet> recoTracklets12;
      // vector<Tracklet> recoTracklets13;
      // vector<Tracklet> recoTracklets23;

      // if (mimicPixelCounting) {
      //    protoTracklets12 = recoProtoTracklets(layer1, layer1);
      //    protoTracklets13 = recoProtoTracklets(layer2, layer2);
      //    protoTracklets23 = recoProtoTracklets(layer3, layer3);
      // } else {
      //    protoTracklets12 = recoProtoTracklets(layer1, layer2);
      //    protoTracklets13 = recoProtoTracklets(layer1, layer3);
      //    protoTracklets23 = recoProtoTracklets(layer2, layer3);
      // }

      // recoTracklets12 = cleanTracklets(protoTracklets12, 0, cuts);
      // recoTracklets13 = cleanTracklets(protoTracklets13, 0, cuts);
      // recoTracklets23 = cleanTracklets(protoTracklets23, 0, cuts);

      // Move the Vertex back
      // if (smearVertex!=0) {
      //    tdata12.vz[1] = trackletVertex - smear;
      //    tdata13.vz[1] = trackletVertex - smear;
      //    tdata23.vz[1] = trackletVertex - smear;
      // }

      // Vertex Compatibility information

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
      tdata12.nL1ABit    = par.nL1ABit;
      tdata12.nL1TBit    = par.nL1TBit;
      tdata12.xi         = par.xi;
      tdata12.passDS     = par.passDS;
      tdata12.passSingleTrack = par.passSingleTrack;
      tdata12.ntrks      = par.ntrks;
      tdata12.ntrksCut   = par.ntrksCut;
      tdata12.nPU        = nPileUp;
      tdata12.recoPU     = recoPU;

      for (int j=0; j<(int)par.nHltBit; j++)
         tdata12.hltBit[j] = par.hltBit[j];
      for (int j=0; j<(int)par.nL1ABit; j++)
         tdata12.l1ABit[j] = par.l1ABit[j];
      for (int j=0; j<(int)par.nL1TBit; j++)
         tdata12.l1TBit[j] = par.l1TBit[j];

      int ntracklet12s = 0;
      int ntracklet12b = 0;
      for (int j=0; j<(int)tdata12.nTracklet; j++) {
         tdata12.eta1[j] = recoTracklets12[j].eta1();
         tdata12.eta2[j] = recoTracklets12[j].eta2();
         tdata12.r1[j]   = recoTracklets12[j].r1();
         tdata12.r2[j]   = recoTracklets12[j].r2();
         tdata12.cs1[j]  = recoTracklets12[j].cs1();
         tdata12.cs2[j]  = recoTracklets12[j].cs2();
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

      tdata13.nTracklet  = recoTracklets13.size();
      tdata13.nhit1      = layer1.size();
      tdata13.nhit2      = layer3.size();
      tdata13.nRun       = par.nRun;
      tdata13.nEv        = par.nEv;
      tdata13.nLumi      = par.nLumi;
      tdata13.nBX        = par.nBX;
      tdata13.nHFn       = par.nHFp;
      tdata13.nHFp       = par.nHFn;
      tdata13.nHits      = layer1.size() + layer2.size() + layer3.size();
      tdata13.nL1ABit    = par.nL1ABit;
      tdata13.nL1TBit    = par.nL1TBit;
      tdata13.xi         = par.xi;
      tdata13.passDS     = par.passDS;
      tdata13.passSingleTrack = par.passSingleTrack;
      tdata13.ntrks      = par.ntrks;
      tdata13.ntrksCut   = par.ntrksCut;
      tdata13.nPU        = nPileUp;
      tdata13.recoPU     = recoPU;

      for (int j=0; j<(int)par.nHltBit; j++)
         tdata13.hltBit[j] = par.hltBit[j];
      for (int j=0; j<(int)par.nL1ABit; j++)
         tdata13.l1ABit[j] = par.l1ABit[j];
      for (int j=0; j<(int)par.nL1TBit; j++)
         tdata13.l1TBit[j] = par.l1TBit[j];

      int ntracklet13s = 0;
      int ntracklet13b = 0;
      for (int j=0; j<(int)tdata13.nTracklet; j++) {
         tdata13.eta1[j] = recoTracklets13[j].eta1();
         tdata13.eta2[j] = recoTracklets13[j].eta2();
         tdata13.r1[j]   = recoTracklets13[j].r1();
         tdata13.r2[j]   = recoTracklets13[j].r2();
         tdata13.cs1[j]  = recoTracklets13[j].cs1();
         tdata13.cs2[j]  = recoTracklets13[j].cs2();
         tdata13.phi1[j] = recoTracklets13[j].phi1();
         tdata13.phi2[j] = recoTracklets13[j].phi2();
         tdata13.deta[j] = recoTracklets13[j].deta();
         tdata13.dphi[j] = recoTracklets13[j].dphi();
         if (fabs(tdata13.deta[j])<0.1) {
            if (fabs(tdata13.dphi[j])<1.0)
               ntracklet13s++;
            if (fabs(tdata13.dphi[j])>1.0 && fabs(tdata13.dphi[j])<2.0)
               ntracklet13b++;
         }
      }
      tdata13.mult = ntracklet13s - ntracklet13b;
      tdata13.mult2 = layer1Cut.size();
      tdata13.npart = 0;
      for (int j=0; j<12; j++) tdata13.nhad[j] = 0;
      for (int j=0; j<par.npart; j++) {
         if (fabs(par.eta[j])>3 || par.chg[j]==0 || fabs(par.pdg[j])==11 || fabs(par.pdg[j])==13) continue;
         tdata13.eta[tdata13.npart] = par.eta[j];
         tdata13.phi[tdata13.npart] = par.phi[j];
         tdata13.chg[tdata13.npart] = par.chg[j];
         tdata13.pdg[tdata13.npart] = par.pdg[j];
         tdata13.pt[tdata13.npart] = par.pt[j];
         tdata13.npart++;
         int bin = (int)((par.eta[j]+3)*2);
         int pdg = (int)abs(par.pdg[j]);
         if (pdg==211 || pdg==321 || pdg==2212 || pdg==3122) tdata13.nhad[bin]++;
      }

      tdata13.evtType = par.evtType;
      tdata13.pro2 = pro2;
      for (int j=0; j<par.nv; j++)
         tdata13.vz[j] += vzShift;

      trackletTree13->Fill();

      tdata23.nTracklet  = recoTracklets23.size();
      tdata23.nhit1      = layer2.size();
      tdata23.nhit2      = layer3.size();
      tdata23.nRun       = par.nRun;
      tdata23.nEv        = par.nEv;
      tdata23.nLumi      = par.nLumi;
      tdata23.nBX        = par.nBX;
      tdata23.nHFn       = par.nHFp;
      tdata23.nHFp       = par.nHFn;
      tdata23.nHits      = layer1.size() + layer2.size() + layer3.size();
      tdata23.nL1ABit    = par.nL1ABit;
      tdata23.nL1TBit    = par.nL1TBit;
      tdata23.xi         = par.xi;
      tdata23.passDS     = par.passDS;
      tdata23.passSingleTrack = par.passSingleTrack;
      tdata23.ntrks      = par.ntrks;
      tdata23.ntrksCut   = par.ntrksCut;
      tdata23.nPU        = nPileUp;
      tdata23.recoPU     = recoPU;

      for (int j=0; j<(int)par.nHltBit; j++)
         tdata23.hltBit[j] = par.hltBit[j];
      for (int j=0; j<(int)par.nL1ABit; j++)
         tdata23.l1ABit[j] = par.l1ABit[j];
      for (int j=0; j<(int)par.nL1TBit; j++)
         tdata23.l1TBit[j] = par.l1TBit[j];

      int ntracklet23s = 0;
      int ntracklet23b = 0;
      for (int j=0; j<(int)tdata23.nTracklet; j++) {
         tdata23.eta1[j] = recoTracklets23[j].eta1();
         tdata23.eta2[j] = recoTracklets23[j].eta2();
         tdata23.r1[j]   = recoTracklets23[j].r1();
         tdata23.r2[j]   = recoTracklets23[j].r2();
         tdata23.cs1[j]  = recoTracklets23[j].cs1();
         tdata23.cs2[j]  = recoTracklets23[j].cs2();
         tdata23.phi1[j] = recoTracklets23[j].phi1();
         tdata23.phi2[j] = recoTracklets23[j].phi2();
         tdata23.deta[j] = recoTracklets23[j].deta();
         tdata23.dphi[j] = recoTracklets23[j].dphi();
         if (fabs(tdata23.deta[j])<0.1) {
            if (fabs(tdata23.dphi[j])<1.0)
               ntracklet23s++;
            if (fabs(tdata23.dphi[j])>1.0 && fabs(tdata23.dphi[j])<2.0)
               ntracklet23b++;
         }
      }
      tdata23.mult = ntracklet23s - ntracklet23b;
      tdata23.mult2 = layer1Cut.size();
      tdata23.npart = 0;
      for (int j=0; j<12; j++) tdata23.nhad[j] = 0;
      for (int j=0; j<par.npart; j++) {
         if (fabs(par.eta[j])>3 || par.chg[j]==0 || fabs(par.pdg[j])==11 || abs(par.pdg[j])==13) continue;
         tdata23.eta[tdata23.npart] = par.eta[j];
         tdata23.phi[tdata23.npart] = par.phi[j];
         tdata23.chg[tdata23.npart] = par.chg[j];
         tdata23.pdg[tdata23.npart] = par.pdg[j];
         tdata23.pt[tdata23.npart] = par.pt[j];
         tdata23.npart++;
         int bin = (int)((par.eta[j]+3)*2);
         int pdg = (int)abs(par.pdg[j]);
         if (pdg==211 || pdg==321 || pdg==2212 || pdg==3122) tdata23.nhad[bin]++;
      }

      tdata23.evtType = par.evtType;
      tdata23.pro2 = pro2;
      for (int j=0; j<par.nv; j++)
         tdata23.vz[j] += vzShift;

      trackletTree23->Fill();

      tdata14.nTracklet  = recoTracklets14.size();
      tdata14.nhit1      = layer1.size();
      tdata14.nhit2      = layer4.size();
      tdata14.nRun       = par.nRun;
      tdata14.nEv        = par.nEv;
      tdata14.nLumi      = par.nLumi;
      tdata14.nBX        = par.nBX;
      tdata14.nHFn       = par.nHFp;
      tdata14.nHFp       = par.nHFn;
      tdata14.nHits      = 0;
      tdata14.nL1ABit    = par.nL1ABit;
      tdata14.nL1TBit    = par.nL1TBit;
      tdata14.xi         = par.xi;
      tdata14.passDS     = par.passDS;
      tdata14.passSingleTrack = par.passSingleTrack;
      tdata14.ntrks      = par.ntrks;
      tdata14.ntrksCut   = par.ntrksCut;
      tdata14.nPU        = nPileUp;
      tdata14.recoPU     = recoPU;

      for (int j=0; j<(int)par.nHltBit; j++)
         tdata14.hltBit[j] = par.hltBit[j];
      for (int j=0; j<(int)par.nL1ABit; j++)
         tdata14.l1ABit[j] = par.l1ABit[j];
      for (int j=0; j<(int)par.nL1TBit; j++)
         tdata14.l1TBit[j] = par.l1TBit[j];

      int ntracklet14s = 0;
      int ntracklet14b = 0;
      for (int j=0; j<(int)tdata14.nTracklet; j++) {
         tdata14.eta1[j] = recoTracklets14[j].eta1();
         tdata14.eta2[j] = recoTracklets14[j].eta2();
         tdata14.r1[j]   = recoTracklets14[j].r1();
         tdata14.r2[j]   = recoTracklets14[j].r2();
         tdata14.cs1[j]  = recoTracklets14[j].cs1();
         tdata14.cs2[j]  = recoTracklets14[j].cs2();
         tdata14.phi1[j] = recoTracklets14[j].phi1();
         tdata14.phi2[j] = recoTracklets14[j].phi2();
         tdata14.deta[j] = recoTracklets14[j].deta();
         tdata14.dphi[j] = recoTracklets14[j].dphi();
         if (fabs(tdata14.deta[j])<0.1) {
            if (fabs(tdata14.dphi[j])<1.0)
               ntracklet14s++;
            if (fabs(tdata14.dphi[j])>1.0 && fabs(tdata14.dphi[j])<2.0)
               ntracklet14b++;
         }
      }
      tdata14.mult = ntracklet14s - ntracklet14b;
      tdata14.mult2 = layer1Cut.size();
      tdata14.npart = 0;
      for (int j=0; j<12; j++) tdata14.nhad[j] = 0;
      for (int j=0; j<par.npart; j++) {
         if (fabs(par.eta[j])>3 || par.chg[j]==0 || fabs(par.pdg[j])==11 || abs(par.pdg[j])==13) continue;
         tdata14.eta[tdata14.npart] = par.eta[j];
         tdata14.phi[tdata14.npart] = par.phi[j];
         tdata14.chg[tdata14.npart] = par.chg[j];
         tdata14.pdg[tdata14.npart] = par.pdg[j];
         tdata14.pt[tdata14.npart] = par.pt[j];
         tdata14.npart++;
         int bin = (int)((par.eta[j]+3)*2);
         int pdg = (int)abs(par.pdg[j]);
         if (pdg==211 || pdg==321 || pdg==2212 || pdg==3122) tdata14.nhad[bin]++;
      }

      tdata14.evtType = par.evtType;
      tdata14.pro2 = pro2;
      for (int j=0; j<par.nv; j++)
         tdata14.vz[j] += vzShift;

      trackletTree14->Fill();

      tdata15.nTracklet  = recoTracklets15.size();
      tdata15.nhit1      = layer1.size();
      tdata15.nhit2      = layer5.size();
      tdata15.nRun       = par.nRun;
      tdata15.nEv        = par.nEv;
      tdata15.nLumi      = par.nLumi;
      tdata15.nBX        = par.nBX;
      tdata15.nHFn       = par.nHFp;
      tdata15.nHFp       = par.nHFn;
      tdata15.nHits      = 0;
      tdata15.nL1ABit    = par.nL1ABit;
      tdata15.nL1TBit    = par.nL1TBit;
      tdata15.xi         = par.xi;
      tdata15.passDS     = par.passDS;
      tdata15.passSingleTrack = par.passSingleTrack;
      tdata15.ntrks      = par.ntrks;
      tdata15.ntrksCut   = par.ntrksCut;
      tdata15.nPU        = nPileUp;
      tdata15.recoPU     = recoPU;

      for (int j=0; j<(int)par.nHltBit; j++)
         tdata15.hltBit[j] = par.hltBit[j];
      for (int j=0; j<(int)par.nL1ABit; j++)
         tdata15.l1ABit[j] = par.l1ABit[j];
      for (int j=0; j<(int)par.nL1TBit; j++)
         tdata15.l1TBit[j] = par.l1TBit[j];

      int ntracklet15s = 0;
      int ntracklet15b = 0;
      for (int j=0; j<(int)tdata15.nTracklet; j++) {
         tdata15.eta1[j] = recoTracklets15[j].eta1();
         tdata15.eta2[j] = recoTracklets15[j].eta2();
         tdata15.r1[j]   = recoTracklets15[j].r1();
         tdata15.r2[j]   = recoTracklets15[j].r2();
         tdata15.cs1[j]  = recoTracklets15[j].cs1();
         tdata15.cs2[j]  = recoTracklets15[j].cs2();
         tdata15.phi1[j] = recoTracklets15[j].phi1();
         tdata15.phi2[j] = recoTracklets15[j].phi2();
         tdata15.deta[j] = recoTracklets15[j].deta();
         tdata15.dphi[j] = recoTracklets15[j].dphi();
         if (fabs(tdata15.deta[j])<0.1) {
            if (fabs(tdata15.dphi[j])<1.0)
               ntracklet15s++;
            if (fabs(tdata15.dphi[j])>1.0 && fabs(tdata15.dphi[j])<2.0)
               ntracklet15b++;
         }
      }
      tdata15.mult = ntracklet15s - ntracklet15b;
      tdata15.mult2 = layer1Cut.size();
      tdata15.npart = 0;
      for (int j=0; j<12; j++) tdata15.nhad[j] = 0;
      for (int j=0; j<par.npart; j++) {
         if (fabs(par.eta[j])>3 || par.chg[j]==0 || fabs(par.pdg[j])==11 || abs(par.pdg[j])==13) continue;
         tdata15.eta[tdata15.npart] = par.eta[j];
         tdata15.phi[tdata15.npart] = par.phi[j];
         tdata15.chg[tdata15.npart] = par.chg[j];
         tdata15.pdg[tdata15.npart] = par.pdg[j];
         tdata15.pt[tdata15.npart] = par.pt[j];
         tdata15.npart++;
         int bin = (int)((par.eta[j]+3)*2);
         int pdg = (int)abs(par.pdg[j]);
         if (pdg==211 || pdg==321 || pdg==2212 || pdg==3122) tdata15.nhad[bin]++;
      }

      tdata15.evtType = par.evtType;
      tdata15.pro2 = pro2;
      for (int j=0; j<par.nv; j++)
         tdata15.vz[j] += vzShift;

      trackletTree15->Fill();

      tdata45.nTracklet  = recoTracklets45.size();
      tdata45.nhit1      = layer4.size();
      tdata45.nhit2      = layer5.size();
      tdata45.nRun       = par.nRun;
      tdata45.nEv        = par.nEv;
      tdata45.nLumi      = par.nLumi;
      tdata45.nBX        = par.nBX;
      tdata45.nHFn       = par.nHFp;
      tdata45.nHFp       = par.nHFn;
      tdata45.nHits      = 0;
      tdata45.nL1ABit    = par.nL1ABit;
      tdata45.nL1TBit    = par.nL1TBit;
      tdata45.xi         = par.xi;
      tdata45.passDS     = par.passDS;
      tdata45.passSingleTrack = par.passSingleTrack;
      tdata45.ntrks      = par.ntrks;
      tdata45.ntrksCut   = par.ntrksCut;
      tdata45.nPU        = nPileUp;
      tdata45.recoPU     = recoPU;

      for (int j=0; j<(int)par.nHltBit; j++)
         tdata45.hltBit[j] = par.hltBit[j];
      for (int j=0; j<(int)par.nL1ABit; j++)
         tdata45.l1ABit[j] = par.l1ABit[j];
      for (int j=0; j<(int)par.nL1TBit; j++)
         tdata45.l1TBit[j] = par.l1TBit[j];

      int ntracklet45s = 0;
      int ntracklet45b = 0;
      for (int j=0; j<(int)tdata45.nTracklet; j++) {
         tdata45.eta1[j] = recoTracklets45[j].eta1();
         tdata45.eta2[j] = recoTracklets45[j].eta2();
         tdata45.r1[j]   = recoTracklets45[j].r1();
         tdata45.r2[j]   = recoTracklets45[j].r2();
         tdata45.cs1[j]  = recoTracklets45[j].cs1();
         tdata45.cs2[j]  = recoTracklets45[j].cs2();
         tdata45.phi1[j] = recoTracklets45[j].phi1();
         tdata45.phi2[j] = recoTracklets45[j].phi2();
         tdata45.deta[j] = recoTracklets45[j].deta();
         tdata45.dphi[j] = recoTracklets45[j].dphi();
         if (fabs(tdata45.deta[j])<0.1) {
            if (fabs(tdata45.dphi[j])<1.0)
               ntracklet45s++;
            if (fabs(tdata45.dphi[j])>1.0 && fabs(tdata45.dphi[j])<2.0)
               ntracklet45b++;
         }
      }
      tdata45.mult = ntracklet45s - ntracklet45b;
      tdata45.mult2 = layer1Cut.size();
      tdata45.npart = 0;
      for (int j=0; j<12; j++) tdata45.nhad[j] = 0;
      for (int j=0; j<par.npart; j++) {
         if (fabs(par.eta[j])>3 || par.chg[j]==0 || fabs(par.pdg[j])==11 || abs(par.pdg[j])==13) continue;
         tdata45.eta[tdata45.npart] = par.eta[j];
         tdata45.phi[tdata45.npart] = par.phi[j];
         tdata45.chg[tdata45.npart] = par.chg[j];
         tdata45.pdg[tdata45.npart] = par.pdg[j];
         tdata45.pt[tdata45.npart] = par.pt[j];
         tdata45.npart++;
         int bin = (int)((par.eta[j]+3)*2);
         int pdg = (int)abs(par.pdg[j]);
         if (pdg==211 || pdg==321 || pdg==2212 || pdg==3122) tdata45.nhad[bin]++;
      }

      tdata45.evtType = par.evtType;
      tdata45.pro2 = pro2;
      for (int j=0; j<par.nv; j++)
         tdata45.vz[j] += vzShift;

      trackletTree45->Fill();
   }

   outf->Write();
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

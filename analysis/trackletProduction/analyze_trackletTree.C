#define _DZ 0.12
#define _DPHI 0.05

#define _PI 3.14159265358979

#include <deque>

#include "TFile.h"
#include "TNtuple.h"
#include "TTimeStamp.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TGraph.h"
#include "TFitResult.h"

#include "pdfs.h"
#include "tracklet.h"
#include "ClusterVertexCompat.h"

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

int analyze_trackletTree(const char* infile = "PixelTree.root", // Input PixelTree file
                         const char* outfile = "output.root",   // Ouptut Tracklet Tree
                         uint64_t start_entry = 0,              // Starting entry number in the PixelTree
                         uint64_t end_entry = 1000000000,       // Ending entry number in the PixelTree
                         bool use_random_vertex = 0,            // Use random vertex (instead of the reco one)
                         bool skip_cvqual = 0,                  // Skip calculation of the cluster-vertex compatibility
                         double pileup = 0,                     // Artifically overlap event to mimic pile-up
                         bool reweight_vertex = 1,              // Reweight vertex distribution to match data
                         bool reweight_mult = 1,                // Reweight multiplicity distribution
                         float add_bkg_l1 = 0,                  // Add random background to first pixel layer
                         float add_bkg_l2 = 0,                  // Add random background to second pixel layer
                         float add_bkg_l3 = 0,                  // Add random background to third pixel layer
                         double split_prob = 0,                 // Splitting probability of the pixel hit
                         double drop_prob = 0,                  // Emulate efficiency loss
                         bool smear_pixels = 0,                 // Smear pixel hits
                         bool cut_cluster_size = 0,             // Cut on cluster size to reduce background
                         bool do_two_hit_corr = 0,              // Create a pixel counting tree instead of tracklet tree
                         bool use_kk_vertex = 0,                // Use vertex from other reco vertex collection
                         double smear_vertex = 0)               // Add additional smearing to vertex position
{
   TTimeStamp myTime;
   gRandom->SetSeed(myTime.GetNanoSec());
   printf(" # init random: %f\n", gRandom->Rndm());

   TFile* inf = TFile::Open(infile);
   TTree* t = dynamic_cast<TTree*>(inf->Get("ana/PixelTree"));
   TTree* hltTree = (TTree*)inf->Get("hltanalysis/HltTree");

   TFile* outf = new TFile(outfile, "recreate");
   TTree* trackletTree12 = new TTree("TrackletTree12", "Tree of Reconstructed Tracklets");
   TTree* trackletTree13 = new TTree("TrackletTree13", "Tree of Reconstructed Tracklets");
   TTree* trackletTree23 = new TTree("TrackletTree23", "Tree of Reconstructed Tracklets");

   TrackletData tdata12;
   TrackletData tdata13;
   TrackletData tdata23;

   setTrackletTreeBranch(trackletTree12, tdata12);
   setTrackletTreeBranch(trackletTree13, tdata13);
   setTrackletTreeBranch(trackletTree23, tdata23);

   bool isMC = 0;
   double vzShift = 0;

   if (t->GetEntries("nRun<10")!=0) {
      isMC = true;
      vzShift = -0.551191;
#if defined(_EPOS_5TEV)
      pileup = 0.003;
#elif defined(_EPOS_8TEV) || defined(_HIJING_8TEV)
      pileup = 0.006;
#endif
      printf("$ Monte Carlo analysis\n");
      printf(" # vzShift: %f\n", vzShift);
   } else {
      pileup = 0;
      reweight_vertex = 0;
      reweight_mult = 0;
      add_bkg_l1 = 0;
      add_bkg_l2 = 0;
      add_bkg_l3 = 0;
      smear_vertex = 0;
      printf("$ data analysis\n");
   }

   TH3D* hLayer1Hit = new TH3D("hLayer1Hit", "", 75, 0, 15, 60, -3, 3, 64, -3.2, 3.2);
   TH3D* hLayer2Hit = new TH3D("hLayer2Hit", "", 75, 0, 15, 60, -3, 3, 64, -3.2, 3.2);
   TH3D* hLayer3Hit = new TH3D("hLayer3Hit", "", 75, 0, 15, 60, -3, 3, 64, -3.2, 3.2);

   // Prepare hit spectra for random hit
   if (add_bkg_l1) {
      printf(" # projecting background to layer 1\n");
      t->Project("hLayer1Hit", "phi1:eta1:r1");
   }
   if (add_bkg_l2) {
      printf(" # projecting background to layer 2\n");
      t->Project("hLayer2Hit", "phi2:eta2:r2");
   }
   if (add_bkg_l3) {
      printf(" # projecting background to layer 3\n");
      t->Project("hLayer3Hit", "phi3:eta3:r3");
   }

   Parameters par;
   getPixelTreeBranch(t, par);

   int HLT_MB_path = 0;
   hltTree->SetBranchStatus("*", 0);
   hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1", 1);
   hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1", &HLT_MB_path);

   TH1F* hmult_weights = 0;
#ifdef _EPOS_8TEV
   hmult_weights = get_8tev_epos_mult_weights();
#endif
#ifdef _HIJING_8TEV
   hmult_weights = get_8tev_hijing_mult_weights();
#endif
#ifdef _EPOS_5TEV
   hmult_weights = get_5tev_epos_mult_weights();
#endif

   printf(" # number of events: %lli\n", t->GetEntries());

   int n_pu_entry = 1;
   if (pileup!=0)
      printf("$ pileup mixing with mean: %f\n", pileup);
   if (add_bkg_l1)
      printf("$ adding noise\n");
   if (reweight_vertex)
      printf("$ reweighting vertex\n");
   if (reweight_mult)
      printf("$ reweighting mult\n");

   if (use_random_vertex)
      printf("$ using random vertex\n");
   else if (use_kk_vertex)
      printf("$ using reco vertex\n");
   else
      printf("$ using tracklet vertex\n");

#ifdef _POKE_HOLES
   if (isMC) {
      printf("$ poking holes into the pixels\n");
   } else {
      printf("  ! poking holes in data\n");
      printf("  ! comment out '#define _POKE_HOLES' in recohit.h and try again\n");
      return 1;
   }
#endif
   printf("................................................................\n");

   TF1* csfitf = new TF1("csfit", csfit, -4, 4, 1);
   csfitf->SetParameters(1,8882, 0);
   csfitf->SetParLimits(0, 1.8882, 1.8882);

   uint64_t nentries = t->GetEntries();
   // Main loop ===============================================================
   for (uint64_t i=start_entry; i<nentries && i<end_entry; i+=n_pu_entry) {
      t->GetEntry(i);
      if (i % 1000 == 0) {
         printf("   run: %i, event: %lu | %lli %lli %lli\n", par.nRun, i,
                trackletTree12->GetEntries(),
                trackletTree13->GetEntries(),
                trackletTree23->GetEntries());
      }

      if (par.nRun == 285832 && par.nLumi <= 164)
         continue;

      hltTree->GetEntry(i);
      if (!HLT_MB_path && !isMC)
         continue;

      float event_weight = 1.;
      // Reweight MC vertex distribution to match data
      if (reweight_vertex) {
         double myVz = par.vz[1];
         if (myVz < -90) {
            // TF1* f = new TF1("f", "1", -20, 20);
            // TF1* f = new TF1("f", "gaus", -30, 30);
            // f->SetParameters(1, -1.56743, 6.43514);
            // myVz = f->GetRandom();
            // delete f;
            myVz = par.vz[0];
         }

         double data_pdf = 0;
#ifdef _EPOS_5TEV
         // 5 TeV pPb Run 285090
         data_pdf = TMath::Gaus(myVz, 1.09219, 6.27013 - vzShift, 1);
#endif
         // 8 TeV pPb Run 285517
         // data_pdf = TMath::Gaus(myVz, -0.3164 - vzShift, 4.7283, 1);
#if defined(_EPOS_8TEV) || defined(_HIJING_8TEV)
         // 8 TeV pPb Run 285832
         data_pdf = TMath::Gaus(myVz, 0.97423 - vzShift, 4.60961, 1);
#endif

         double mc_pdf = 1;
#ifdef _EPOS_5TEV
         // PixelTree-EPOS-5TeV-HLT.root
         mc_pdf = TMath::Gaus(myVz, -1.79326, 6.50467, 1);
#endif
#ifdef _EPOS_8TEV
         // PixelTree-EPOS-8TeV-HLT.root
         mc_pdf = TMath::Gaus(myVz, -1.56743, 6.43514, 1);
#endif
#ifdef _HIJING_8TEV
         // PixelTree-HIJING-8TeV-HLT.root
         mc_pdf = TMath::Gaus(myVz, -1.61986, 6.46244, 1);
#endif

         event_weight = event_weight * data_pdf / mc_pdf;
      }

      // Fill reco vertex information
      tdata12.nv = par.nv+1;
      tdata13.nv = par.nv+1;
      tdata23.nv = par.nv+1;
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

      // Add background hits
      int bckHits1 = par.nhits1 * add_bkg_l1;
      int bckHits2 = par.nhits2 * add_bkg_l2;
      int bckHits3 = par.nhits3 * add_bkg_l3;

      if (add_bkg_l1!=0) {
         for (int i=par.nhits1; i<par.nhits1+bckHits1; i++) {
            double eta, phi, r;
            hLayer1Hit->GetRandom3(r, eta, phi);
            par.eta1[i] = eta;
            par.phi1[i] = phi;
            par.r1[i] = r;
         }
         par.nhits1 += bckHits1;
      }

      if (add_bkg_l2!=0) {
         for (int i=par.nhits2; i<par.nhits2+bckHits2; i++) {
            double eta, phi, r;
            hLayer2Hit->GetRandom3(r, eta, phi);
            par.eta2[i] = eta;
            par.phi2[i] = phi;
            par.r2[i] = r;
         }
         par.nhits2 += bckHits2;
      }

      if (add_bkg_l3!=0) {
         for (int i=par.nhits3; i<par.nhits3+bckHits3; i++) {
            double eta, phi, r;
            hLayer3Hit->GetRandom3(r, eta, phi);
            par.eta3[i] = eta;
            par.phi3[i] = phi;
            par.r3[i] = r;
         }
         par.nhits3 += bckHits3;
      }

      Parameters pu_par[12];
      pu_par[0] = par;
      double vzPileUp[12];
      vzPileUp[0] = par.vz[0];
      int recoPU = 0;

      std::vector<RecoHit> layer1raw, layer2raw, layer3raw;
      prepareHits(layer1raw, par, 1, 0, 0, 0, split_prob, drop_prob, smear_pixels, cut_cluster_size);
      prepareHits(layer2raw, par, 2, 0, 0, 0, split_prob, drop_prob, smear_pixels, cut_cluster_size);
      // prepareHits(layer3raw, par, 3, 0, 0, 0, split_prob, drop_prob, smear_pixels, cut_cluster_size);

      std::vector<RecoHit> allhits, fallhits;
      prepareHits(allhits, par, 1, 0, 0, 0, split_prob, drop_prob, smear_pixels, cut_cluster_size);
      prepareHits(allhits, par, 2, 0, 0, 0, split_prob, drop_prob, smear_pixels, cut_cluster_size);
      // prepareHits(allhits, par, 3, 0, 0, 0, split_prob, drop_prob, smear_pixels, cut_cluster_size);

      if (pileup!=0) {
         n_pu_entry = 0;
         do {
            n_pu_entry = gRandom->Poisson(pileup);
         } while (n_pu_entry == 0);

         for (int p=1; p<n_pu_entry; p++) {
            if (i + p == nentries) {
               n_pu_entry = p;
               break;
            }
            t->GetEntry(i+p);
            vzPileUp[p] = par.vz[0];

            bckHits1 = par.nhits1 * add_bkg_l1;
            bckHits2 = par.nhits2 * add_bkg_l2;
            bckHits3 = par.nhits3 * add_bkg_l3;

            if (add_bkg_l1!=0) {
               for (int i=par.nhits1; i<par.nhits1+bckHits1; i++) {
                  double eta, phi, r;
                  hLayer1Hit->GetRandom3(r, eta, phi);
                  par.eta1[i] = eta;
                  par.phi1[i] = phi;
                  par.r1[i] = r;
               }
               par.nhits1 += bckHits1;
            }

            if (add_bkg_l2!=0) {
               for (int i=par.nhits2; i<par.nhits2+bckHits2; i++) {
                  double eta, phi, r;
                  hLayer2Hit->GetRandom3(r, eta, phi);
                  par.eta2[i] = eta;
                  par.phi2[i] = phi;
                  par.r2[i] = r;
               }
               par.nhits2 += bckHits2;
            }

            if (add_bkg_l3!=0) {
               for (int i=par.nhits3; i<par.nhits3+bckHits3; i++) {
                  double eta, phi, r;
                  hLayer3Hit->GetRandom3(r, eta, phi);
                  par.eta3[i] = eta;
                  par.phi3[i] = phi;
                  par.r3[i] = r;
               }
               par.nhits3 += bckHits3;
            }
            pu_par[p] = par;

            prepareHits(layer1raw, par, 1, 0, 0, 0, split_prob, drop_prob, smear_pixels, cut_cluster_size);
            prepareHits(layer2raw, par, 2, 0, 0, 0, split_prob, drop_prob, smear_pixels, cut_cluster_size);
            // prepareHits(layer3raw, par, 3, 0, 0, 0, split_prob, drop_prob, smear_pixels, cut_cluster_size);

            prepareHits(allhits, par, 1, 0, 0, 0, split_prob, drop_prob, smear_pixels, cut_cluster_size);
            prepareHits(allhits, par, 2, 0, 0, 0, split_prob, drop_prob, smear_pixels, cut_cluster_size);
         }
         t->GetEntry(i);
         par = pu_par[0];
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

      // double smear = 0;
      // if (smear_vertex != 0) {
      //    while (smear == 0) {
      //       double x = gRandom->Rndm()*2-1;
      //       if (gRandom->Rndm()<TMath::Gaus(x, 0, smear_vertex, 1))
      //          smear = x;
      //    }
      //    trackletVertex += smear;
      // }

      if (use_random_vertex) {
         tdata12.vz[1] = gRandom->Rndm() * 28 - 13 - vzShift;
         tdata13.vz[1] = tdata12.vz[1];
         tdata23.vz[1] = tdata12.vz[1];
      } else if (use_kk_vertex) {
         tdata12.vx[1] = par.vx[1];
         tdata12.vy[1] = par.vy[1];
         tdata12.vz[1] = par.vz[1];
         tdata13.vx[1] = par.vx[1];
         tdata13.vy[1] = par.vy[1];
         tdata13.vz[1] = par.vz[1];
         tdata23.vx[1] = par.vx[1];
         tdata23.vy[1] = par.vy[1];
         tdata23.vz[1] = par.vz[1];
      } else {
         tdata12.vz[1] = trackletVertex;
         tdata13.vz[1] = trackletVertex;
         tdata23.vz[1] = trackletVertex;
      }

      // Process hits with Vz constraint:
      std::vector<RecoHit> layer1Cut;
      prepareHits(layer1Cut, par, 1, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], split_prob, drop_prob, 0, 1);

      std::vector<RecoHit> layer1, layer2, layer3, layer4, layer5;
      prepareHits(layer1, par, 1, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], split_prob, drop_prob, smear_pixels, cut_cluster_size);
      prepareHits(layer2, par, 2, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], split_prob, drop_prob, smear_pixels, cut_cluster_size);
      prepareHits(layer3, par, 3, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], split_prob, drop_prob, smear_pixels, cut_cluster_size);

      if (n_pu_entry>0) {
         for (int p=1; p<n_pu_entry; p++) {
            t->GetEntry(i+p);
            par = pu_par[p];
            prepareHits(layer1, par, 1, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], split_prob, drop_prob, smear_pixels, cut_cluster_size);
            prepareHits(layer2, par, 2, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], split_prob, drop_prob, smear_pixels, cut_cluster_size);
            prepareHits(layer3, par, 3, tdata12.vx[1], tdata12.vy[1], tdata12.vz[1], split_prob, drop_prob, smear_pixels, cut_cluster_size);
         }
         t->GetEntry(i);
         par = pu_par[0];
      }

      std::vector<Tracklet> recoTracklets12, recoTracklets13, recoTracklets23;

      std::vector<RecoHit> combinedhits;
      if (do_two_hit_corr) {
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
      }

      // Cluster Vertex Compatibility =========================================
      double clusVtxQual1 = skip_cvqual ? 0 : getClusVtxCompat(layer1, 1);;
      double clusVtxQual2 = skip_cvqual ? 0 : getClusVtxCompat(layer2, 2);;
      double clusVtxQual3 = skip_cvqual ? 0 : getClusVtxCompat(layer3, 3);;

      // Move the Vertex back
      // if (smear_vertex!=0) {
      //    tdata12.vz[1] = trackletVertex - smear;
      //    tdata13.vz[1] = trackletVertex - smear;
      //    tdata23.vz[1] = trackletVertex - smear;
      // }

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
   tdata##q##w.passHLT    = HLT_MB_path;                       \
   tdata##q##w.clusVtxQual1 = clusVtxQual1;                    \
   tdata##q##w.clusVtxQual2 = clusVtxQual2;                    \
   tdata##q##w.clusVtxQual3 = clusVtxQual3;                    \
                                                               \
   int ntracklet##q##w##s = 0;                                 \
   int ntracklet##q##w##b = 0;                                 \
   for (int j=0; j<(int)tdata##q##w.nTracklet; j++) {          \
      tdata##q##w.eta1[j] = recoTracklets##q##w[j].eta1;       \
      tdata##q##w.eta2[j] = recoTracklets##q##w[j].eta2;       \
      tdata##q##w.r1[j]   = recoTracklets##q##w[j].r1;         \
      tdata##q##w.r2[j]   = recoTracklets##q##w[j].r2;         \
      tdata##q##w.cs1[j]  = recoTracklets##q##w[j].cs1;        \
      tdata##q##w.cs2[j]  = recoTracklets##q##w[j].cs2;        \
      tdata##q##w.ch1[j]  = recoTracklets##q##w[j].ch1;        \
      tdata##q##w.ch2[j]  = recoTracklets##q##w[j].ch2;        \
      tdata##q##w.phi1[j] = recoTracklets##q##w[j].phi1;       \
      tdata##q##w.phi2[j] = recoTracklets##q##w[j].phi2;       \
      tdata##q##w.deta[j] = recoTracklets##q##w[j].deta;       \
      tdata##q##w.dphi[j] = recoTracklets##q##w[j].dphi;       \
      if (fabs(tdata##q##w.deta[j])<0.1) {                     \
         if (fabs(tdata##q##w.dphi[j])<1.0)                    \
            ntracklet##q##w##s++;                              \
         if (fabs(tdata##q##w.dphi[j])>1.0 &&                  \
             fabs(tdata##q##w.dphi[j])<2.0)                    \
            ntracklet##q##w##b++;                              \
      }                                                        \
   }                                                           \
   tdata##q##w.mult = ntracklet##q##w##s - ntracklet##q##w##b; \
   tdata##q##w.nhit1_cut = layer1Cut.size();                   \
                                                               \
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
   for (int j=0; j<par.nv; j++)                                \
      tdata##q##w.vz[j] += vzShift;                            \
}
      fillTrackletTree(1, 2);
      fillTrackletTree(1, 3);
      fillTrackletTree(2, 3);

      if (reweight_mult && hmult_weights)
         event_weight = event_weight * hmult_weights->GetBinContent(hmult_weights->FindBin(tdata23.mult));

      tdata12.weight = event_weight;
      tdata13.weight = event_weight;
      tdata23.weight = event_weight;

      trackletTree12->Fill();
      trackletTree13->Fill();
      trackletTree23->Fill();
   }

   inf->Close();

   trackletTree12->Write("", TObject::kOverwrite);
   trackletTree13->Write("", TObject::kOverwrite);
   trackletTree23->Write("", TObject::kOverwrite);

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

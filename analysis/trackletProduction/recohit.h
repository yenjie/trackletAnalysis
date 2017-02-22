#ifndef _RECOHIT_H
#define _RECOHIT_H

#define _POKE_HOLES

#define _MAX_ENTRY 4000

#include <vector>

#ifdef _POKE_HOLES
#include "pdfs.h"
#endif

#include "TTree.h"
#include "TRandom.h"
#include "TMath.h"
#include "Math/Vector3D.h"

class RecoHit {
   public:
      RecoHit(double eta, double phi, double r, double cs, double ch) :
         eta(eta), phi(phi), r(r), cs(cs), ch(ch)
      {}
      RecoHit(double eta, double phi, double r, double cs, double ch, int layer) :
         eta(eta), phi(phi), r(r), cs(cs), ch(ch), layer(layer)
      {}

      ~RecoHit() {};

      double eta;
      double phi;
      double r;
      double cs;   // cluster size
      double ch;   // cluster charge
      int layer;
};

class Parameters {
   public:
      int nRun, nEv, nLumi, nBX, nHFn, nHFp, nHits;
      float beamSpotX, beamSpotY, beamSpotZ;
      float vx[8], vy[8], vz[8];
      float eta1[_MAX_ENTRY], phi1[_MAX_ENTRY], r1[_MAX_ENTRY], cs1[_MAX_ENTRY], ch1[_MAX_ENTRY];
      float eta2[_MAX_ENTRY], phi2[_MAX_ENTRY], r2[_MAX_ENTRY], cs2[_MAX_ENTRY], ch2[_MAX_ENTRY];
      float eta3[_MAX_ENTRY], phi3[_MAX_ENTRY], r3[_MAX_ENTRY], cs3[_MAX_ENTRY], ch3[_MAX_ENTRY];
      float eta[_MAX_ENTRY], phi[_MAX_ENTRY], pt[_MAX_ENTRY];
      int nhits1, nhits2, nhits3, mult, nv, npart, evtType, chg[_MAX_ENTRY], pdg[_MAX_ENTRY];
};

bool comp_phi(RecoHit a, RecoHit b) {
   return a.phi < b.phi;
}

void prepareHits(std::vector<RecoHit>& cleaned_hits, Parameters par, Int_t layer,
                 double vx, double vy, double vz,
                 double split_prob = 0, double drop_prob = 0, bool smear_pixels = 0,
                 bool cut_on_cluster_size = 0) {
   double x0, y0; // beamspot for each run
   if (par.nRun == 285090) {
      x0 = 0.06149;
      y0 = 0.1056;
   } else if (par.nRun == 285517) {
      x0 = 0.05746;
      y0 = 0.1051;
   } else if (par.nRun == 285832) {
      x0 = 0.05774;
      y0 = 0.10872;
   } else {
      x0 = 0.1048;
      y0 = 0.1687;
   }

#ifdef _POKE_HOLES
   TH2D* hholes = get_holes();
#endif

   std::vector<RecoHit> hits;
   if (layer == 1) {
      for (int ihit = 0; ihit < par.nhits1; ++ihit) {
         if (par.nRun == 285517) {
            if ((fabs(par.phi1[ihit] - 1.78408) < 0.00001 && fabs(par.eta1[ihit] + 1.16674) < 0.00001))
               continue;
         } else if (par.nRun == 285832) {
            if ((fabs(par.phi1[ihit] - 1.78421) < 0.00001 && fabs(par.eta1[ihit] + 1.16703) < 0.00001) ||
                (fabs(par.phi1[ihit] - 1.43093) < 0.00001 && fabs(par.eta1[ihit] + 0.66483) < 0.00001) ||
                (fabs(par.phi1[ihit] - 1.42661) < 0.00001 && fabs(par.eta1[ihit] + 0.66445) < 0.00001) ||
                (fabs(par.phi1[ihit] - 1.40936) < 0.00001 && fabs(par.eta1[ihit] + 0.66284) < 0.00001) ||
                (fabs(par.phi1[ihit] - 1.40721) < 0.00001 && fabs(par.eta1[ihit] + 0.66262) < 0.00001) ||
                (fabs(par.phi1[ihit] - 2.62996) < 0.00001 && fabs(par.eta1[ihit] + 2.10101) < 0.00001) ||
                (fabs(par.phi1[ihit] - 1.35537) < 0.00001 && fabs(par.eta1[ihit] + 2.48903) < 0.00001) ||
                (fabs(par.phi1[ihit] - 1.35055) < 0.00001 && fabs(par.eta1[ihit] + 2.48966) < 0.00001) ||
                (fabs(par.phi1[ihit] - 1.40721) < 0.00001 && fabs(par.eta1[ihit] + 0.66262) < 0.00001) ||
                (fabs(par.phi1[ihit] - 1.40936) < 0.00001 && fabs(par.eta1[ihit] + 0.66284) < 0.00001) ||
                (fabs(par.phi1[ihit] - 1.41582) < 0.00001 && fabs(par.eta1[ihit] + 0.66346) < 0.00001) ||
                (fabs(par.phi1[ihit] - 1.42661) < 0.00001 && fabs(par.eta1[ihit] + 0.66445) < 0.00001) ||
                (fabs(par.phi1[ihit] - 1.43093) < 0.00001 && fabs(par.eta1[ihit] + 0.66483) < 0.00001) ||
                (fabs(par.phi1[ihit] - 3.09539) < 0.00001 && fabs(par.eta1[ihit] - 2.08781) < 0.00001) ||
                (fabs(par.phi1[ihit] - 1.15823) < 0.00001 && fabs(par.eta1[ihit] + 2.09717) < 0.00001) ||
                (fabs(par.phi1[ihit] - 1.16068) < 0.00001 && fabs(par.eta1[ihit] + 2.09732) < 0.00001))
               continue;
         }
#ifdef _POKE_HOLES
         // drop hits in MC to match hit distributions in data
         if (par.nRun < 10) {
            int hit_bin = hholes->FindBin(par.phi1[ihit], par.eta1[ihit]);
            float hit_ratio = hholes->GetBinContent(hit_bin);
            if (hit_ratio != 0 && gRandom->Rndm() > hit_ratio)
               continue;
         }
#endif
         if (gRandom->Rndm() < drop_prob)
            continue;
         RecoHit tmp(par.eta1[ihit], par.phi1[ihit], par.r1[ihit], par.cs1[ihit], par.ch1[ihit], 1);
         hits.push_back(tmp);
         // put artifical split hits
         if (gRandom->Rndm() < split_prob)
            hits.push_back(tmp);
      }
   } else if (layer == 2) {
      for (int ihit = 0; ihit < par.nhits2; ++ihit) {
         // Run 285517
         if (par.nRun == 285517) {
            if ((fabs(par.phi2[ihit] + 2.02799) < 0.00001 && fabs(par.eta2[ihit] - 0.34798) < 0.00001) ||
                (fabs(par.phi2[ihit] + 2.03061) < 0.00001 && fabs(par.eta2[ihit] - 0.34795) < 0.00001) ||
                (fabs(par.phi2[ihit] + 2.29179) < 0.00001 && fabs(par.eta2[ihit] - 0.70377) < 0.00001) ||
                (fabs(par.phi2[ihit] - 1.50481) < 0.00001 && fabs(par.eta2[ihit] - 1.41267) < 0.00001) ||
                (fabs(par.phi2[ihit] + 1.91605) < 0.00001 && fabs(par.eta2[ihit] - 1.12815) < 0.00001) ||
                (fabs(par.phi2[ihit] - 0.11765) < 0.00001 && fabs(par.eta2[ihit] - 0.43056) < 0.00001))
               continue;
         } else if (par.nRun == 285832) {
            if ((fabs(par.phi2[ihit] + 2.02809) < 0.00001 && fabs(par.eta2[ihit] - 0.34774) < 0.00001) ||
                (fabs(par.phi2[ihit] + 2.03070) < 0.00001 && fabs(par.eta2[ihit] - 0.34771) < 0.00001) ||
                (fabs(par.phi2[ihit] + 2.29325) < 0.00001 && fabs(par.eta2[ihit] - 0.70344) < 0.00001) ||
                (fabs(par.phi2[ihit] + 2.29048) < 0.00001 && fabs(par.eta2[ihit] - 0.70358) < 0.00001) ||
                (fabs(par.phi2[ihit] - 1.50475) < 0.00001 && fabs(par.eta2[ihit] - 1.41259) < 0.00001) ||
                (fabs(par.phi2[ihit] + 1.91616) < 0.00001 && fabs(par.eta2[ihit] - 1.12794) < 0.00001) ||
                (fabs(par.phi2[ihit] + 0.96467) < 0.00001 && fabs(par.eta2[ihit] - 0.65614) < 0.00001) ||
                (fabs(par.phi2[ihit] + 0.96730) < 0.00001 && fabs(par.eta2[ihit] - 0.65611) < 0.00001) ||
                (fabs(par.phi2[ihit] + 0.34368) < 0.00001 && fabs(par.eta2[ihit] - 1.84190) < 0.00001) ||
                (fabs(par.phi2[ihit] + 2.50245) < 0.00001 && fabs(par.eta2[ihit] - 0.28293) < 0.00001) ||
                (fabs(par.phi2[ihit] + 2.50521) < 0.00001 && fabs(par.eta2[ihit] - 0.28301) < 0.00001) ||
                (fabs(par.phi2[ihit] + 0.88593) < 0.00001 && fabs(par.eta2[ihit] - 0.68248) < 0.00001) ||
                (fabs(par.phi2[ihit] + 0.96467) < 0.00001 && fabs(par.eta2[ihit] - 0.65614) < 0.00001) ||
                (fabs(par.phi2[ihit] + 0.96730) < 0.00001 && fabs(par.eta2[ihit] - 0.65611) < 0.00001) ||
                (fabs(par.phi2[ihit] + 2.12876) < 0.00001 && fabs(par.eta2[ihit] + 0.34704) < 0.00001) ||
                (fabs(par.phi2[ihit] + 2.13015) < 0.00001 && fabs(par.eta2[ihit] + 0.34905) < 0.00001) ||
                (fabs(par.phi2[ihit] + 2.13154) < 0.00001 && fabs(par.eta2[ihit] + 0.34711) < 0.00001) ||
                (fabs(par.phi2[ihit] + 0.96534) < 0.00001 && fabs(par.eta2[ihit] + 0.24780) < 0.00001))
               continue;
         }
         if (gRandom->Rndm() < drop_prob)
            continue;
         RecoHit tmp(par.eta2[ihit], par.phi2[ihit], par.r2[ihit], par.cs2[ihit], par.ch2[ihit], 2);
         hits.push_back(tmp);
         // put artifical split hits
         if (gRandom->Rndm() < split_prob)
            hits.push_back(tmp);
      }
   } else if (layer == 3) {
      for (int ihit = 0; ihit < par.nhits3; ++ihit) {
         if (par.nRun == 285517) {
            if ((fabs(par.phi3[ihit] - 0.01456) < 0.00001 && fabs(par.eta3[ihit] + 0.06655) < 0.00001))
               continue;
         } else if (par.nRun == 285832) {
            if ((fabs(par.phi3[ihit] - 0.01456) < 0.00001 && fabs(par.eta3[ihit] + 0.06675) < 0.00001) ||
                (fabs(par.phi3[ihit] + 0.72997) < 0.00001 && fabs(par.eta3[ihit] + 0.62795) < 0.00001) ||
                (fabs(par.phi3[ihit] + 1.97283) < 0.00001 && fabs(par.eta3[ihit] - 0.73149) < 0.00001))
               continue;
         }
         if (gRandom->Rndm() < drop_prob)
            continue;
         RecoHit tmp(par.eta3[ihit], par.phi3[ihit], par.r3[ihit], par.cs3[ihit], par.ch3[ihit], 3);
         hits.push_back(tmp);
         // put artifical split hits
         if (gRandom->Rndm() < split_prob)
            hits.push_back(tmp);
      }
   }

#ifdef _POKE_HOLES
   hholes->Delete();
#endif

   sort(hits.begin(), hits.end(), comp_phi);

   for (std::size_t ihit = 0; ihit < hits.size(); ++ihit) {
      double x = hits[ihit].r * cos(hits[ihit].phi);
      double y = hits[ihit].r * sin(hits[ihit].phi);
      double z = hits[ihit].r / tan(atan(exp(-hits[ihit].eta)) * 2);

      if (smear_pixels) {
         x += gRandom->Gaus(0, 0.005);
         y += gRandom->Gaus(0, 0.005);
         z += gRandom->Gaus(0, 0.005);
      }

      ROOT::Math::XYZVector rel_vector(x - x0, y - y0, z - vz);
      RecoHit rel_hit(rel_vector.eta(), rel_vector.phi(), rel_vector.rho(), hits[ihit].cs, hits[ihit].ch, hits[ihit].layer);

      double eta = rel_vector.eta();
      if (cut_on_cluster_size) {
         if (fabs(eta) <= 0.5 &&                    hits[ihit].cs < 1) continue;
         if (fabs(eta) <= 1.0 && fabs(eta) > 0.5 && hits[ihit].cs < 2) continue;
         if (fabs(eta) <= 1.5 && fabs(eta) > 1.0 && hits[ihit].cs < 3) continue;
         if (fabs(eta) <= 2.0 && fabs(eta) > 1.5 && hits[ihit].cs < 4) continue;
         if (fabs(eta) <= 2.5 && fabs(eta) > 2.0 && hits[ihit].cs < 6) continue;
         if (fabs(eta) <= 5.0 && fabs(eta) > 2.5 && hits[ihit].cs < 9) continue;
      }
      cleaned_hits.push_back(rel_hit);
   }
}

void getPixelTreeBranch(TTree* t, Parameters& par) {
   t->SetBranchAddress("nRun", &par.nRun);
   t->SetBranchAddress("nEv", &par.nEv);
   t->SetBranchAddress("nLumi", &par.nLumi);
   t->SetBranchAddress("nBX", &par.nBX);
   t->SetBranchAddress("nHFn", &par.nHFn);
   t->SetBranchAddress("nHFp", &par.nHFp);

   t->SetBranchAddress("eta1", par.eta1);
   t->SetBranchAddress("phi1", par.phi1);
   t->SetBranchAddress("r1", par.r1);
   t->SetBranchAddress("cs1", par.cs1);
   t->SetBranchAddress("ch1", par.ch1);
   t->SetBranchAddress("eta2", par.eta2);
   t->SetBranchAddress("phi2", par.phi2);
   t->SetBranchAddress("r2", par.r2);
   t->SetBranchAddress("cs2", par.cs2);
   t->SetBranchAddress("ch2", par.ch2);
   t->SetBranchAddress("eta3", par.eta3);
   t->SetBranchAddress("phi3", par.phi3);
   t->SetBranchAddress("r3", par.r3);
   t->SetBranchAddress("cs3", par.cs3);
   t->SetBranchAddress("ch3", par.ch3);
   t->SetBranchAddress("nhits1", &par.nhits1);
   t->SetBranchAddress("nhits2", &par.nhits2);
   t->SetBranchAddress("nhits3", &par.nhits3);

   t->SetBranchAddress("vx", par.vx);
   t->SetBranchAddress("vy", par.vy);
   t->SetBranchAddress("vz", par.vz);
   t->SetBranchAddress("beamSpotX", &par.beamSpotX);
   t->SetBranchAddress("beamSpotY", &par.beamSpotY);
   t->SetBranchAddress("beamSpotZ", &par.beamSpotZ);
   t->SetBranchAddress("nv", &par.nv);
   t->SetBranchAddress("npart", &par.npart);
   t->SetBranchAddress("eta", &par.eta);
   t->SetBranchAddress("phi", &par.phi);
   t->SetBranchAddress("pt", &par.pt);
   t->SetBranchAddress("chg", &par.chg);
   t->SetBranchAddress("pdg", &par.pdg);
   t->SetBranchAddress("evtType", &par.evtType);
}

#endif

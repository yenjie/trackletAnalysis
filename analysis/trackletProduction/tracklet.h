#ifndef _TRACKLET_H
#define _TRACKLET_H

#define PI 3.141592653589

#include "recohit.h"

double calc_dphi(double phi1, double phi2) {
   double dphi = phi1 - phi2;
   if (dphi > PI)
      dphi = 2 * PI - dphi;
   else if (dphi < -PI)
      dphi = -2 * PI - dphi;

   return dphi;
}

class Tracklet {
   public:
      Tracklet(double eta1, double eta2, double phi1, double phi2, double r1, double r2, double cs1, double cs2, double ch1, double ch2) :
         eta1(eta1), eta2(eta2), phi1(phi1), phi2(phi2), r1(r1), r2(r2), cs1(cs1), cs2(cs2), ch1(ch1), ch2(ch2),
         deta(eta1 - eta2), dphi(calc_dphi(phi1, phi2))
      {}

      ~Tracklet() {};

      double eta1;
      double eta2;
      double phi1;
      double phi2;
      double r1;
      double r2;
      double cs1;
      double cs2;
      double ch1;
      double ch2;

      double deta;
      double dphi;
};

typedef struct phit {
   int index;
   bool paired;
   int lindex;
   int pindex;
   int cindex;
} phit;

std::vector<Tracklet> recoTracklets(std::vector<RecoHit> allhits, int l1, int l2) {
   std::vector<phit> hits[5];
   std::vector<phit> cands;
   std::vector<bool> valid;
   std::vector<Tracklet> tracklets;

   hits[--l1].reserve(allhits.size());
   hits[--l2].reserve(allhits.size());

   int tl[5];
   for (std::size_t a = 0; a < allhits.size(); a++) {
      int l = allhits[a].layer - 1;
      int sl = l1 + l2 - l;

      int tl1 = l1;
      int tl2 = l2;
      int tsl = sl;
      if (l1 == l2) {
         tl1 = l;
         tl2 = 3 - l;
         tsl = tl2;
      }

      phit nexthit;
      nexthit.index = a;
      nexthit.paired = 0;
      nexthit.lindex = hits[sl].size();
      hits[l].push_back(nexthit);

      std::vector<phit>::iterator thisl = hits[l].end() - 1;
      std::vector<phit>::iterator thatl = hits[sl].end() - 1;
      for (; !hits[sl].empty() && !(*thisl).paired && thatl + 1 != hits[sl].begin(); --thatl) {
         if ((*thatl).index == (*thisl).index) continue;
         if ((*thatl).paired) {
            valid[(*thatl).cindex] = 0;
            if (abs(allhits[(*thatl).index].eta - allhits[(*thatl).pindex].eta) < abs(allhits[(*thatl).index].eta - allhits[(*thisl).index].eta)) {
               tl[l] = (*thatl).pindex;
               tl[tsl] = (*thatl).index;
               Tracklet tracklet(allhits[tl[tl1]].eta, allhits[tl[tl2]].eta, allhits[tl[tl1]].phi, allhits[tl[tl2]].phi, allhits[tl[tl1]].r, allhits[tl[tl2]].r, allhits[tl[tl1]].cs, allhits[tl[tl2]].cs, allhits[tl[tl1]].ch, allhits[tl[tl2]].ch);
               tracklets.push_back(tracklet);

               hits[sl].erase(thatl);

               std::vector<phit>::iterator it = hits[l].begin();
               for (; it != hits[l].end() && (*it).index < tl[l]; ++it);
               for (; it != hits[l].end(); ++it)
                  (*it).lindex --;
            } else {
               phit unpaired = cands[(*thatl).cindex];
               unpaired.paired = 0;

               (*thisl).paired = 1;
               (*thisl).pindex = (*thatl).index;
               (*thisl).cindex = cands.size();
               (*thatl).paired = 1;
               (*thatl).pindex = (*thisl).index;
               (*thatl).cindex = cands.size();
               cands.push_back(*thatl);
               hits[sl].erase(thatl);
               valid.push_back(1);

               for (thisl = hits[l].begin(); thisl != hits[l].end() && (*thisl).index < unpaired.index; ++thisl);
               hits[l].insert(thisl, unpaired);
               thatl = hits[sl].begin() + unpaired.lindex;
            }
         } else {
            (*thisl).paired = 1;
            (*thisl).pindex = (*thatl).index;
            (*thisl).cindex = cands.size();
            (*thatl).paired = 1;
            (*thatl).pindex = (*thisl).index;
            (*thatl).cindex = cands.size();
            cands.push_back(*thatl);
            hits[sl].erase(thatl);
            valid.push_back(1);
         }
      }
   }

   for (std::size_t t = 0; t < cands.size(); t++) {
      if (valid[t]) {
         int l = allhits[cands[t].index].layer - 1;
         int sl = l1 + l2 - l;

         int tl1 = l1;
         int tl2 = l2;
         int tsl = sl;
         if (l1 == l2) {
            tl1 = l;
            tl2 = 3 - l;
            tsl = tl2;
         }

         tl[l] = cands[t].index;
         tl[tsl] = cands[t].pindex;
         Tracklet tracklet(allhits[tl[tl1]].eta, allhits[tl[tl2]].eta, allhits[tl[tl1]].phi, allhits[tl[tl2]].phi, allhits[tl[tl1]].r, allhits[tl[tl2]].r, allhits[tl[tl1]].cs, allhits[tl[tl2]].cs, allhits[tl[tl1]].ch, allhits[tl[tl2]].ch);
         tracklets.push_back(tracklet);
      }
   }

   return tracklets;
}

class TrackletData {
   public:
      int nRun, nEv, nLumi, nBX, nHFn, nHFp, nHits;
      int passHLT;
      float eta1[_MAX_ENTRY], phi1[_MAX_ENTRY], r1[_MAX_ENTRY], cs1[_MAX_ENTRY], ch1[_MAX_ENTRY];
      float eta2[_MAX_ENTRY], phi2[_MAX_ENTRY], r2[_MAX_ENTRY], cs2[_MAX_ENTRY], ch2[_MAX_ENTRY];
      float deta[_MAX_ENTRY], dphi[_MAX_ENTRY];
      float vx[8], vy[8], vz[8];
      float eta[_MAX_ENTRY], phi[_MAX_ENTRY], pt[_MAX_ENTRY], nhad[12];
      int chg[_MAX_ENTRY], pdg[_MAX_ENTRY];
      int nTracklet, nhit1, nhit2, mult, nhit1_cut, nv, npart, evtType, trackletType;
      float weight;
      float clusVtxQual1, clusVtxQual2, clusVtxQual3;
};

void setTrackletTreeBranch(TTree* trackletTree, TrackletData& tdata) {
   trackletTree->Branch("nRun", &tdata.nRun, "nRun/I");
   trackletTree->Branch("nEv", &tdata.nEv, "nEv/I");
   trackletTree->Branch("nLumi", &tdata.nLumi, "nLumi/I");
   trackletTree->Branch("nBX", &tdata.nBX, "nBX/I");
   trackletTree->Branch("nHFn", &tdata.nHFn, "nHFn/I");
   trackletTree->Branch("nHFp", &tdata.nHFp, "nHFp/I");
   trackletTree->Branch("nHits", &tdata.nHits, "nHits/I");

   trackletTree->Branch("passHLT", &tdata.passHLT, "passHLT/I");

   trackletTree->Branch("nTracklets", &tdata.nTracklet, "nTracklets/I");
   trackletTree->Branch("nhit1", &tdata.nhit1, "nhit1/I");
   trackletTree->Branch("nhit2", &tdata.nhit2, "nhit2/I");
   trackletTree->Branch("mult", &tdata.mult, "mult/I");
   trackletTree->Branch("nhit1_cut", &tdata.nhit1_cut, "nhit1_cut/I");
   trackletTree->Branch("nv", &tdata.nv, "nv/I");
   trackletTree->Branch("vx", tdata.vx, "vx[nv]/F");
   trackletTree->Branch("vy", tdata.vy, "vy[nv]/F");
   trackletTree->Branch("vz", tdata.vz, "vz[nv]/F");
   trackletTree->Branch("eta1", tdata.eta1, "eta1[nTracklets]/F");
   trackletTree->Branch("phi1", tdata.phi1, "phi1[nTracklets]/F");
   trackletTree->Branch("r1", tdata.r1, "r1[nTracklets]/F");
   trackletTree->Branch("cs1", tdata.cs1, "cs1[nTracklets]/F");
   trackletTree->Branch("ch1", tdata.ch1, "ch1[nTracklets]/F");
   trackletTree->Branch("eta2", tdata.eta2, "eta2[nTracklets]/F");
   trackletTree->Branch("phi2", tdata.phi2, "phi2[nTracklets]/F");
   trackletTree->Branch("r2", tdata.r2, "r2[nTracklets]/F");
   trackletTree->Branch("cs2", tdata.cs2, "cs2[nTracklets]/F");
   trackletTree->Branch("ch2", tdata.ch1, "ch2[nTracklets]/F");
   trackletTree->Branch("deta", tdata.deta, "deta[nTracklets]/F");
   trackletTree->Branch("dphi", tdata.dphi, "dphi[nTracklets]/F");

   trackletTree->Branch("weight", &tdata.weight, "weight/F");

   trackletTree->Branch("clusVtxQual1", &tdata.clusVtxQual1, "clusVtxQual1/F");
   trackletTree->Branch("clusVtxQual2", &tdata.clusVtxQual2, "clusVtxQual2/F");
   trackletTree->Branch("clusVtxQual3", &tdata.clusVtxQual3, "clusVtxQual3/F");

   trackletTree->Branch("npart", &tdata.npart, "npart/I");
   trackletTree->Branch("eta", tdata.eta, "eta[npart]/F");
   trackletTree->Branch("phi", tdata.phi, "phi[npart]/F");
   trackletTree->Branch("pdg", tdata.pdg, "pdg[npart]/I");
   trackletTree->Branch("chg", tdata.chg, "chg[npart]/I");
   trackletTree->Branch("nhad", tdata.nhad, "nhad[12]/F");
   trackletTree->Branch("pt", tdata.pt, "pt[npart]/F");
   trackletTree->Branch("evtType", &tdata.evtType, "evtType/I");

   trackletTree->SetAlias("dR", "sqrt(deta*deta+dphi*dphi)");
   trackletTree->SetAlias("dRR", "sqrt(deta*deta+dphi*dphi+(r1-r2)*(r1-r2))") ;
   trackletTree->SetAlias("z1", "r1/tan(atan(exp(-eta1))*2)+vz[3]");
   trackletTree->SetAlias("z2", "r2/tan(atan(exp(-eta2))*2)+vz[3]");
   trackletTree->SetAlias("x1", "r1*cos(phi1)");
   trackletTree->SetAlias("x2", "r2*cos(phi2)");
   trackletTree->SetAlias("y1", "r1*sin(phi1)");
   trackletTree->SetAlias("y2", "r2*sin(phi2)");
}

#endif

#ifndef RECOHIT
#define RECOHIT

#define maxEntry 4000
#define maxEntry2 4000

#include <vector>
#include <algorithm>
#include <Math/Vector3D.h>
#include <TRandom.h>
#include <TTree.h>
#include <TMath.h>

using namespace std;

class RecoHit {
   public:
      RecoHit(double _eta, double _phi, double _r, double _cs, double _ch) {
         eta = _eta;
         phi = _phi;
         r = _r;
         cs = _cs;
         ch = _ch;
      };
      RecoHit(double _eta, double _phi, double _r, double _cs, double _ch, int _l) {
         eta = _eta;
         phi = _phi;
         r = _r;
         cs = _cs;
         ch = _ch;
         layer = _l;
      };

      ~RecoHit() {};

      double eta;
      double phi;
      double r;
      double cs;   // cluster size
      double ch;   // cluster charge
      int layer;
};

class SelectionCriteria {
   public:
      double drCut   ;       // to remove double hit
      double dPhiCut ;       // to remove double hit
      double dEtaCut ;       // to remove double hit
      double vzCut   ;       // vertex cut

      bool verbose_ ;
      bool useDeltaPhi_;
      bool useDeltaRho_;
      bool checkSecondLayer_;
};

class Parameters {
   public:
      int nRun, nEv, nLumi, nHltBit, nL1ABit, nL1TBit, nBX, nHFn, nHFp, nHits;
      bool hltBit[500], l1ABit[500], l1TBit[500];
      int L1_BPTX_AND, L1_BPTX_OR, L1_BPTX_plus, L1_BPTX_minus;
      float beamSpotX, beamSpotY, beamSpotZ;
      float vx[8], vy[8], vz[8];
      float eta1[maxEntry], phi1[maxEntry], r1[maxEntry], cs1[maxEntry], ch1[maxEntry];
      float eta2[maxEntry], phi2[maxEntry], r2[maxEntry], cs2[maxEntry], ch2[maxEntry];
      float eta3[maxEntry], phi3[maxEntry], r3[maxEntry], cs3[maxEntry], ch3[maxEntry];
      float eta4[maxEntry], phi4[maxEntry], r4[maxEntry], cs4[maxEntry], ch4[maxEntry];
      float eta5[maxEntry], phi5[maxEntry], r5[maxEntry], cs5[maxEntry], ch5[maxEntry];
      float eta[maxEntry], phi[maxEntry], pt[maxEntry];
      int nhits1, nhits2, nhits3, nhits4, nhits5, mult, nv, npart, evtType, chg[maxEntry], pdg[maxEntry];
      float xi;
      bool passDS, passSingleTrack;
      int ntrks, ntrksCut;
};

class TrackletData {
   public:
      int nRun, nEv, nLumi, nHltBit, nL1ABit, nL1TBit, nBX, nHFn, nHFp, nHits;
      bool hltBit[500], l1ABit[500], l1TBit[500];
      int L1_BPTX_AND, L1_BPTX_OR, L1_BPTX_plus, L1_BPTX_minus;
      float eta1[maxEntry], phi1[maxEntry], r1[maxEntry], cs1[maxEntry], ch1[maxEntry];
      float eta2[maxEntry], phi2[maxEntry], r2[maxEntry], cs2[maxEntry], ch2[maxEntry];
      float deta[maxEntry], dphi[maxEntry];
      float vx[8], vy[8], vz[8];
      float eta[maxEntry2], phi[maxEntry2], pt[maxEntry2], nhad[12];
      int chg[maxEntry2], pdg[maxEntry2];
      float pro2;
      int nTracklet, nhit1, nhit2, mult, mult2, nv, npart, evtType, trackletType;
      float xi;
      bool passDS, passSingleTrack;
      int ntrks, ntrksCut;
      int nPU, recoPU;
      float vzPU[12];
      float clusVtxQual1, clusVtxQual2, clusVtxQual3;
};

bool compareEta(RecoHit a, RecoHit b) {
   return a.eta < b.eta;
}
bool comparePhi(RecoHit a, RecoHit b) {
   return a.phi < b.phi;
}
bool compareAbsEta(RecoHit a, RecoHit b) {
   return fabs(a.eta) < fabs(b.eta);
}

double calcDphi(double phi1, double phi2);

void prepareHits(vector<RecoHit> &cleanedHits, Parameters par, SelectionCriteria cuts, Int_t layer,
                  double vx, double vy, double vz, double splitProb = 0, double dropProb = 0,
                  bool cutOnClusterSize = 0, double runNum = 0, double nLumi = 0, bool smearPixels = 0)
{
   vector<RecoHit> hits;
   static bool firstCall = 0;

   double x0, y0;
   // The beamspot for each run
   if (par.nRun == 247324) {
      x0 = 0.052;
      y0 = 0.160;
   } else {
      x0 = 0.032;
      y0 = 0;
   }

   if (vz!=0 && firstCall==0) {
      cout << "Beamspot X0 = " << x0 << " Y0 = " << y0 << endl;
      firstCall = 1;
   }

   if (layer == 1) {
      for (int ihit = 0; ihit < par.nhits1; ++ihit) {
         if (par.phi1[ihit]>-0.3904153 && par.phi1[ihit]<-0.3904152 && par.eta1[ihit]>-1.45722 && par.eta1[ihit]<-1.45721)
            continue;
         if (par.phi1[ihit]>1.3385 && par.phi1[ihit]<1.3386 && par.eta1[ihit]>-2.5014 && par.eta1[ihit]<-2.5013)
            continue;
         if (par.phi1[ihit]>2.0982 && par.phi1[ihit]<2.0984 && par.eta1[ihit]>-0.9114 && par.eta1[ihit]<-0.9112)
            continue;
         RecoHit tmp(par.eta1[ihit], par.phi1[ihit], par.r1[ihit], par.cs1[ihit], par.ch1[ihit], 1);
         if (gRandom->Rndm()<dropProb) continue;
         hits.push_back(tmp);
         // put artifical split hits
         if (gRandom->Rndm()<splitProb) hits.push_back(tmp);
      }
   } else if (layer == 2) {
      for (int ihit = 0; ihit < par.nhits2; ++ihit) {
         if (par.eta2[ihit]>-1.85779 && par.eta2[ihit]<-1.85778 && par.phi2[ihit]>2.18935 && par.phi2[ihit]<2.18936)
            continue;
         if (par.eta2[ihit]>0.12535 && par.eta2[ihit]<0.12536 && par.phi2[ihit]>1.03344 && par.phi2[ihit]<1.03345)
            continue;
         if (par.eta2[ihit]>0.12626 && par.eta2[ihit]<0.12631 && par.phi2[ihit]>1.034 && par.phi2[ihit]<1.03475)
            continue;
         if (par.eta2[ihit]>-1.2623 && par.eta2[ihit]<-1.26229 && par.phi2[ihit]>2.8906 && par.phi2[ihit]<2.8907)
            continue;
         if (par.phi2[ihit]>2.1157 && par.phi2[ihit]<2.1158 && par.eta2[ihit]>-1.9016 && par.eta2[ihit]<-1.9015)
            continue;
         RecoHit tmp(par.eta2[ihit], par.phi2[ihit], par.r2[ihit], par.cs2[ihit], par.ch2[ihit], 2);
         if (gRandom->Rndm()<dropProb) continue;
         hits.push_back(tmp);
         // put artifical split hits
         if (gRandom->Rndm()<splitProb) hits.push_back(tmp);
      }
   } else if (layer == 3) {
      for (int ihit = 0; ihit < par.nhits3; ++ihit) {
         if (par.eta3[ihit]>0.76085 && par.eta3[ihit]<0.76087 && par.phi3[ihit]>-1.22261 && par.phi3[ihit]<-1.2226)
            continue;
         if (par.eta3[ihit]>-1.39076 && par.eta3[ihit]<-1.39074 && par.phi3[ihit]>1.116954 && par.phi3[ihit]<1.116955)
            continue;
         if (par.eta3[ihit]>-1.07629 && par.eta3[ihit]<-1.07628 && par.phi3[ihit]>-2.59858 && par.phi3[ihit]<-2.59857)
            continue;
         if (par.eta3[ihit]>0.99493 && par.eta3[ihit]<0.99494 && par.phi3[ihit]>0.75637 && par.phi3[ihit]<0.7563704)
            continue;
         if (par.eta3[ihit]>1.30654 && par.eta3[ihit]<1.30655 && par.phi3[ihit]>-0.49127 && par.phi3[ihit]<-0.49126)
            continue;
         if (par.eta3[ihit]>-0.17039 && par.eta3[ihit]<-0.17038 && par.phi3[ihit]>-1.44136 && par.phi3[ihit]<-1.441353)
            continue;
         if (par.phi3[ihit]>1.11 && par.phi3[ihit]<1.17 && par.eta3[ihit]>-1.4 && par.eta3[ihit]<-1.34)
            continue;
         if (par.phi3[ihit]>-0.522 && par.phi3[ihit]<-0.509 && par.eta3[ihit]>0.8862 && par.eta3[ihit]<0.8864)
            continue;
         if (par.phi3[ihit]>0.5434 && par.phi3[ihit]<0.5446 && par.eta3[ihit]>0.9554 && par.eta3[ihit]<0.9561)
            continue;
         if (par.phi3[ihit]>1.269 && par.phi3[ihit]<1.2694 && par.eta3[ihit]>0.64 && par.eta3[ihit]<0.72)
            continue;
         if (par.phi3[ihit]>-0.2274 && par.phi3[ihit]<-0.2271 && par.eta3[ihit]>1.43 && par.eta3[ihit]<1.475)
            continue;
         if (par.phi3[ihit]>-1.6671 && par.phi3[ihit]<-1.6669 && par.eta3[ihit]>1.408 && par.eta3[ihit]<1.4084)
            continue;
         if (par.phi3[ihit]>0.7193 && par.phi3[ihit]<0.7194 && par.eta3[ihit]>-0.2824 && par.eta3[ihit]<-0.2823)
            continue;
         if (par.phi3[ihit]>-1.2439 && par.phi3[ihit]<-1.2438 && par.eta3[ihit]>-0.693 && par.eta3[ihit]<-0.692)
            continue;
         if (par.phi3[ihit]>0.064 && par.phi3[ihit]<0.073 && par.eta3[ihit]>-1.04 && par.eta3[ihit]<-0.98)
            continue;
         if (par.phi3[ihit]>1.024 && par.phi3[ihit]<1.0242 && par.eta3[ihit]>0.7992 && par.eta3[ihit]<0.7995)
            continue;
         if (par.phi3[ihit]>0.07 && par.phi3[ihit]<0.073 && par.eta3[ihit]>-0.63 && par.eta3[ihit]<-0.54)
            continue;
         RecoHit tmp(par.eta3[ihit], par.phi3[ihit], par.r3[ihit], par.cs3[ihit], par.ch3[ihit], 3);
         if (gRandom->Rndm()<dropProb) continue;
         hits.push_back(tmp);
         // put artifical split hits
         if (gRandom->Rndm()<splitProb) hits.push_back(tmp);
      }
   } else if (layer%10 == 4) {
      for (int ihit = 0; ihit < par.nhits4; ++ihit) {
         if ((layer/10==2 && par.r4[ihit]<9.0) || (layer/10==3 && par.r4[ihit]>9.0)) {
            RecoHit tmp(par.eta4[ihit], par.phi4[ihit], par.r4[ihit], par.cs4[ihit], par.ch4[ihit], layer/10);
            hits.push_back(tmp);
         }
      }
   } else if (layer%10 == 5) {
      for (int ihit = 0; ihit < par.nhits5; ++ihit) {
         if ((layer/10==2 && par.r5[ihit]<9.0) || (layer/10==3 && par.r5[ihit]>9.0)) {
            RecoHit tmp(par.eta5[ihit], par.phi5[ihit], par.r5[ihit], par.cs5[ihit], par.ch5[ihit], layer/10);
            hits.push_back(tmp);
         }
      }
   }
   sort (hits.begin(), hits.end(), comparePhi);

   for (int ihit = 0; ihit<(int)hits.size(); ++ihit) {
      // double dr = 0;
      // double dphi = 10;
      // double deta = 10;
      // int flag = 0;
      // if (ihit != 0) {
      //    for (int k=ihit-1; k<ihit; k++) {
      //       dphi = fabs(calcDphi(hits[k].phi, hits[ihit].phi));
      //       deta = fabs(hits[k].eta - hits[ihit].eta);
      //       dr   = fabs(hits[k].r - hits[ihit].r);
      //       // no double hit removal...
      //       // if (dr>cuts.drCut && dphi<cuts.dPhiCut) flag=1;
      //       // if (dphi > cuts.dPhiCut) k=0;
      //    }
      // }
      // if (flag==1) continue;

      double x = hits[ihit].r*cos(hits[ihit].phi);
      double y = hits[ihit].r*sin(hits[ihit].phi);
      double z = hits[ihit].r/tan(atan(exp(-hits[ihit].eta))*2);

      if (smearPixels) {
         x += gRandom->Gaus(0, 0.005);
         y += gRandom->Gaus(0, 0.005);
         z += gRandom->Gaus(0, 0.005);
      }

      ROOT::Math::XYZVector tmpVector(x-x0, y-y0, z-vz);
      RecoHit tmpHit(tmpVector.eta(), tmpVector.phi(), tmpVector.rho(), hits[ihit].cs, hits[ihit].ch, hits[ihit].layer);
      double eta = tmpVector.eta();

      if (cutOnClusterSize && fabs(eta)<=0.5 &&                  hits[ihit].cs<1) continue;
      if (cutOnClusterSize && fabs(eta)<=1.0 && fabs(eta)>0.5 && hits[ihit].cs<2) continue;
      if (cutOnClusterSize && fabs(eta)<=1.5 && fabs(eta)>1.0 && hits[ihit].cs<3) continue;
      if (cutOnClusterSize && fabs(eta)<=2.0 && fabs(eta)>1.5 && hits[ihit].cs<4) continue;
      if (cutOnClusterSize && fabs(eta)<=2.5 && fabs(eta)>2.0 && hits[ihit].cs<6) continue;
      if (cutOnClusterSize && fabs(eta)<=5.0 && fabs(eta)>2.5 && hits[ihit].cs<9) continue;
      cleanedHits.push_back(tmpHit);
   }
}

// RecoHit RandomHit(double etaMin, double etaMax, double phiMin, double phiMax) {
//    double eta = etaMin+(etaMax-etaMin)*gRandom->Rndm();
//    double phi = phiMin+(phiMax-phiMin)*gRandom->Rndm();
//    RecoHit myRandomHit(eta, phi, 0, 100);
//    return myRandomHit;
// }

double calcDphi(double phi1_, double phi2_) {
   double pi = 3.14159265358979;
   double dphi = phi1_ - phi2_;

   if (dphi>0){
      while (dphi>2*pi) dphi -= 2*pi;
      if (dphi>pi) dphi = 2*pi-dphi;
   } else {
      while (dphi<-2*pi) dphi += 2*pi;
      if (dphi<-pi) dphi = -2*pi-dphi;
   }

   return dphi;
}

// void combineRecHit(vector<RecoHit> &c, vector<RecoHit> a, vector<RecoHit> b) {
//    for (unsigned int i=0; i<a.size(); i++)
//       c.push_back(a[i]);
//    for (unsigned int i=0; i<b.size(); i++)
//       c.push_back(b[i]);
// }

void getPixelTreeBranch(TTree* t, Parameters &par) {
   t->SetBranchAddress("nRun", &par.nRun);
   t->SetBranchAddress("nEv", &par.nEv);
   t->SetBranchAddress("nLumi", &par.nLumi);
   t->SetBranchAddress("nBX", &par.nBX);
   t->SetBranchAddress("nHFn", &par.nHFn);
   t->SetBranchAddress("nHFp", &par.nHFp);

   t->SetBranchAddress("nHLTBit", &par.nHltBit);
   t->SetBranchAddress("hltBit", par.hltBit);
   t->SetBranchAddress("nL1A", &par.nL1ABit);
   t->SetBranchAddress("L1A", par.l1ABit);
   t->SetBranchAddress("nL1T", &par.nL1TBit);
   t->SetBranchAddress("L1T", par.l1TBit);

   t->SetBranchAddress("L1Tech_BPTX_plus_AND_minus.v0", &par.L1_BPTX_AND);
   t->SetBranchAddress("L1Tech_BPTX_plus_OR_minus.v0", &par.L1_BPTX_OR);
   t->SetBranchAddress("L1Tech_BPTX_plus.v0", &par.L1_BPTX_plus);
   t->SetBranchAddress("L1Tech_BPTX_minus.v0", &par.L1_BPTX_minus);

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
   t->SetBranchAddress("eta4", par.eta4);
   t->SetBranchAddress("phi4", par.phi4);
   t->SetBranchAddress("r4", par.r4);
   t->SetBranchAddress("cs4", par.cs4);
   t->SetBranchAddress("ch4", par.ch4);
   t->SetBranchAddress("eta5", par.eta5);
   t->SetBranchAddress("phi5", par.phi5);
   t->SetBranchAddress("r5", par.r5);
   t->SetBranchAddress("cs5", par.cs5);
   t->SetBranchAddress("ch5", par.ch5);
   t->SetBranchAddress("nhits1", &par.nhits1);
   t->SetBranchAddress("nhits2", &par.nhits2);
   t->SetBranchAddress("nhits3", &par.nhits3);
   t->SetBranchAddress("nhits4", &par.nhits4);
   t->SetBranchAddress("nhits5", &par.nhits5);

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

   t->SetBranchAddress("xi", &par.xi);
   t->SetBranchAddress("passDS", &par.passDS);
   t->SetBranchAddress("passSingleTrack", &par.passSingleTrack);
   t->SetBranchAddress("ntrks", &par.ntrks);
   t->SetBranchAddress("ntrksCut", &par.ntrksCut);
}

#endif /* RECHOHIT */

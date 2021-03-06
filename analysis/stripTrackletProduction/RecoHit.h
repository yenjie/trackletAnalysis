#define maxEntry 10000
#define maxEntry2 10000

#include <vector>
#include <algorithm>
#include <Math/Vector3D.h>
#include <TRandom.h>
#include <TTree.h>
#include <TMath.h>

using namespace std;

class RecoHit {
   public:
      RecoHit(double _eta, double _phi, double _r, double _cs) {
         eta = _eta;
         phi = _phi;
         r = _r;
         cs = _cs;
      };
      RecoHit(double _eta, double _phi, double _r, double _cs, int _l) {
         eta = _eta;
         phi = _phi;
         r = _r;
         cs = _cs;
         layer = _l;
      };

      ~RecoHit() {};

      double eta;
      double phi;
      double r;
      int layer;
      double cs;   // cluster size
      double ch;   // cluster charge
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
      float beamSpotX, beamSpotY, beamSpotZ;
      float vx[8], vy[8], vz[8];
      float eta1[maxEntry], phi1[maxEntry], r1[maxEntry], cs1[maxEntry], ch1[maxEntry];
      float eta2[maxEntry], phi2[maxEntry], r2[maxEntry], cs2[maxEntry], ch2[maxEntry];
      float eta3[maxEntry], phi3[maxEntry], r3[maxEntry], cs3[maxEntry], ch3[maxEntry];
      float eta4[maxEntry], phi4[maxEntry], r4[maxEntry], cs4[maxEntry], ch4[maxEntry];
      float eta[maxEntry], phi[maxEntry], pt[maxEntry];
      int nhits1, nhits2, nhits3, nhits4, mult, nv, npart, evtType, chg[maxEntry], pdg[maxEntry];
      float xi;
      bool passDS, passSingleTrack;
      int ntrks, ntrksCut;
};

class TrackletData {
   public:
      int nRun, nEv, nLumi, nHltBit, nL1ABit, nL1TBit, nBX, nHFn, nHFp, nHits;
      bool hltBit[500], l1ABit[500], l1TBit[500];
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
   x0 = 0.032;
   y0 = 0.0;

   if (layer == 1) {
      for (int ihit = 0; ihit < par.nhits1; ++ihit) {
         if (par.eta1[ihit]>-1.7 && par.eta1[ihit]<0.1 && par.phi1[ihit]>-1 && par.phi1[ihit]<0)
            continue;
         if (par.eta1[ihit]>0.06 && par.eta1[ihit]<0.45 && par.phi1[ihit]>-0.04 && par.phi1[ihit]<0.22)
            continue;
         if (par.eta1[ihit]>0.84 && par.eta1[ihit]<1.1 && par.phi1[ihit]>-0.24 && par.phi1[ihit]<0.22)
            continue;
         RecoHit tmp(par.eta1[ihit], par.phi1[ihit], par.r1[ihit], par.cs1[ihit], 1);
         if (gRandom->Rndm()<dropProb) continue;
         hits.push_back(tmp);
         // put artifical split hits
         if (gRandom->Rndm()<splitProb) hits.push_back(tmp);
      }
   } else if (layer == 2) {
      for (int ihit = 0; ihit < par.nhits2; ++ihit) {
         RecoHit tmp(par.eta2[ihit], par.phi2[ihit], par.r2[ihit], par.cs2[ihit], 2);
         if (gRandom->Rndm()<dropProb) continue;
         hits.push_back(tmp);
         // put artifical split hits
         if (gRandom->Rndm()<splitProb) hits.push_back(tmp);
      }
   } else if (layer == 3) {
      for (int ihit = 0; ihit < par.nhits3; ++ihit) {
         RecoHit tmp(par.eta3[ihit], par.phi3[ihit], par.r3[ihit], par.cs3[ihit], 3);
         if (gRandom->Rndm()<dropProb) continue;
         hits.push_back(tmp);
         // put artifical split hits
         if (gRandom->Rndm()<splitProb) hits.push_back(tmp);
      }
   } else if (layer == 4) {
      for (int ihit = 0; ihit < par.nhits4; ++ihit) {
         RecoHit tmp(par.eta4[ihit], par.phi4[ihit], par.r4[ihit], par.cs4[ihit], 4);
         if (gRandom->Rndm()<dropProb) continue;
         hits.push_back(tmp);
         // put artifical split hits
         if (gRandom->Rndm()<splitProb) hits.push_back(tmp);
      }
   }
   sort (hits.begin(), hits.end(), comparePhi);

   for (int ihit = 0; ihit < (int)hits.size(); ++ihit) {
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

      if (vz!=0 && firstCall==0) {
         cout << "Beamspot X0 = " << x0 << " Y0 = " << y0 << endl;
         firstCall = 1;
      }

      ROOT::Math::XYZVector tmpVector(x-x0, y-y0, z-vz);
      RecoHit tmpHit(tmpVector.eta(), tmpVector.phi(), tmpVector.rho(), hits[ihit].cs);
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
   t->SetBranchAddress("nhits1", &par.nhits1);
   t->SetBranchAddress("nhits2", &par.nhits2);
   t->SetBranchAddress("nhits3", &par.nhits3);
   t->SetBranchAddress("nhits4", &par.nhits4);

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

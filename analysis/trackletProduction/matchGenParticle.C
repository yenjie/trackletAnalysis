#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include "Tracklet.h"
#include <iostream>


double angularRangeReduce(double x)
{
  const double cody_waite_x_max = 1608.4954386379741381;
  const double two_pi_0 = 6.2831853071795649157;
  const double two_pi_1 = 2.1561211432631314669e-14;
  const double two_pi_2 = 1.1615423895917441336e-27;
  double ret = 0;

  if(x >= -cody_waite_x_max && x <= cody_waite_x_max) {
    const double inverse_two_pi =
      0.15915494309189534197;
    const double k = rint(x * inverse_two_pi);
    ret = ((x - (k * two_pi_0)) - k * two_pi_1) -
      k * two_pi_2;
  }

  return ret;
}

double deltaPhi(double phi1, double phi2) {
  return angularRangeReduce(phi1 - phi2);
}

double deltaEta(double eta1, double eta2) {
  return eta1 - eta2;
}
 
double deltaR2(double eta1, double phi1, double eta2, double phi2) {
    double deta = deltaEta(eta1,eta2);
    double dphi = deltaPhi(phi1,phi2);
    return (deta*deta+dphi*dphi);
}

double deltaR(double eta1, double phi1, double eta2, double phi2) {
    double r2 = (deltaR2(eta1,phi1,eta2,phi2));
    return TMath::Sqrt(r2);
}

void setBranchAddress(TTree *trackletTree, TrackletData &tdata)
{
   trackletTree->SetBranchAddress("nRun", &tdata.nRun);
   trackletTree->SetBranchAddress("nEv", &tdata.nEv);
   trackletTree->SetBranchAddress("nLumi", &tdata.nLumi);
   trackletTree->SetBranchAddress("nBX", &tdata.nBX);
   trackletTree->SetBranchAddress("nHFn", &tdata.nHFn);
   trackletTree->SetBranchAddress("nHFp", &tdata.nHFp);
   trackletTree->SetBranchAddress("nHits", &tdata.nHits);

   trackletTree->SetBranchAddress("nHltBit", &tdata.nHltBit);
   trackletTree->SetBranchAddress("hltBit", tdata.hltBit);
   trackletTree->SetBranchAddress("nL1ABit", &tdata.nL1ABit);
   trackletTree->SetBranchAddress("l1ABit", tdata.l1ABit);
   trackletTree->SetBranchAddress("nL1TBit", &tdata.nL1TBit);
   trackletTree->SetBranchAddress("l1TBit", tdata.l1TBit);

   trackletTree->SetBranchAddress("nTracklets", &tdata.nTracklet);
   trackletTree->SetBranchAddress("nhit1", &tdata.nhit1);
   trackletTree->SetBranchAddress("nhit2", &tdata.nhit2);
   trackletTree->SetBranchAddress("mult", &tdata.mult);
   trackletTree->SetBranchAddress("mult2", &tdata.mult2);
   trackletTree->SetBranchAddress("nv", &tdata.nv);
   trackletTree->SetBranchAddress("vx", tdata.vx);
   trackletTree->SetBranchAddress("vy", tdata.vy);
   trackletTree->SetBranchAddress("vz", tdata.vz);
   trackletTree->SetBranchAddress("eta1", tdata.eta1);
   trackletTree->SetBranchAddress("phi1", tdata.phi1);
   trackletTree->SetBranchAddress("r1", tdata.r1);
   trackletTree->SetBranchAddress("cs1", tdata.cs1);
   trackletTree->SetBranchAddress("eta2", tdata.eta2);
   trackletTree->SetBranchAddress("phi2", tdata.phi2);
   trackletTree->SetBranchAddress("r2", tdata.r2);
   trackletTree->SetBranchAddress("cs2", tdata.cs2);
   trackletTree->SetBranchAddress("deta", tdata.deta);
   trackletTree->SetBranchAddress("dphi", tdata.dphi);
   trackletTree->SetBranchAddress("recoPU", &tdata.recoPU);

   trackletTree->SetBranchAddress("npart", &tdata.npart);
   trackletTree->SetBranchAddress("eta", tdata.eta);
   trackletTree->SetBranchAddress("phi", tdata.phi);
   trackletTree->SetBranchAddress("pdg", tdata.pdg);
   trackletTree->SetBranchAddress("chg", tdata.chg);
   trackletTree->SetBranchAddress("nhad", tdata.nhad);
   trackletTree->SetBranchAddress("pt", tdata.pt);
   trackletTree->SetBranchAddress("evtType", &tdata.evtType);
   trackletTree->SetBranchAddress("pro2", &tdata.pro2);

   trackletTree->SetBranchAddress("nPU", &tdata.nPU);

   trackletTree->SetBranchAddress("xi", &tdata.xi);
   trackletTree->SetBranchAddress("passDS", &tdata.passDS);
   trackletTree->SetBranchAddress("passSingleTrack", &tdata.passSingleTrack);
   trackletTree->SetBranchAddress("ntrks", &tdata.ntrks);
   trackletTree->SetBranchAddress("ntrksCut", &tdata.ntrksCut);

/*
   trackletTree->SetAlias("dR", "sqrt(deta*deta+dphi*dphi)");
   trackletTree->SetAlias("dRR", "sqrt(deta*deta+dphi*dphi+(r1-r2)*(r1-r2))") ;
   trackletTree->SetAlias("z1", "r1/tan(atan(exp(-eta1))*2)+vz[3]");
   trackletTree->SetAlias("z2", "r2/tan(atan(exp(-eta2))*2)+vz[3]");
   trackletTree->SetAlias("x1", "r1*cos(phi1)");
   trackletTree->SetAlias("x2", "r2*cos(phi2)");
   trackletTree->SetAlias("y1", "r1*sin(phi1)");
   trackletTree->SetAlias("y2", "r2*sin(phi2)");
*/
}

void matchGenParticle(char *infname, char *outfname)
{
   TFile *inf = new TFile(infname);
   TTree *t = (TTree*)inf->Get("TrackletTree14");
   

   TFile *outf = new TFile(outfname,"recreate");
   TNtuple *nt = new TNtuple("nt","","pt:eta:phi:eta1:phi1:dphi:deta:dR");
   TrackletData tdata;
   
   setBranchAddress(t,tdata);
   cout <<t->GetEntries()<<endl;
   for (int i=0;i<t->GetEntries();i++)
   {
      t->GetEntry(i);
      for (int j=0;j<tdata.npart;j++)
      {
         double minDR = 1000;
	 int index = -10;
	 for (int k=0;k<tdata.nTracklet;k++)
	 {
	    double dR = deltaR(tdata.eta1[k],tdata.phi1[k],tdata.eta[j],tdata.phi[j]);
	    if (dR<minDR) {
	       minDR=dR;
	       index = k;
	    }
	 }
	 nt->Fill(tdata.pt[j],tdata.eta[j],tdata.phi[j],tdata.eta1[index],tdata.phi1[index],tdata.dphi[index],tdata.deta[index],minDR);
      }
   }
   
   nt->Write();
   
}

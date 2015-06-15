#ifndef CLUSTERVERTEXCOMPAT
#define CLUSTERVERTEXCOMPAT

#include <TFile.h>
#include <TNtuple.h>
#include <TTimeStamp.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TString.h>

#include "RecoHit.h"

Int_t getContainedHits(const std::vector<RecoHit> &hits, double z0, double &diff, double slope = 0.43);

Double_t getClusVtxCompat(const std::vector<RecoHit> &hits, const Int_t layer = 1) {

  const Double_t   minZ_ = -30.;          // beginning z-vertex position
  const Double_t   maxZ_ = 30.;           // end z-vertex position
  const Double_t   zStep_ = 0.01;         // size of steps in z-vertex test

  Double_t slope = 0.43;
  if(layer==1)      slope = 0.43;
  else if(layer==2) slope = 0.25;
  else if(layer==3) slope = 0.18;

  // estimate z-position from cluster lengths
  Double_t zest = 0.;
  Int_t nhits = 0, nhits_max = 0;
  Double_t diff = 0, diff_max = 1e+9;
  for(Double_t z0 = minZ_; z0 <= maxZ_; z0 += zStep_) {
    nhits = getContainedHits(hits, z0, diff, slope);
    if(nhits == 0)
      continue;
    if(nhits > nhits_max) {
      diff_max = 1e+9;
      nhits_max = nhits;
    }
    if(nhits >= nhits_max && diff < diff_max) {
      diff_max = diff;
      zest = z0;
    }
  }

  diff = 0;
  int nbest=0, nminus=0, nplus=0;
  nbest = getContainedHits(hits,zest,diff,slope);
  nminus = getContainedHits(hits,zest-10.,diff,slope);
  nplus = getContainedHits(hits,zest+10.,diff,slope);
  
  Double_t clusVtxQual=0.0;
  if ((nminus+nplus)> 0)
    clusVtxQual = (2.0*nbest)/(nminus+nplus);  // A/B
  else if(nbest>0)
    clusVtxQual = 1e6; //some very large number
  
  return clusVtxQual;
}

Int_t getContainedHits(const std::vector<RecoHit> &hits, double z0, double &diff, double slope) {
  // Calculate number of hits contained in v-shaped window in cluster y-width vs. z-position.
  Int_t n = 0;
  diff   = 0.;
     
  for(std::vector<RecoHit>::const_iterator hit = hits.begin(); hit!= hits.end(); ++hit) {
    Double_t z = hit->r/tan(2.*atan(exp(-hit->eta)));
    Double_t exp = 1.5+slope*fabs(z-z0);
    if(fabs(exp - hit->cs) <= 1.) { //v-shape with 1 cm width
      diff += fabs(exp - hit->cs);
      n++;
    }
  }
  return n;
}

#endif

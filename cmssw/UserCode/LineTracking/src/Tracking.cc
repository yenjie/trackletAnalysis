#include "../interface/Tracking.h"

#include <fstream>
#include <iostream>
#include <cmath>

#include "TMatrixD.h"

#define sqr(x) ((x) * (x))

using namespace std;

/*****************************************************************************/
class Cluster
{
 public:
  vector<double> eta;  // z position
  double phi;  // sigma_z^2
  int    n;    // number of attached points
  vector<Cluster>::iterator minCluster;  // closest cluster
  double                    minDistance; // smallest distance
  unsigned int              minVertex;   // primary vertex favored
  
  vector<int> list; // list of points

  bool use;

  pair<unsigned long int, unsigned long int> detids;

//  double r,z;
};

/*****************************************************************************/
Tracking::Tracking
  (const double & dMax_, const vector<double> & z0_) : dMax(dMax_), z0(z0_)
{
}

/*****************************************************************************/
bool Tracking::findClosest
  (vector<Cluster> & clusters, vector<Cluster>::iterator c1)
{
  bool isFirst = true;

  for(vector<Cluster>::iterator c2 = clusters.begin();
                                c2!= clusters.end(); c2++)
  if(c2->use)
  if(c1 != c2)
  for(unsigned int j = 0; j < z0.size(); j++) // try all primary vertices
  if(c1->n == 1 || c1->minVertex == j)
  if(c2->n == 1 || c2->minVertex == j)
  if(c1->detids != c2->detids) // fIXME
  {
    double dphi = fabs(c1->phi - c2->phi);
    if(dphi > M_PI) dphi -= 2*M_PI;

    double dist = sqr(c1->eta[j]  - c2->eta[j]) + sqr(dphi);

    if(dist < c1->minDistance || isFirst)
    {
      c1->minCluster  = c2;
      c1->minDistance = dist;
      c1->minVertex   = j;

      isFirst = false;
    }
  }

  return !isFirst;
}

/*****************************************************************************/
void Tracking::run
  (const vector<TVector3> & points,
   const vector<pair<unsigned long int, unsigned long int> > & detids,
   unsigned int & nOptimal,
   vector<vector<int> > & lists,
   unsigned int maxResult)
{
  // Setup
  nOptimal = points.size();

  vector<Cluster> clusters;

  // Initialize clusters
  for(vector<TVector3>::const_iterator p = points.begin();
                                       p!= points.end(); p++)
  {
    Cluster cluster;

    for(vector<double>::const_iterator z = z0.begin();
                                       z!= z0.end(); z++)
    {
      TVector3 V(0,0,*z);
      cluster.eta.push_back((*p - V).Eta());
    }

    cluster.phi = p->Phi();
    cluster.n   = 1;
    cluster.use = true;

    cluster.detids = detids[int(p - points.begin())];

    vector<int> l;
    l.push_back(p - points.begin());

    cluster.list = l;

    clusters.push_back(cluster);
  }

  lists.clear();

  for(vector<Cluster>::iterator c1 = clusters.begin();
                                c1!= clusters.end(); c1++)
  if(c1->use)
    lists.push_back(c1->list);

  unsigned int nUse = clusters.size() - 1;

  // Find nearest neighbors
  bool ok = false;
  for(vector<Cluster>::iterator c1 = clusters.begin();
                                c1!= clusters.end(); c1++)
    if(findClosest(clusters, c1)) ok = true;

  // return if there are no pairs to look at
  if(!ok) return;

  int nSteps  = 0;
  int nUpdate = 0;

  if(nUse > 0)
  while(nUse > 0)
  {
    vector<Cluster>::iterator c[2];
    double minDistance = 0;
    bool isFirst = true;
    unsigned int cv = 0;

    // Find smallest distance
    for(vector<Cluster>::iterator c1 = clusters.begin();
                                  c1!= clusters.end(); c1++)
    if(c1->use)
    if(c1->minDistance < minDistance || isFirst)
    {
      minDistance = c1->minDistance;

      c[0] = c1 ;
      c[1] = c1->minCluster;
      cv   = c1->minVertex;

      isFirst = false;
    }

    // Join 
    double eta = (c[0]->n*c[0]->eta[cv] + c[1]->n*c[1]->eta[cv])/
                 (c[0]->n + c[1]->n);

    double phi = atan2(c[0]->n*sin(c[0]->phi) + c[1]->n*sin(c[1]->phi),
                       c[0]->n*cos(c[0]->phi) + c[1]->n*cos(c[1]->phi));

    // Update
    c[0]->eta[cv] = eta;
    c[0]->phi = phi;     
    c[0]->n   = c[0]->n + c[1]->n;
    c[0]->minVertex = cv;

    for(vector<int>::iterator il = c[1]->list.begin();
                              il!= c[1]->list.end(); il++) 
      c[0]->list.push_back(*il);

    // Remove c[1]
    c[1]->use  = false;

    //
    if(minDistance < sqr(dMax))
    {
      nOptimal = nUse;

      lists.clear();

      for(vector<Cluster>::iterator c1 = clusters.begin();
                                    c1!= clusters.end(); c1++)
      if(c1->use)
        lists.push_back(c1->list);
    }

    for(vector<Cluster>::iterator c1 = clusters.begin();
                                  c1!= clusters.end(); c1++)
    if(c1->use && (c1->minCluster == c[0] ||
                   c1->minCluster == c[1]))
    {
      nUpdate++;
      findClosest(clusters, c1);
    }

    nSteps++;

    nUse--;
  }
}

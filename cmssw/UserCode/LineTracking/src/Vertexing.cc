#include "../interface/Vertexing.h"

#include "TMatrixD.h"
#include <iostream>
#include <cmath>

#define sqr(x) ((x) * (x))

using namespace std;

/*****************************************************************************/
class VCluster
{
 public:
  double pos;   // z position
  double sig2;  // sigma_z^2
  int    n;     // number of attached tracks
  vector<VCluster>::iterator minCluster;  // closest cluster
  double                     minDistance; // smallest distance
  
  vector<int> list; // list of tracks

  bool use;
};

/*****************************************************************************/
Vertexing::Vertexing(double dMax_)
{
  dMax = dMax_;
}

/*****************************************************************************/
void Vertexing::findClosest
  (vector<VCluster> & clusters, vector<VCluster>::iterator c1)
{
  bool isFirst = true;

  for(vector<VCluster>::iterator c2 = clusters.begin();
                                 c2!= clusters.end(); c2++)
  if(c2->use)
  if(c1 != c2)
  {
    double dist = sqr(c1->pos  - c2->pos ) /
                     (c1->sig2 + c2->sig2);

    if(dist < c1->minDistance || isFirst)
    {
      c1->minCluster  = c2;
      c1->minDistance = dist;

      isFirst = false;
    }
  }
}

/*****************************************************************************/
void Vertexing::run
  (const vector<pair<double,double> > & points,
   vector<pair<TVectorD, TVectorD> > & result,
   unsigned int & nOptimal,
   vector<vector<int> > & lists,
   unsigned int maxResult)
{
  // Setup
  nOptimal = points.size();

  vector<VCluster> clusters;

  // Initialize clusters
  for(vector<pair<double,double> >::const_iterator p = points.begin();
                                                   p!= points.end(); p++)
  {
    VCluster cluster;

    cluster.pos   = p->first;
    cluster.sig2  = p->second;
    cluster.n     = 1;
    cluster.use   = true;

    vector<int> l;
    l.push_back(p - points.begin());

    cluster.list = l;

    clusters.push_back(cluster);
  }

  lists.clear();

  for(vector<VCluster>::iterator c1 = clusters.begin();
                                c1!= clusters.end(); c1++)
  if(c1->use)
    lists.push_back(c1->list);

  unsigned int nUse = clusters.size() - 1;

  // Find nearest neighbors
  for(vector<VCluster>::iterator c1 = clusters.begin();
                                c1!= clusters.end(); c1++)
    findClosest(clusters, c1);

  int nSteps  = 0;
  int nUpdate = 0;

  if(nUse > 0)
  while(nUse > 0)
  {
    vector<VCluster>::iterator c[2];
    double minDistance = 0;
    bool isFirst = true;

    // Find smallest distance
    for(vector<VCluster>::iterator c1 = clusters.begin();
                                  c1!= clusters.end(); c1++)
    if(c1->use)
    if(c1->minDistance < minDistance || isFirst)
    {
      minDistance = c1->minDistance;

      c[0] = c1 ;
      c[1] = c1->minCluster ;

      isFirst = false;
    }

    // Join 
    double sig2 = 1 /
                  (        1 / c[0]->sig2 +         1 / c[1]->sig2);
    double pos  = (c[0]->pos / c[0]->sig2 + c[1]->pos / c[1]->sig2) * sig2;

    // Update
    c[0]->pos  = pos;     
    c[0]->sig2 = sig2;     
    c[0]->n    = c[0]->n + c[1]->n;

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

      for(vector<VCluster>::iterator c1 = clusters.begin();
                                     c1!= clusters.end(); c1++)
      if(c1->use)
        lists.push_back(c1->list);
    }
    else
      cerr
        << "  fPNN " << nUse << " stopped at : " << sqrt(minDistance)
        << " with " << c[0]->n << " tracks" << endl;

    for(vector<VCluster>::iterator c1 = clusters.begin();
                                   c1!= clusters.end(); c1++)
    if(c1->use && (c1->minCluster == c[0] ||
                   c1->minCluster == c[1]))
    {
      nUpdate++;
      findClosest(clusters, c1);
    }

    nSteps++;

    // Save to clusters
    if(nUse <= maxResult)
    {
      int k = 0;

      TVectorD mu(nUse);
      TVectorD P(nUse);

      for(unsigned int i = 0  ; i < clusters.size(); i++)
      if(clusters[i].use)
      {
        // average position
        mu(k) = clusters[i].pos;

        // probability, with weight
        P(k) = float(clusters[i].n)/nUse;

        k++;
      }

      result[nUse] = pair<TVectorD,TVectorD>(mu,P);
    }
    nUse--;
  }
  else
  {
    // Single track
    nUse = 1;

    int k = 0;

    TVectorD mu(nUse);
    TVectorD P(nUse);

    for(unsigned int i = 0  ; i < clusters.size(); i++)
    if(clusters[i].use)
    {
      // average position
      mu(k) = clusters[i].pos;

      // probability, with weight
      P(k) = float(clusters[i].n)/nUse;

      k++;
    }

    result[nUse] = pair<TVectorD,TVectorD>(mu,P);
  }
}

#ifndef _Vertexing_
#define _Vertexing_

#include <utility>
#include <vector>

#include "TVectorD.h"

class VCluster;

class Vertexing
{
 public:
  Vertexing(double dMax);

  void run(const std::vector<std::pair<double,double> > & points,
                 std::vector<std::pair<TVectorD, TVectorD> > & clusters,
                 unsigned int & nOptimal,
                 std::vector<std::vector<int> > & lists,
                 unsigned int maxClusters); 

 private:
  void findClosest(std::vector<VCluster> & clusters,
                   std::vector<VCluster>::iterator c1);
  double dMax;
};

#endif

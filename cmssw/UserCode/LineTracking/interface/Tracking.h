#ifndef _Tracking_
#define _Tracking_

#include <utility>
#include <vector>

#include "TVectorD.h"
#include "TVector3.h"

class Cluster;

class Tracking
{
 public:
  Tracking(const double & dMax_, const std::vector<double> & z0_);

  void run(const std::vector<TVector3> & points,
           const std::vector<std::pair<unsigned long int, unsigned long int> > & detids,
           unsigned int & nOptimal,
           std::vector<std::vector<int> > & lists,
           unsigned int maxClusters); 

 private:
  bool findClosest(std::vector<Cluster> & clusters,
                   std::vector<Cluster>::iterator c1);
  double dMax;
  std::vector<double> z0;
};

#endif

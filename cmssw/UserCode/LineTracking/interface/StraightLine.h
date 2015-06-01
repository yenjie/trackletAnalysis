#ifndef _StraightLine_
#define _StraightLine_

#include <vector>

#include "TVector3.h"

class StraightLine
{
 public:
  StraightLine() {}

  double fit(const std::vector<TVector3> & points,
                   std::vector<double> & x);

 private:
  double value(const std::vector<TVector3> & points,
               const std::vector<double> & x);

  int nHitsForSum;
};

#endif

#include "../interface/StraightLine.h"

#include <iostream>
#include <cmath>

#include "TVectorD.h"
#include "TMatrixD.h"

#include "TDecompLU.h"

using namespace std;

#define sqr(x) ((x) * (x))

/*****************************************************************************/
double StraightLine::value(const vector<TVector3> & points,
                           const vector<double> & x)
{
  const double & z0    = x[0];
  const double & d0    = x[1];
  const double & theta = x[2];
  const double & phi   = x[3];

  double sumd2 = 0.;

  for(vector<TVector3>::const_iterator p = points.begin();
                                       p!= points.end(); p++)
  if(p < points.begin() + nHitsForSum)
  {
    TVector3 V(-d0*sin(theta), d0*cos(phi), z0);
    TVector3 u(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));

    TVector3 x = u.Cross(*p - V);

    sumd2 += x.Mag2();
  }

  return sumd2;
}


/*****************************************************************************/
double StraightLine::fit(const vector<TVector3> & points,
                               vector<double> & x)
{
  nHitsForSum = 3;

  const TVector3 & a = points.front();
  const TVector3 & b = points.back();

  double ra = a.Perp();
  double rb = b.Perp();
  
  // Initial guess
  /* z0    */ x[0] = a.z() - ra/(rb-ra) * (b.z() - a.z());
  /* d0    */ x[1] = 0.;
  /* theta */ x[2] = (b-a).Theta();
  /* phi   */ x[3] = (b-a).Phi();

  for(int k = 0; k < 5; k++)
  { 
  const double de = 1e-3;

  TVectorD F(4);
  for(int i = 0; i < 4; i++)
  {
    vector<double> v(2);
    x[i] +=   de; v[1] = value(points, x);
    x[i] -= 2*de; v[0] = value(points, x);
    x[i] +=   de;

    if(i != 1) // zero for d0
      F(i) = (v[1] - v[0])/(2*de);
  }

  TMatrixD J(4,4);
  for(int i = 0; i < 4; i++)
  for(int j = 0; j < 4; j++)
  if(i != 1 && j != 1)
  {
    double v00,v10,v01,v11;

    x[i] +=   de; x[j] +=   de; v11 = value(points, x);
                  x[j] -= 2*de; v10 = value(points, x);
    x[i] -= 2*de;               v00 = value(points, x);
                  x[j] += 2*de; v01 = value(points, x);
    x[i] +=   de; x[j] -=   de;

    J(i,j) = (v11 + v00 - v10 - v01) / (2*de) / (2*de);
  }
  else
    J(i,j) = (i==j ? 1 : 0);

  if(J(3,3) == 0.) J(3,3) = 1.;

  TDecompLU lu; lu.SetMatrix(J);

  if(lu.Decompose())
  {
    TVectorD dx(F); lu.Solve(dx); dx *= -1.;
 
    vector<double> x_new(4);
    for(int i = 0; i < 4; i++) x_new[i] = x[i] + dx(i);
  
    if(value(points, x_new) < value(points, x))
      for(int i = 0; i < 4; i++) x[i] += dx(i);
    else
      for(int i = 0; i < 4; i++) x[i] -= 1e-4*F(i);
  }
  }

  nHitsForSum = 4;
  return sqrt(value(points, x) / min(nHitsForSum, int(points.size())));
}


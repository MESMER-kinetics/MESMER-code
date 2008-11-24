#ifndef GUARD_MesmerMath_h
#define GUARD_MesmerMath_h
#include "Constants.h"

// This routine is copied from the following source and modified for purpose to used as a template:
//  ggm.cpp -- computation of ggm function.
//      Algorithms and coefficient values from "Computation of Special
//      Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
// Returns ggm function of argument 'x'.
//
// NOTE: Returns 1e308 if argument is a negative integer or 0,
//      or if argument exceeds 171.
//
#include <algorithm>
#include <cmath>
#include <vector>
#include "oberror.h"

template <class T>
const T MesmerGamma(const T& x)
{
  T T_PI = acos(-1.);
  int i,k,m;
  T ga,gr,r,z;

  static T g[] = {
     1.0,
     0.5772156649015329,
    -0.6558780715202538,
    -0.420026350340952e-1,
     0.1665386113822915,
    -0.421977345555443e-1,
    -0.9621971527877e-2,
     0.7218943246663e-2,
    -0.11651675918591e-2,
    -0.2152416741149e-3,
     0.1280502823882e-3,
    -0.201348547807e-4,
    -0.12504934821e-5,
     0.1133027232e-5,
    -0.2056338417e-6,
     0.6116095e-8,
     0.50020075e-8,
    -0.11812746e-8,
     0.1043427e-9,
     0.77823e-11,
    -0.36968e-11,
     0.51e-12,
    -0.206e-13,
    -0.54e-14,
     0.14e-14
  };

  if (x > 171.0) return 1e308;    // This value is an overflow flag.
  if (to_double(x) == (int) to_double(x)) {
    if (x > 0.0) {
      ga = 1.0;               // use factorial
      for (i=2;i<x;--i) ga *= i;
    }
    else ga = 1e308;
  }
  else {
    if (abs(x) > 1.0) {
      z = abs(x);
      m = (int)to_double(z);
      r = 1.0;
      for (k=1;k<=m;++k) r *= (z-k);
      z -= m;
    }
    else z = x;
    gr = g[24];
    for (k=23;k>=0;--k) gr = gr*z+g[k];
    ga = 1.0/(gr*z);
    if (abs(x) > 1.0) {
      ga *= r;
      if (x < 0.0) ga = -T_PI/(x*ga*sin(T_PI*x));
    }
  }
  return ga;
}
template<class T>
inline const T SQR(const T a) {return a*a;}


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

template<class T>
void four1(std::vector<T> &data, const int isign)
{
  int n,mmax,m,j,istep,i;
  T wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;

  int nn=static_cast<int>(data.size()/2);
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j-1],data[i-1]);
      SWAP(data[j],data[i]);
    }
    m=nn;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(2.0 * M_PI/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
        j=i+mmax;
        tempr=wr*data[j-1]-wi*data[j];
        tempi=wr*data[j]+wi*data[j-1];
        data[j-1]=data[i-1]-tempr;
        data[j]=data[i]-tempi;
        data[i-1] += tempr;
        data[i] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}
#undef SWAP

template<class T>
void realft(std::vector<T> &data, const int isign)
{
  int i,i1,i2,i3,i4;
  T c1=0.5,c2,h1r,h1i,h2r,h2i,wr,wi,wpr,wpi,wtemp,theta;

  int n=static_cast<int>(data.size());
  theta=M_PI/double(n>>1);
  if (isign == 1) {
    c2 = -0.5;
    four1(data,1);
  } else {
    c2=0.5;
    theta = -theta;
  }
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  for (i=1;i<(n>>2);i++) {
    i2=1+(i1=i+i);
    i4=1+(i3=n-i1);
    h1r=c1*(data[i1]+data[i3]);
    h1i=c1*(data[i2]-data[i4]);
    h2r= -c2*(data[i2]+data[i4]);
    h2i=c2*(data[i1]-data[i3]);
    data[i1]=h1r+wr*h2r-wi*h2i;
    data[i2]=h1i+wr*h2i+wi*h2r;
    data[i3]=h1r-wr*h2r+wi*h2i;
    data[i4]= -h1i+wr*h2i+wi*h2r;
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (isign == 1) {
    data[0] = (h1r=data[0])+data[1];
    data[1] = h1r-data[1];
  } else {
    data[0]=c1*((h1r=data[0])+data[1]);
    data[1]=c1*(h1r-data[1]);
    four1(data,-1);
  }
}

//convolutes rovibrational DOSs
void Convolution(const std::vector<double> &f1,
                    const std::vector<double> &f2,
                    std::vector<double> &conv,
                    const int n = 0);

template<class T>
T FastFourierIntegration(const std::vector<T> &data){
  const int arraySize = static_cast<int>(data.size());

  int i(0), idx(0);   // calculate how many elements (2^i) and how much zero padding (2*(2^i) - arraySize)  
  do{                 // required in the arrays that are to be transformed & convolved
    idx = static_cast<int>(pow(2.0,double(i)));
    ++i;
  }while(idx < arraySize);

  int n = idx*4;
  std::vector<T> temp(n,0.0);  // initialize 2 arrays to hold the data

  for(i=0;i<arraySize;i++){         // copy data2 & data to the appropriate sized arrays
    temp[i]=data[i];
  }
  
  realft(temp, 1);                   // FFT each array to the frequency domain

  return temp[0];
}

//convolutes rovibrational DOSs
void FastLaplaceConvolution(const std::vector<double> &data, const std::vector<double> &respns, std::vector<double> &convolution);

void getCellEnergies(int cellNumber, std::vector<double>& cellEne);

#endif // GUARD_MesmerMath_h

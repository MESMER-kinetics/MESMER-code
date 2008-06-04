#include "fftw++.h"

#ifndef __convolution_h__
#define __convolution_h__ 1

class convolution {
protected:
  unsigned int n;
  unsigned int m;
  unsigned int n2;
  rcfft1d *rc;
  crfft1d *cr;
public:  
  convolution(unsigned int n, unsigned int m) : n(n), m(m), n2(n/2) {
    rc=new rcfft1d(n);
    cr=new crfft1d(n);
  }
  
// Compute H = F (*) G, where F and G are the non-negative Fourier
// components of real functions f and g, respectively. Dealiasing via
// zero-padding is implemented automatically.
//
// Arrays F[n/2+1], g[n/2+1] must be distinct.
// Input F[i] (0 <= i < m), where 3*m <= n, g[i] (0 <= i < n/2).
// Output H[i] = F (*) G  (0 <= i < m), F[i]=f[i], g[i] (0 <= i < n/2).
//
// Array H[n/2+1] can coincide with either F or g, in which case the output H
// subsumes F or g, respectively.

  void fft0(Complex *H, Complex *F, Complex *g) {
    for(unsigned int i=m; i <= n2; i++) F[i]=0.0;
    cr->fft(F);
	
    double ninv=1.0/n;
    for(unsigned int i=0; i < n2; i++)
      H[i]=Complex(F[i].real()*g[i].real()*ninv,F[i].imag()*g[i].imag()*ninv);
	
    rc->fft(H);
  }

// Compute H = F (*) G, where F and G are the non-negative Fourier
// components of real functions f and g, respectively. Dealiasing via
// zero-padding is implemented automatically.
//
// Arrays F[n/2+1], G[n/2+1] must be distinct.
// Input F[i], G[i] (0 <= i < m), where 3*m <= n.
// Output H[i] = F (*) G  (0 <= i < m), F[i]=f[i], G[i]=g[i] (0 <= i < n/2).
//
// Array H[n/2+1] can coincide with either F or G, in which case the output H
// subsumes f or g, respectively.

  void fft(Complex *H, Complex *F, Complex *G) {
    for(unsigned int i=m; i <= n2; i++) G[i]=0.0;
    cr->fft(G);
	
    fft0(H,F,G);
  }	

// Compute H = F (*) G, where F and G contain the non-negative Fourier
// components of real functions f and g, respectively, via direct convolution
// instead of a Fast Fourier Transform technique.
//
// Input F[i], G[i] (0 <= i < m).
// Output H[i] = F (*) G  (0 <= i < m), F and G unchanged.
//
// Array H[m] must be distinct from F[m] and G[m].

  void direct(Complex *H, Complex *F, Complex *G) {
    for(unsigned int i=0; i < m; i++) {
      Complex sum=0.0;
      for(unsigned int j=0; j <= i; j++) sum += F[j]*G[i-j];
      for(unsigned int j=i+1; j < m; j++) sum += F[j]*conj(G[j-i]);
      for(unsigned int j=1; j < m-i; j++) sum += conj(F[j])*G[i+j];
      H[i]=sum;
    }
  }	

};

#endif

#include "fftw++.h"

// Compilation: g++ example0.cc fftw++.cc -lfftw3

using std::cout;

int main()
{
  unsigned int n=4; 
  Complex *f=FFTWComplex(n);
  
  fft1d Forward(n,-1);
  
  for(unsigned int i=0; i < n; i++) f[i]=i;
	
  Forward.fft(f);
	
  for(unsigned int i=0; i < n; i++) cout << f[i] << endl;
  
  FFTWdelete(f);
}

#include "Array.h"
#include "fftw++.h"

// Compilation: g++ example1r.cc fftw++.cc -lfftw3

using std::cout;
using Array::array1;

int main()
{
  unsigned int n=4;
  unsigned int np=n/2+1;
  size_t align=sizeof(Complex);
  
  array1<double> f(n,align);
  array1<Complex> g(np,align);
  
  rcfft1d Forward(f,g);
  crfft1d Backward(g,f);
  
  for(unsigned int i=0; i < n; i++) f[i]=i;
	
  Forward.fft(f,g);
  
  cout << g << endl;
  
  Backward.fftNormalized(g,f);
  
  cout << f << endl;
}

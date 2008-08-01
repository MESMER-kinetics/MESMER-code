#include "MesmerMath.h"
#include "Constants.h"

using namespace std;
using namespace mesmer;
using namespace Constants;

//convolutes DOSs
void Convolution(const std::vector<double> &f1,
                    const std::vector<double> &f2,
                    std::vector<double> &conv,
                    const int n)
{
  std::vector<double>::size_type vSize = f1.size();
  if (n) vSize = static_cast<std::vector<double>::size_type>(n);
  for (std::vector<double>::size_type i = 0; i < vSize; ++i){
    conv[i] = 0.;
    for (std::vector<double>::size_type j = 0; j <= i; ++j){
      conv[i] += f1[j] * f2[i - j];
    }
  }
}

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(std::vector<double> &data, const int isign)
{
  int n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;

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


void realft(std::vector<double> &data, const int isign)
{
  int i,i1,i2,i3,i4;
  double c1=0.5,c2,h1r,h1i,h2r,h2i,wr,wi,wpr,wpi,wtemp,theta;

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



void FastLaplaceConvolution(std::vector<double> &data1, std::vector<double> 
                            &data2, std::vector<double> &convolution)
{ /* this routine takes as input three identically sized vectors, data1, data2, and convolution.
   It takes the FFT of data1 & data2, multiplies them in the frequency domain, and performs an inverse
   FFT of their product to the time domain.  the time convolution is returned in the convolution
   input vector.  note: the zero padding calculated at the beginning of this function ensures that the 
   convolution returned is a Laplace type convolution, i.e.:
                      k=t
                     -----
                      \
   (data1*data2)(t)=   )   data1(t-k) x data2(k)dt   
                      /          
                     -----
                      k=0                                                                             */

  int i, n, idx, no2;
  double tmp;
  int arraySize = static_cast<int>(data1.size()); // check the size of data1 and data2, which should be identical

  i=0;                // calculate how many elements (2^i) and how much zero padding (2*(2^i) - arraySize)  
  do{                 // required in the arrays that are to be transformed & convolved
    idx = static_cast<int>(pow(2.0,double(i)));
    ++i;
  }while(idx < arraySize);

  n = idx*4;

  std::vector<double> temp(n,0.0);  // initialize 2 arrays to hold the data
  std::vector<double> ans(n,0.0);

  for(i=0;i<arraySize;i++){         // copy data2 & data to the appropriate sized arrays
    ans[i]=data1[i];
    temp[i]=data2[i];
  }

  realft(ans, 1);                   // FFT each array to the frequency domain
  realft(temp, 1);

  no2 = n/2;                        // normalization for inverse FFT

  for(i=2;i<n;i+=2){                // perform convolution by multiplying the FFTS in the frequency domain 
    tmp=ans[i];
    ans[i]=(ans[i]*temp[i]-ans[i+1]*temp[i+1])/no2;
    ans[i+1]=(ans[i+1]*temp[i]+tmp*temp[i+1])/no2;
  }
  ans[0]=ans[0]*temp[0]/no2;
  ans[1]=ans[1]*temp[1]/no2;

  realft(ans, -1);                   // inverse FFT back to the time/energy domain

  double convolution_value,initial_value;     // replace those values in the convolution vector that are 
  int replace(0);                             // unreliable because of precision issues
  do{
    convolution_value = 0;
    initial_value = ans[replace];
    for (int j(0); j <= replace; ++j)
      convolution_value += data2[replace-j] * data1[j];
    ans[replace] = convolution_value;
    ++replace;
  }while(convolution_value/initial_value < 0.99 || convolution_value/initial_value > 1.01);

  cinfo << endl << replace << " values in the FFT convolution routine were replaced by standard convolution" << endl;

  for(i=0;i<arraySize;i++)           // copy the inverse FT into the convolution vector for output
    convolution[i]=ans[i];

}

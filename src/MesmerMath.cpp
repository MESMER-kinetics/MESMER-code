#include "MesmerMath.h"

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

void FastLaplaceConvolution(const std::vector<double> &data1, const std::vector<double> &data2, 
                            std::vector<double> &convolution)
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
  const int arraySize = static_cast<int>(data1.size()); // check the size of data1 and data2, which should be identical
  convolution.clear();
  convolution.resize(arraySize,0.0);

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

  int chk(20),replace(0);    // chk specifies how many values are required to demonstrate that
  int ctr(0);                // the FFT convolution gives results within 1% of the slow convolution  
  std::vector<double> conv_value(chk,0.0),initial_value(chk,0.0); 
                             
  do{                                        // replace those values in the convolution vector that are 
    ctr = 0;                                 // unreliable because of precision issues - especially significant
    for(int jj(0); jj<chk; ++jj){            // for convolution of data sets including very small values, e.g.,
      conv_value[jj] = 0;                    // convolution of tunneling probabilities with TS Sum of States
      initial_value[jj] = ans[replace];
      for(int j(0); j <= replace; ++j)
        conv_value[jj] += data2[replace-j] * data1[j];
      ans[replace] = conv_value[jj];
      if(conv_value[jj]/initial_value[jj]<1.01 && conv_value[jj]/initial_value[jj]>0.99 || initial_value[jj]==0)
        ctr += 1;                            // increment ctr by one if it the convergence criteria are specified
      ++replace;                                
    }
  }while(ctr!=chk);

  cinfo << endl << replace << " values in the FFT convolution routine were replaced by standard convolution" << endl;

  for(i=0;i<arraySize;i++)           // copy the inverse FT into the convolution vector for output
    convolution[i]=ans[i];

}


void getCellEnergies(int cellNumber, std::vector<double>& cellEne)
{
  cellEne.clear();
  for (int i = 0 ; i < cellNumber ; ++i ) {
    cellEne.push_back(double(i) + 0.5);
  }
}


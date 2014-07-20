#include "MesmerMath.h"
#include "Constants.h"
#include <stdexcept>

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



void FastLaplaceConvolution(const vector<double> &data1original,
  const vector<double> &data2original,
  vector<double> &convolution)
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

  ErrorContext e(__FUNCTION__);

  // Check the size of data1 and data2, which should be identical.
  const size_t arraySize = data1original.size(); 
  if(data2original.size() != arraySize){
    cerr << "Fast Laplace Convolution: data arrays are not of identical size" << endl;
  }

  // Calculate how many elements (2^i) and how much zero padding (2*(2^i) - arraySize)
  // required in the arrays that are to be transformed & convolved.
  size_t n(2), no2(0);
  do {
    n *= 2 ;
  } while(n < arraySize);
  n  *= 4;
  no2 = n/2;                        // normalization for inverse FFT

  std::vector<double> temp(n,0.0), ans(n,0.0);  // initialize 2 arrays to hold the data
  std::vector<double> data1(arraySize,0.0), data2(arraySize,0.0);
  convolution.clear();
  convolution.resize(arraySize,0.0);

  // Check if there are any zeros at the beginning of the arrays, and count them.
  // If there are leading zero elements, then we can speed up the convolution by
  // removing them, and then placing them back in after we finish the convolution.

  size_t i, NumberZeroElements1(0), NumberZeroElements2(0);
  while(data1original[NumberZeroElements1]==0.0 && NumberZeroElements1 < arraySize) {++NumberZeroElements1;}
  while(data2original[NumberZeroElements2]==0.0 && NumberZeroElements2 < arraySize) {++NumberZeroElements2;}

  if (NumberZeroElements1 == arraySize || NumberZeroElements2 == arraySize ) {
    cinfo << endl << "A null array has been detected in the convolution method." << endl;
  }

  for(i=0;i<(arraySize-NumberZeroElements1);++i)
    ans[i]=data1[i]=data1original[NumberZeroElements1+i];

  for(i=0;i<(arraySize-NumberZeroElements2);++i)
    temp[i]=data2[i]=data2original[NumberZeroElements2+i];

  // keep track of total Number of Zero Elements for shifting following Convolution
  size_t TotalNumberLeadingZeros = NumberZeroElements1 + NumberZeroElements2;

  realft(ans, 1);                   // FFT each array to the frequency domain
  realft(temp, 1);

  for(i=2;i<n;i+=2){                // perform convolution by multiplying the FFTS in the frequency domain
    double tmp=ans[i];
    ans[i]=(ans[i]*temp[i]-ans[i+1]*temp[i+1])/no2;
    ans[i+1]=(ans[i+1]*temp[i]+tmp*temp[i+1])/no2;
  }
  ans[0]=ans[0]*temp[0]/no2;
  ans[1]=ans[1]*temp[1]/no2;

  realft(ans, -1);                   // inverse FFT back to the time/energy domain


  // Replace those values in the convolution vector that are unreliable because 
  // of precision issues - especially significant for convolution of data sets
  // including very small values, e.g., convolution of tunneling probabilities
  // with TS Sum of States.
  // chk specifies how many values are required to demonstrate that the FFT
  // convolution gives results within 1% of the slow convolution
  int chk(20), replace(0), ctr(0);
  const double tolerance(0.01) ;
  do {
    ctr = 0;
    for(int jj(0); jj<chk; ++jj){
      double conv_value = 0.0;
      double initial_value = ans[replace];
      for(int j(0); j <= replace; ++j)
        conv_value += data2[replace-j] * data1[j];
      ans[replace] = conv_value ;
      double tmp = (initial_value != 0.0) ? conv_value/initial_value - 1.0 : 0.0 ;
      if(fabs(tmp) < tolerance || initial_value==0)
        ctr++ ;
      ++replace;
    }
  } while(ctr!=chk);

  cinfo << replace << " values in the FFT convolution routine were replaced by standard convolution" << once << endl;

  // Copy the inverse FT into the convolution vector for output.
  // If there were Leading Zeros, put them back in.
  for(i=0;i<(convolution.size() - TotalNumberLeadingZeros);++i)
    convolution[TotalNumberLeadingZeros+i]=ans[i];

}


void getCellEnergies(const size_t &cellNumber, const double &cellSize, std::vector<double>& cellEne)
{
  cellEne.clear();
  for (size_t i(0) ; i < cellNumber ; ++i ) {
    cellEne.push_back(cellSize*(double(i) + 0.5));
  }
}

bool convertToFourierCoefficients(const size_t expansion, vector<double>& ak, vector<double>& bk, double& a0, vector<double> pesEnes){
  a0 = 0.0;
  size_t n = pesEnes.size();
  for(size_t i(0); i < n; ++i) {
    a0 += pesEnes[i] / double(n);
  }

  ak.resize(expansion, 0.0);
  bk.resize(expansion, 0.0);

  for(size_t k(0); k < expansion; ++k) {
    for(size_t i(0); i < n; ++i) {
      double theta = (2.0 * double(i) / double(n) - 0.95) * M_PI;
      double nTheta = (double(k)+1.0) * theta;
      ak[k] += 2.0 * pesEnes[i] * cos(nTheta) / double(n);
      bk[k] += 2.0 * pesEnes[i] * sin(nTheta) / double(n);
    }
  }

  return true;
}

double getEnergyFromFourierCoefficients(double theta, const vector<double> ak, const vector<double> bk, const double a0){
  size_t expansionak = ak.size();
  size_t expansionbk = bk.size();
  double returnValue = a0;
  for(size_t k(0); k < expansionak; ++k) {
    returnValue += ak[k] * cos((double(k)+1.0) * theta);
  }
  for(size_t k(0); k < expansionbk; ++k) {
    returnValue += bk[k] * sin((double(k)+1.0) * theta);
  }
  return returnValue;
}

void airy(double x, double& ai, double& aip, double& bi, double& bip)

{
  double z;
  double zz;
  double t;
  double f;
  double g;
  double uf;
  double ug;
  double k;
  double zeta;
  double theta;
  int domflg;
  double c1;
  double c2;
  double sqrt3;
  double sqpii;
  double afn;
  double afd;
  double agn;
  double agd;
  double apfn;
  double apfd;
  double apgn;
  double apgd;
  double an;
  double ad;
  double apn;
  double apd;
  double bn16;
  double bd16;
  double bppn;
  double bppd;

  sqpii = 5.64189583547756286948E-1;
  c1 = 0.35502805388781723926;
  c2 = 0.258819403792806798405;
  sqrt3 = 1.732050807568877293527;
  domflg = 0;

  if( x>25.77 )
  {
    ai = 0;
    aip = 0;
    bi = 1.0e+300;
    bip = 1.0e+300;
    return;
  }
  if( x<-2.09 )
  {
    domflg = 15;
    t = sqrt(-x);
    zeta = -2.0*x*t/3.0;
    t = sqrt(t);
    k = sqpii/t;
    z = 1.0/zeta;
    zz = z*z;
    afn = -1.31696323418331795333E-1;
    afn = afn*zz-6.26456544431912369773E-1;
    afn = afn*zz-6.93158036036933542233E-1;
    afn = afn*zz-2.79779981545119124951E-1;
    afn = afn*zz-4.91900132609500318020E-2;
    afn = afn*zz-4.06265923594885404393E-3;
    afn = afn*zz-1.59276496239262096340E-4;
    afn = afn*zz-2.77649108155232920844E-6;
    afn = afn*zz-1.67787698489114633780E-8;
    afd = 1.00000000000000000000E0;
    afd = afd*zz+1.33560420706553243746E1;
    afd = afd*zz+3.26825032795224613948E1;
    afd = afd*zz+2.67367040941499554804E1;
    afd = afd*zz+9.18707402907259625840E0;
    afd = afd*zz+1.47529146771666414581E0;
    afd = afd*zz+1.15687173795188044134E-1;
    afd = afd*zz+4.40291641615211203805E-3;
    afd = afd*zz+7.54720348287414296618E-5;
    afd = afd*zz+4.51850092970580378464E-7;
    uf = 1.0+zz*afn/afd;
    agn = 1.97339932091685679179E-2;
    agn = agn*zz+3.91103029615688277255E-1;
    agn = agn*zz+1.06579897599595591108E0;
    agn = agn*zz+9.39169229816650230044E-1;
    agn = agn*zz+3.51465656105547619242E-1;
    agn = agn*zz+6.33888919628925490927E-2;
    agn = agn*zz+5.85804113048388458567E-3;
    agn = agn*zz+2.82851600836737019778E-4;
    agn = agn*zz+6.98793669997260967291E-6;
    agn = agn*zz+8.11789239554389293311E-8;
    agn = agn*zz+3.41551784765923618484E-10;
    agd = 1.00000000000000000000E0;
    agd = agd*zz+9.30892908077441974853E0;
    agd = agd*zz+1.98352928718312140417E1;
    agd = agd*zz+1.55646628932864612953E1;
    agd = agd*zz+5.47686069422975497931E0;
    agd = agd*zz+9.54293611618961883998E-1;
    agd = agd*zz+8.64580826352392193095E-2;
    agd = agd*zz+4.12656523824222607191E-3;
    agd = agd*zz+1.01259085116509135510E-4;
    agd = agd*zz+1.17166733214413521882E-6;
    agd = agd*zz+4.91834570062930015649E-9;
    ug = z*agn/agd;
    theta = zeta+0.25*M_PI;
    f = sin(theta);
    g = cos(theta);
    ai = k*(f*uf-g*ug);
    bi = k*(g*uf+f*ug);
    apfn = 1.85365624022535566142E-1;
    apfn = apfn*zz+8.86712188052584095637E-1;
    apfn = apfn*zz+9.87391981747398547272E-1;
    apfn = apfn*zz+4.01241082318003734092E-1;
    apfn = apfn*zz+7.10304926289631174579E-2;
    apfn = apfn*zz+5.90618657995661810071E-3;
    apfn = apfn*zz+2.33051409401776799569E-4;
    apfn = apfn*zz+4.08718778289035454598E-6;
    apfn = apfn*zz+2.48379932900442457853E-8;
    apfd = 1.00000000000000000000E0;
    apfd = apfd*zz+1.47345854687502542552E1;
    apfd = apfd*zz+3.75423933435489594466E1;
    apfd = apfd*zz+3.14657751203046424330E1;
    apfd = apfd*zz+1.09969125207298778536E1;
    apfd = apfd*zz+1.78885054766999417817E0;
    apfd = apfd*zz+1.41733275753662636873E-1;
    apfd = apfd*zz+5.44066067017226003627E-3;
    apfd = apfd*zz+9.39421290654511171663E-5;
    apfd = apfd*zz+5.65978713036027009243E-7;
    uf = 1.0+zz*apfn/apfd;
    apgn = -3.55615429033082288335E-2;
    apgn = apgn*zz-6.37311518129435504426E-1;
    apgn = apgn*zz-1.70856738884312371053E0;
    apgn = apgn*zz-1.50221872117316635393E0;
    apgn = apgn*zz-5.63606665822102676611E-1;
    apgn = apgn*zz-1.02101031120216891789E-1;
    apgn = apgn*zz-9.48396695961445269093E-3;
    apgn = apgn*zz-4.60325307486780994357E-4;
    apgn = apgn*zz-1.14300836484517375919E-5;
    apgn = apgn*zz-1.33415518685547420648E-7;
    apgn = apgn*zz-5.63803833958893494476E-10;
    apgd = 1.00000000000000000000E0;
    apgd = apgd*zz+9.85865801696130355144E0;
    apgd = apgd*zz+2.16401867356585941885E1;
    apgd = apgd*zz+1.73130776389749389525E1;
    apgd = apgd*zz+6.17872175280828766327E0;
    apgd = apgd*zz+1.08848694396321495475E0;
    apgd = apgd*zz+9.95005543440888479402E-2;
    apgd = apgd*zz+4.78468199683886610842E-3;
    apgd = apgd*zz+1.18159633322838625562E-4;
    apgd = apgd*zz+1.37480673554219441465E-6;
    apgd = apgd*zz+5.79912514929147598821E-9;
    ug = z*apgn/apgd;
    k = sqpii*t;
    aip = -k*(g*uf+f*ug);
    bip = k*(f*uf-g*ug);
    return;
  }
  if( x>=2.09 )
  {
    domflg = 5;
    t = sqrt(x);
    zeta = 2.0*x*t/3.0;
    g = exp(zeta);
    t = sqrt(t);
    k = 2.0*t*g;
    z = 1.0/zeta;
    an = 3.46538101525629032477E-1;
    an = an*z+1.20075952739645805542E1;
    an = an*z+7.62796053615234516538E1;
    an = an*z+1.68089224934630576269E2;
    an = an*z+1.59756391350164413639E2;
    an = an*z+7.05360906840444183113E1;
    an = an*z+1.40264691163389668864E1;
    an = an*z+9.99999999999999995305E-1;
    ad = 5.67594532638770212846E-1;
    ad = ad*z+1.47562562584847203173E1;
    ad = ad*z+8.45138970141474626562E1;
    ad = ad*z+1.77318088145400459522E2;
    ad = ad*z+1.64234692871529701831E2;
    ad = ad*z+7.14778400825575695274E1;
    ad = ad*z+1.40959135607834029598E1;
    ad = ad*z+1.00000000000000000470E0;
    f = an/ad;
    ai = sqpii*f/k;
    k = -0.5*sqpii*t/g;
    apn = 6.13759184814035759225E-1;
    apn = apn*z+1.47454670787755323881E1;
    apn = apn*z+8.20584123476060982430E1;
    apn = apn*z+1.71184781360976385540E2;
    apn = apn*z+1.59317847137141783523E2;
    apn = apn*z+6.99778599330103016170E1;
    apn = apn*z+1.39470856980481566958E1;
    apn = apn*z+1.00000000000000000550E0;
    apd = 3.34203677749736953049E-1;
    apd = apd*z+1.11810297306158156705E1;
    apd = apd*z+7.11727352147859965283E1;
    apd = apd*z+1.58778084372838313640E2;
    apd = apd*z+1.53206427475809220834E2;
    apd = apd*z+6.86752304592780337944E1;
    apd = apd*z+1.38498634758259442477E1;
    apd = apd*z+9.99999999999999994502E-1;
    f = apn/apd;
    aip = f*k;
    if( x>8.3203353 )
    {
      bn16 = -2.53240795869364152689E-1;
      bn16 = bn16*z+5.75285167332467384228E-1;
      bn16 = bn16*z-3.29907036873225371650E-1;
      bn16 = bn16*z+6.44404068948199951727E-2;
      bn16 = bn16*z-3.82519546641336734394E-3;
      bd16 = 1.00000000000000000000E0;
      bd16 = bd16*z-7.15685095054035237902E0;
      bd16 = bd16*z+1.06039580715664694291E1;
      bd16 = bd16*z-5.23246636471251500874E0;
      bd16 = bd16*z+9.57395864378383833152E-1;
      bd16 = bd16*z-5.50828147163549611107E-2;
      f = z*bn16/bd16;
      k = sqpii*g;
      bi = k*(1.0+f)/t;
      bppn = 4.65461162774651610328E-1;
      bppn = bppn*z-1.08992173800493920734E0;
      bppn = bppn*z+6.38800117371827987759E-1;
      bppn = bppn*z-1.26844349553102907034E-1;
      bppn = bppn*z+7.62487844342109852105E-3;
      bppd = 1.00000000000000000000E0;
      bppd = bppd*z-8.70622787633159124240E0;
      bppd = bppd*z+1.38993162704553213172E1;
      bppd = bppd*z-7.14116144616431159572E0;
      bppd = bppd*z+1.34008595960680518666E0;
      bppd = bppd*z-7.84273211323341930448E-2;
      f = z*bppn/bppd;
      bip = k*t*(1.0+f);
      return;
    }
  }
  f = 1.0;
  g = x;
  t = 1.0;
  uf = 1.0;
  ug = x;
  k = 1.0;
  z = x*x*x;
  while(t>5.0e-16)
  {
    uf = uf*z;
    k = k+1.0;
    uf = uf/k;
    ug = ug*z;
    k = k+1.0;
    ug = ug/k;
    uf = uf/k;
    f = f+uf;
    k = k+1.0;
    ug = ug/k;
    g = g+ug;
    t = fabs(uf/f);
  }
  uf = c1*f;
  ug = c2*g;
  if( domflg%2==0 )
  {
    ai = uf-ug;
  }
  if( domflg/2%2==0 )
  {
    bi = sqrt3*(uf+ug);
  }
  k = 4.0;
  uf = x*x/2.0;
  ug = z/3.0;
  f = uf;
  g = 1.0+ug;
  uf = uf/3.0;
  t = 1.0;
  while(t>5.0e-16)
  {
    uf = uf*z;
    ug = ug/k;
    k = k+1.0;
    ug = ug*z;
    uf = uf/k;
    f = f+uf;
    k = k+1.0;
    ug = ug/k;
    uf = uf/k;
    g = g+ug;
    k = k+1.0;
    t = fabs(ug/g);
  }
  uf = c1*f;
  ug = c2*g;
  if( domflg/4%2==0 )
  {
    aip = uf-ug;
  }
  if( domflg/8%2==0 )
  {
    bip = sqrt3*(uf+ug);
  }
}

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void airy2(const double x, double &ai)
{
  const double PI=3.141592653589793238, ONOVRT=0.577350269189626;
  const double THIRD=(1.0/3.0), TWOTHR=2.0*THIRD;
  double absx,ri,rip,rj,rjp,rk,rkp,rootx,ry,ryp,z;

  absx=fabs(x);
  rootx=sqrt(absx);
  z=TWOTHR*absx*rootx;
  if (x > 0.0) {
    bessik(z,THIRD,ri,rk,rip,rkp);
    ai=rootx*ONOVRT*rk/PI;
  } else if (x < 0.0) {
    bessjy(z,THIRD,rj,ry,rjp,ryp);
    ai=0.5*rootx*(rj-ONOVRT*ry);
  } else {
    ai=0.355028053887817;
  }
}

void bessik(const double x, const double xnu, double &ri, double &rk, double &rip, double &rkp)
{
  const int MAXIT=10000;
  const double EPS=numeric_limits<double>::epsilon();
  const double FPMIN=numeric_limits<double>::min()/EPS;
  const double XMIN=2.0, PI=3.141592653589793;
  double a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,gam1,gam2,
    gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,ripl,
    ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2;
  int i,l,nl;

  if (x <= 0.0 || xnu < 0.0) nrerror("bad arguments in bessik");
  nl=int(xnu+0.5);
  xmu=xnu-nl;
  xmu2=xmu*xmu;
  xi=1.0/x;
  xi2=2.0*xi;
  h=xnu*xi;
  if (h < FPMIN) h=FPMIN;
  b=xi2*xnu;
  d=0.0;
  c=h;
  for (i=0;i<MAXIT;i++) {
    b += xi2;
    d=1.0/(b+d);
    c=b+1.0/c;
    del=c*d;
    h=del*h;
    if (fabs(del-1.0) <= EPS) break;
  }
  if (i >= MAXIT)
    nrerror("x too large in bessik; try asymptotic expansion");
  ril=FPMIN;
  ripl=h*ril;
  ril1=ril;
  rip1=ripl;
  fact=xnu*xi;
  for (l=nl-1;l >= 0;l--) {
    ritemp=fact*ril+ripl;
    fact -= xi;
    ripl=fact*ritemp+ril;
    ril=ritemp;
  }
  f=ripl/ril;
  if (x < XMIN) {
    x2=0.5*x;
    pimu=PI*xmu;
    fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
    d = -log(x2);
    e=xmu*d;
    fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
    beschb(xmu,gam1,gam2,gampl,gammi);
    ff=fact*(gam1*cosh(e)+gam2*fact2*d);
    sum=ff;
    e=exp(e);
    p=0.5*e/gampl;
    q=0.5/(e*gammi);
    c=1.0;
    d=x2*x2;
    sum1=p;
    for (i=1;i<=MAXIT;i++) {
      ff=(i*ff+p+q)/(i*i-xmu2);
      c *= (d/i);
      p /= (i-xmu);
      q /= (i+xmu);
      del=c*ff;
      sum += del;
      del1=c*(p-i*ff);
      sum1 += del1;
      if (fabs(del) < fabs(sum)*EPS) break;
    }
    if (i > MAXIT) nrerror("bessk series failed to converge");
    rkmu=sum;
    rk1=sum1*xi2;
  } else {
    b=2.0*(1.0+x);
    d=1.0/b;
    h=delh=d;
    q1=0.0;
    q2=1.0;
    a1=0.25-xmu2;
    q=c=a1;
    a = -a1;
    s=1.0+q*delh;
    for (i=1;i<MAXIT;i++) {
      a -= 2*i;
      c = -a*c/(i+1.0);
      qnew=(q1-b*q2)/a;
      q1=q2;
      q2=qnew;
      q += c*qnew;
      b += 2.0;
      d=1.0/(b+a*d);
      delh=(b*d-1.0)*delh;
      h += delh;
      dels=q*delh;
      s += dels;
      if (fabs(dels/s) <= EPS) break;
    }
    if (i >= MAXIT) nrerror("bessik: failure to converge in cf2");
    h=a1*h;
    rkmu=sqrt(PI/(2.0*x))*exp(-x)/s;
    rk1=rkmu*(xmu+x+0.5-h)*xi;
  }
  rkmup=xmu*xi*rkmu-rk1;
  rimu=xi/(f*rkmu-rkmup);
  ri=(rimu*ril1)/ril;
  rip=(rimu*rip1)/ril;
  for (i=1;i <= nl;i++) {
    rktemp=(xmu+i)*xi2*rk1+rkmu;
    rkmu=rk1;
    rk1=rktemp;
  }
  rk=rkmu;
  rkp=xnu*xi*rkmu-rk1;
}

void bessjy(const double x, const double xnu, double &rj, double &ry, double &rjp, double &ryp)
{
  const int MAXIT=10000;
  const double EPS=numeric_limits<double>::epsilon();
  const double FPMIN=numeric_limits<double>::min()/EPS;
  const double XMIN=2.0, PI=3.141592653589793;
  double a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e,f,fact,fact2,
    fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q,r,rjl,
    rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,sum1,
    temp,w,x2,xi,xi2,xmu,xmu2;
  int i,isign,l,nl;

  if (x <= 0.0 || xnu < 0.0)
    nrerror("bad arguments in bessjy");
  nl=(x < XMIN ? int(xnu+0.5) : max(0,int(xnu-x+1.5)));
  xmu=xnu-nl;
  xmu2=xmu*xmu;
  xi=1.0/x;
  xi2=2.0*xi;
  w=xi2/PI;
  isign=1;
  h=xnu*xi;
  if (h < FPMIN) h=FPMIN;
  b=xi2*xnu;
  d=0.0;
  c=h;
  for (i=0;i<MAXIT;i++) {
    b += xi2;
    d=b-d;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=b-1.0/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=c*d;
    h=del*h;
    if (d < 0.0) isign = -isign;
    if (fabs(del-1.0) <= EPS) break;
  }
  if (i >= MAXIT)
    nrerror("x too large in bessjy; try asymptotic expansion");
  rjl=isign*FPMIN;
  rjpl=h*rjl;
  rjl1=rjl;
  rjp1=rjpl;
  fact=xnu*xi;
  for (l=nl-1;l>=0;l--) {
    rjtemp=fact*rjl+rjpl;
    fact -= xi;
    rjpl=fact*rjtemp-rjl;
    rjl=rjtemp;
  }
  if (rjl == 0.0) rjl=EPS;
  f=rjpl/rjl;
  if (x < XMIN) {
    x2=0.5*x;
    pimu=PI*xmu;
    fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
    d = -log(x2);
    e=xmu*d;
    fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
    beschb(xmu,gam1,gam2,gampl,gammi);
    ff=2.0/PI*fact*(gam1*cosh(e)+gam2*fact2*d);
    e=exp(e);
    p=e/(gampl*PI);
    q=1.0/(e*PI*gammi);
    pimu2=0.5*pimu;
    fact3 = (fabs(pimu2) < EPS ? 1.0 : sin(pimu2)/pimu2);
    r=PI*pimu2*fact3*fact3;
    c=1.0;
    d = -x2*x2;
    sum=ff+r*q;
    sum1=p;
    for (i=1;i<=MAXIT;i++) {
      ff=(i*ff+p+q)/(i*i-xmu2);
      c *= (d/i);
      p /= (i-xmu);
      q /= (i+xmu);
      del=c*(ff+r*q);
      sum += del;
      del1=c*p-i*del;
      sum1 += del1;
      if (fabs(del) < (1.0+fabs(sum))*EPS) break;
    }
    if (i > MAXIT)
      nrerror("bessy series failed to converge");
    rymu = -sum;
    ry1 = -sum1*xi2;
    rymup=xmu*xi*rymu-ry1;
    rjmu=w/(rymup-f*rymu);
  } else {
    a=0.25-xmu2;
    p = -0.5*xi;
    q=1.0;
    br=2.0*x;
    bi=2.0;
    fact=a*xi/(p*p+q*q);
    cr=br+q*fact;
    ci=bi+p*fact;
    den=br*br+bi*bi;
    dr=br/den;
    di = -bi/den;
    dlr=cr*dr-ci*di;
    dli=cr*di+ci*dr;
    temp=p*dlr-q*dli;
    q=p*dli+q*dlr;
    p=temp;
    for (i=1;i<MAXIT;i++) {
      a += 2*i;
      bi += 2.0;
      dr=a*dr+br;
      di=a*di+bi;
      if (fabs(dr)+fabs(di) < FPMIN) dr=FPMIN;
      fact=a/(cr*cr+ci*ci);
      cr=br+cr*fact;
      ci=bi-ci*fact;
      if (fabs(cr)+fabs(ci) < FPMIN) cr=FPMIN;
      den=dr*dr+di*di;
      dr /= den;
      di /= -den;
      dlr=cr*dr-ci*di;
      dli=cr*di+ci*dr;
      temp=p*dlr-q*dli;
      q=p*dli+q*dlr;
      p=temp;
      if (fabs(dlr-1.0)+fabs(dli) <= EPS) break;
    }
    if (i >= MAXIT) nrerror("cf2 failed in bessjy");
    gam=(p-f)/q;
    rjmu=sqrt(w/((p-f)*gam+q));
    rjmu=SIGN(rjmu,rjl);
    rymu=rjmu*gam;
    rymup=rymu*(p+q/gam);
    ry1=xmu*xi*rymu-rymup;
  }
  fact=rjmu/rjl;
  rj=rjl1*fact;
  rjp=rjp1*fact;
  for (i=1;i<=nl;i++) {
    rytemp=(xmu+i)*xi2*ry1-rymu;
    rymu=ry1;
    ry1=rytemp;
  }
  ry=rymu;
  ryp=xnu*xi*rymu-ry1;
}

void beschb(const double x, double &gam1, double &gam2, double &gampl, double &gammi)
{
  const int NUSE1=7, NUSE2=8;
  static const double c1_d[] = {
    -1.142022680371168e0,6.5165112670737e-3,
    3.087090173086e-4,-3.4706269649e-6,6.9437664e-9,
    3.67795e-11,-1.356e-13};
    static const double c2_d[] = {
      1.843740587300905e0,-7.68528408447867e-2,
      1.2719271366546e-3,-4.9717367042e-6,-3.31261198e-8,
      2.423096e-10,-1.702e-13,-1.49e-15};

      double xx;

      static std::vector<double> c1( c1_d,c1_d + sizeof( c1_d ) / sizeof( c1_d[0] ) );
      static std::vector<double> c2( c2_d,c2_d + sizeof( c2_d ) / sizeof( c2_d[0] ) );

      xx=8.0*x*x-1.0;
      gam1=chebev(-1.0,1.0,c1,NUSE1,xx);
      gam2=chebev(-1.0,1.0,c2,NUSE2,xx);
      gampl= gam2-x*gam1;
      gammi= gam2+x*gam1;
}

double chebev(const double a, const double b, std::vector<double> &c, const int m, const double x)
{
  double d=0.0,dd=0.0,sv,y,y2;
  int j;

  if ((x-a)*(x-b) > 0.0)
    nrerror("x not in range in routine chebev");
  y2=2.0*(y=(2.0*x-a-b)/(b-a));
  for (j=m-1;j>0;j--) {
    sv=d;
    d=y2*d-dd+c[j];
    dd=sv;
  }
  return y*d-dd+0.5*c[0];
}

void nrerror(std::string message)
{
  cerr << message << endl;
}

#undef SIGN

//
// Calculation of the modifed Bessel function, Io(x), for real x.
//
double ModifiedBessalFunction(const double x)  
{
  static double p1(1.0),          p2(3.5156229),     p3(3.0899424) ;
  static double p4(1.2067492),    p5(0.2659732),     p6(0.360768e-1) ;
  static double p7(0.45813e-2) ;

  static double q1(0.39894228),   q2(0.1328592e-1),  q3(0.225319e-2) ;
  static double q4(-0.157565e-2), q5(0.916281e-2),   q6(-0.2057706e-1) ;
  static double q7(0.2635537e-1), q8(-0.1647633e-1), q9(0.392377e-2) ;

  //  Bessel index

  if (abs(x) < 3.75) {

    double y1 = (x/3.75) ;
    y1 *= y1 ;

    return (p1+y1* (p2+y1* (p3+y1* (p4+y1* (p5+y1* (p6+ y1*p7)))))) ;

  } else {

    double ax = abs(x) ;
    double y1 = 3.75/ax ;
    double b1 = (exp(ax)/sqrt(ax)) ;
    double fd = (q6+y1* (q7+y1* (q8+y1*q9))) ;
    double b2 = (q1+y1* (q2+y1* (q3+y1* (q4+y1* (q5+y1*fd))))) ;

    return b1*b2 ;

  }
}

//
// Calculation of the Chi^2 probability function.
//
double ChiSquaredPrbFn(double x, double a) {

  if (x < 0.0 || a < 0.0)
    throw(std::runtime_error("Incorrect parameters to __FUNCTION__.")) ;

  static const size_t maxItr(200) ; 
  static const double eps(2.22e-16) ;
  static const double fpMin(2.23e-308) ;

  double chiSquaredPrb(0.0) ;
  if (x < (a + 1.0)) {
    double ap(a) ;
    double sum(1.0/a) ;
    double del(sum) ; 
    size_t i(0) ;
    for ( ; i < maxItr && (abs(del) >= abs(sum)*eps) ; i++) {
      ap  += 1.0 ;
      del *= x/ap ;
      sum += del ;
    }
    if (i >= maxItr)
      cinfo << "Warning: Maximum iterations exceed in calculation of Gamma function." << endl ;
    chiSquaredPrb = 1.0 - sum*exp(-x + a*log(x))/MesmerGamma(a) ;
  } else {
    double b(x + 1.0 - a) ; 
    double c(1.0/fpMin) ;
    double d(1.0/b) ;
    double h(d) ;
    double del(0.0) ;
    size_t i(1) ;
    for ( ; i < maxItr && (abs(del - 1.0) >= eps) ; i++) {
      double an = -double(i)*(double(i) - a) ;
      b += 2.0 ; 
      d = an*d + b ;
      if (abs(d) < fpMin) d = fpMin ;
      c = b + an/c ;
      if (abs(c) < fpMin) c = fpMin ;
      d = 1.0/d ;
      del = d*c ;
      h *= del ;
    }
    if (i >= maxItr)
      cinfo << "Warning: Maximum iterations exceed in calculation of Gamma function." << endl ;
    chiSquaredPrb = h*exp(-x + a*log(x))/MesmerGamma(a) ;
  }

  return chiSquaredPrb ;

}

// Calculate the Cosine Fourier coefficients for general data.
void FourierCosCoeffs(vector<double> &angle, vector<double> &value, vector<double> &cosCoeff, size_t expansion)
{
  size_t ndata = angle.size() ;
  cosCoeff.clear() ;

  for(size_t k(0); k < expansion ; ++k) {
    double sum(0.0) ;
    for(size_t i(0); i < ndata; ++i) {
      double nTheta = double(k) * angle[i];
      sum += value[i] * cos(nTheta);
    }
    cosCoeff.push_back(2.0*sum/double(ndata)) ;
  }
  cosCoeff[0] /= 2.0 ;
}

// Calculate the Sine Fourier coefficients for general data.
void FourierSinCoeffs(vector<double> &angle, vector<double> &value, vector<double> &sinCoeff, size_t expansion)
{
  size_t ndata = angle.size() ;
  sinCoeff.clear() ;

  for(size_t k(0); k < expansion; ++k) {
    double sum(0.0) ;
    for(size_t i(0); i < ndata; ++i) {
      double nTheta = double(k) * angle[i];
      sum += value[i] * sin(nTheta);
    }
    sinCoeff.push_back(2.0*sum/double(ndata)) ;
  }
  sinCoeff[0] = 0.0 ;

  return ;
}


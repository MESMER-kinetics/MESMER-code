//-------------------------------------------------------------------------------------------
//
// AnalyticalRepresentation.cpp
//
// Author: Struan Robertson
// Date:   21/Jul/2013
//
// This class implements the methods to that determine and analytical representation of a 
// master equation representation.
//
//-------------------------------------------------------------------------------------------

#include <fstream>
#include <iomanip>

#include "../System.h"
#include "../calcmethod.h"
#include "../dMatrix.h"
#include "FittingUtils.h"

namespace mesmer
{
  class AnalyticalRepresentation : public CalcMethod, private FittingUtils
  {
  public:

    AnalyticalRepresentation(const char* id) : FittingUtils(), 
      m_id(id),
      m_format(),
      m_reactionRef1(),
      m_reactionRef2(),
      m_TMax(0.0),
      m_TMin(0.0),
      m_CMax(0.0),
      m_CMin(0.0),
      m_NTpt(0),
      m_NCpt(0),
      RpTMin(0),
      RpTMax(0),
      lgCMin(0),
      lgCMax(0),
      m_ExpanSizeT(0),
      m_ExpanSizeC(0)
    { Register(); }

    virtual ~AnalyticalRepresentation() {}
    virtual const char* getID()  { return m_id; }

    //Read in data for this method from XML
    virtual bool ParseData(PersistPtr pp);

    //Function to do the work
    virtual bool DoCalculation(System* pSys);

  protected:

  private:

    //Funtion to calculation a chebyshev polynomial of a given order
    virtual double Cheb_poly(int order, double argument);

    //Function to give constant coefficient of chebyshev expansion
    virtual double Coefficient(int i, int j);

    // Convererts the Chebyshev gridpoints to back to Temperature/Concetrations
    virtual vector<double> Transform(vector<double> gridpoints, double max, double min );

    const char* m_id;
    string m_format ;
    vector<string> m_reactionRef1 ;
    vector<string> m_reactionRef2 ;

    double m_TMax ;
    double m_TMin ;
    double m_CMax ;
    double m_CMin ;

    double RpTMin ;
    double RpTMax ;
    double lgCMin ;
    double lgCMax ;

    int    m_NTpt ;
    int    m_NCpt ;

    int m_ExpanSizeT;
    int m_ExpanSizeC;

  };

  ////////////////////////////////////////////////
  //Global instance
  AnalyticalRepresentation theAnalyticalRepresentation("analyticalRepresentation");
  ///////////////////////////////////////////////

  bool AnalyticalRepresentation::ParseData(PersistPtr pp)
  {
    /* Typical data
    <me:control>
    ...
    <me:calcMethod xsi:type="me:analyticalRepresentation">
    <me:fittingTolerance>0.1</me:fittingTolerance>
    <me:fittingIterations>5</me:fittingIterations>
    </me:calcMethod>
    ...
    </me:control>
    */
    //Read in fitting parameters, or use values from defaults.xml.

    m_NTpt = pp->XmlReadInteger("me:chebNumTemp");
    m_NCpt = pp->XmlReadInteger("me:chebNumConc");

    m_TMax = pp->XmlReadDouble("me:chebMaxTemp");
    m_TMin = pp->XmlReadDouble("me:chebMinTemp");
    m_CMax = pp->XmlReadDouble("me:chebMaxConc");
    m_CMin = pp->XmlReadDouble("me:chebMinConc");

    m_ExpanSizeT = pp->XmlReadInteger("me:chebTExSize");
    m_ExpanSizeC = pp->XmlReadInteger("me:chebPExSize");

    PersistPtr ppARepRef = pp->XmlMoveTo("me:analyticalRepRef");
    ppARepRef = ppARepRef->XmlMoveTo("me:reaction");
    while(ppARepRef){
      m_reactionRef1.push_back(ppARepRef->XmlReadValue("ref1"));
      m_reactionRef2.push_back(ppARepRef->XmlReadValue("ref2"));
      ppARepRef = ppARepRef->XmlMoveTo("me:reaction");
    }


    RpTMin = 1.0 / m_TMin ;
    RpTMax = 1.0 / m_TMax ;
    lgCMin = log10(m_CMin);
    lgCMax = log10(m_CMax);


    return true;
  }

  bool AnalyticalRepresentation::DoCalculation(System* pSys)
  {
    //Do not output all the intermediate results to XML
    pSys->m_Flags.overwriteXmlAnalysis = true;

    // Use the same grain numbers for for all calcuations regardless of 
    // temperature (i.e. reduce the number of times micro-rates are caluclated).
    pSys->m_Flags.useTheSameCellNumber = true;

    // Warnings and less not sent to console.
    ChangeErrorLevel e(obError); 

    // First gets some points. NEED TO GET CHEBYSHEV GRID POINTS AND TRANSFORM BACK TO TP CONDITION SET

    vector<double> TGrid(m_NTpt);
    for (int i = 0; i< m_NTpt; ++i){
      TGrid[i] = cos(((2.0*(i+1.0))-1.0)  / (2.0 * double(m_NTpt)) *M_PI) ;
    }

    vector<double> CGrid(m_NCpt);
    for (int i = 0; i< m_NCpt; ++i){
      CGrid[i] = cos(((2.0*(i+1.0))-1.0)  / (2.0 * double(m_NCpt)) *M_PI);
    }


    vector<double> Temperature=Transform(TGrid, RpTMin, RpTMax);
    for (int i = 0; i< m_NTpt; ++i){
      Temperature[i] = 1/Temperature[i];
    }

    vector<double> Concentration=Transform(CGrid, lgCMin, lgCMax);
    for (int i = 0; i< m_NCpt; ++i){
      Concentration[i] = pow(10,(Concentration[i]));
    }
    vector<double> RateCoefficients(m_reactionRef1.size());

    // Get rate coefficients
    vector<double> v(m_reactionRef1.size());
    vector<vector<double> >v2(m_NCpt,v);
    vector<vector<vector<double> > > PTgrid(m_NTpt, v2) ;
    for (int m=0; m < m_NTpt; ++m ){
      for(int n=0; n < m_NCpt; ++n){
        for(size_t k=0; k < m_reactionRef1.size(); ++k){
          pSys->calculate(Temperature[m], Concentration[n], m_reactionRef1, m_reactionRef2, RateCoefficients) ;
          PTgrid[m][n][k] = RateCoefficients[k];
        } 
      }
    }

    //Calculate chebyshev coefficients. Three indicies are required in order calculate Chebyshev coefficients for each specified BW rate.
    vector<double> v3(m_reactionRef1.size());
    vector<vector<double> >v4(m_ExpanSizeC,v3);	
    vector<vector<vector<double> > > ChebyshevCoeff(m_ExpanSizeT,v4);
    for (int i=0; i < m_ExpanSizeT; ++i ){
      for(int j=0; j < m_ExpanSizeC ; ++j ){
        for(size_t k=0; k < m_reactionRef1.size(); ++k){
          for(int m = 0; m <m_NTpt ; ++m){
            for(int n = 0; n <m_NCpt ; ++n){
              ChebyshevCoeff[i][j][k] += log10(PTgrid[m][n][k])*Cheb_poly(i, TGrid[m])*Cheb_poly(j, CGrid[n]);
            }
          }				

          ChebyshevCoeff[i][j][k] *= Coefficient(i, j) / (double(m_NTpt) * double(m_NCpt)) ;
        }			
      }
    }

    //Print out table of Chebyshev coefficients for each BW rate specified
    for (size_t k=0; k < m_reactionRef1.size(); ++k ){
      ctest << setw(16) << "Chebeyshev Coefficients for reaction" << k << endl; 
      for (int i=0; i < m_ExpanSizeT; ++i ){
        for(int j=0; j < m_ExpanSizeC ; ++j ){
          ctest << setw(16) << ChebyshevCoeff[i][j][k];
        }
        ctest << endl;
      }
      ctest << endl;
    }

    for (size_t k=0; k < m_reactionRef1.size(); ++k ){
      ctest << setw(16) << "Comparison of fitted rate coefficients for reaction" << k << endl;

      for(int n = 0; n <m_NCpt ; ++n){
        ctest << setw(16) << Concentration[n];
      }
      ctest <<  endl;

      for(int m = 0; m <m_NTpt ; ++m){
        for(int n = 0; n <m_NCpt ; ++n){
          double ChebRate(0.);
          for (int i=0; i < m_ExpanSizeT; ++i ){
            for(int j=0; j < m_ExpanSizeC ; ++j ){
              ChebRate += ChebyshevCoeff[i][j][k]*Cheb_poly(i, TGrid[m])*Cheb_poly(j, CGrid[n]);
            }
          }	
          ChebRate = pow(10.0,ChebRate);
          ctest << setw(16) << Temperature[m];
          ctest << setw(16) << ChebRate <<"/" << PTgrid[m][n][k];

        }
        ctest << endl;
      }
      ctest << endl;
    }




    return true;
  }

  vector<double> AnalyticalRepresentation::Transform(vector<double> gridpoints, double max, double min ){
    vector<double> conditions(gridpoints.size());
    for(size_t i=0; i < gridpoints.size(); ++i ){
      conditions[i] = (gridpoints[i] * (max - min) + max + min) / 2.0;
    }
    return conditions;

  }

  double AnalyticalRepresentation::Coefficient( int i , int j ){
    double coeff=0.;

    if(i != 0 && j != 0){
      coeff = 4.0;
    }

    else if((i==0 && j!=0) || ( j == 0 && i != 0)){
      coeff = 2.0;
    }
    else if( i == 0 && j == 0 ){
      coeff = 1.0;
    }
    else{
      // this branch should never be executed
    }

    return coeff;
  }

  //returns the value of a chebeyshev polynomil of order i for a given argument value
  double AnalyticalRepresentation::Cheb_poly(int order, double arg){

    double value = cos(((order+1.0) - 1.0) * acos(arg));
    return value;	
  }

}//namespace


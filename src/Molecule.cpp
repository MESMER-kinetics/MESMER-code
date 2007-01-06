//-------------------------------------------------------------------------------------------
//
// Molecule.cpp 
//
// Author: Struan Robertson 
// Date:   5/Jan/2003
//
// This file contains the implementation of the Molecule class and its derived classes.
//
//-------------------------------------------------------------------------------------------

#include <iostream>
#include <sstream>
#include <math.h>
#include "Molecule.h"
#include "Constants.h"

using namespace std ;
using namespace Constants ;
namespace mesmer
{
    Molecule::Molecule():
m_Mass(0.0),
m_Sigma(0.0), 
m_Epsilon(0.0)
{ }

ModelledMolecule::ModelledMolecule():
m_MmtIntA(0.0), 
m_MmtIntB(0.0), 
m_MmtIntC(0.0), 
m_Sym(0.0), 
m_ZPE(0.0),
m_MaxCell(MAXCELL),
m_MaxGrn(MAXGRN),
m_GrnSz(0),
m_cdos(NULL), 
m_ecll(NULL) 
{ }

CollidingMolecule::CollidingMolecule():
m_egme(NULL), 
m_DeltaEdown(0.0) 
{ }

ModelledMolecule::~ModelledMolecule()
{
    // Free any memory assigned for calculating densities of states. 
    if (m_cdos != NULL) m_alloc.deallocate(m_cdos,MAXCELL) ; 
    if (m_ecll != NULL) m_alloc.deallocate(m_ecll,MAXCELL) ; 
}

CollidingMolecule::~CollidingMolecule()
{
    if (m_egme != NULL) delete m_egme ; 
}

TransitionState::TransitionState()
{ }

/* Will need Clone() functions
Molecule::Molecule(const Molecule& molecule) {
// Copy constructor - define later SHR 23/Feb/2003
}

Molecule& Molecule::operator=(const Molecule&) {
// Assignment operator - define later SHR 23/Feb/2003

return *this ;
}
*/

//
// Read the Molecular data from I/O stream.
//
bool Molecule::ReadFromXML(TiXmlElement* pnMol) 
{
    m_Name = pnMol->Attribute("id");
    if(m_Name.empty())
    {
        cerr << "A <cml:molecule> element has no valid id attribute";
        return false;
    }
    return true;
}

bool BathGasMolecule::ReadFromXML(TiXmlElement* pnMol) 
{
    //Read base class parameters first
    if(!Molecule::ReadFromXML(pnMol))
        return false;

    TiXmlElement* pnPropList = pnMol->FirstChildElement("propertyList");
    if(pnPropList)
    {
        const char* txt;

        txt= ReadProperty(pnPropList,"me:epsilon");
        if(!txt)
            return false;
        {
            istringstream idata(txt);
            idata >> m_Epsilon;
        }

        txt= ReadProperty(pnPropList,"me:sigma");
        if(!txt)
            return false;
        {
            istringstream idata(txt);
            idata >> m_Sigma;
        }

        txt= ReadProperty(pnPropList,"me:MW");
        if(!txt)
            return false;
        {
            istringstream idata(txt);
            idata >> m_Mass;
        }
    }
    return true;
}

bool ModelledMolecule::ReadFromXML(TiXmlElement* pnMol) 
{
    //Read base class parameters first
    if(!Molecule::ReadFromXML(pnMol))
        return false;

    TiXmlElement* pnPropList = pnMol->FirstChildElement("propertyList");
    if(pnPropList)
    {
        const char* txt;

        txt= ReadProperty(pnPropList,"me:ZPE");
        if(!txt)
            return false;
        {
            istringstream idata(txt);
            idata >> m_ZPE ;
        }
        txt= ReadProperty(pnPropList,"me:rotConsts");
        if(!txt)
            return false;
        {
            istringstream idata(txt);
            idata >> m_MmtIntA >> m_MmtIntB >> m_MmtIntC ;
        }
        txt= ReadProperty(pnPropList,"me:symmetryNumber");
        if(!txt)
            return false;
        {
            istringstream idata(txt);
            idata >> m_Sym;
        }

        txt= ReadProperty(pnPropList,"me:vibFreqs");
        if(!txt)
            return false;
        {
            istringstream idata(txt);
            double x ;
            while (idata >> x)
                m_VibFreq.push_back(x) ;
        }
    }

    //
    // Now the Molecular data is available, calculate density of states.
    //
    calcDensityOfStates() ;

    return true;
}

bool CollidingMolecule::ReadFromXML(TiXmlElement* pnMol) 
{
    //Read base class parameters first
    if(!ModelledMolecule::ReadFromXML(pnMol))
        return false;

    TiXmlElement* pnPropList = pnMol->FirstChildElement("propertyList");
    if(pnPropList)
    {
        const char* txt;

        txt= ReadProperty(pnPropList,"me:epsilon");
        if(!txt)
            return false;
        {
            istringstream idata(txt);
            idata >> m_Epsilon;
        }

        txt= ReadProperty(pnPropList,"me:sigma");
        if(!txt)
            return false;
        {
            istringstream idata(txt);
            idata >> m_Sigma;
        }

        txt= ReadProperty(pnPropList,"me:MW");
        if(!txt)
            return false;
        {
            istringstream idata(txt);
            idata >> m_Mass;
        }

        txt= ReadProperty(pnPropList,"me:deltaEDown");
        if(!txt)
            return false;
        {
            istringstream idata(txt);
            idata >> m_DeltaEdown;
        }
    }
    return true;
}
/*
m_VibFreq.clear() ;

int bracket_count = 1 ;                       // When zero end of molecule.

while (bracket_count > 0) {                   // Main input loop begins.

//
// Find first delinater character.
//
char c ;
while (in.get(c) && c != '@' && c != '}' ) ;

if (c == '}') {                                       // End of data. 
bracket_count-- ;
continue ;
}

string keyword ;
while (in.get(c) && c != '{' )
keyword += c ;

if        (keyword == "Name") {                    // Molecule name.

istringstream idata(collectData(in)) ;

idata >> m_Name ;

continue ;

} else if (keyword == "Rotational constants") {

istringstream idata(collectData(in)) ;

idata >> m_MmtIntA >> m_MmtIntB >> m_MmtIntC ;

continue ;

} else if (keyword == "Symmetry number") {

istringstream idata(collectData(in)) ;

idata >> m_Sym ;

continue ;

} else if (keyword == "Vibrational frequencies") {

istringstream idata(collectData(in)) ;

double x ;

while (idata >> x)
m_VibFreq.push_back(x) ;

} else if (keyword == "Delta Edown") {

istringstream idata(collectData(in)) ;

idata >> m_DeltaEdown ;

continue ;

} else if (keyword == "Epsilon") {

istringstream idata(collectData(in)) ;

idata >> m_Epsilon ;

continue ;

} else if (keyword == "Sigma") {

istringstream idata(collectData(in)) ;

idata >> m_Sigma ;

continue ;

} else if (keyword == "Mass") {

istringstream idata(collectData(in)) ;

idata >> m_Mass ;

continue ;

}

}                                             // Main input loop ends.
*/

//
// Get density of states.
//
void ModelledMolecule::densityOfStates(double *ddos) {

    // If density of states have not already been calcualted then do so.

    if (!m_cdos)
        calcDensityOfStates() ;

    for (int i = 0 ; i < MAXCELL ; ++i )
        ddos[i] = m_cdos[i] ;
}

//
// Get cell energies.
//
void ModelledMolecule::cellEnergies(double *decll) {

    // If energies have not already been calcualted then do so.

    if (!m_ecll)
        calcDensityOfStates() ;

    for (int i = 0 ; i < MAXCELL ; ++i )
        decll[i] = m_ecll[i] ;
}

//-------------------------------------------------------------------------------------------
//
// Methods related to the Master Equation.

//
// Initialize the Collision Operator.
//
void CollidingMolecule::initCollisionOperator(double beta) {

    //
    //     i) Determine Probabilities of Energy Transfer.
    //    ii) Normalisation of Probability matrix.
    //   iii) Symmetrise Collision Matrix.
    //
    int i, j ;

    double alpha = 1.0/m_DeltaEdown ;
    //
    // Allocate memmory.
    //
    if (m_egme)                                 // Delete any existing matrix.
        delete m_egme ;

    m_egme = new dMatrix(MAXGRN) ;              // Collision operator matrix.

    double *work = m_alloc.allocate(MAXGRN) ;   // Work space.
    //
    // Initialisation and error checking.
    //
    for ( i = 0 ; i < MAXGRN ; i++ ) {
        work[i] = 0.0 ;
        if (m_gdos[i] <= 0.0) {
            cout << "     ********* Error: Data indicates that      ************" << endl
                << "     ********* there is a grain with no states ************" << endl ;
            exit(1) ;
        }
    }
    //
    // The collision operator.
    //
    for ( i = 0 ; i < MAXGRN ; i++ ) {
        double ei = m_egrn[i] ;
        double ni = m_gdos[i] ;
        for ( j = i ; j < MAXGRN ; j++ ) {
            //
            // Transfer to lower Energy -
            //
            (*m_egme)[i][j] = exp(-alpha*(m_egrn[j] - ei)) ;
            //
            // Transfer to higher Energy (via detailed balance) -
            //
            (*m_egme)[j][i] = (*m_egme)[i][j] * (m_gdos[j]/ni)*exp(-beta*(m_egrn[j] - ei)) ;
        }
    }
    //
    // Normalization of Probability matrix.
    // Normalising coefficients are found by using the fact that column sums
    // are unity. The procedure leads to a matrix that is of upper triangular
    // form and the normalisation constants are found by back substitution.
    //
    for ( i = MAXGRN - 1 ; i > -1 ; i-- ) {

        double sm1(0.0) ;
        for ( j = 0 ; j <= i ; j++ )
            sm1 += (*m_egme)[j][i] ;

        double sm2(0.0) ;
        for ( j = i + 1 ; j < MAXGRN ; j++ )
            sm2 += (*m_egme)[j][i] * work[j] ;

        if (sm1 <= 0.0) {
            cout << "     ********* Error: Normalization Coeff.     ************" << endl
                << "     ********* in EGME is < or equal to zero.  ************" << endl  
                << "     i = " << i << " j = " << j                              << endl ;
            exit(1) ;
        }
        work[i] = (1.0-sm2)/sm1 ;
    }
    //
    // Apply normalization coefficients and account for collisional
    // loss by subrtacting unity from the leading diagonal.
    //
    cout << endl << "Normalization constants" << endl << endl ;
    for ( i = 0 ; i < MAXGRN ; i++ ) {
        cout <<  m_egrn[i] << " " << work[i] << endl ;
        (*m_egme)[i][i] *= work[i] ;
        (*m_egme)[i][i] -= 1.0 ;
        for ( j = i + 1 ; j < MAXGRN ; j++ ) {
            (*m_egme)[j][i] *= work[j] ;
            (*m_egme)[i][j] *= work[j] ;
        }
    }

    cout << endl << "Column Sums" << endl << endl ;
    for ( i = 0 ; i < MAXGRN ; i++ ) {
        double sum(0.0) ;
        for ( j = 0 ; j < MAXGRN ; j++ )  
            sum += (*m_egme)[j][i] ;

        cout << sum << endl ;
    }
    //
    // Determine the equilibrium vector for symmetrization. Note the work
    // array is over written and now contains square root of the Boltzman
    // distribution.
    //
    for ( i = 0 ; i < MAXGRN ; i++ ) {
        double ei = log(m_gdos[i]) - beta*m_egrn[i] + 10.0 ;
        ei = exp(ei) ;
        work[i] = sqrt(ei) ;
    }
    //
    // Symmeterisation of the collision matrix.
    //
    for ( i = 1 ; i < MAXGRN ; i++ ) {
        for ( j = 0 ; j < i ; j++ ) {
            (*m_egme)[j][i] *= (work[i]/work[j]) ;
            (*m_egme)[i][j]  = (*m_egme)[j][i] ;
        }
    }
    //
    // Deallocate memory.
    //
    if (work != NULL) m_alloc.deallocate(work,MAXGRN ) ; 

}

//
// Diagonalize the Collision Operator.
//
void CollidingMolecule::diagCollisionOperator() {

    // Allocate space for eigenvalues.

	int msize = m_egme->size() ;
    vector<double> rr(msize, 0.0) ;

    m_egme->diagonalize(&rr[0]) ;

    cout << endl ;
    for (int i(0) ; i < msize ; i++) 
        cout << rr[i] << endl; 

    cout << endl ;
    for (int i(0) ; i < msize ; i++)  {
        formatFloat(cout, m_egrn[i],              6, 15) ;
        formatFloat(cout, (*m_egme)[i][MAXGRN-1], 6, 15) ;
        formatFloat(cout, (*m_egme)[i][MAXGRN-2], 6, 15) ;
        cout << endl ;
    }

}

//
// Calculate collision frequency.
//
double CollidingMolecule::collisionFrequency(double temp, double conc, double bthMass, double bthSigma, double bthEpsilon) {

	//
	// Lennard-Jones Collision frequency. The collision integral is calculated 
	// using the formula of Neufeld et al., J.C.P. Vol. 57, Page 1100 (1972).
	// CONCentration is in molec/cm^3.
	//

	double A = 1.16145 ;
	double B = 0.14874 ;
	double C = 0.52487 ;
	double D = 0.77320 ;
	double E = 2.16178 ;
	double F = 2.43787 ;
	double amu = 1.6606E-27 ;

	double pi = acos(-1.0) ;
	//
	// Calculate collision parameter averages.
	//
	double mu   = amu*m_Mass*bthMass/(m_Mass + bthMass) ;
	double eam  = sqrt(m_Epsilon*bthEpsilon) ;
	double sam  = (m_Sigma + bthSigma)*0.5 ;
	double tstr = temp/eam ;
	//
	// Calculate collision integral.
	//
	double collFrq = A*exp(-log(tstr)*B) + C*exp(-D*tstr) + E*exp(-F*tstr) ;
	//
	// Calculate molecular collision frequency.
	//
	collFrq *= pi*sam*sam*1.e-20*sqrt(8.*1.381e-23*temp/(pi*mu)) ;
	//
	// Calculate overall collision frequency.
	//
	collFrq *= conc*1.e+06 ;

	return collFrq ;
}

//
// Calculate a reaction matrix element.
//
double CollidingMolecule::matrixElement(int eigveci, int eigvecj, vector<double> &k, int ndim) {

    double sum = 0.0 ;
    for (int i = 0 ; i < ndim ; i++)
        sum +=  k[i]*(*m_egme)[i][eigveci]*(*m_egme)[i][eigvecj] ;

    return sum ;
}

//
// Copy collision operator to diagonal block of system matrix.
//
void CollidingMolecule::copyCollisionOperator(dMatrix *CollOptr, const int size, const int locate) const {

    // Find size of system matrix.

    int smsize = CollOptr->size() ;

    // Check there is enough space in system matrix.

    if (locate + size > smsize) {
        cout << "Error in the size of the system matrix" << endl ;
        exit(1) ;
    }

    // Copy collision operator to the diagonal block indicated by "locate".

    for (int i(0) ; i < size ; i++) {
        int ii(locate + i) ;
        for (int j(0) ; j < size ; j++) {
            int jj(locate + j) ;
            (*CollOptr)[ii][jj] = (*m_egme)[i][j] ;
        }
    }

}



//-------------------------------------------------------------------------------------------
//
// Private methods.

//
// Calculate the rovibrational density of states.
//
void ModelledMolecule::calcDensityOfStates() {

    cout << endl << "The number of Frequencies is " << m_VibFreq.size() << endl ;

    m_cdos = m_alloc.allocate(MAXCELL) ; 
    m_ecll = m_alloc.allocate(MAXCELL) ; 

    //
    // Initialize density of states array using calculated rotational
    // density of state.
    //
    double cnt = sqrt(4./(m_MmtIntA * m_MmtIntB * m_MmtIntC))/m_Sym ;
    int i ;
    for ( i = 0 ; i < MAXCELL ; ++i ) {
        m_ecll[i] = static_cast<double>(i) + 0.5 ;
        m_cdos[i]  = cnt*sqrt(m_ecll[i]) ;
    }

    // Implementation of the Bayer-Swinehart algorithm.

    for ( vector<double>::size_type j = 0 ; j < m_VibFreq.size() ; ++j ) { 
        int iFreq = static_cast<int>(m_VibFreq[j]) ;
        for ( i = 0 ; i < MAXCELL - iFreq ; ++i ) 
            m_cdos[i + iFreq] += m_cdos[i] ;
    }

    // Calculate grain averages.

    calcGrainAverages() ;

    // Need to replace the following with XML output.
    //    if ( get_verbosity() ) 
    testDensityOfStates() ;

    return ;
}

//
// Test the rovibrational density of states.
//
void ModelledMolecule::testDensityOfStates() {

    cout << endl << "Test density of states: " << m_Name << endl << endl ;

    double temp ;
    double beta ;
    double pi = acos(-1.0) ;

    for ( int n = 0 ; n < 29 ; ++n ) {

        temp = 100.0*static_cast<double>(n + 2) ;
        beta = 1.0/(boltzmann*temp) ;

        // Calculate partition functions based on cells.

        double sumc  = 0.0 ;
        for ( int i = 0 ; i < MAXCELL ; ++i ) {
            sumc += m_cdos[i]*exp(-beta*m_ecll[i]) ;
        }

        // Calculate partition functions based on grains.

        double sumg  = 0.0 ;
        for ( int i = 0 ; i < MAXGRN ; ++i ) {
            sumg += m_gdos[i]*exp(-beta*m_egrn[i]) ;
        }

        // Calculate partition functions using analytical formula.

        double qtot = 1.0 ;
        for ( vector<double>::size_type j = 0 ; j < m_VibFreq.size() ; ++j ) { 
            qtot /= (1.0 - exp(-beta*m_VibFreq[j])) ;
        }
        qtot *= sqrt(pi/(m_MmtIntA * m_MmtIntB * m_MmtIntC))*(pow(beta,-1.5))/m_Sym ;
        formatFloat(cout, temp,  6,  7) ;
        formatFloat(cout, qtot,  6, 15) ;
        formatFloat(cout, sumc,  6, 15) ;
        formatFloat(cout, sumg,  6, 15) ;
        cout << endl ;

    }

}

//
// Calculate the average grain energy and then number of states per grain.
//
void ModelledMolecule::calcGrainAverages() {

    m_egrn.resize(MAXGRN) ;
    m_gdos.resize(MAXGRN) ;

    int igsz = MAXCELL/MAXGRN ;

    // Check that there are enough cells.

    if (igsz < 1) {
        cout << "     ********* Not enought Cells to produce ************" << endl
            << "     ********* requested number of Grains.  ************" << endl ;
        exit(1) ;
    }

    int idx1 = 0 ;
    int idx2 = 0 ;

    for (int i = 0 ; i < MAXGRN ; i++ ) {

        int idx3 = idx1 ;

        // Calculate the number of states in a grain.

        double smt = 0.0 ;
        for (int j = 0 ; j < igsz ; j++, idx1++ ) 
            smt += m_cdos[idx1] ;

        // Calculate average energy of the grain if it contains sum states.

        if ( smt > 0.0 ) {

            double smat = 0.0 ;
            for (int j = 0 ; j < igsz ; j++, idx3++ ) 
                smat += m_ecll[idx3] * m_cdos[idx3] ;

            m_gdos[idx2] = smt ;
            m_egrn[idx2] = smat/smt ;
            idx2++ ;
        }
    }

    // Issue warning if number of grains produced is less that requested.

    if ( idx2 < MAXGRN ) {
        cout <<  endl
            <<  "     WARNING: Number of grains produced is less than requested" << endl
            <<  "     Number of grains requested: " << MAXGRN << endl
            <<  "     Number of grains produced : " << idx2 << endl ;
    }
}
}//namespace

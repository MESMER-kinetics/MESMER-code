#ifndef GUARD_ThermodynamicUtils_h
#define GUARD_ThermodynamicUtils_h

//-------------------------------------------------------------------------------------------
//
// ThermodynamicUtils.h
//
// Author: Struan Robertson
// Date:   25/Nov/2018
//
// Definition of a utility class that is inherited by ThermodynamicTable and
// AnalyticalRrepresentation methods. 
//
//-------------------------------------------------------------------------------------------

namespace mesmer
{
  // Forward class declarations.
  class Molecule;

	class ThermodynamicUtils
	{
	public:

		~ThermodynamicUtils() {};

	protected:

    double SdivR(vector<double>::iterator i, double T) const;

    string WriteChemKinNASAPoly(Molecule* pmol, vector<double> coeffs, double TLo, double TMid, double THi);

    string WriteCanteraNASAPoly(Molecule* pmol, vector<double> coeffs, double TLo, double TMid, double THi);

  };

}  //namespace

#endif // GUARD_ThermodynamicUtils_h

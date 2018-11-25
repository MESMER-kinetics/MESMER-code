//-------------------------------------------------------------------------------------------
//
// ThermodynamicUtils.cpp
//
// Author: Struan Robertson
// Date:   25/Nov/2018
//
// This class implements the methods to be inherited by ThermodynamicTable and
// AnalyticalRrepresentation methods. 
//
//-------------------------------------------------------------------------------------------

#include "../System.h"
#include "../Molecule.h"
#include "../gStructure.h"
#include "ThermodynamicUtils.h"

namespace mesmer
{

  double ThermodynamicUtils::SdivR(vector<double>::iterator i, double T) const
  {
    //S/R  = a1 lnT + a2 T + a3 T^2 /2 + a4 T^3 /3 + a5 T^4 /4 + a7
    //return *i*log(T) + *(i+1)*T + *(i+2)*T*T / 2 + *(i+3)*T*T*T / 3 + *(i+4)*T*T*T*T / 4 + *(i+6);
    return *i*log(T) + T * (*(i + 1) + T * (*(i + 2) / 2 + T * (*(i + 3) / 3 + T * (*(i + 4) / 4)))) + *(i + 6);
  }

  string ThermodynamicUtils::WriteChemKinNASAPoly(Molecule* pmol, vector<double> coeffs, double TLo, double TMid, double THi)
  {
    stringstream ss;
    unsigned int i;

#if _MSC_VER && _MSC_VER<1900
    unsigned oldf = _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    ss << '\n';
    ss << left << setw(24) << pmol->getName().substr(0, 24);
    map<string, int> Comp = pmol->getStruc().GetElementalComposition();
    int npad = 4 - Comp.size();
    map<string, int>::const_iterator itr = Comp.begin();
    for (; itr != Comp.end(); itr++)
      ss << left << setw(2) << itr->first << right << setw(3) << itr->second;
    for (; npad; --npad)
      ss << "     ";
    ss << right << 'G' << fixed << setprecision(3) << setw(10) << TLo;
    ss << setw(10) << THi << setw(9) << TMid << "    01" << '\n';

    ss << scientific << setprecision(7);
    for (i = 0; i < 5; ++i)
      ss << setw(15) << coeffs[i];
    ss << "    2\n";
    for (i = 5; i < 10; ++i)
      ss << setw(15) << coeffs[i];
    ss << "    3\n";
    for (i = 10; i < 15; ++i)
      ss << setw(15) << coeffs[i];
    ss << "    4" << endl;

#if _MSC_VER && _MSC_VER<1900
    _set_output_format(oldf);
#endif

    return ss.str();
  }

  string ThermodynamicUtils::WriteCanteraNASAPoly(Molecule* pmol, vector<double> coeffs, double TLo, double TMid, double THi)
  {
    stringstream ss;
    unsigned int i;

#if _MSC_VER && _MSC_VER<1900
    unsigned oldf = _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    ss << endl;
    ss << "species(name = \"" << pmol->getName() << "\"," << endl;

    ss << "        atoms = \"";
    map<string, int> Comp = pmol->getStruc().GetElementalComposition();
    map<string, int>::const_iterator itr = Comp.begin();
    for (; itr != Comp.end(); itr++)
      ss << right << setw(2) << itr->first << ":" << left << setw(3) << itr->second;
    ss << "\"," << endl;

    ss << "        thermo = (" << endl;
    ss << "            NASA( [" << fixed << setprecision(3) << right << setw(10) << TLo << ", " << TMid << "], [";

    ss << scientific << setprecision(7);
    for (i = 0; i < 2; ++i)
      ss << setw(15) << coeffs[i] << ", ";
    ss << endl << setw(21) << " ";
    for (i = 2; i < 5; ++i)
      ss << setw(15) << coeffs[i] << ", ";
    ss << endl << setw(21) << " ";
    ss << setw(15) << coeffs[5] << ", " << setw(15) << coeffs[6] << "] )," << endl;

    ss << "            NASA( [" << fixed << setprecision(3) << right << setw(10) << TMid << ", " << THi << "], [";

    ss << scientific << setprecision(7);
    for (i = 7; i < 9; ++i)
      ss << setw(15) << coeffs[i] << ", ";
    ss << endl << setw(21) << " ";
    for (i = 9; i < 12; ++i)
      ss << setw(15) << coeffs[i] << ", ";
    ss << endl << setw(21) << " ";
    ss << setw(15) << coeffs[12] << ", " << setw(15) << coeffs[13] << "] ) ) )" << endl;

#if _MSC_VER && _MSC_VER<1900
    _set_output_format(oldf);
#endif

    return ss.str();
  }

}//namespace


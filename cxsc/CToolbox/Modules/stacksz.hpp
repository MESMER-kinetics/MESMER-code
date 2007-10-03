//============================================================================
//
//                              Program/Module
//                                   from
//                 C++ TOOLBOX FOR VERIFIED COMPUTING I
//                         Basic Numerical Problems
//
//      Copyright (c) 1995   Rolf Hammer, Matthias Hocks, Dietmar Ratz
//
// This program/module is free software for non-commercial use. For details
// on theory, algorithms, and programs, see the book
//
//  R. Hammer, M. Hocks, U. Kulisch, D. Ratz:  C++ Toolbox for
//  Verified Computing I - Basic Numerical Problems. Springer-Verlag,
//  Heidelberg, New York, 1995.
//
// This program/module is distributed WITHOUT ANY WARRANTY. For details,
// see the "Disclaimer / Legal Matters" of the book (page iv).
//
//============================================================================
//----------------------------------------------------------------------------
// File: stacksz (header)
// Purpose: Some C++ compilers do not allow to change the stack size for a
//    program by a compiler option. This module can be used instead. It must
//    be included with an #include-directive.
//----------------------------------------------------------------------------
#ifndef __STACKSZ_HPP
#define __STACKSZ_HPP


using namespace cxsc;
using namespace std;


// Borland C++ V3.1
//-----------------
// extern unsigned _stklen = 16384U;  // To get 32 KB stack size (default size
                                   // 4KB). Should be increased if a program
                                   // terminates abnormally without any error
                                   // message. Maximum stack size 64 KB!
#endif





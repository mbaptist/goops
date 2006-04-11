 // -*- C++ -*-

/////////////////////////////////////////////////////////////////////
// Interface to FFTW, using slate++ for scalar and vector fields.  //
// Manuel Baptista, 23th February 2004.                            //
// Last Modified 1st March 2004.                                   //
/////////////////////////////////////////////////////////////////////

#ifndef FFTI_H
#define FFTI_H


//scalar and  vector
//complex to complex and real to complex
//transforms
#include "fft.h"

//scalar real/complex transforms, heterogeneous along z
//Only defined for D=3
//in x,y s_rfft is used, in z sin or cos is used according
//to template flag 
#include "s_rfft_h_z.h"

//vector real/complex transforms, heterogeneous along z
//Only defined for D=3, N=3
#include "v_rfft_h_z_ccs.h"
#include "v_rfft_h_z_ssc.h"

#endif


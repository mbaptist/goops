 // -*- C++ -*-

/////////////////////////////////////////////////////////////////////
// Interface to FFTW, using slate++ for scalar and vector fields.  //
// Manuel Baptista, 23th February 2004.                            //
// Last Modified 1st March 2004.                                   //
/////////////////////////////////////////////////////////////////////

#ifndef FFT_H
#define FFT_H


//scalar complex/complex transforms
#include "s_fft.h"

//scalar complex/complex transforms, heterogeneous along z
//Only defined for D=3
//in x,y s_rfft is used, in z sin or cos is used according
//to template flag 
#include "s_fft_h_z.h"

//scalar real/complex transforms
#include "s_rfft.h"

//vector real/complex transforms
#include "v_rfft.h"

//scalar cosine transforms in one dimension
#include "s_cosfft_1d.h"

//scalar sine transforms in one dimension
#include "s_sinfft_1d.h"

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

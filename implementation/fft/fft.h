/*

Copyright 2004,2005,2006 Manuel Baptista

This file is part of GOOPS

GOOPS is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

GOOPS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

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
#include "fftclass.h"

#include "v_rfft.h"

//scalar real/complex transforms, heterogeneous along z
//Only defined for D=3
//in x,y s_rfft is used, in z sin or cos is used according
//to template flag 
#include "s_rfft_h_z.h"

//vector real/complex transforms, heterogeneous along z
//Only defined for D=3, N=3
#include "v_rfft_h_z.h"

#endif


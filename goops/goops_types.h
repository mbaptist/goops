// -*- C++ -*-
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



#ifndef GOOPS_TYPES_H
#define GOOPS_TYPES_H

#include <cat.h>

#if 1

//typedef float Real;//Real scalars
typedef double Real;//Real scalars
//typedef long double Real;//Real scalars
typedef std::complex<Real> Complex;//Complex scalars
typedef cat::Tvector<Real,3> RV;//Real vectors
typedef cat::Tvector<Complex,3> CV;//Complex Vectors
typedef cat::Array<Real,3> RSF;//Real Scalar Fields
typedef cat::Array<Complex,3> CSF;//Complex Scalar Fields
typedef cat::Array<RV,3> RVF;//Real Vector Fields
typedef cat::Array<CV,3> CVF;//Complex Vector Fields

#endif

#if 0
template <class Real,class ArrayRank,class VectorSize>
struct goops_types
{
typedef Real RT;///Real scalars
typedef std::complex<RT> CT;//Complex scalars
typedef cat::Tvector<RT,VectorSize> RV;//Real vectors
typedef cat::Tvector<CT,VectorSize> CV;//Complex Vectors
typedef cat::Array<RT,ArrayRank> RSF;//Real Scalar Fields
typedef cat::Array<CT,ArrayRank> CSF;//Complex Scalar Fields
typedef cat::Array<RV,ArrayRank> RVF;//Real Vector Fields
typedef cat::Array<CV,ArrayRank> CVF;//Complex Vector Fields
};
#endif

#endif


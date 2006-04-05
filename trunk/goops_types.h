// -*- C++ -*-

#ifndef GOOPS_TYPES_H
#define GOOPS_TYPES_H

#include <cat.h>

#if 1

//typedef float Real;//Real scalars
typedef double Real;//Real scalars
//typedef long double Real;//Real scalars
typedef std::complex<Real> Complex;//Complex scalars
typedef cat::tvector<Real,3> RV;//Real vectors
typedef cat::tvector<Complex,3> CV;//Complex Vectors
typedef cat::array<Real,3> RSF;//Real Scalar Fields
typedef cat::array<Complex,3> CSF;//Complex Scalar Fields
typedef cat::array<RV,3> RVF;//Real Vector Fields
typedef cat::array<CV,3> CVF;//Complex Vector Fields

#endif

#if 0
template <class Real,class ArrayRank,class VectorSize>
struct goops_types
{
typedef Real RT;///Real scalars
typedef std::complex<RT> CT;//Complex scalars
typedef cat::tvector<RT,VectorSize> RV;//Real vectors
typedef cat::tvector<CT,VectorSize> CV;//Complex Vectors
typedef cat::array<RT,ArrayRank> RSF;//Real Scalar Fields
typedef cat::array<CT,ArrayRank> CSF;//Complex Scalar Fields
typedef cat::array<RV,ArrayRank> RVF;//Real Vector Fields
typedef cat::array<CV,ArrayRank> CVF;//Complex Vector Fields
};
#endif

#endif

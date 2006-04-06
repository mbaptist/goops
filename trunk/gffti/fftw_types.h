// -*- C++ -*-

////////////////////////////////////////////////////
// fftw_types
//
// Provides a C++ interface to fftw types
//
///////////////////////////////////////////////////
#ifndef FFTW_TYPES_H
#define FFTW_TYPES_H

#include <complex>






template <class REAL>
struct FFTW_Precision
{
	typedef REAL Real;
	typedef	std::complex<Real> Complex;
	typedef	fftw_plan Plan;
}
template <>
struct FFTW_Precision<float>
{
	typedef float Real;
	typedef	std::complex<Real> Complex;
	typedef	fftw_plan Plan;
}




#endif


// -*- C++ -*-

#ifndef GFFTI_TYPES_H
#define GFFTI_TYPES_H

#include <cat.h>

template <class Real,class ArrayRank,class VectorSize>
struct GFFTI
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


template <class Real>
struct FFTW
{
	typedef Real RT;///Real scalars
	typedef std::complex<RT> CT;//Complex scalars
	typedef fftw_plan Plan;
};

template <class Real>
struct FFTW<float>
{
	typedef float RT;///Real scalars
	typedef std::complex<RT> CT;//Complex scalars
	typedef fftwf_plan Plan;
};

template <class Real>
struct FFTW<long double>
{
	typedef long double RT;///Real scalars
	typedef std::complex<RT> CT;//Complex scalars
	typedef fftwl_plan Plan;
};




#endif

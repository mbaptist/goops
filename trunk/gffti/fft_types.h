// -*- C++ -*-

////////////////////////////////////////////////////
// fft_types
//
// Provides a C++ interface to fftw types
//
///////////////////////////////////////////////////
#ifndef FFT_TYPES_H
#define FFT_TYPES_H


#include <cat.h>

#include <complex>

#include <fftw3.h>

namespace goops
{

typedef double RS;
typedef std::complex<RS> CS;
typedef fftw_plan Plan;

template <class RealType,class FourierType>
	struct FFT_TYPES;

//Complex to Complex
//Scalars Fields
template <int D>
	struct FFT_TYPES<cat::array<CS,D>,cat::array<CS,D> >
{
	typedef CS RDT;
	typedef CS FDT;
};
//Vector Fields
template <int D,int N>
	struct FFT_TYPES<cat::array<cat::tvector<CS,N>,D>,cat::array<cat::tvector<CS,N>,D> >
{
	typedef CS RDT;
	typedef CS FDT;
};

//Real to Complex / Complex to Real
//Scalars Fields
template <int D>
	struct FFT_TYPES<cat::array<RS,D>,cat::array<CS,D> >
{
	typedef RS RDT;
	typedef CS FDT;
};
//Vector Fields
template <int D,int N>
	struct FFT_TYPES<cat::array<cat::tvector<RS,N>,D>,cat::array<cat::tvector<CS,N>,D> >
{
	typedef RS RDT;
	typedef CS FDT;
};

//Real to Real
//Scalars Fields
template <int D>
	struct FFT_TYPES<cat::array<RS,D>,cat::array<RS,D> >
{
	typedef RS RDT;
	typedef RS FDT;
};
//Vector Fields
template <int D,int N>
	struct FFT_TYPES<cat::array<cat::tvector<RS,N>,D>,cat::array<cat::tvector<RS,N>,D> >
{
	typedef RS RDT;
	typedef RS FDT;
};

}

#endif


#if 0
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
	typedef	std::complex<float> Complex;
	typedef	fftw_plan Plan;
}
#endif
#if 0
#include <cat.h>
template <class Real,class ArrayRank,class VectorSize>
struct GFFTI
{
typedef Real RS;///Real scalars
typedef std::complex<RT> CV;//Complex scalars
typedef cat::tvector<RT,VectorSize> RV;//Real vectors
typedef cat::tvector<CT,VectorSize> CV;//Complex Vectors
typedef cat::array<RT,ArrayRank> RSF;//Real Scalar Fields
typedef cat::array<CT,ArrayRank> CSF;//Complex Scalar Fields
typedef cat::array<RV,ArrayRank> RVF;//Real Vector Fields
typedef cat::array<CV,ArrayRank> CVF;//Complex Vector Fields
};
#endif

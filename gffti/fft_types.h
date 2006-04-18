// -*- C++ -*-

////////////////////////////////////////////////////
// fft_types
//
// Provides a C++ interface to fftw types
//
///////////////////////////////////////////////////
#ifndef FFT_TYPES_H
#define FFT_TYPES_H
// namespace goops
// {

#include <cat.h>

#include <complex>

#include <fftw3.h>

typedef double RS;
typedef std::complex<RS> CS;
typedef fftw_plan fftPlan;

//Generic definition
template <class Type>
struct fftTypes;

//Specialization for constant types
template <class Type>
struct fftTypes<const Type>
{
	typedef  const typename fftTypes<Type>::NumericType NumericType;
	typedef const typename fftTypes<Type>::ElementType ElementType;
	typedef typename fftTypes<Type>::NumericType ccNumericType;
	static const int Rank=fftTypes<Type>::Rank;
	static const int VectorRank=fftTypes<Type>::VectorRank;
	typedef typename fftTypes<Type>::fftNumericType fftNumericType;
};

//Partial Specializations for non-constant types
//Constant versions are automatically defined by the code above

//Complex Scalars Fields
template <int D>
struct fftTypes<cat::array<CS,D> >
{
	typedef CS NumericType;
	typedef CS ElementType;
	static const int Rank=D;
	static const int VectorRank=1;
	typedef fftw_complex fftNumericType;
};

//Complex Vector Fields
template <int D,int N>
	struct fftTypes<cat::array<cat::tvector<CS,N>,D> >
{
	typedef CS NumericType;
	typedef cat::tvector<CS,N> ElementType;
	static const int Rank=D;
	static const int VectorRank=N;
	typedef fftw_complex fftNumericType;
};

//Real Scalars Fields
template <int D>
struct fftTypes<cat::array<RS,D> >
{
	typedef RS NumericType;
	typedef RS ElementType;
	static const int Rank=D;
	static const int VectorRank=1;
	typedef RS fftNumericType;
};

//Real Vector Fields
template <int D,int N>
struct fftTypes<cat::array<cat::tvector<RS,N>,D> >
{
	typedef RS NumericType;
	typedef cat::tvector<RS,N> ElementType;
	static const int Rank=D;
	static const int VectorRank=N;
	typedef RS fftNumericType;
};

//}

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



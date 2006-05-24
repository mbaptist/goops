//
// C++ Interface: fft
//
// Description:
//
//
// Author: Manuel Baptista <mbaptist@fc.up.pt>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef FFT_H
#define FFT_H

// namespace goops
// {

#include "fft_types.h"

#include "plan.h"

#include <cmath>

//Generic FFT class
template <class RealType,class FourierType>
class FFT
{
public:
	Plan<FourierType,const RealType> direct_plan;
	Plan<RealType,const FourierType> inverse_plan;
	//Ctors
	FFT();//Default Ctor
	//Dtor
	~FFT();
private:
	//Forbidden Ctors
	FFT(const FFT &);//Copy ctor
public:
	//Public methods
	void direct_transform(FourierType & fourierfield,const RealType & realfield);
	void inverse_transform(RealType & realfield,const FourierType & fourierfield);
	//FourierType direct_transform(const RealType & realfield);
	//RealType inverse_transform(const FourierType & fourierfield);
};

template <>
template <int D>
class FFT<cat::array<CS,D>,cat::array<CS,D> >
{
	typedef cat::array<CS,D> RealType;
	typedef cat::array<CS,D> FourierType;
public:
	Plan<FourierType,const RealType> direct_plan;
	Plan<RealType,const FourierType> inverse_plan;
	//Ctors
	FFT();//Default Ctor
	//Dtor
	~FFT();
private:
	//Forbidden Ctors
	FFT(const FFT &);//Copy ctor
public:
	//Public methods
	void direct_transform(FourierType & fourierfield,const RealType & realfield);
	void inverse_transform(RealType & realfield,const FourierType & fourierfield);
	//FourierType direct_transform(const RealType & realfield);
	//RealType inverse_transform(const FourierType & fourierfield);
};

template <>
template <int D>
class FFT<cat::array<RS,D>,cat::array<RS,D> >
{
	typedef cat::array<RS,D> RealType;
	typedef cat::array<RS,D> FourierType;
	typedef void (FFT::* DTT)(FourierType & fourierfield,const RealType & realfield);
	typedef void (FFT::* ITT)(RealType & realfield,const FourierType & fourierfield);
public:
	Plan<FourierType,const RealType> direct_plan;
	Plan<RealType,const FourierType> inverse_plan;
	//Ctors
	FFT(const string & subtype__);//Default Ctor
	//Dtor
	~FFT();
private:
	//Forbidden Ctors
	FFT();//Default Ctor
	FFT(const FFT &);//Copy ctor
public:
	//Public methods
	void direct_transform(FourierType & fourierfield,const RealType & realfield);
	void inverse_transform(RealType & realfield,const FourierType & fourierfield);
	//FourierType direct_transform(const RealType & realfield);
	//RealType inverse_transform(const FourierType & fourierfield);
private:
	const string subtype;
	fftw_r2r_kind r2r_direct_kind(const string & subtype__);
	fftw_r2r_kind r2r_inverse_kind(const string & subtype__);
	DTT direct_transform_imp;
	ITT inverse_transform_imp;
	void sin_direct_transform(FourierType & fourierfield,const RealType & realfield);
	void sin_inverse_transform(RealType & realfield,const FourierType & fourierfield);
	void cos_direct_transform(FourierType & fourierfield,const RealType & realfield);
	void cos_inverse_transform(RealType & realfield,const FourierType & fourierfield);
};

#include "fft.C"

#endif

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

template <class RealType,class FourierType>
	class FFT_BASE
{
protected:
	Plan <class FourierType,class RealType> direct_plan;
	Plan <class RealType,class FourierType> inverse_plan;
	//Ctors
	FFT_BASE();//Default Ctor
  FFT_BASE(const string & direction);
	FFT_BASE(const string & subtype,const string & direction);
	//Dtor
	virtual ~FFT_BASE();
private:
	//Forbidden Ctors
	FFT_BASE(const FFT_BASE &);//Copy ctor
public:
	//Public methods
	void direct_transform(FourierType & fourierfield,const RealType & realfield);
	void inverse_transform(RealType & realfield,const FourierType & fourierfield);
	//FourierType direct_transform(const RealType & realfield);
	//RealType inverse_transform(const FourierType & fourierfield);
};

template <class RealType,class FourierType>
class FFT:public FFT_BASE<RealType,FourierType>
{
	using FFT_BASE<RealType,FourierType>;
};

template <>
template<int D>
class FFT<cat::array<CS,D>,cat::array<CS,D> >:public FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >
{
	using FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >;
public:
	FFT(const string & direction);	
private:
	FFT();
	FFT(const FFT &);
};

template <>
template<int D>
class FFT<cat::array<RS,D>,cat::array<RS,D> >:public FFT_BASE<cat::array<RS,D>,cat::array<RS,D> >
{
	using FFT_BASE<cat::array<RS,D>,cat::array<RS,D> >;
public:
	FFT(const string & subtype,const string & direction);
private:
	FFT();
	FFT(const FFT &);
};

#include "fft.C"

#endif

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
class FFT
{
public:
	Plan<FourierType,const RealType> direct_plan;
	Plan<RealType,const FourierType> inverse_plan;
	//Ctors
	FFT();//Default Ctor
	FFT(const string &);
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






// template <class RealType,class FourierType>
// class FFT:public FFT_BASE<RealType,FourierType>
// {
// 	using FFT_BASE<RealType,FourierType>::direct_plan;
// 	using FFT_BASE<RealType,FourierType>::inverse_plan;
// };
// 
// template <>
// template<int D>
// class FFT<cat::array<CS,D>,cat::array<CS,D> >:public FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >
// {
// 	using FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >::direct_plan;
// 	using FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >::inverse_plan;
// public:
// 	FFT();	
// private:
// 	FFT(const FFT &);
// };
// 
// template <>
// class FFT<cat::array<RS,1>,cat::array<RS,1> >:public FFT_BASE<cat::array<RS,1>,cat::array<RS,1> >
// {
// 	using FFT_BASE<cat::array<RS,1>,cat::array<RS,1> >::direct_plan;
// 	using FFT_BASE<cat::array<RS,1>,cat::array<RS,1> >::inverse_plan;
// public:
// 	FFT(const string & subtype);
// private:
// 	FFT();
// 	FFT(const FFT &);
// public:
// 	void direct_transform(cat::array<RS,1> & fourierfield, const cat::array<RS,1> & realfield);
// 	void inverse_transform(cat::array<RS,1> & realfield,const cat::array<RS,1> & fourierfield);
// };

#include "fft.C"

#endif

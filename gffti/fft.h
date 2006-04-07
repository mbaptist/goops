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

namespace goops
{

#include "fft_types.h"

template <class RealType,class FourierType>
	class FFT
{
private:
	//Types
	typedef typename FFT_TYPES<RealType,FourierType>::RDT RDT;
	typedef typename FFT_TYPES<RealType,FourierType>::FDT FDT;
	typename FFT_TYPES<RealType,FourierType>::Plan Plan;
	//Members
	RDT * realdata;
	FDT * fourierdata;
	Plan direct_plan,inverse_plan;
	//Sizes
	size_t realsize;
	size_t fouriersize;
	
public:
	//Ctors
	FFT();//Default Ctor
	//Dtor
	~FFT();
	
private:
	//Forbidden Ctors
	FFT(const & FFT);//Copy ctor
	
private:
	//Private methods
	void create_plans();
	void destroy_plans();
	
public:
	//Public methods
	void direct_transform(FourierType & fourierfield,RealType & realfield);
	void inverse_transform(RealType & realfield,FourierType & fourierfield);
	//FourierType direct_transform(const RealType & realfield);
	//RealType inverse_transform(const FourierType & fourierfield);
};



template <int D>
	template <class RealType,class FourierType>
	class FFT<cat::array<RS,D>,cat::array<RS,D> > : public FFT<RealType,FourierType>
{
public:
	FFT(const string & subtype);
};

template <int D>
	template <class RealType,class FourierType>
	class FFT<cat::array<cat::tvector<RS,N>,D>,cat::array<cat::tvector<RS,N>,D> > : public FFT<RealType,FourierType>
{
public:
	FFT(const string & subtype);
};


}




#include "fft.C"

#endif
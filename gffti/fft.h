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

#include <fftw3.h>

#include <cmath>

template <class RealType,class FourierType>
	class FFT_BASE
{
protected:
	static const int DIM=1;	
	//Types
	typedef typename FFT_TYPES<RealType,FourierType>::RDT RDT;
	typedef typename FFT_TYPES<RealType,FourierType>::FDT FDT;
	//Members
	RDT * realdata;
	FDT * fourierdata;
	Plan direct_plan;
	Plan inverse_plan;
	//Size
	int realsize;
	int fouriersize;
	int * realshape;
	int * fouriershape;
	bool plans_created;
	
protected:
	//Ctors
	FFT_BASE();//Default Ctor
	//Dtor
	~FFT_BASE();
	
private:
	//Forbidden Ctors
	FFT_BASE(const FFT_BASE &);//Copy ctor
	
protected:
	//Private methods
	void create_plans();
	virtual void do_create_plans()=0;
	void destroy_plans();
	void switch_data(RealType & realfield,FourierType & fourierfield);
	int evalsize(int * shape);
public:
	//Public methods
	void direct_transform(FourierType & fourierfield,const RealType & realfield);
	void inverse_transform(RealType & realfield,const FourierType & fourierfield);
	//FourierType direct_transform(const RealType & realfield);
	//RealType inverse_transform(const FourierType & fourierfield);
};


template <class RealType,class FourierType>
	class FFT: public FFT_BASE<RealType,FourierType>
{
};

template <>
	template <int D>
	class FFT<cat::array<CS,D>,cat::array<CS,D> >: public FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >
{
	static const int DIM=D;
	using FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >::realdata;
	using FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >::fourierdata;
	using FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >::direct_plan;
	using FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >::inverse_plan;
	using FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >::realshape;
	using FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >::fouriershape;
	using FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >::realsize;
	using FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >::fouriersize;
	using FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >::plans_created;
	void do_create_plans();
};

template <>
	template <int D>
	class FFT<cat::array<RS,D>,cat::array<CS,D> >: public FFT_BASE<cat::array<RS,D>,cat::array<CS,D> >
{
	static const int DIM=D;
	using FFT_BASE<cat::array<RS,D>,cat::array<CS,D> >::realdata;
	using FFT_BASE<cat::array<RS,D>,cat::array<CS,D> >::fourierdata;
	using FFT_BASE<cat::array<RS,D>,cat::array<CS,D> >::direct_plan;
	using FFT_BASE<cat::array<RS,D>,cat::array<CS,D> >::inverse_plan;
	using FFT_BASE<cat::array<RS,D>,cat::array<CS,D> >::realshape;
	using FFT_BASE<cat::array<RS,D>,cat::array<CS,D> >::fouriershape;
	using FFT_BASE<cat::array<RS,D>,cat::array<CS,D> >::realsize;
	using FFT_BASE<cat::array<RS,D>,cat::array<CS,D> >::fouriersize;
	using FFT_BASE<cat::array<RS,D>,cat::array<CS,D> >::plans_created;
	void do_create_plans();
 };


template <>
	template <int D>
	class FFT<cat::array<RS,D>,cat::array<RS,D> >: public FFT_BASE<cat::array<RS,D>,cat::array<RS,D> >
{
private:
	static const int DIM=D;
	using FFT_BASE<cat::array<RS,D>,cat::array<RS,D> >::realdata;
	using FFT_BASE<cat::array<RS,D>,cat::array<RS,D> >::fourierdata;
	using FFT_BASE<cat::array<RS,D>,cat::array<RS,D> >::direct_plan;
	using FFT_BASE<cat::array<RS,D>,cat::array<RS,D> >::inverse_plan;
	using FFT_BASE<cat::array<RS,D>,cat::array<RS,D> >::realshape;
	using FFT_BASE<cat::array<RS,D>,cat::array<RS,D> >::fouriershape;
	using FFT_BASE<cat::array<RS,D>,cat::array<RS,D> >::realsize;
	using FFT_BASE<cat::array<RS,D>,cat::array<RS,D> >::fouriersize;
	using FFT_BASE<cat::array<RS,D>,cat::array<RS,D> >::plans_created;
	void do_create_plans();
	fftw_r2r_kind r2r_kind_direct;
	fftw_r2r_kind r2r_kind_inverse;
private:
	FFT();
public:
	FFT(const string & subtype);
	void direct_transform(cat::array<RS,D> & fourierfield,const cat::array<RS,D> & realfield);
	void inverse_transform(cat::array<RS,D> & realfield,const cat::array<RS,D> & fourierfield);
};


// }

#include "fft.C"

#endif

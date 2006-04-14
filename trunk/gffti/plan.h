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

#ifndef PLAN_H
#define PLAN_H

// namespace goops
// {
#include "fft_types.h"

#include <fftw3.h>

#include <cmath>

template <class TypeIn,class TypeOut>
	class PLAN_BASE
{
protected:
	//Types
	typedef typename FFT_TYPES<RealType>::RDT IDT;
	typedef typename FFT_TYPES<RealType,FourierType>::FDT ODT;
	//Members
	TypeInElement * datain;
	TypeOutElement * dataout;
	Plan plan;
	//Size
	int realsize;
	int fouriersize;
	int * realshape_ptr;
	int * fouriershapr_ptr;
	bool plan_exists;
	
protected:
	//Ctors
	PLAN_BASE();//Default Ctor
	//Dtor
	~PLAN_BASE();
	
private:
	//Forbidden Ctors
	PLAN_BASE(const PLAN_BASE &);//Copy ctor
	
protected:
	//Private methods
	void create_plans();
	virtual void do_create_plans()=0;
	void destroy_plans();
	void switch_data(RealType & realfield,FourierType & fourierfield);
	
public:
	//Public methods
	void switch_data(RealType & realfield,FourierType & fourierfield);
	void execute();
};


template <class RealType,class FourierType>
	class FFT: public FFT_BASE<RealType,FourierType>
{
};

template <>
	template <int D>
	class FFT<cat::array<CS,D>,cat::array<CS,D> >: public FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >
{
	using FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >::realdata;
	using FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >::fourierdata;
	using FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >::direct_plan;
	using FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >::inverse_plan;
	using FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >::realsize;
	using FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >::fouriersize;
	using FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >::plans_created;
	void do_create_plans();
};

template <>
	template <int D>
	class FFT<cat::array<RS,D>,cat::array<CS,D> >: public FFT_BASE<cat::array<RS,D>,cat::array<CS,D> >
{
	using FFT_BASE<cat::array<RS,D>,cat::array<CS,D> >::realdata;
	using FFT_BASE<cat::array<RS,D>,cat::array<CS,D> >::fourierdata;
	using FFT_BASE<cat::array<RS,D>,cat::array<CS,D> >::direct_plan;
	using FFT_BASE<cat::array<RS,D>,cat::array<CS,D> >::inverse_plan;
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
	using FFT_BASE<cat::array<RS,D>,cat::array<RS,D> >::realdata;
	using FFT_BASE<cat::array<RS,D>,cat::array<RS,D> >::fourierdata;
	using FFT_BASE<cat::array<RS,D>,cat::array<RS,D> >::direct_plan;
	using FFT_BASE<cat::array<RS,D>,cat::array<RS,D> >::inverse_plan;
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
	void direct_transform(cat::array<RS,D> & fourierfield, const cat::array<RS,D> & realfield);
	void inverse_transform(cat::array<RS,D> & realfield, const cat::array<RS,D> & fourierfield);
};


// }

#include "fft.C"

#endif

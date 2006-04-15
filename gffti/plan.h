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

template <class TypeOut,class TypeIn>
	class PlanBase
{
protected:
	//Types
	typedef typename FFT_TYPES<TypeIn>::TypeElement TypeInElement;
	typedef typename FFT_TYPES<TypeOut>::TypeElement TypeOutElement;
	//Members
	TypeInElement * datain;
	TypeOutElement * dataout;
	fftwPlan plan;
	//Size
	int datainsize;
	int dataoutsize;
	int * datainshape;
	int * dataoutshape;
	bool plan_exists;
	
protected:
	//Ctors
	PlanBase();//Default Ctor
	//Dtor
	~PlanBase();
	
private:
	//Forbidden Ctors
	PlanBase(const PlanBase &);//Copy ctor
	
protected:
	//Private methods
	void create_plan();
	virtual void do_create_plans()=0;
	void destroy_plans();
	void switch_data(RealType & realfield,FourierType & fourierfield);
	
public:
	//Public methods
	void switch_data(RealType & realfield,FourierType & fourierfield);
	void execute();
};


template <class RealType,class FourierType>
	class Plan: public FFT_BASE<RealType,FourierType>
{
};

template <>
	template <int D>
	class FFT<cat::array<CS,D>,cat::array<CS,D> >: public FFT_BASE<cat::array<CS,D>,cat::array<CS,D> >
{
	void do_create_plans();
};

template <>
	template <int D>
	class FFT<cat::array<RS,D>,cat::array<CS,D> >: public FFT_BASE<cat::array<RS,D>,cat::array<CS,D> >
{
	using FFT_BASE<cat::array<RS,D>,cat::array<CS,D> >::realdata;
	void do_create_plans();
 };


template <>
	template <int D>
	class FFT<cat::array<RS,D>,cat::array<RS,D> >: public FFT_BASE<cat::array<RS,D>,cat::array<RS,D> >
{
private:
	using FFT_BASE<cat::array<RS,D>,cat::array<RS,D> >::realdata;
	void do_create_plans();
	fftw_r2r_kind r2r_kind_direct;
	fftw_r2r_kind r2r_kind_inverse;
private:
	FFT();
public:
	FFT(const string & subtype);
};


// }

#include "fft.C"

#endif

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

//Class PlanBase
template <class TypeOut,class TypeIn>
	class PlanBase
{
protected:
	//Types
	typedef typename fftTypes<TypeIn>::fftNumericType fftTypeIn;
	typedef typename fftTypes<TypeOut>::fftNumericType fftTypeOut;
	//Members
	//output data
	fftTypeOut * dataout;
	const int * dataoutshape;
	int dataoutsize;
	//input data
	fftTypeIn * datain;
	const int * datainshape;
	int datainsize;
	//plan variables
	bool plan_exists;//Plan existence flag
	fftPlan plan;//Plan variable
protected:
	//Ctors
	PlanBase();//Default Ctor
	//Dtor
	virtual ~PlanBase();
private:
	//Forbidden Ctors
	PlanBase(const PlanBase &);//Copy ctor
protected:
	//Private methods
	void create_plan();
	virtual void do_create_plan()=0;
	void destroy_plan();
	int evalsize(const int * shape,int rank);
public:
	//Public methods
	void switch_data(TypeOut & fieldout,TypeIn & fieldin);
	void execute();
};

//Generic class Plan (inherits from class PlanBase)
template <class TypeOut,class TypeIn>
	class Plan: public PlanBase<TypeOut,TypeIn>
{
	typedef PlanBase<TypeOut,TypeIn> BaseClass;
	using BaseClass::dataout;
	using BaseClass::dataoutshape;
	using BaseClass::dataoutsize;
	using BaseClass::datain;
	using BaseClass::datainshape;
	using BaseClass::datainsize;
	using BaseClass::plan_exists;
	using BaseClass::plan;
	Plan();
};

//Specialisations of class Plan

//Complex to Complex transforms (both directions)
template <>
	template <int D>
	class Plan<cat::array<CS,D>,const cat::array<CS,D> >: public PlanBase<cat::array<CS,D>,const cat::array<CS,D> >
{
	typedef PlanBase<cat::array<CS,D>,const cat::array<CS,D> > BaseClass;
	using BaseClass::dataout;
	using BaseClass::dataoutshape;
	using BaseClass::dataoutsize;
	using BaseClass::datain;
	using BaseClass::datainshape;
	using BaseClass::datainsize;
	using BaseClass::plan_exists;
	using BaseClass::plan;
	const int direction;
public:
	Plan(const int & direction__);
private:
	Plan();//Forbids default ctor, since the direction of the transform must be specified
	void do_create_plan();
};

//Real to Complex transforms (direct transform)
template <>
template <int D>
class Plan<cat::array<CS,D>,const cat::array<RS,D> >: public PlanBase<cat::array<CS,D>,const cat::array<RS,D> >
{
private:
	typedef PlanBase<cat::array<CS,D>,const cat::array<RS,D> > BaseClass;
	using BaseClass::dataout;
	using BaseClass::dataoutshape;
	using BaseClass::dataoutsize;
	using BaseClass::datain;
	using BaseClass::datainshape;
	using BaseClass::datainsize;
	using BaseClass::plan_exists;
	using BaseClass::plan;
	void do_create_plan();
};

//Complex to Real transforms (inverse transform)
template <>
	template <int D>
	class Plan<cat::array<RS,D>,const cat::array<CS,D> >: public PlanBase<cat::array<RS,D>,const cat::array<CS,D> >
{
private:
	typedef PlanBase<cat::array<RS,D>,const cat::array<CS,D> > BaseClass;
	using BaseClass::dataout;
	using BaseClass::dataoutshape;
	using BaseClass::dataoutsize;
	using BaseClass::datain;
	using BaseClass::datainshape;
	using BaseClass::datainsize;
	using BaseClass::plan_exists;
	using BaseClass::plan;
	void do_create_plan();
 };

//1D Real to Real transforms (both directions for sine and co-sine transforms)
template <>
	class Plan<cat::array<RS,1>,const cat::array<RS,1> >: public PlanBase<cat::array<RS,1>,const cat::array<RS,1> >
{
private:
	typedef PlanBase<cat::array<RS,1>,const cat::array<RS,1> > BaseClass;
	using BaseClass::dataout;
	using BaseClass::dataoutshape;
	using BaseClass::dataoutsize;
	using BaseClass::datain;
	using BaseClass::datainshape;
	using BaseClass::datainsize;
	using BaseClass::plan_exists;
	using BaseClass::plan;
	fftw_r2r_kind r2r_kind;
public:
	Plan(const fftw_r2r_kind & r2r_kind__);
	~Plan();
private:
	Plan();
	void do_create_plan();
public:
	//void switch_data(cat::array<RS,1> & fieldout,const cat::array<RS,1> & fieldin);
};


// }

#include "plan.C"

#endif


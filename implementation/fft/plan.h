/*

Copyright 2004,2005,2006 Manuel Baptista

This file is part of GOOPS

GOOPS is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

GOOPS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

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
	void execute(TypeOut & fieldout,TypeIn & fieldin);
	virtual void guru_execute()=0;
};

//Generic class Plan (inherits from class PlanBase)
template <class TypeOut,class TypeIn>
	class Plan: public PlanBase<TypeOut,TypeIn>
{
};

//Specialisations of class Plan

//Complex to Complex transforms (both directions)
template <>
	template <int D>
	class Plan<cat::Array<CS,D>,const cat::Array<CS,D> >: public PlanBase<cat::Array<CS,D>,const cat::Array<CS,D> >
{
	typedef cat::Array<CS,D> TypeOut;
	typedef const cat::Array<CS,D> TypeIn;
	typedef PlanBase<cat::Array<CS,D>,const cat::Array<CS,D> > BaseClass;
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
public:
	void guru_execute();
};

//Real to Complex transforms (direct transform)
template <>
template <int D>
class Plan<cat::Array<CS,D>,const cat::Array<RS,D> >: public PlanBase<cat::Array<CS,D>,const cat::Array<RS,D> >
{
private:
	typedef cat::Array<CS,D> TypeOut;
	typedef const cat::Array<RS,D> TypeIn;
	typedef PlanBase<cat::Array<CS,D>,const cat::Array<RS,D> > BaseClass;
	using BaseClass::dataout;
	using BaseClass::dataoutshape;
	using BaseClass::dataoutsize;
	using BaseClass::datain;
	using BaseClass::datainshape;
	using BaseClass::datainsize;
	using BaseClass::plan_exists;
	using BaseClass::plan;
	void do_create_plan();
public:
	void guru_execute();
};

//Complex to Real transforms (inverse transform)
template <>
	template <int D>
	class Plan<cat::Array<RS,D>,const cat::Array<CS,D> >: public PlanBase<cat::Array<RS,D>,const cat::Array<CS,D> >
{
private:
	typedef cat::Array<RS,D> TypeOut;
	typedef const cat::Array<CS,D> TypeIn;
	typedef PlanBase<cat::Array<RS,D>,const cat::Array<CS,D> > BaseClass;
	using BaseClass::dataout;
	using BaseClass::dataoutshape;
	using BaseClass::dataoutsize;
	using BaseClass::datain;
	using BaseClass::datainshape;
	using BaseClass::datainsize;
	using BaseClass::plan_exists;
	using BaseClass::plan;
	void do_create_plan();
public:
	void guru_execute();
 };

//1D Real to Real transforms (both directions for sine and co-sine transforms)
template <>
template <int D>
	class Plan<cat::Array<RS,D>,const cat::Array<RS,D> >: public PlanBase<cat::Array<RS,D>,const cat::Array<RS,D> >
{
private:
	typedef cat::Array<RS,D> TypeOut;
	typedef const cat::Array<RS,D> TypeIn;
	typedef PlanBase<cat::Array<RS,D>,const cat::Array<RS,D> > BaseClass;
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
private:
	Plan();
	void do_create_plan();
public:
	void guru_execute();
};


// }

#include "plan.C"

#endif


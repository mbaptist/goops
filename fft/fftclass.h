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
class FFT<cat::Array<CS,D>,cat::Array<CS,D> >
{
	typedef cat::Array<CS,D> RealType;
	typedef cat::Array<CS,D> FourierType;
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
class FFT<cat::Array<RS,D>,cat::Array<RS,D> >
{
	typedef cat::Array<RS,D> RealType;
	typedef cat::Array<RS,D> FourierType;
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

#include "fftclass.C"

#endif

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
	//Ctors
	FFT_BASE();//Default Ctor
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

#include "fft.C"

#endif

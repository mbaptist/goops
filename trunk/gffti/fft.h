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

#include "fft_types.h"

template <class RealType,class FourierType>
class FFT
{
private:
	//Types	
	typedef typename FFT_TYPES<RealType,FourierType>::RDT RDT;
	typedef typename FFT_TYPES<RealType,FourierType>::FDT FDT;
	typename FFT_TYPES<RealType,FourierType>::Plan Plan;
	typedef typename FFT_TYPES<RealType,FourierType>::SizeT SizeT;
	//Members
	RDT * realdata;
	FDT * fourierdata;
	Plan direct_plan,inverse_plan;
	//Sizes
	ST realshape;
	ST fouriershape;

public:
	//Ctors
	FFT(typename FFT_TYPES<RealType,FourierType>::SizeT shape);//Ctor from size
	FFT(RealType & realfield,FourierType	& fourierfield);//Ctor from fields
	//Dtor
	~FFT();

private:
	//Forbidden Ctors
	FFT();//Default ctor
	FFT(const & FFT);//Copy ctor
	
private:
	//Private methods 	
	void set_data(RealType & realfield,FourierType	& fourierfield);
	
	void allocate_data();
	void free_data();	

	void create_plans();
	void destroy_plans();

public:
	//Public methods
	//void switchmode();
	
	
};
